#!/usr/bin/env python3
"""Validation and isolated-run support for the editable relaxation .com scripts.

Python standard library only. No spectrum processing occurs on import.
"""
import argparse
from collections import Counter
from contextlib import contextmanager
import csv
from decimal import Decimal, InvalidOperation
import fcntl
import glob
import math
import os
from pathlib import Path
import re
import shutil
import signal
import struct
import subprocess
import sys
import tempfile
import textwrap

NUMBER = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'
TIME = re.compile(rf'({NUMBER})\s*(s|m|ms|u|us|µs|μs)?\Z')


def finite_number(value, description, minimum=None):
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f'{description}: expected a number, got {value!r}') from error
    if not math.isfinite(result) or (minimum is not None and result < minimum):
        raise ValueError(f'{description}: invalid value {value!r}')
    return result


def meaningful_lines(path):
    return [line.strip() for line in Path(path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith('#')]


def parse_delay(value):
    match = TIME.fullmatch(value.strip())
    if not match:
        raise ValueError(f'Invalid delay {value!r}; use seconds, s, m/ms or u/us/µs')
    scale = {None: 1.0, 's': 1.0, 'm': .001, 'ms': .001,
             'u': .000001, 'us': .000001, 'µs': .000001, 'μs': .000001}[match[2]]
    return finite_number(float(match[1]) * scale, 'delay in seconds', 0)


def create_delays(base, vc, vd, ncyc, source='auto', block=.01, constant=.01,
                  cycle_seconds='', cycle_constant=0):
    choices = {}
    for kind, pattern in (('vc', vc), ('vd', vd), ('ncyc', ncyc)):
        path = Path(pattern)
        matches = sorted(Path(p) for p in glob.glob(str(path if path.is_absolute() else Path(base) / path))
                         if Path(p).is_file())
        if source == kind or source == 'auto':
            if len(matches) > 1:
                raise ValueError(f'Multiple {kind} lists match {pattern!r}; configure one exact path')
            if matches:
                choices[kind] = matches[0]
    if len(choices) != 1:
        raise ValueError('Select exactly one existing delay source (vc, vd or ncyc); configure delay_source and its path')
    kind, path = next(iter(choices.items()))
    values = meaningful_lines(path)
    if not values:
        raise ValueError(f'{path}: empty delay list')
    if kind == 'vd':
        delays = [parse_delay(value) for value in values]
    else:
        step = finite_number(block if kind == 'vc' else cycle_seconds, 'cycle/block duration', 0)
        offset = finite_number(constant if kind == 'vc' else cycle_constant, 'constant delay', 0)
        if step == 0:
            raise ValueError('Cycle/block duration must be positive')
        counts = [finite_number(value, 'cycle count', 0) for value in values]
        if any(not count.is_integer() for count in counts):
            raise ValueError('Cycle counts must be integers')
        delays = [finite_number(offset + count * step, 'converted delay', 0) for count in counts]
    Path('delaylist').write_text(''.join(f'{v:.17g}\n' for v in delays))
    Path('relax_header.txt').write_text(f'{len(delays)} ' + ' '.join(f'{v:.17g}' for v in delays) + '\n')
    return delays


def counted_vector(path, positive=False, nonnegative=True):
    lines = meaningful_lines(path)
    if len(lines) != 1:
        raise ValueError(f'{path}: expected one counted vector')
    fields = lines[0].split()
    if not fields or not fields[0].isdigit() or int(fields[0]) < 1:
        raise ValueError(f'{path}: invalid count')
    n = int(fields[0])
    if len(fields) != n + 1:
        raise ValueError(f'{path}: count {n} differs from {len(fields) - 1} values')
    result = [finite_number(v, str(path), 0 if nonnegative else None) for v in fields[1:]]
    if positive and any(v <= 0 for v in result):
        raise ValueError(f'{path}: noise SDs must be positive')
    return result


def nonempty(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f'Missing or empty output: {path}')


def pipe_header(path):
    """Read just the standard 512-float NMRPipe header, in either byte order."""
    nonempty(path)
    if Path(path).stat().st_size <= 2048:
        raise ValueError(f'{path}: missing spectrum payload')
    with Path(path).open('rb') as handle:
        raw = handle.read(2048)
    if len(raw) != 2048:
        raise ValueError(f'{path}: incomplete NMRPipe header')
    for endian in ('<', '>'):
        values = struct.unpack(endian + '512f', raw)
        if math.isclose(values[2], 2.345, abs_tol=1e-5):
            return values
    raise ValueError(f'{path}: invalid NMRPipe byte-order marker')


def integer_header(value, name):
    if not math.isfinite(value) or value < 1 or not value.is_integer():
        raise ValueError(f'Invalid NMRPipe {name}: {value}')
    return int(value)


def plane_count(path):
    header = pipe_header(path)
    if header[9] != 3 or header[26] != 3:
        raise ValueError(f'{path}: expected a 3D NMRPipe series with the plane axis third')
    return integer_header(header[15], 'FDF3SIZE')


def check_fid(path, delays):
    header = pipe_header(path)
    if header[9] != 2:
        raise ValueError(f'{path}: expected a 2D interleaved FID')
    vectors = integer_header(header[219], 'FDSPECNUM')
    n = len(meaningful_lines(delays))
    if not n or vectors % n:
        raise ValueError('Indirect vector count is not divisible by the number of delays')


def check_conversion(delays):
    """Accept the original combined-FID or already-split conversion routes."""
    combined = Path('test.fid')
    planes = list(Path('fid').glob('test*.fid'))
    if combined.exists() and planes:
        raise ValueError('Conversion produced both test.fid and FID planes; select one output layout')
    if combined.exists():
        check_fid(combined, delays)
    elif planes:
        check_fid_series('fid', delays)
    else:
        raise ValueError('Conversion produced neither test.fid nor fid/test*.fid')


def check_fid_series(directory, delays):
    """Validate converted H/N planes with a real third (delay) axis.

    Conversion must already arrange the axes and handle gradient shuffling.
    This function only checks headers, file sizes and the expected plane names.
    """
    n = len(meaningful_lines(delays))
    if not n:
        raise ValueError('Empty delay list')
    directory = Path(directory)
    paths = sorted(directory.glob('test*.fid'))
    expected = [directory / f'test{i:03d}.fid' for i in range(1, n + 1)]
    if paths != expected:
        raise ValueError(f'{directory}: converted FID planes do not match the delays')
    grid = None
    for path in paths:
        header = pipe_header(path)
        if header[9] != 3 or tuple(header[24:28]) != (2, 1, 3, 4) or header[57] != 0:
            raise ValueError(f'{path}: expected an XYZ plane series with the delay axis third')
        if header[51] != 1 or header[15] != n or header[442] != n:
            raise ValueError(f'{path}: expected a real delay axis with {n} planes')
        if header[220] != 0 or header[222] != 0:
            raise ValueError(f'{path}: expected untransformed FID axes')
        if header[56] not in (0, 1) or header[55] not in (0, 1):
            raise ValueError(f'{path}: invalid quadrature flags')
        # FDSPECNUM counts the stored rows; FDSIZE counts X real or imaginary
        # points, so only the X complex flag adds a payload factor here.
        shape = (integer_header(header[99], 'FDSIZE'), integer_header(header[219], 'FDSPECNUM'))
        signature = shape + (header[56], header[55])
        if grid is not None and signature != grid:
            raise ValueError(f'{directory}: inconsistent converted FID grids')
        grid = signature
        expected_size = 2048 + 4 * shape[0] * shape[1] * (2 if header[56] == 0 else 1)
        if path.stat().st_size != expected_size:
            raise ValueError(f'{path}: FID payload size differs from its header')


def noise_value(path):
    text = Path(path).read_text()
    matches = re.findall(rf'\bvEstNoise\b\s*(?:=|:)?\s*({NUMBER})(?![\w.])', text)
    if len(matches) != 1:
        raise ValueError(f'{path}: expected one vEstNoise value')
    noise = finite_number(matches[0], 'estimated noise SD', 0)
    if noise <= 0:
        raise ValueError('Estimated noise SD must be positive')
    return noise


def check_processing():
    delays = [finite_number(v, 'delay', 0) for v in meaningful_lines('delaylist')]
    if not delays or counted_vector('relax_header.txt') != delays:
        raise ValueError('Delay header does not match delaylist')
    if len(counted_vector('vEstNoise.txt', positive=True)) != len(delays):
        raise ValueError('Noise and delay counts differ')
    if plane_count('m.ft3') != len(delays):
        raise ValueError('Spectrum plane and delay counts differ')
    for directory, suffix in (('fid', 'fid'), ('ft', 'ft2')):
        paths = sorted(Path(directory).glob(f'test*.{suffix}'))
        expected = [Path(directory) / f'test{i:03d}.{suffix}' for i in range(1, len(delays) + 1)]
        if paths != expected:
            raise ValueError(f'{directory}: plane files do not match the delays')
        grid = None
        for path in paths:
            header = pipe_header(path)
            shape = (integer_header(header[99], 'FDSIZE'), integer_header(header[219], 'FDSPECNUM'))
            if grid is not None and shape != grid:
                raise ValueError(f'{directory}: inconsistent plane grids')
            grid = shape
            expected_dimension = 3 if directory == 'fid' else 2
            expected_planes = len(delays) if directory == 'fid' else 1
            if (header[9] != expected_dimension or header[15] != expected_planes or
                    header[442] != expected_planes):
                raise ValueError(f'{path}: inconsistent plane/series header')
    nonempty('proc.ucsf')


@contextmanager
def run_directory(source, kind):
    source = Path(source).resolve()
    # One lock for both scripts also prevents analysis during publication of new spectra.
    with (source / '.relax.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise ValueError('Another relaxation run is active in this folder') from error
        with tempfile.TemporaryDirectory(prefix=f'.relax-{kind}-', dir=source) as temporary:
            yield Path(temporary)


def publish(stage, destination, names):
    """Replace owned outputs after validation; restore the previous set on failure."""
    destination = Path(destination)
    backup = Path(stage) / '.previous'
    backup.mkdir()
    names = list(names)
    if len(set(names)) != len(names) or any(Path(n).name != n for n in names):
        raise ValueError('Invalid output names')
    for name in names:
        path = Path(stage) / name
        if not path.exists() or path.is_symlink():
            raise ValueError(f'Missing generated output: {name}')
        target = destination / name
        if target.is_symlink():
            raise ValueError(f'Refusing to replace a symlink output: {target}')
    replaced, saved = [], []
    try:
        for name in names:
            target = destination / name
            if target.exists():
                target.rename(backup / name)
                saved.append(name)
            (Path(stage) / name).rename(target)
            replaced.append(name)
    except BaseException:
        for name in reversed(replaced):
            (destination / name).rename(Path(stage) / name)
        for name in reversed(saved):
            (backup / name).rename(destination / name)
        raise


def run_processing(script, source):
    source, script = Path(source).resolve(), Path(script).resolve()
    nonempty(source / 'ser')
    with run_directory(source, 'proc') as stage:
        (stage / 'ser').symlink_to(source / 'ser')
        shell = shutil.which('tcsh')
        if shell is None:
            raise ValueError('tcsh is required')
        worker_env = dict(os.environ, FTNMR_RELAX_STAGE=str(stage))
        subprocess.run([shell, '-f', str(script), '--worker', str(source)],
                       cwd=stage, env=worker_env, check=True)
        # Validate again in the controller before publishing anything.
        old = Path.cwd()
        try:
            os.chdir(stage)
            check_processing()
        finally:
            os.chdir(old)
        names = ['delaylist', 'relax_header.txt', 'fid', 'ft', 'm.ft3', 'proc.ucsf', 'vEstNoise.txt']
        if (stage / 'series.tab').is_file():
            names.append('series.tab')
        publish(stage, source, names)



def set_series_metadata(paths):
    """The installed series.com settings, with every sethdr exit status checked."""
    n = len(paths)
    if not n:
        raise ValueError('No series files')
    for path in paths:
        nonempty(path)
        subprocess.run(['sethdr', str(path), '-pipeFlag', '0', '-cubeFlag', '0',
                        '-fileCount', str(n), '-zN', str(n), '-zT', str(n), '-zTD', str(n),
                        '-zMODE', 'Real', '-zFT', 'Freq', '-zSW', str(n), '-zOBS', '1.0',
                        '-zCAR', '0.0', '-zLAB', 'ID', '-ndim', '3'], check=True)


def read_peaks(path, peak_format):
    lines = meaningful_lines(path)
    if not lines:
        raise ValueError('Empty peak list')
    header = lines[0].split()
    columns = {name.lower(): i for i, name in enumerate(header)}
    dimensions = 3 if peak_format == 'sparky3d' else 2
    required = ['assignment', 'overlap'] + [f'w{i}' for i in range(1, dimensions + 1)]
    if not set(required).issubset(columns):
        raise ValueError(f'{path}: missing assignment, coordinate or overlap columns')
    peaks, seen = [], set()
    for index, line in enumerate(lines[1:], 1):
        fields = line.split()
        if len(fields) <= max(columns[name] for name in required):
            raise ValueError(f'{path}: incomplete peak row {index}')
        assignment = fields[columns['assignment']]
        coords = [fields[columns[f'w{i}']] for i in range(1, dimensions + 1)]
        try:
            numeric = tuple(Decimal(value) for value in coords)
        except InvalidOperation as error:
            raise ValueError(f'{path}: invalid coordinates at row {index}') from error
        if not all(v.is_finite() and v >= 0 for v in numeric):
            raise ValueError(f'{path}: coordinates must be finite and nonnegative for data_xtr')
        key = (assignment, numeric)
        if key in seen:
            raise ValueError(f'{path}: duplicate assignment AND coordinates at row {index}')
        seen.add(key)
        overlap = fields[columns['overlap']]
        if overlap not in ('True', 'False'):
            raise ValueError(f'{path}: invalid overlap flag at row {index}')
        residues = re.findall(r'(?:^|-)[A-Za-z]{1,3}(\d+)[A-Za-z0-9]*', assignment)
        if not residues or len(set(residues)) != 1:
            raise ValueError(f'Cannot identify one residue in assignment {assignment!r}')
        # These tags obey the legacy sparky2d/3d regular expressions.
        tag = f'FT{index:06d}-N-H' if dimensions == 3 else f'FT{index:06d}-H'
        peaks.append({'peak_id': tag, 'assignment': assignment, 'coordinates': coords,
                      'residue': int(residues[0]), 'overlap': overlap == 'True'})
    if not peaks:
        raise ValueError('Peak list has no peaks')
    counts = Counter(p['assignment'] for p in peaks)
    for peak in peaks:
        peak['duplicate_assignment'] = counts[peak['assignment']] > 1
    return peaks


def parse_extracted(text, peaks, n):
    expected = {p['peak_id'] for p in peaks}
    result = {}
    for line in text.splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        fields = line.split()
        tag = fields[0]
        if tag not in expected or tag in result or len(fields) != n + 1:
            raise ValueError('Extracted peak identifiers or intensity counts do not match the inputs')
        result[tag] = [finite_number(v, f'intensity for {tag}') for v in fields[1:]]
    if set(result) != expected:
        raise ValueError('Extracted peak identifiers are missing')
    return result


def parse_fit(text, mode):
    records, warnings, global_warnings = {}, {}, []
    width = 2 if mode == 'noe' else 6
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith('#'):
            if 'warning' in stripped.lower():
                site = re.search(r'\bsite\s+(\S+)\s*$', stripped)
                if site:
                    warnings.setdefault(site[1], []).append(line)
                else:
                    global_warnings.append(line)
            continue
        fields = stripped.split()
        if len(fields) != width + 1 or fields[0] in records:
            raise ValueError('Malformed or duplicate fitted peak record')
        try:
            values = list(map(float, fields[1:]))
        except ValueError as error:
            raise ValueError('Malformed fitted numeric values') from error
        records[fields[0]] = values
    if set(warnings) - set(records):
        raise ValueError('Fitter warnings name unexpected peak identifiers')
    return records, {tag: warnings.get(tag, []) + global_warnings for tag in records}


def gnuplot_string(text):
    return '"' + str(text).replace('\\', '\\\\').replace('"', '\\"').replace('\n', '\\n').replace('\r', ' ') + '"'


def prepare_reports(stage, peaks, extracted, fit_text, times, noises, mode, prefix, pdf_name, run_tag):
    records, warnings = parse_fit(fit_text, mode)
    if set(records) != {p['peak_id'] for p in peaks}:
        raise ValueError('Missing or unexpected fitted peak identifiers')
    rows = []
    for peak in peaks:
        values = records[peak['peak_id']]
        flags = []
        valid = all(math.isfinite(v) for v in values) and values[1] >= 0
        if mode != 'noe':
            valid = valid and values[3] >= 0 and values[4] >= 0 and values[5] >= 0 and values[2] != 0
            if math.isfinite(values[0]) and values[0] <= 0:
                flags.append('nonpositive_rate')
            if values[4] > values[5]:
                flags.append('chi2_exceeds_threshold')
        if not valid:
            flags.append('invalid_fit')
        if peak['overlap']:
            flags.append('overlap')
        if peak['duplicate_assignment']:
            flags.append('duplicate_assignment')
        warning_text = '\n'.join(warnings[peak['peak_id']])
        if warning_text:
            flags.append('fitter_warning')
            if 'uncertainties scaled' in warning_text.lower():
                flags.append('uncertainty_rescaled')
        rows.append(dict(peak, values=values, valid=valid, flags=flags, warnings=warning_text))
    if not any(row['valid'] for row in rows):
        raise ValueError('No valid fitted peaks; previous results retained')
    tag_header = f'# {run_tag}\n'
    names = []
    for suffix, include in (('results.out', lambda r: True),
                            ('results_peak_overlap.out', lambda r: r['overlap']),
                            ('results_filtered.out', lambda r: not r['overlap'] and r['valid'])):
        name = f'{prefix}_{suffix}'
        text = tag_header + ''.join(f"{r['residue']} {r['values'][0]:.12g} {r['values'][1]:.12g}\n" for r in rows if include(r))
        (stage / name).write_text(text)
        names.append(name)
    diagnostics = f'{prefix}_fit_diagnostics.tsv'
    fields = ['peak_id', 'assignment', 'w1', 'w2', 'w3', 'residue', 'value', 'error', 'i0', 'i0_error',
              'chi2', 'chi2_threshold', 'overlap', 'duplicate_assignment', 'valid_fit', 'flags', 'warnings']
    with (stage / diagnostics).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for row in rows:
            v = row['values']
            info = {key: row[key] for key in ('peak_id', 'assignment', 'residue', 'overlap', 'duplicate_assignment', 'warnings')}
            info.update({f'w{i}': c for i, c in enumerate(row['coordinates'], 1)})
            info.update(value=v[0], error=v[1], i0=v[2] if mode != 'noe' else '',
                        i0_error=v[3] if mode != 'noe' else '', chi2=v[4] if mode != 'noe' else '',
                        chi2_threshold=v[5] if mode != 'noe' else '', valid_fit=row['valid'], flags=';'.join(row['flags']))
            writer.writerow(info)
    names.append(diagnostics)
    plot = [tag_header.rstrip(), 'set terminal pdfcairo noenhanced color size 6,5 font "Sans,9"',
            'set output ' + gnuplot_string(pdf_name), 'set grid', 'set key off', 'set autoscale']
    blocks = []
    for index, row in enumerate(rows):
        values = extracted[row['peak_id']]
        block = [f"# {row['assignment']} peak_id={row['peak_id']} index {index}"]
        if mode == 'noe':
            block += [f'{i} {y:.17g} {sd:.17g}' for i, (y, sd) in enumerate(zip(values, noises), 1)]
        else:
            rate, _, i0 = row['values'][:3]
            normalized = math.isfinite(i0) and i0 != 0 and row['valid']
            for t, y, sd in zip(times, values, noises):
                if normalized:
                    try:
                        prediction = i0 * math.exp(-rate * t)
                    except OverflowError as error:
                        raise ValueError('Fitted prediction overflows') from error
                    if not math.isfinite(prediction):
                        raise ValueError('Nonfinite fitted prediction')
                    block.append(f'{t:.17g} {y/i0:.17g} {sd/abs(i0):.17g} '
                                 f'{prediction/i0:.17g} {(y-prediction)/i0:.17g} {(y-prediction)/sd:.17g}')
                else:
                    block.append(f'{t:.17g} {y:.17g} {sd:.17g} NaN NaN NaN')
            label = row['assignment'] + ' [' + row['peak_id'] + ']'
            if row['flags']:
                label += '\n' + textwrap.fill(', '.join(row['flags']), width=76)
            span = max(times) - min(times) or max(abs(max(times)), .01)
            lo, hi = min(times) - .03 * span, max(times) + .03 * span
            plot += [f'set xrange [{lo:.17g}:{hi:.17g}]', 'set multiplot layout 2,1 title ' + gnuplot_string(label),
                     'set xlabel "Delay (s)"', 'set ylabel ' + gnuplot_string('Intensity / fitted I0' if normalized else 'Raw intensity (invalid fit/I0)')]
            command = f'plot "intensity.out" index {index} using 1:2:3 with yerrorbars pt 7'
            if normalized:
                command += f', exp(-({rate:.17g})*x) with lines'
            plot += [command, 'set ylabel "Residual / input noise SD"']
            if normalized:
                plot.append(f'plot "intensity.out" index {index} using 1:6 with points pt 7, 0 with lines')
            else:
                plot += ['set label 999 "Residual unavailable: invalid fit/I0" at graph 0.1,0.5',
                         'plot [0:1] [-1:1] 0 with lines', 'unset label 999']
            plot.append('unset multiplot')
        blocks.append('\n'.join(block))
    (stage / 'intensity.out').write_text(tag_header + '\n\n\n'.join(blocks) + '\n')
    plot += ['set autoscale x', 'set title ' + gnuplot_string(mode.upper() + ' | ' + run_tag), 'set xlabel "Residue"',
             'set ylabel ' + gnuplot_string('hetNOE' if mode == 'noe' else mode.upper() + ' (1/s)'),
             'set key top right', 'set datafile missing "NaN"']
    groups = [(False, 'Isolated peaks', '#285879'), (True, 'Overlap flagged', '#c17b19')]
    groups = [(flag, label, color) for flag, label, color in groups
              if any(row['valid'] and row['overlap'] == flag for row in rows)]
    plot.append('plot ' + ', '.join(
        '\'-\' using 1:2:3 with yerrorbars pt 7 lc rgb ' + gnuplot_string(color) + ' title ' + gnuplot_string(label)
        for flag, label, color in groups))
    for flag, label, color in groups:
        plot.extend(f"{row['residue']} {row['values'][0]:.17g} {row['values'][1]:.17g}"
                    for row in rows if row['valid'] and row['overlap'] == flag)
        plot.append('e')
    # NOE receives only this ratio overview, never an exponential fit or decay residuals.
    (stage / 'plot.gnu').write_text('\n'.join(plot) + '\n')
    names += ['intensity.out', 'plot.gnu']
    return names


def checked_run(command, cwd, stdout=None):
    print('+ ' + ' '.join(map(str, command)), flush=True)
    subprocess.run(list(map(str, command)), cwd=cwd, stdout=stdout, check=True)


def analyze(args):
    source = Path.cwd().resolve()
    matches = sorted(Path(p).resolve() for p in glob.glob(args.peaks) if Path(p).is_file())
    if len(matches) != 1:
        raise ValueError('Select exactly one peak list; set SPARKY_PEAKS to its filename')
    peak_path = matches[0]
    data_path, header_path, noise_path = (Path(p).resolve() for p in (args.data, args.header, args.noise))
    times = counted_vector(header_path, nonnegative=args.mode != 'noe')
    noises = counted_vector(noise_path, positive=True)
    if len(times) != len(noises):
        raise ValueError('Delay and noise counts differ')
    if args.format == 'nmrpipe' and plane_count(data_path) != len(times):
        raise ValueError('Spectrum plane, delay and noise counts differ')
    if args.mode != 'noe' and len(set(times)) < 2:
        raise ValueError('Exponential fitting needs at least two distinct time points')
    peaks = read_peaks(peak_path, args.peak_format)
    prefix = args.prefix or peak_path.name.split('.')[0]
    if Path(prefix).name != prefix or prefix in ('', '.', '..'):
        raise ValueError('Output prefix must be a filename without directories')
    pdf_name = peak_path.name + '_intensity.pdf'
    run_tag = source.parent.name + '/' + source.name
    with run_directory(source, 'analyze') as stage:
        clean = stage / 'peaks.list'
        clean.write_text('assignment ' + ' '.join(f'w{i}' for i in range(1, len(peaks[0]['coordinates']) + 1)) + '\n' +
                         ''.join(p['peak_id'] + ' ' + ' '.join(p['coordinates']) + '\n' for p in peaks))
        raw = stage / 'extracted_ids.out'
        with raw.open('w') as handle:
            checked_run(['data_xtr', '-data', data_path, '-dfmt', args.format,
                         '-peak', clean, '-pfmt', args.peak_format], stage, handle)
        extracted = parse_extracted(raw.read_text(), peaks, len(times))
        # Use original order deliberately; identification is independent of fitter order.
        fit_input = header_path.read_text().rstrip() + '\n' + ''.join(
            p['peak_id'] + ' ' + ' '.join(f'{v:.17g}' for v in extracted[p['peak_id']]) + '\n' for p in peaks)
        (stage / 'rlx_dfit.inp').write_text(fit_input)
        shutil.copyfile(noise_path, stage / 'vEstNoise.txt')
        command = ['rlx_dfit', '-mode', 'noe' if args.mode == 'noe' else 'exp2', '-rms_noise',
                   '-data', 'rlx_dfit.inp', '-d_esd', 'vEstNoise.txt']
        if args.mode != 'noe':
            command += ['-r_init', '1.0' if args.mode == 'r1' else '10.0']
        checked_run(command + ['-out', 'rlx_dfit.out'], stage)
        nonempty(stage / 'rlx_dfit.out')
        names = prepare_reports(stage, peaks, extracted, (stage / 'rlx_dfit.out').read_text(),
                                times, noises, args.mode, prefix, pdf_name, run_tag)
        restored = ''.join(p['assignment'] + ' ' + ' '.join(f'{v:.17g}' for v in extracted[p['peak_id']]) + '\n' for p in peaks)
        (stage / 'data_xtr.out').write_text(restored)
        checked_run(['gnuplot', 'plot.gnu'], stage)
        nonempty(stage / pdf_name)
        names += ['data_xtr.out', 'rlx_dfit.inp', 'rlx_dfit.out', pdf_name]
        publish(stage, source, names)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    p = commands.add_parser('delays')
    p.add_argument('--base', required=True)
    p.add_argument('--vc', default='vclist'); p.add_argument('--vd', default='vdlist')
    p.add_argument('--ncyc', default='ncyc_list'); p.add_argument('--source', choices=('auto', 'vc', 'vd', 'ncyc'), default='auto')
    p.add_argument('--block', default='.01'); p.add_argument('--constant', default='.01')
    p.add_argument('--cycle-seconds', default=''); p.add_argument('--cycle-constant', default='0')
    p = commands.add_parser('check-fid'); p.add_argument('path'); p.add_argument('delays')
    p = commands.add_parser('check-conversion'); p.add_argument('delays')
    p = commands.add_parser('check-fid-series'); p.add_argument('directory'); p.add_argument('delays')
    p = commands.add_parser('noise-value'); p.add_argument('path')
    commands.add_parser('validate-proc')
    p = commands.add_parser('series'); p.add_argument('paths', nargs='+')
    p = commands.add_parser('run-proc'); p.add_argument('--script', required=True); p.add_argument('--source', required=True)
    p = commands.add_parser('analyze')
    p.add_argument('--peaks', required=True); p.add_argument('--data', required=True)
    p.add_argument('--format', default='nmrpipe'); p.add_argument('--mode', choices=('r1', 'r2', 'noe'), required=True)
    p.add_argument('--peak-format', choices=('sparky2d', 'sparky3d'), default='sparky3d')
    p.add_argument('--header', required=True); p.add_argument('--noise', required=True); p.add_argument('--prefix', default='')
    args = parser.parse_args(argv)
    try:
        if args.command == 'delays':
            create_delays(args.base, args.vc, args.vd, args.ncyc, args.source, args.block,
                          args.constant, args.cycle_seconds, args.cycle_constant)
        elif args.command == 'check-fid':
            check_fid(args.path, args.delays)
        elif args.command == 'check-conversion':
            check_conversion(args.delays)
        elif args.command == 'check-fid-series':
            check_fid_series(args.directory, args.delays)
        elif args.command == 'noise-value':
            print(f'{noise_value(args.path):.17g}')
        elif args.command == 'validate-proc':
            check_processing()
        elif args.command == 'series':
            set_series_metadata(args.paths)
        elif args.command == 'run-proc':
            run_processing(args.script, args.source)
        elif args.command == 'analyze':
            analyze(args)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        print(f'Relaxation run failed: {error}', file=sys.stderr)
        return 1
    return 0


def handle_termination(signum, frame):
    raise InterruptedError('Run interrupted by SIGTERM')


if __name__ == '__main__':
    signal.signal(signal.SIGTERM, handle_termination)
    sys.exit(main())
