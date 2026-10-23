#!/usr/bin/env python3
"""Validation and isolated-run support for the editable relaxation .com scripts.

Python standard library only. No spectrum processing occurs on import.
"""
import argparse
from collections import Counter
from contextlib import contextmanager
import csv
from decimal import Decimal, InvalidOperation
from datetime import datetime, timezone
import fcntl
import glob
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import shlex
import statistics
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


def decay_panel_layout(residual=False, top=.90, bottom=.12):
    """Align panels; usable plot-area heights have an exact 80:20 ratio."""
    lower = (top - bottom - .03) * .2
    top, bottom = (bottom + lower, bottom) if residual else (top, bottom + lower + .03)
    return ['set lmargin at screen 0.14', 'set rmargin at screen 0.96',
            f'set tmargin at screen {top:.6g}', f'set bmargin at screen {bottom:.6g}',
            'set format x "%g"' if residual else 'set format x ""',
            'set xlabel "Delay (s)"' if residual else 'unset xlabel']


AUTO_PLOT_MARGINS = ('set lmargin -1', 'set rmargin -1', 'set tmargin -1', 'set bmargin -1')


def gnuplot_string(text):
    return '"' + str(text).replace('\\', '\\\\').replace('"', '\\"').replace('\n', '\\n').replace('\r', ' ') + '"'


# Versioned, portable reporting. All calculations below use saved measurements
# and reported fit parameters; none of these functions estimates fit parameters.
REPORT_SCHEMA = '1'
ASSESSMENT = dict(alpha=.05, coherence_min_peaks=10, coherence_sign_fraction=.8,
                  coherence_median_sd=2.)
CONCLUSIONS = {
    'departures_detected': 'Departures detected at supplied noise levels',
    'no_evidence_against': 'No evidence against single-exponential decay at supplied noise levels',
    'inconclusive': 'Inconclusive',
    'not_applicable': 'Not applicable: NOE ratio analysis',
}
TRACE_FIELDS = ['schema_version', 'analysis_id', 'source_collection', 'source_experiment',
                'peak_list', 'output_prefix', 'pdf_name', 'experiment_type', 'fit_model',
                'fit_origin', 'fit_revision', 'report_revision', 'fit_input_sha256',
                'fit_output_sha256', 'noise_sha256']
CSV_FIELDS = TRACE_FIELDS + [
    'peak_index', 'peak_id', 'assignment', 'residue', 'w1', 'w2', 'w3', 'overlap',
    'duplicate_assignment', 'valid_fit', 'plane_index', 'delay_s', 'plane_code',
    'intensity', 'noise_sd', 'prediction', 'residual', 'residual_noise_sd',
    'value', 'value_error', 'value_units', 'i0', 'i0_error', 'fitter_chi2', 'fitter_threshold',
    'uncertainty_rescaled', 'warnings', 'flags', 'chi2_input', 'reduced_chi2_input', 'dof',
    'p_input', 'p_holm', 'primary_eligible', 'exclusion_reason', 'repeat_chi2', 'repeat_dof',
    'repeat_pair_z', 'repeat_delay_chi2', 'repeat_delay_dof', 'repeat_delay_pair_z', 'assessment_method', 'dataset_conclusion', 'primary_peak_count', 'holm_rejection_count',
    'coherent_planes', 'alpha', 'coherence_min_peaks', 'coherence_sign_fraction',
    'coherence_median_sd', 'plane_median_sd', 'plane_q1_sd', 'plane_q3_sd',
    'plane_sign_fraction', 'plane_bias', 'delay_median_sd', 'delay_q1_sd', 'delay_q3_sd']


def file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def code_revision():
    # A copied helper must not accidentally report the surrounding data repo's
    # Git revision. The content fingerprint identifies the actual executed code.
    return 'sha256:' + file_sha256(__file__)


def chi_square_sf(value, dof):
    """Upper regularized gamma Q(dof/2, value/2), Numerical Recipes expansion.

    Goodness-of-fit probabilities are approximate for the nonlinear model and
    conditional on the supplied independent Gaussian measurement uncertainties.
    """
    if not math.isfinite(value) or value < 0 or dof <= 0:
        raise ValueError('Invalid chi-square statistic or degrees of freedom')
    if value == 0:
        return 1.
    a, x = dof / 2., value / 2.
    factor = math.exp(a * math.log(x) - x - math.lgamma(a))
    if x < a + 1:
        term = total = 1. / a
        for n in range(1, 10001):
            term *= x / (a + n)
            total += term
            if abs(term) <= abs(total) * 1e-14:
                return min(1., max(0., 1. - factor * total))
    else:
        tiny = 1e-300
        b = x + 1. - a
        c, d = 1. / tiny, 1. / max(b, tiny)
        h = d
        for n in range(1, 10001):
            an = -n * (n - a)
            b += 2.
            d = an * d + b
            c = b + an / c
            if abs(d) < tiny:
                d = tiny
            if abs(c) < tiny:
                c = tiny
            d = 1. / d
            delta = d * c
            h *= delta
            if abs(delta - 1.) <= 1e-14:
                return min(1., max(0., factor * h))
    raise ValueError('Chi-square probability did not converge')


def holm_adjust(pvalues):
    ordered = sorted(pvalues, key=pvalues.get)
    result, previous = {}, 0.
    for rank, key in enumerate(ordered):
        previous = max(previous, min(1., (len(ordered) - rank) * pvalues[key]))
        result[key] = previous
    return result


def quantile(values, fraction):
    values = sorted(values)
    if not values:
        return None
    pos = (len(values) - 1) * fraction
    lo, hi = math.floor(pos), math.ceil(pos)
    return values[lo] + (values[hi] - values[lo]) * (pos - lo)


def residual_summary(values):
    if not values:
        return dict(median=None, q1=None, q3=None, sign_fraction=None, bias=False, n=0)
    median = statistics.median(values)
    signs = max(sum(v > 0 for v in values), sum(v < 0 for v in values)) / len(values)
    bias = (len(values) >= ASSESSMENT['coherence_min_peaks'] and
            signs >= ASSESSMENT['coherence_sign_fraction'] and
            abs(median) > ASSESSMENT['coherence_median_sd'])
    return dict(median=median, q1=quantile(values, .25), q3=quantile(values, .75),
                sign_fraction=signs, bias=bias, n=len(values))


def report_context(mode, prefix, pdf_name, run_tag, **extra):
    parts = run_tag.rsplit('/', 1)
    context = dict(schema_version=REPORT_SCHEMA, source_collection=parts[0],
                   source_experiment=parts[1] if len(parts) == 2 else '',
                   peak_list=prefix + '.list', output_prefix=prefix, pdf_name=pdf_name,
                   experiment_type=mode, fit_model='noe' if mode == 'noe' else 'exp2',
                   fit_origin='saved-results', fit_revision='unknown',
                   fit_input_sha256='', fit_output_sha256='', noise_sha256='')
    context.update(extra)
    context['report_revision'] = code_revision()
    return context


def build_assessment(rows, extracted, times, noises, context):
    n = len(times)
    if not n or len(noises) != n or any(not math.isfinite(v) or v <= 0 for v in noises):
        raise ValueError('Invalid measurement/noise counts')
    mode = context['experiment_type']
    if mode not in ('r1', 'r1rho', 'r2', 'noe') or any(not math.isfinite(t) for t in times):
        raise ValueError('Invalid report experiment type or plane values')
    if mode != 'noe' and any(t < 0 for t in times):
        raise ValueError('Negative decay delay')
    if set(extracted) != {r['peak_id'] for r in rows}:
        raise ValueError('Report measurements and peak identities differ')
    groups = {}
    for i, time in enumerate(times):
        groups.setdefault(time, []).append(i)
    primary = []
    for row in rows:
        y = extracted[row['peak_id']]
        if len(y) != n or any(not math.isfinite(v) for v in y):
            raise ValueError('Invalid saved intensity count or value')
        row.update(points=[], chi2_input=None, reduced_chi2_input=None, p_input=None,
                   p_holm=None, dof=n-2 if mode != 'noe' else None, repeat_chi2=None,
                   repeat_dof=0 if mode != 'noe' else None, repeat_pair_z=None, repeat_groups={})
        for t, intensity, sd in zip(times, y, noises):
            prediction = residual = z = None
            if mode != 'noe' and row['valid']:
                try:
                    prediction = row['values'][2] * math.exp(-row['values'][0] * t)
                except OverflowError as error:
                    raise ValueError('Fitted prediction overflows') from error
                if not math.isfinite(prediction):
                    raise ValueError('Nonfinite fitted prediction')
                residual = intensity - prediction
                z = residual / sd
            row['points'].append(dict(intensity=intensity, noise=sd, prediction=prediction,
                                      residual=residual, z=z))
        if mode != 'noe':
            repeat_ss, repeat_df, pairs = 0., 0, []
            for delay, indices in groups.items():
                if len(indices) < 2:
                    continue
                weights = [1. / noises[i]**2 for i in indices]
                mean = math.fsum(y[i]*w for i, w in zip(indices, weights)) / math.fsum(weights)
                scatter = math.fsum((y[i]-mean)**2/noises[i]**2 for i in indices)
                repeat_ss += scatter
                repeat_df += len(indices)-1
                pair_z = None
                if len(indices) == 2:
                    i, j = indices
                    pair_z = (y[i]-y[j]) / math.hypot(noises[i], noises[j])
                    pairs.append(pair_z)
                row['repeat_groups'][delay] = dict(chi2=scatter, dof=len(indices)-1, pair_z=pair_z)
            row.update(repeat_chi2=repeat_ss if repeat_df else None, repeat_dof=repeat_df,
                       repeat_pair_z=pairs[0] if len(pairs) == 1 else None)
            if row['valid']:
                # Explicit accurate summation is reproducible across Python 3.10
                # and 3.12, whose built-in floating-point sum algorithms differ.
                chi2 = math.fsum(point['z']**2 for point in row['points'])
                if not math.isfinite(chi2):
                    raise ValueError('Nonfinite original-noise chi-square')
                row['chi2_input'] = chi2
                if n > 2:
                    row['reduced_chi2_input'] = chi2 / (n-2)
                    row['p_input'] = chi_square_sf(chi2, n-2)
        reasons = []
        if mode == 'noe': reasons.append('not_decay_analysis')
        if not row['valid']: reasons.append('invalid_fit')
        if row['overlap']: reasons.append('overlap')
        if mode != 'noe' and n <= 2: reasons.append('no_positive_residual_dof')
        row['primary_eligible'] = not reasons
        row['exclusion_reason'] = ';'.join(reasons)
        if not reasons: primary.append(row)
    adjusted = holm_adjust({r['peak_id']: r['p_input'] for r in primary})
    for row in primary:
        row['p_holm'] = adjusted[row['peak_id']]
        if row['p_holm'] < ASSESSMENT['alpha']:
            row['flags'].append('original_noise_lack_of_fit')
    planes = [residual_summary([r['points'][i]['z'] for r in primary]) for i in range(n)]
    delays = {t: residual_summary([r['points'][i]['z'] for r in primary for i in indices])
              for t, indices in groups.items()}
    coherent = [i+1 for i, summary in enumerate(planes) if summary['bias']]
    rejected = sum(r['p_holm'] < ASSESSMENT['alpha'] for r in primary)
    conclusion = ('not_applicable' if mode == 'noe' else 'inconclusive' if not primary else
                  'departures_detected' if rejected or coherent else 'no_evidence_against')
    # Result identity excludes reporting-code revision and experiment-label changes,
    # but binds all saved numerical inputs, fitted parameters, mappings and warnings.
    identity = dict(source=[context['source_collection'], context['source_experiment']],
                    model=context['fit_model'], times=list(map(float,times)), noises=list(map(float,noises)),
                    peaks=[dict(peak_id=r['peak_id'], assignment=r['assignment'],
                                coordinates=r['coordinates'], residue=r['residue'],
                                overlap=r['overlap'], values=list(map(float,r['values'])), warnings=r['warnings'],
                                intensities=list(map(float,extracted[r['peak_id']]))) for r in rows])
    context = dict(context, analysis_id='a1-' + hashlib.sha256(
        json.dumps(identity, sort_keys=True, separators=(',', ':')).encode()).hexdigest()[:20])
    return dict(context=context, rows=rows, times=list(map(float,times)), noises=list(map(float,noises)), planes=planes,
                delays=delays, conclusion=conclusion, primary_count=len(primary),
                rejected=rejected, coherent=coherent)


def analysis_csv_rows(report):
    context = report['context']
    mode = context['experiment_type']
    global_fields = dict(dataset_conclusion=report['conclusion'], primary_peak_count=report['primary_count'],
                         holm_rejection_count=report['rejected'],
                         coherent_planes=';'.join(map(str, report['coherent'])), **ASSESSMENT)
    for peak_index, row in enumerate(report['rows'], 1):
        v = row['values']
        fields = dict(context, **global_fields, peak_index=peak_index, peak_id=row['peak_id'],
                      assignment=row['assignment'], residue=row['residue'], overlap=row['overlap'],
                      duplicate_assignment=row['duplicate_assignment'], valid_fit=row['valid'],
                      value=v[0], value_error=v[1], value_units='ratio' if mode == 'noe' else 's^-1',
                      i0=v[2] if mode != 'noe' else None, i0_error=v[3] if mode != 'noe' else None,
                      fitter_chi2=v[4] if mode != 'noe' else None,
                      fitter_threshold=v[5] if mode != 'noe' else None,
                      uncertainty_rescaled='uncertainty_rescaled' in row['flags'],
                      warnings=row['warnings'], flags=';'.join(row['flags']),
                      **{key: row[key] for key in ('chi2_input','reduced_chi2_input','dof','p_input',
                         'p_holm','primary_eligible','exclusion_reason','repeat_chi2','repeat_dof','repeat_pair_z')})
        fields.update({f'w{i}': c for i, c in enumerate(row['coordinates'], 1)})
        for i, (time, point) in enumerate(zip(report['times'], row['points'])):
            plane, delay = report['planes'][i], report['delays'][time]
            repeat = row['repeat_groups'].get(time,{})
            yield dict(fields, plane_index=i+1, delay_s=time if mode != 'noe' else None,
                       plane_code=time if mode == 'noe' else None, intensity=point['intensity'],
                       noise_sd=point['noise'], prediction=point['prediction'], residual=point['residual'],
                       residual_noise_sd=point['z'], plane_median_sd=plane['median'],
                       plane_q1_sd=plane['q1'], plane_q3_sd=plane['q3'],
                       plane_sign_fraction=plane['sign_fraction'], plane_bias=plane['bias'],
                       repeat_delay_chi2=repeat.get('chi2'), repeat_delay_dof=repeat.get('dof'),
                       repeat_delay_pair_z=repeat.get('pair_z'),
                       assessment_method=('Exponential suitability not applicable; NOE ratio analysis; supplied intensities, noise and ratio uncertainties unchanged' if mode=='noe' else
                           'I(t)=I0 exp(-R t); independent Gaussian input noise; chi2 df=N-2; approximate nonlinear GOF; Holm alpha=0.05; descriptive shared-plane screen >=10 peaks, >=80% same sign, |median residual/SD|>2; few repeats limit inference; rejection or uncertainty rescaling does not establish biexponentiality'),
                       delay_median_sd=delay['median'], delay_q1_sd=delay['q1'], delay_q3_sd=delay['q3'])


def write_analysis_csv(path, report):
    with Path(path).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in analysis_csv_rows(report):
            writer.writerow({field: row.get(field) for field in CSV_FIELDS})


def strict_bool(value):
    if value not in ('True', 'False'):
        raise ValueError(f'Invalid saved boolean: {value!r}')
    return value == 'True'


def safe_output_name(name):
    if not name or Path(name).name != name or name in ('.', '..'):
        raise ValueError('Report output names must be filenames without directories')
    return name


def load_analysis_csv(path):
    with Path(path).open(newline='') as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != CSV_FIELDS:
            raise ValueError('Unsupported or incomplete analysis CSV schema')
        records = list(reader)
    if not records or any(None in row or None in row.values() for row in records):
        raise ValueError('Empty or incomplete analysis CSV')
    context = {key: records[0][key] for key in TRACE_FIELDS}
    if context['schema_version'] != REPORT_SCHEMA:
        raise ValueError('Unsupported analysis CSV schema version')
    safe_output_name(context['output_prefix']); safe_output_name(context['pdf_name'])
    if context['fit_model'] != ('noe' if context['experiment_type'] == 'noe' else 'exp2'):
        raise ValueError('CSV model/experiment type mismatch')
    grouped = {}
    for record in records:
        if any(record[field] != context[field] for field in TRACE_FIELDS):
            raise ValueError('Inconsistent repeated CSV provenance')
        index = int(record['peak_index']); plane = int(record['plane_index'])
        if index < 1 or plane < 1:
            raise ValueError('Invalid peak or plane index')
        slots = grouped.setdefault(index, {})
        if plane in slots:
            raise ValueError('Duplicate peak/plane CSV row')
        slots[plane] = record
    if sorted(grouped) != list(range(1, len(grouped)+1)):
        raise ValueError('Missing CSV peak index')
    peaks, extracted, fit_lines, times, noises, seen_ids, seen_entries = [], {}, [], None, None, set(), set()
    warning_map = {}
    varying = {'plane_index','delay_s','plane_code','intensity','noise_sd','prediction','residual',
               'residual_noise_sd','plane_median_sd','plane_q1_sd','plane_q3_sd',
               'plane_sign_fraction','plane_bias','delay_median_sd','delay_q1_sd','delay_q3_sd',
               'repeat_delay_chi2','repeat_delay_dof','repeat_delay_pair_z'}
    for index, slots in sorted(grouped.items()):
        if sorted(slots) != list(range(1,len(slots)+1)):
            raise ValueError('Missing CSV measurement plane')
        first = slots[1]
        for record in slots.values():
            if any(record[key] != first[key] for key in CSV_FIELDS if key not in varying):
                raise ValueError('Inconsistent repeated CSV fit or peak information')
        tag = first['peak_id']
        coords = [first[f'w{i}'] for i in (1,2,3) if first[f'w{i}'] != '']
        if tag in seen_ids or len(coords) not in (2,3):
            raise ValueError('Duplicate peak ID or invalid saved coordinates')
        numeric_coords = tuple(finite_number(c,'coordinate',0) for c in coords)
        entry = (first['assignment'],numeric_coords)
        if entry in seen_entries:
            raise ValueError('Duplicate saved assignment and coordinates')
        seen_ids.add(tag);seen_entries.add(entry)
        peaks.append(dict(peak_id=tag, assignment=first['assignment'], residue=int(first['residue']),
                          coordinates=coords, overlap=strict_bool(first['overlap']),
                          duplicate_assignment=strict_bool(first['duplicate_assignment'])))
        field = 'plane_code' if context['experiment_type'] == 'noe' else 'delay_s'
        local_times = [finite_number(slots[i][field],'plane value') for i in sorted(slots)]
        local_noises = [finite_number(slots[i]['noise_sd'],'noise SD',0) for i in sorted(slots)]
        if times is not None and (times != local_times or noises != local_noises):
            raise ValueError('CSV peaks disagree on delay/noise vectors')
        times, noises = local_times, local_noises
        extracted[tag] = [finite_number(slots[i]['intensity'],'intensity') for i in sorted(slots)]
        fields = ['value','value_error'] + ([] if context['experiment_type']=='noe' else
                                          ['i0','i0_error','fitter_chi2','fitter_threshold'])
        fitted = [float(first[key]) for key in fields]
        warning_map[tag] = first['warnings']
        fit_lines.append(tag+' '+' '.join(map(str,fitted)))
    counts = Counter(p['assignment'] for p in peaks)
    if any(p['duplicate_assignment'] != (counts[p['assignment']]>1) for p in peaks):
        raise ValueError('CSV duplicate-assignment flags disagree with identities')
    bundle = dict(peaks=peaks, extracted=extracted, fit_text='\n'.join(fit_lines)+'\n',
                  times=list(map(float,times)), noises=list(map(float,noises)), context=context, warning_map=warning_map)
    # Recompute all diagnostics and verify the result identity and every exported
    # value. This prevents silently accepting edited or inconsistent derived fields.
    rows = fitted_rows(peaks,bundle['fit_text'],context['experiment_type'],warning_map)
    computed = build_assessment(rows,extracted,times,noises,context)
    if computed['context']['analysis_id'] != context['analysis_id']:
        raise ValueError('CSV analysis identity does not match its contents')
    expected = {(str(r['peak_index']),str(r['plane_index'])): r for r in analysis_csv_rows(computed)}
    for record in records:
        wanted = expected[(record['peak_index'],record['plane_index'])]
        for key in CSV_FIELDS:
            if key == 'report_revision': continue
            actual = record[key]; value = wanted.get(key)
            canonical = '' if value is None else str(value)
            if actual == canonical: continue
            try:
                a, b = float(actual), float(canonical)
            except ValueError:
                raise ValueError(f'CSV field {key} disagrees with saved measurements/fit')
            if not ((math.isnan(a) and math.isnan(b)) or math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12)):
                raise ValueError(f'CSV field {key} disagrees with saved measurements/fit')
    return bundle


def close_saved(a,b):
    return (math.isnan(a) and math.isnan(b)) or a == b or math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12)


def load_saved_analysis(source, prefix, mode, peak_list_name=None):
    source = Path(source)
    safe_output_name(prefix)
    if mode not in ('r1','r1rho','r2','noe'):
        raise ValueError('Saved-result refresh requires --experiment-type')
    diagnostic = source / (prefix+'_fit_diagnostics.tsv')
    with diagnostic.open(newline='') as handle:
        mappings = list(csv.DictReader(handle,delimiter='\t'))
    if not mappings: raise ValueError('Empty saved diagnostics mapping')
    peaks, seen, counts = [], set(), Counter(r['assignment'] for r in mappings)
    entries = set()
    for row in mappings:
        if row['peak_id'] in seen: raise ValueError('Duplicate saved peak ID')
        seen.add(row['peak_id'])
        coords = [row[f'w{i}'] for i in (1,2,3) if row[f'w{i}']]
        if len(coords) not in (2,3) or any(not math.isfinite(float(c)) or float(c)<0 for c in coords):
            raise ValueError('Invalid saved peak coordinates')
        entry=(row['assignment'],tuple(map(float,coords)))
        if entry in entries or strict_bool(row['duplicate_assignment']) != (counts[row['assignment']]>1):
            raise ValueError('Duplicate saved assignment/coordinates or inconsistent duplicate flag')
        entries.add(entry)
        peaks.append(dict(peak_id=row['peak_id'],assignment=row['assignment'],residue=int(row['residue']),
                          coordinates=coords,overlap=strict_bool(row['overlap']),
                          duplicate_assignment=counts[row['assignment']]>1))
    fit_text=(source/'rlx_dfit.out').read_text()
    fitted, _=parse_fit(fit_text,mode)
    if set(fitted) != seen: raise ValueError('Saved fit identifiers differ from diagnostics')
    for row in mappings:
        fields=['value','error']+([] if mode=='noe' else ['i0','i0_error','chi2','chi2_threshold'])
        if any(not ((math.isnan(float(row[key])) and math.isnan(value)) or float(row[key])==value)
               for key,value in zip(fields,fitted[row['peak_id']])):
            raise ValueError('Saved fit parameters disagree with diagnostics')
    lines=meaningful_lines(source/'rlx_dfit.inp')
    header=lines[0].split();n=int(header[0]);times=[finite_number(t,'saved delay') for t in header[1:]]
    if n!=len(times) or n<1: raise ValueError('Saved fitter delay header count mismatch')
    extracted=parse_extracted('\n'.join(lines[1:]),peaks,n)
    noises=counted_vector(source/'vEstNoise.txt',positive=True)
    if len(noises)!=n or counted_vector(source/'relax_header.txt',nonnegative=mode!='noe')!=times:
        raise ValueError('Saved input, noise or delay vectors disagree')
    # Check against the saved plotting table to detect mixed fit/noise generations.
    blocks={};tag=None
    for line in (source/'intensity.out').read_text().splitlines():
        match=re.search(r'^# (.*?) peak_id=(\S+) index (\d+)\s*$',line)
        if match:
            tag=match[2]
            if tag in blocks:raise ValueError('Duplicate saved intensity block')
            blocks[tag]=dict(assignment=match[1],index=int(match[3]),values=[])
        elif line.strip() and not line.lstrip().startswith('#'):
            if tag is None:raise ValueError('Saved intensity data lacks identity mapping')
            blocks[tag]['values'].append(list(map(float,line.split())))
    if set(blocks)!=seen:raise ValueError('Saved intensity blocks differ from fit identifiers')
    rows=fitted_rows(peaks,fit_text,mode)
    for saved, row in zip(mappings, rows):
        if strict_bool(saved['valid_fit']) != row['valid'] or saved['warnings'] != row['warnings']:
            raise ValueError('Saved diagnostics validity/warnings disagree with native fit')
    for index,row in enumerate(rows):
        block=blocks[row['peak_id']]
        if block['assignment']!=row['assignment'] or block['index']!=index or len(block['values'])!=n:
            raise ValueError('Saved intensity identity/order mismatch')
        v=row['values'];normalized=mode!='noe' and row['valid']
        for i,(t,y,sd,actual) in enumerate(zip(times,extracted[row['peak_id']],noises,block['values'])):
            if mode=='noe': expected=[i+1,y,sd]
            elif normalized:
                prediction=v[2]*math.exp(-v[0]*t)
                expected=[t,y/v[2],sd/abs(v[2]),prediction/v[2],(y-prediction)/v[2],(y-prediction)/sd]
            else:expected=[t,y,sd,math.nan,math.nan,math.nan]
            if len(actual)!=len(expected) or any(not close_saved(a,b) for a,b in zip(actual,expected)):
                raise ValueError('Saved intensity/noise table disagrees with fitter inputs/parameters')
    if peak_list_name is None:
        outputs=[line for line in (source/'plot.gnu').read_text().splitlines() if line.startswith('set output ')]
        if len(outputs)!=1:raise ValueError('Provide --peak-list-name for saved report refresh')
        pdf_name=shlex.split(outputs[0])[2]
        suffix='_intensity.pdf'
        if not pdf_name.endswith(suffix):raise ValueError('Provide --peak-list-name for saved report refresh')
        peak_list_name=pdf_name[:-len(suffix)]
    safe_output_name(peak_list_name)
    context=report_context(mode,prefix,peak_list_name+'_intensity.pdf',source.parent.name+'/'+source.name,
                           peak_list=peak_list_name,fit_input_sha256=file_sha256(source/'rlx_dfit.inp'),
                           fit_output_sha256=file_sha256(source/'rlx_dfit.out'),
                           noise_sha256=file_sha256(source/'vEstNoise.txt'))
    # Reports produced by this schema retain the original fit provenance on later
    # legacy-prefix refreshes. Unknown old fit-producing code stays unknown.
    if all(key in mappings[0] for key in TRACE_FIELDS):
        previous={key:mappings[0][key] for key in TRACE_FIELDS}
        if any(any(row[key]!=previous[key] for key in TRACE_FIELDS) for row in mappings):
            raise ValueError('Inconsistent saved diagnostics provenance')
        for key in ('fit_input_sha256','fit_output_sha256','noise_sha256'):
            if previous[key] and previous[key]!=context[key]:
                raise ValueError('Saved diagnostics fingerprints disagree with native inputs/output')
        if previous['fit_model']!=context['fit_model'] or previous['output_prefix']!=prefix:
            raise ValueError('Saved diagnostics model/prefix mismatch')
        context.update({key:previous[key] for key in ('source_collection','source_experiment',
                                                     'fit_origin','fit_revision')})
    return dict(peaks=peaks,extracted=extracted,fit_text=fit_text,times=times,noises=noises,context=context)


def refresh_report(args):
    if args.csv:
        input_path=Path(args.csv).resolve()
        bundle=load_analysis_csv(input_path)
        source=input_path.parent
        watched=[input_path]
        if args.experiment_type and args.experiment_type!=bundle['context']['experiment_type']:
            raise ValueError('CSV experiment type is fixed; use saved-result refresh to relabel legacy fits')
    else:
        source=Path.cwd().resolve()
        bundle=load_saved_analysis(source,args.saved_prefix,args.experiment_type,args.peak_list_name)
        watched=[source/name for name in ('rlx_dfit.inp','rlx_dfit.out','vEstNoise.txt','relax_header.txt',
                  'intensity.out',args.saved_prefix+'_fit_diagnostics.tsv')]
    destination=Path(args.output_dir).resolve() if args.output_dir else source
    if not destination.is_dir():raise ValueError('Report output directory must already exist')
    # Re-read under the same lock used by fitting/publication to avoid mixed runs.
    with run_directory(destination,'report') as stage:
        if args.csv: bundle=load_analysis_csv(input_path)
        elif destination==source:bundle=load_saved_analysis(source,args.saved_prefix,args.experiment_type,args.peak_list_name)
        else:raise ValueError('Legacy refresh writes beside saved results; use CSV to export elsewhere')
        before={p:file_sha256(p) for p in watched}
        c=bundle['context'];mode=c['experiment_type']
        names=prepare_reports(stage,bundle['peaks'],bundle['extracted'],bundle['fit_text'],
                              bundle['times'],bundle['noises'],mode,c['output_prefix'],c['pdf_name'],
                              c['source_collection']+'/'+c['source_experiment'],context=c, warning_map=bundle.get('warning_map'))
        checked_run(['gnuplot','plot.gnu'],stage)
        nonempty(stage/c['pdf_name']);names.append(c['pdf_name'])
        if any(file_sha256(p)!=digest for p,digest in before.items()):
            raise ValueError('Saved report inputs changed during rendering')
        with (stage/(c['output_prefix']+'_analysis.csv')).open(newline='') as handle:
            report_id=next(csv.DictReader(handle))['analysis_id']
        existing=[name for name in names if (destination/name).is_file()]
        if existing:
            stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
            backup=destination/'.report-backups'/report_id/stamp
            backup.mkdir(parents=True)
            for name in existing:shutil.copy2(destination/name,backup/name)
        publish(stage,destination,names)
        print(f'Report {report_id}: {destination/(c["output_prefix"]+"_analysis.csv")}',flush=True)



def display_number(value, precision=4):
    return 'NA' if value is None or not math.isfinite(value) else f'{value:.{precision}g}'


def fit_annotation(row, mode):
    v = row['values']
    unit = 'ratio' if mode == 'noe' else 's^-1'
    result = f"{'hetNOE' if mode == 'noe' else 'R1rho' if mode == 'r1rho' else mode.upper()} = {display_number(v[0])} +/- {display_number(v[1])} {unit}"
    if mode != 'noe':
        result += f"; I0 = {display_number(v[2])} +/- {display_number(v[3])}"
        result += f"\nOriginal-noise chi2/df = {display_number(row['chi2_input'])}/{row['dof']}; Holm p = {display_number(row['p_holm'],3)}"
    return result


def text_page(title, lines):
    """One PDF page; label 9000 is a persistent trace footer on every page."""
    commands = list(AUTO_PLOT_MARGINS) + ['unset border', 'unset xtics', 'unset ytics', 'unset grid',
                'unset xlabel', 'unset ylabel', 'set key off', 'set title ' + gnuplot_string(title),
                'set xrange [0:1]', 'set yrange [0:1]']
    for i, line in enumerate(lines, 1):
        commands.append(f'set label {i} ' + gnuplot_string(line) +
                        f' at screen .045,{.88-(i-1)*.034:.5f} left font "Sans,7" front')
    commands += ['plot 0 lc rgb "white" notitle']
    commands += [f'unset label {i}' for i in range(1,len(lines)+1)]
    return commands + ['set border', 'set xtics', 'set ytics', 'set grid', 'set autoscale']


def report_intro(report):
    c = report['context']; mode = c['experiment_type']; rows = report['rows']
    name = 'T1rho / R1rho' if mode == 'r1rho' else 'hetNOE' if mode == 'noe' else mode.upper()
    footer = c['analysis_id'] + ' | ' + c['source_collection'] + '/' + c['source_experiment']
    commands = ['set label 9000 ' + gnuplot_string(footer) + ' at screen .02,.022 font "Sans,5.5" front',
                'set label 9001 ' + gnuplot_string('Report code: ' + c['report_revision']) + ' at screen .02,.008 font "Sans,4.5" front']
    lines = [CONCLUSIONS[report['conclusion']],
             f"{len(rows)} peaks; {len(report['times'])} acquisition planes; {len(set(report['times']))} distinct plane values",
             f"Primary assessment: {report['primary_count']} eligible; {len(rows)-report['primary_count']} excluded",
             f"Holm rejections: {report['rejected']}; descriptive shared-plane flags: {report['coherent'] or 'none'}",
             'Source: ' + c['source_collection'] + '/' + c['source_experiment'],
             'Peak list: ' + c['peak_list'], 'Model: ' + c['fit_model'] + '; saved-fit origin: ' + c['fit_origin'],
             'Analysis ID: ' + c['analysis_id'], 'Fit-producing code: ' + c['fit_revision']]
    lines += textwrap.wrap('Reporting code: ' + c['report_revision'],100)
    for key in ('fit_input_sha256','fit_output_sha256','noise_sha256'):
        lines += textwrap.wrap(key + ': ' + (c[key] or 'not recorded'),100)
    lines += ['No extraction, processing or refitting is performed by report regeneration.',
              'CSV carries full precision, measurements, uncertainties, identities, warnings and methods settings.']
    commands += text_page(name + ' | analysis summary',lines)
    if mode != 'noe':
        methods = [
            'Saved model: I(t) = I0 exp(-R t). Parameters and reported uncertainties are unchanged.',
            'Original-noise chi2 = sum(((I - prediction)/input SD)^2); residual df = N - 2.',
            'Approximate upper-tail chi2 p-values; Holm family adjustment at alpha = 0.05.',
            'Primary coverage: valid, non-overlapping fits with positive residual df.',
            'Overlap, invalid fits and other exclusions are listed on results pages and in the CSV.',
            'Quality flags do not change compatibility filtering: only overlap and invalid fits are excluded.',
            'Assumes independent Gaussian measurement errors with correctly supplied per-plane SDs.',
            'Probabilities approximate nonlinear-fit goodness of fit; supplied noise is not re-estimated.',
            'Cross-peak independence is not required by Holm; shared patterns are screened descriptively.',
            'Shared-plane screen: >=10 eligible peaks, >=80% same residual sign, |median residual/SD| > 2.',
            'This screen is descriptive, not an independent hypothesis test. Delay summaries pool repeats.',
            'Repeat scatter uses inverse-variance weighted means and sum((I - mean)^2/SD^2), df = k - 1.',
            'Pair contrast = (I_first - I_second)/sqrt(SD_first^2 + SD_second^2).',
            'Few repeats give limited power and cannot establish noise calibration or drift on their own.',
            'Neither rejection nor uncertainty rescaling establishes a biexponential mechanism.',
            'Residual inspection complements statistical tests; supplied noise levels condition conclusions.',
            'Decay panels use fitted-I0 normalization; zero/nonfinite I0 receives a flagged raw panel.',
            'Main plot and residual plot areas use an 80:20 height ratio.',
            'References: NIST Engineering Statistics Handbook, model validation (4.4).',
            'https://itl.nist.gov/div898/handbook/pmd/section4/pmd44.htm',
            'Holm method: R stats p.adjust documentation (valid under arbitrary dependence).',
            'https://www.stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html']
        commands += text_page('Methods and assessment coverage',methods)
        eligible = [r for r in rows if r['primary_eligible']]
        if eligible:
            commands += list(AUTO_PLOT_MARGINS) + ['set title "Original-noise goodness of fit (eligible peaks)"',
                        'set xlabel "Residue (eligible peaks)"', 'set ylabel "Original-noise reduced chi2"',
                        'set autoscale', 'set key off', 'plot "-" using 1:2 with points pt 7 lc rgb "#285879"']
            commands += [f"{r['residue']} {r['reduced_chi2_input']:.17g}" for r in eligible] + ['e']
            for kind in ('plane','delay'):
                summaries = ([(i+1,s) for i,s in enumerate(report['planes'])] if kind=='plane' else sorted(report['delays'].items()))
                commands += ['set title ' + gnuplot_string('Shared residual summaries by acquisition plane' if kind=='plane' else 'Shared residual summaries by delay'),
                             'set xlabel ' + gnuplot_string('Acquisition plane (original order)' if kind=='plane' else 'Delay (s); repeated delays pooled'),
                             'set ylabel "Residual / input SD (median and IQR)"', 'set autoscale', 'set key top right',
                             'plot "-" using 1:2:3:4 with yerrorbars pt 7 title "Median, IQR", 2 dt 2 title "+2 SD", -2 dt 2 title "-2 SD"']
                # Pad the axis without conflating original-order planes and delays.
                # Insert the range before this plot's command and inline data.
                span = summaries[-1][0] - summaries[0][0] or 1.
                commands.insert(len(commands)-1, f'set xrange [{summaries[0][0]-.04*span:.17g}:{summaries[-1][0]+.04*span:.17g}]')
                commands += [f"{x:.17g} {s['median']:.17g} {s['q1']:.17g} {s['q3']:.17g}" for x,s in summaries]+['e']
            repeated = [r for r in eligible if r['repeat_dof']]
            if repeated:
                commands += text_page('Repeated-delay scatter (eligible peaks)',[
                    f"Repeat degrees of freedom per peak: {repeated[0]['repeat_dof']}",
                    f"Median weighted within-repeat chi2: {statistics.median(r['repeat_chi2'] for r in repeated):.5g}",
                    'Per-peak totals and each repeated-delay group are included in the CSV.',
                    'Each repeated measurement contributes separately to the decay goodness-of-fit statistic.',
                    'Few repeats limit interpretation; noise estimates and fitted uncertainties remain unchanged.'])
    else:
        commands += text_page('NOE ratio methods',[
            'Exponential suitability: not applicable. Ratios and their fitter uncertainties are retained.',
            'Per-plane intensities, supplied noise SDs, peak identity and trace fields are included in the CSV.',
            'Compatibility filtering excludes overlap and invalid fits. No decay curves are calculated.'])
    for offset in range(0,len(rows),16):
        lines = []
        for row in rows[offset:offset+16]:
            eligible = 'primary' if row['primary_eligible'] else 'excluded: '+row['exclusion_reason']
            lines.append(f"{row['peak_id']} | {row['assignment']} | {display_number(row['values'][0])} +/- {display_number(row['values'][1])} | {eligible}")
            detail = ('; '.join(row['flags']) or 'no flags')
            if mode != 'noe': detail = f"original reduced chi2={display_number(row['reduced_chi2_input'])}; Holm p={display_number(row['p_holm'],3)}; " + detail
            lines.append(detail)
        # 32 lines need a denser page; split into eight peaks per text page.
        for page in range(0,len(lines),16):
            commands += text_page('Peak results | values +/- fitter errors',lines[page:page+16])
    return commands



def fitted_rows(peaks, fit_text, mode, warning_map=None):
    records, warnings = parse_fit(fit_text, mode)
    if warning_map is not None:
        if set(warning_map) != set(records):
            raise ValueError('Warning and fit identities differ')
        warnings = {key: text.splitlines() for key, text in warning_map.items()}
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
    return rows


def prepare_reports(stage, peaks, extracted, fit_text, times, noises, mode, prefix, pdf_name, run_tag, context=None, warning_map=None):
    rows = fitted_rows(peaks, fit_text, mode, warning_map)
    context = dict(context) if context is not None else report_context(mode, prefix, pdf_name, run_tag)
    context['report_revision'] = code_revision()
    report = build_assessment(rows, extracted, times, noises, context)
    context = report['context']
    tag_header = f'# {run_tag}\n' + ''.join(f'# {key}: {context[key]}\n' for key in TRACE_FIELDS if key in context)
    csv_name = prefix + '_analysis.csv'
    write_analysis_csv(stage / csv_name, report)
    names = [csv_name]
    for suffix, include in (('results.out', lambda r: True),
                            ('results_peak_overlap.out', lambda r: r['overlap']),
                            ('results_filtered.out', lambda r: not r['overlap'] and r['valid'])):
        name = f'{prefix}_{suffix}'
        text = tag_header + ''.join(f"{r['residue']} {r['values'][0]:.12g} {r['values'][1]:.12g}\n" for r in rows if include(r))
        (stage / name).write_text(text)
        names.append(name)
    diagnostics = f'{prefix}_fit_diagnostics.tsv'
    fields = ['peak_id', 'assignment', 'w1', 'w2', 'w3', 'residue', 'value', 'error', 'i0', 'i0_error',
              'chi2', 'chi2_threshold', 'overlap', 'duplicate_assignment', 'valid_fit', 'flags', 'warnings'] + TRACE_FIELDS
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
            info.update(context)
            writer.writerow(info)
    names.append(diagnostics)
    plot = [tag_header.rstrip(), 'set terminal pdfcairo noenhanced color size 6,5 font "Sans,9"',
            'set output ' + gnuplot_string(pdf_name), 'set grid', 'set key off', 'set autoscale']
    plot += report_intro(report)
    plot += ['unset title', 'set key off', 'set autoscale']
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
            label += '\n' + fit_annotation(row, mode)
            if row['flags']:
                label += '\n' + textwrap.fill(', '.join(row['flags']), width=76)
            span = max(times) - min(times) or max(abs(max(times)), .01)
            lo, hi = min(times) - .03 * span, max(times) + .03 * span
            plot += [f'set xrange [{lo:.17g}:{hi:.17g}]', 'set multiplot title ' + gnuplot_string(label) + ' font "Sans,8"']
            plot += ['set autoscale y'] + decay_panel_layout(top=.82, bottom=.11)
            plot.append('set ylabel ' + gnuplot_string('Intensity / fitted I0' if normalized else 'Raw intensity (invalid fit/I0)'))
            command = f'plot "intensity.out" index {index} using 1:2:3 with yerrorbars pt 7'
            if normalized:
                command += f', exp(-({rate:.17g})*x) with lines'
            plot += [command] + decay_panel_layout(residual=True, top=.82, bottom=.11) + ['set ylabel "Residual / SD"']
            if normalized:
                plot.append(f'plot "intensity.out" index {index} using 1:6 with points pt 7, 0 with lines')
            else:
                plot += ['set label 999 "Residual unavailable: invalid fit/I0" at graph 0.1,0.5',
                         f'plot [{lo:.17g}:{hi:.17g}] [-1:1] 0 with lines', 'unset label 999']
            plot.append('unset multiplot')
        blocks.append('\n'.join(block))
    (stage / 'intensity.out').write_text(tag_header + '\n\n\n'.join(blocks) + '\n')
    plot.extend(AUTO_PLOT_MARGINS)
    plot += ['set autoscale y', 'set format x "%g"']
    plot += ['set autoscale x', 'set title ' + gnuplot_string(mode.upper() + ' | ' + run_tag), 'set xlabel "Residue"',
             'set ylabel ' + gnuplot_string('hetNOE' if mode == 'noe' else mode.upper() + ' (1/s)'),
             'set key top right', 'set datafile missing "NaN"']
    groups = [(False, 'Isolated peaks', '#285879'), (True, 'Overlap flagged', '#c17b19')]
    groups = [(flag, label, color) for flag, label, color in groups
              if any(row['valid'] and row['overlap'] == flag for row in rows)]
    if not groups:
        plot += text_page('No valid fitted values', ['All saved fits are invalid; exponential suitability is not applicable.' if mode=='noe' else 'All saved fits are invalid; assessment is inconclusive.'])
    if groups: plot.append('plot ' + ', '.join(
        '\'-\' using 1:2:3 with yerrorbars pt 7 lc rgb ' + gnuplot_string(color) + ' title ' + gnuplot_string(label)
        for flag, label, color in groups))
    for flag, label, color in groups:
        plot.extend(f"{row['residue']} {row['values'][0]:.17g} {row['values'][1]:.17g}"
                    for row in rows if row['valid'] and row['overlap'] == flag)
        plot.append('e')
    # NOE receives provenance/results pages and a ratio overview, without decay panels.
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
            command += ['-r_init', '1.0' if args.mode in ('r1', 'r1rho') else '10.0']
        checked_run(command + ['-out', 'rlx_dfit.out'], stage)
        nonempty(stage / 'rlx_dfit.out')
        if not any(row['valid'] for row in fitted_rows(peaks,(stage/'rlx_dfit.out').read_text(),args.mode)):
            raise ValueError('No valid fitted peaks; previous results retained')
        names = prepare_reports(stage, peaks, extracted, (stage / 'rlx_dfit.out').read_text(),
                                times, noises, args.mode, prefix, pdf_name, run_tag,
                                context=report_context(args.mode, prefix, pdf_name, run_tag,
                                    peak_list=peak_path.name, fit_origin='analysis', fit_revision=code_revision(),
                                    fit_input_sha256=file_sha256(stage/'rlx_dfit.inp'),
                                    fit_output_sha256=file_sha256(stage/'rlx_dfit.out'),
                                    noise_sha256=file_sha256(stage/'vEstNoise.txt')))
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
    p.add_argument('--format', default='nmrpipe'); p.add_argument('--mode', choices=('r1', 'r1rho', 'r2', 'noe'), required=True)
    p.add_argument('--peak-format', choices=('sparky2d', 'sparky3d'), default='sparky3d')
    p.add_argument('--header', required=True); p.add_argument('--noise', required=True); p.add_argument('--prefix', default='')
    p = commands.add_parser('report', help='Rebuild reports without extraction, processing or fitting')
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument('--csv'); group.add_argument('--saved-prefix')
    p.add_argument('--experiment-type', choices=('r1','r1rho','r2','noe'))
    p.add_argument('--peak-list-name'); p.add_argument('--output-dir')
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
        elif args.command == 'report':
            refresh_report(args)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        print(f'Relaxation run failed: {error}', file=sys.stderr)
        return 1
    return 0


def handle_termination(signum, frame):
    raise InterruptedError('Run interrupted by SIGTERM')


if __name__ == '__main__':
    signal.signal(signal.SIGTERM, handle_termination)
    sys.exit(main())
