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
import uuid

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
REPORT_SCHEMA = '2'
REPORT_VERSION = '2'
EARLY_SCHEMA = REPORT_SCHEMA
EARLY_DEFAULTS = dict(min_peaks=20, min_remaining_delays=5, residual_floor=.005,
                      noise_multiple=3., sign_fraction=.8, score_reduction=.5,
                      max_median_rate_change=.2)
ASSESSMENT = dict(alpha=.05, coherence_min_peaks=10, coherence_sign_fraction=.8,
                  coherence_median_sd=2.)
CONCLUSIONS = {
    'departures_detected': 'Departures detected at supplied noise levels',
    'no_evidence_against': 'No evidence against single-exponential decay at supplied noise levels',
    'inconclusive': 'Inconclusive',
    'not_applicable': 'Not applicable: NOE ratio analysis',
}
LEGACY_TRACE_FIELDS = ['schema_version', 'analysis_id', 'source_collection', 'source_experiment',
                'peak_list', 'output_prefix', 'pdf_name', 'experiment_type', 'fit_model',
                'fit_origin', 'fit_revision', 'report_revision', 'fit_input_sha256',
                'fit_output_sha256', 'noise_sha256']
TRACE_FIELDS = ['schema_version', 'analysis_id', 'source_collection', 'source_experiment',
                'peak_list', 'output_prefix', 'pdf_name', 'experiment_type', 'fit_model',
                'fit_origin', 'report_version']
CSV_FIELDS = TRACE_FIELDS + [
    'peak_index', 'peak_id', 'assignment', 'residue', 'w1', 'w2', 'w3', 'overlap',
    'duplicate_assignment', 'valid_fit', 'plane_index', 'delay_s', 'plane_code',
    'intensity', 'noise_sd', 'prediction', 'residual', 'residual_noise_sd',
    'value', 'value_error', 'value_units', 'i0', 'i0_error', 'fitter_chi2', 'fitter_threshold',
    'uncertainty_rescaled', 'warnings', 'flags', 'chi2_input', 'reduced_chi2_input', 'dof',
    'p_input', 'p_holm', 'primary_eligible', 'exclusion_reason', 'repeat_chi2', 'repeat_dof',
    'repeat_pair_z', 'repeat_delay_chi2', 'repeat_delay_dof', 'repeat_delay_pair_z', 'dataset_conclusion', 'primary_peak_count', 'holm_rejection_count',
    'coherent_planes', 'alpha', 'coherence_min_peaks', 'coherence_sign_fraction',
    'coherence_median_sd', 'plane_median_sd', 'plane_q1_sd', 'plane_q3_sd',
    'plane_sign_fraction', 'plane_bias', 'delay_median_sd', 'delay_q1_sd', 'delay_q3_sd']
EARLY_TRACE_FIELDS = TRACE_FIELDS + ['early_point_audit_json']
FIT_VALUE_FIELDS = ('value', 'error', 'i0', 'i0_error', 'chi2', 'chi2_threshold')
EARLY_CSV_FIELDS = CSV_FIELDS + ['early_point_audit_json', 'used_in_fit', 'fit_basis',
    'selection_reason', 'rate_change_percent'] + [
    group+'_'+field for group in ('all_points', 'candidate') for field in FIT_VALUE_FIELDS] + [
    'all_points_warnings', 'candidate_warnings']


def early_groups(times, noises, extracted, cohort, fitted, normalization):
    """One weighted repeat mean per residue and delay; a fixed all-fit I0 scale."""
    result = []
    for delay in sorted(set(times)):
        indices = [i for i, time in enumerate(times) if time == delay]
        weights = [1. / noises[i]**2 for i in indices]
        total = math.fsum(weights)
        means, errors = [], []
        for tag in cohort:
            amplitude = normalization[tag]
            mean = math.fsum(extracted[tag][i]*w for i,w in zip(indices,weights))/total
            rate, _, i0 = fitted[tag]['values'][:3]
            prediction = i0*math.exp(-rate*delay)
            means.append((mean-prediction)/amplitude)
            errors.append(1./math.sqrt(total)/amplitude)
        median = statistics.median(means)
        result.append(dict(delay_s=delay, planes_1based=[i+1 for i in indices], median=median,
                           normalized_noise=statistics.median(errors),
                           sign_fraction=sum(v*median > 0 for v in means)/len(means)))
    return result


def evaluate_early_point(peaks, extracted, full_rows, times, noises, policy='auto',
                         candidate_rows=None, config=None):
    """Deterministic selection heuristic. Never fits or discards measurements."""
    if policy not in ('auto','report','off'):
        raise ValueError('Invalid early-point policy')
    cfg = dict(EARLY_DEFAULTS if config is None else config)
    if set(cfg) != set(EARLY_DEFAULTS) or any(not math.isfinite(v) for v in cfg.values()):
        raise ValueError('Invalid early-point thresholds')
    if (cfg['min_peaks'] < 1 or int(cfg['min_peaks']) != cfg['min_peaks'] or
        cfg['min_remaining_delays'] < 3 or int(cfg['min_remaining_delays']) != cfg['min_remaining_delays'] or
        cfg['residual_floor'] <= 0 or cfg['noise_multiple'] <= 0 or
        not .5 <= cfg['sign_fraction'] <= 1 or not 0 < cfg['score_reduction'] <= 1 or
        cfg['max_median_rate_change'] <= 0):
        raise ValueError('Invalid early-point thresholds')
    full = {r['peak_id']: r for r in full_rows}
    if set(full) != {p['peak_id'] for p in peaks}:
        raise ValueError('Full-fit identifiers differ from peaks')
    counts = Counter(p['residue'] for p in peaks)
    cohort = [r['peak_id'] for r in full_rows if r['valid'] and not r['overlap'] and
              not r['duplicate_assignment'] and counts[r['residue']] == 1 and
              r['values'][0] > 0 and r['values'][2] > 0]
    shortest = min(times)
    proposed = [i+1 for i,t in enumerate(times) if t == shortest]
    remaining = [i for i,t in enumerate(times) if t != shortest]
    audit = dict(version=1, policy=policy, thresholds=cfg, fit_basis='all_points',
                 reason='disabled' if policy == 'off' else 'pending', candidate_attempted=False,
                 candidate_passed=False, cohort_peak_ids=cohort, cohort_size=len(cohort),
                 candidate_delay_s=shortest, candidate_planes_1based=proposed,
                 excluded_delay_s=None, excluded_planes_1based=[], original_count=len(times),
                 fit_count=len(times), remaining_distinct_delays=len({times[i] for i in remaining}),
                 full_bins=[], candidate_bins=[], before_score=None, after_score=None,
                 score_reduction_fraction=None, median_absolute_rate_change_fraction=None)
    if policy == 'off': return audit
    if len(cohort) < cfg['min_peaks']:
        audit['reason'] = 'insufficient_cohort'; return audit
    if audit['remaining_distinct_delays'] < cfg['min_remaining_delays']:
        audit['reason'] = 'insufficient_remaining_delays'; return audit
    normalization = {tag: full[tag]['values'][2] for tag in cohort}
    bins = early_groups(times,noises,extracted,cohort,full,normalization)
    for group in bins:
        group['limit'] = max(cfg['residual_floor'],cfg['noise_multiple']*group['normalized_noise'])
    audit['full_bins'] = bins
    first = bins[0]
    if abs(first['median']) <= first['limit']:
        audit['reason'] = 'shortest_residual_within_limit'; return audit
    if first['sign_fraction'] < cfg['sign_fraction']:
        audit['reason'] = 'insufficient_sign_agreement'; return audit
    if candidate_rows is None:
        audit['reason'] = 'candidate_required'; return audit
    audit['candidate_attempted'] = True
    candidate = {r['peak_id']:r for r in candidate_rows}
    if set(candidate) != set(full):
        raise ValueError('Candidate-fit identifiers differ from full fit')
    if any(not candidate[tag]['valid'] or candidate[tag]['values'][0] <= 0 or
           candidate[tag]['values'][2] <= 0 for tag in cohort):
        audit['reason'] = 'invalid_candidate_cohort_fit'; return audit
    # Only retained observations enter either score. The original cohort,
    # noise and normalization remain fixed across both comparisons.
    later = early_groups([times[i] for i in remaining], [noises[i] for i in remaining],
                         {tag:[y[i] for i in remaining] for tag,y in extracted.items()},
                         cohort,candidate,normalization)
    limits = {b['delay_s']:b['limit'] for b in bins}
    for group in later:
        group['limit'] = limits[group['delay_s']]
        group['planes_1based'] = [i+1 for i,t in enumerate(times) if t == group['delay_s']]
    before = [b for b in bins if b['delay_s'] != shortest]
    score0 = math.sqrt(math.fsum((b['median']/b['limit'])**2 for b in before)/len(before))
    score1 = math.sqrt(math.fsum((b['median']/b['limit'])**2 for b in later)/len(later))
    change = statistics.median(abs(candidate[tag]['values'][0]/full[tag]['values'][0]-1) for tag in cohort)
    reduction = 1-score1/score0 if score0 else 0.
    audit.update(candidate_bins=later, before_score=score0, after_score=score1,
                 score_reduction_fraction=reduction, median_absolute_rate_change_fraction=change)
    if reduction < cfg['score_reduction']:
        audit['reason'] = 'insufficient_structure_improvement'
    elif any(abs(b['median']) > b['limit'] for b in later):
        audit['reason'] = 'remaining_systematic_residual'
    elif change > cfg['max_median_rate_change']:
        audit['reason'] = 'rate_change_exceeds_limit'
    else:
        audit['candidate_passed'] = True
        audit['reason'] = 'selected_shortest_exclusion' if policy == 'auto' else 'report_only_candidate_passed'
        if policy == 'auto':
            audit.update(fit_basis='exclude_shortest', excluded_delay_s=shortest,
                         excluded_planes_1based=proposed, fit_count=len(remaining))
    return audit


def selection_rows(selection, peaks, mode):
    full = fitted_rows(peaks,selection['all_fit_text'],mode,selection.get('all_warning_map'))
    candidate = (fitted_rows(peaks,selection['candidate_fit_text'],mode,selection.get('candidate_warning_map'))
                 if selection.get('candidate_fit_text') is not None else None)
    return full, candidate


def validate_selection(selection, peaks, extracted, times, noises, mode, chosen_rows):
    full, candidate = selection_rows(selection,peaks,mode)
    saved = selection['audit']
    calculated = evaluate_early_point(peaks,extracted,full,times,noises,saved['policy'],candidate,saved['thresholds'])
    saved = dict(saved)
    saved.pop('archive_sha256', None)  # Obsolete schema-1/early-schema-2 input only.
    selection['audit'] = saved
    if calculated != saved:
        raise ValueError('Early-point decision disagrees with saved measurements/fits')
    expected = candidate if saved['fit_basis'] == 'exclude_shortest' else full
    for actual,wanted in zip(chosen_rows,expected):
        if actual['peak_id'] != wanted['peak_id'] or actual['warnings'] != wanted['warnings'] or any(
                not close_saved(a,b) for a,b in zip(actual['values'],wanted['values'])):
            raise ValueError('Selected fitted parameters disagree with early-point decision')
    for i,row in enumerate(chosen_rows):
        row['all_points_fit'] = full[i]
        row['candidate_fit'] = candidate[i] if candidate is not None else None
    return saved


def write_fit_input(path, peaks, extracted, times, indices):
    Path(path).write_text(str(len(indices))+' '+' '.join(f'{times[i]:.17g}' for i in indices)+'\n'+
        ''.join(p['peak_id']+' '+' '.join(f'{extracted[p["peak_id"]][i]:.17g}' for i in indices)+'\n' for p in peaks))


def fit_with_early_test(stage, peaks, extracted, times, noises, mode, policy='auto'):
    """Run the native fitter once or twice, and retain its exact artifacts."""
    stage = Path(stage)
    def run(stem, indices):
        write_fit_input(stage/(stem+'.inp'),peaks,extracted,times,indices)
        noise_name = stem+'_noise.txt'
        (stage/noise_name).write_text(str(len(indices))+' '+' '.join(f'{noises[i]:.17g}' for i in indices)+'\n')
        command = ['rlx_dfit','-mode','exp2','-rms_noise','-data',stem+'.inp',
                   '-d_esd',noise_name,'-r_init','1.0' if mode in ('r1','r1rho') else '10.0','-out',stem+'.out']
        checked_run(command,stage)
        nonempty(stage/(stem+'.out'))
        text = (stage/(stem+'.out')).read_text()
        return text, fitted_rows(peaks,text,mode)
    full_text, full = run('rlx_dfit_all_points',list(range(len(times))))
    if not any(r['valid'] for r in full):
        raise ValueError('No valid fitted peaks; previous results retained')
    audit = evaluate_early_point(peaks,extracted,full,times,noises,policy)
    candidate_text = None
    stems = ['rlx_dfit_all_points']
    if audit['reason'] == 'candidate_required':
        candidate_text, candidate = run('rlx_dfit_without_shortest',[i for i,t in enumerate(times) if t != min(times)])
        stems.append('rlx_dfit_without_shortest')
        audit = evaluate_early_point(peaks,extracted,full,times,noises,policy,candidate)
    else:
        # Replace any prior generation's candidate artifacts with explicit
        # non-results, so stale candidate files cannot appear current.
        for suffix in ('.inp','.out','_noise.txt'):
            (stage/('rlx_dfit_without_shortest'+suffix)).write_text('# Candidate not attempted: '+audit['reason']+'\n')
        stems.append('rlx_dfit_without_shortest')
    selected = 'rlx_dfit_without_shortest' if audit['fit_basis']=='exclude_shortest' else 'rlx_dfit_all_points'
    for suffix in ('.inp','.out','_noise.txt'):
        shutil.copyfile(stage/(selected+suffix),stage/('rlx_dfit'+suffix))
    names = [stem+suffix for stem in stems for suffix in ('.inp','.out','_noise.txt')]+['rlx_dfit_noise.txt']
    selection = dict(audit=audit,all_fit_text=full_text,candidate_fit_text=candidate_text)
    return selection, names


# Schema 1 is accepted only as an input format. Its obsolete trace fields and
# methods prose are never carried into new reports and are not used as hashes.
LEGACY_CSV_FIELDS = LEGACY_TRACE_FIELDS + CSV_FIELDS[len(TRACE_FIELDS):]
LEGACY_CSV_FIELDS.insert(LEGACY_CSV_FIELDS.index('dataset_conclusion'), 'assessment_method')
LEGACY_EARLY_CSV_FIELDS = LEGACY_CSV_FIELDS + EARLY_CSV_FIELDS[len(CSV_FIELDS):]


def validate_analysis_id(value):
    try:
        parsed = uuid.UUID(value)
    except (ValueError, AttributeError, TypeError) as error:
        raise ValueError('Invalid analysis UUID') from error
    if parsed.version != 4 or str(parsed) != value:
        raise ValueError('Analysis ID must be a canonical UUID4')
    return value



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
                   fit_origin='saved-results', analysis_id=str(uuid.uuid4()),
                   report_version=REPORT_VERSION)
    context.update(extra)
    context['report_version'] = REPORT_VERSION
    return context


def build_assessment(rows, extracted, times, noises, context, selection=None):
    n = len(times)
    if not n or len(noises) != n or any(not math.isfinite(v) or v <= 0 for v in noises):
        raise ValueError('Invalid measurement/noise counts')
    mode = context['experiment_type']
    audit = None
    used = [True]*n
    if selection is not None:
        if mode == 'noe': raise ValueError('NOE cannot have a delay-exclusion decision')
        audit = validate_selection(selection, [dict(r) for r in rows], extracted, times, noises, mode, rows)
        used = [i+1 not in audit['excluded_planes_1based'] for i in range(n)]
        context = dict(context, schema_version=EARLY_SCHEMA,
                       early_point_audit_json=json.dumps(audit,sort_keys=True,separators=(',',':')))
        if audit['fit_basis'] == 'exclude_shortest':
            for row in rows: row['flags'].append('shortest_delay_excluded')
    if mode not in ('r1', 'r1rho', 'r2', 'noe') or any(not math.isfinite(t) for t in times):
        raise ValueError('Invalid report experiment type or plane values')
    if mode != 'noe' and any(t < 0 for t in times):
        raise ValueError('Negative decay delay')
    if set(extracted) != {r['peak_id'] for r in rows}:
        raise ValueError('Report measurements and peak identities differ')
    groups = {}
    for i, time in enumerate(times):
        groups.setdefault(time, [])
        if used[i]: groups[time].append(i)
    primary = []
    for row in rows:
        y = extracted[row['peak_id']]
        if len(y) != n or any(not math.isfinite(v) for v in y):
            raise ValueError('Invalid saved intensity count or value')
        row.update(points=[], chi2_input=None, reduced_chi2_input=None, p_input=None,
                   p_holm=None, dof=sum(used)-2 if mode != 'noe' else None, repeat_chi2=None,
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
                chi2 = math.fsum(point['z']**2 for point,use in zip(row['points'],used) if use)
                if not math.isfinite(chi2):
                    raise ValueError('Nonfinite original-noise chi-square')
                row['chi2_input'] = chi2
                if sum(used) > 2:
                    row['reduced_chi2_input'] = chi2 / (sum(used)-2)
                    row['p_input'] = chi_square_sf(chi2, sum(used)-2)
        reasons = []
        if mode == 'noe': reasons.append('not_decay_analysis')
        if not row['valid']: reasons.append('invalid_fit')
        if row['overlap']: reasons.append('overlap')
        if mode != 'noe' and sum(used) <= 2: reasons.append('no_positive_residual_dof')
        row['primary_eligible'] = not reasons
        row['exclusion_reason'] = ';'.join(reasons)
        if not reasons: primary.append(row)
    adjusted = holm_adjust({r['peak_id']: r['p_input'] for r in primary})
    for row in primary:
        row['p_holm'] = adjusted[row['peak_id']]
        if row['p_holm'] < ASSESSMENT['alpha']:
            row['flags'].append('original_noise_lack_of_fit')
    planes = [residual_summary([r['points'][i]['z'] for r in primary] if used[i] else []) for i in range(n)]
    delays = {t: residual_summary([r['points'][i]['z'] for r in primary for i in indices])
              for t, indices in groups.items()}
    coherent = [i+1 for i, summary in enumerate(planes) if summary['bias']]
    rejected = sum(r['p_holm'] < ASSESSMENT['alpha'] for r in primary)
    conclusion = ('not_applicable' if mode == 'noe' else 'inconclusive' if not primary else
                  'departures_detected' if rejected or coherent else 'no_evidence_against')
    context = dict(context)
    return dict(context=context, rows=rows, times=list(map(float,times)), noises=list(map(float,noises)), planes=planes,
                delays=delays, conclusion=conclusion, primary_count=len(primary),
                rejected=rejected, coherent=coherent, selection=selection, used=used)


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
        if report['selection'] is not None:
            audit=report['selection']['audit']
            fields.update(fit_basis=audit['fit_basis'],selection_reason=audit['reason'])
            for group,key in [('all_points','all_points_fit'),('candidate','candidate_fit')]:
                fit=row[key]
                fields.update({group+'_'+key:fit['values'][j] if fit else None
                               for j,key in enumerate(FIT_VALUE_FIELDS)})
                fields[group+'_warnings']=fit['warnings'] if fit else ''
            full=row['all_points_fit']['values'][0]
            cand=row['candidate_fit']['values'][0] if row['candidate_fit'] else None
            fields['rate_change_percent']=(100*(cand/full-1) if cand is not None and
                math.isfinite(full) and math.isfinite(cand) and full != 0 else None)
        for i, (time, point) in enumerate(zip(report['times'], row['points'])):
            plane, delay = report['planes'][i], report['delays'][time]
            repeat = row['repeat_groups'].get(time,{})
            yield dict(fields, used_in_fit=report['used'][i], plane_index=i+1, delay_s=time if mode != 'noe' else None,
                       plane_code=time if mode == 'noe' else None, intensity=point['intensity'],
                       noise_sd=point['noise'], prediction=point['prediction'], residual=point['residual'],
                       residual_noise_sd=point['z'], plane_median_sd=plane['median'],
                       plane_q1_sd=plane['q1'], plane_q3_sd=plane['q3'],
                       plane_sign_fraction=plane['sign_fraction'], plane_bias=plane['bias'],
                       repeat_delay_chi2=repeat.get('chi2'), repeat_delay_dof=repeat.get('dof'),
                       repeat_delay_pair_z=repeat.get('pair_z'),
                       delay_median_sd=delay['median'], delay_q1_sd=delay['q1'], delay_q3_sd=delay['q3'])


def write_analysis_csv(path, report):
    with Path(path).open('w', newline='') as handle:
        fields = EARLY_CSV_FIELDS if report['selection'] is not None else CSV_FIELDS
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in analysis_csv_rows(report):
            writer.writerow({field: row.get(field) for field in fields})


NUMBERS_BASE_FIELDS = [
    'analysis_id', 'report_version', 'source_collection', 'source_experiment',
    'experiment_type', 'peak_list', 'fit_model', 'fit_origin', 'dataset_conclusion',
    'primary_peak_count', 'holm_rejection_count', 'peak_index', 'peak_id', 'assignment',
    'residue', 'w1', 'w2', 'w3', 'overlap', 'duplicate_assignment', 'valid_fit',
    'primary_eligible', 'exclusion_reason', 'rate_s-1', 'rate_error_s-1', 'i0', 'i0_error',
    'chi2_input', 'reduced_chi2_input', 'dof', 'p_holm', 'uncertainty_rescaled',
    'flags', 'warnings', 'fit_basis', 'selection_reason', 'excluded_shortest_delay_s',
    'excluded_planes_1based', 'all_points_rate_s-1', 'all_points_error_s-1',
    'without_shortest_rate_s-1', 'without_shortest_error_s-1', 'rate_change_percent']


def write_numbers_csv(path, report):
    """Write one complete analysis row per peak, including every measured plane."""
    context = report['context']
    mode = context['experiment_type']
    selection = report['selection']
    audit = selection['audit'] if selection is not None else {}
    value_name='noe_ratio' if mode == 'noe' else 'rate_s-1'
    error_name='noe_ratio_error' if mode == 'noe' else 'rate_error_s-1'
    fields = [value_name if f == 'rate_s-1' else error_name if f == 'rate_error_s-1' else f
              for f in NUMBERS_BASE_FIELDS]
    fields += [f'plane_{i:02d}_{suffix}' for i in range(1, len(report['times']) + 1)
               for suffix in ('delay_s', 'intensity', 'noise_sd', 'used_in_fit')]
    base_rows = {row['peak_id']: row for row in analysis_csv_rows(report)
                 if row['plane_index'] == 1}
    with Path(path).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for peak_index, peak in enumerate(report['rows'], 1):
            source = base_rows[peak['peak_id']]
            row = dict(
                analysis_id=context['analysis_id'], report_version=context['report_version'],
                source_collection=context['source_collection'],
                source_experiment=context['source_experiment'], experiment_type=mode,
                peak_list=context['peak_list'], fit_model=context['fit_model'],
                fit_origin=context['fit_origin'], dataset_conclusion=report['conclusion'],
                primary_peak_count=report['primary_count'], holm_rejection_count=report['rejected'],
                peak_index=peak_index, peak_id=peak['peak_id'], assignment=peak['assignment'],
                residue=peak['residue'], w1=peak['coordinates'][0],
                w2=peak['coordinates'][1] if len(peak['coordinates']) > 1 else '',
                w3=peak['coordinates'][2] if len(peak['coordinates']) > 2 else '',
                overlap=peak['overlap'], duplicate_assignment=peak['duplicate_assignment'],
                valid_fit=peak['valid'], primary_eligible=peak['primary_eligible'],
                exclusion_reason=peak['exclusion_reason'], i0=peak['values'][2] if mode != 'noe' else '',
                i0_error=peak['values'][3] if mode != 'noe' else '',
                chi2_input=peak['chi2_input'], reduced_chi2_input=peak['reduced_chi2_input'],
                dof=peak['dof'], p_holm=peak['p_holm'],
                uncertainty_rescaled='uncertainty_rescaled' in peak['flags'],
                flags=';'.join(peak['flags']), warnings=peak['warnings'],
                fit_basis=audit.get('fit_basis', 'all_points'),
                selection_reason=audit.get('reason', ''),
                excluded_shortest_delay_s=audit.get('excluded_delay_s', ''),
                excluded_planes_1based=';'.join(map(str, audit.get('excluded_planes_1based', []))),
                rate_change_percent=source.get('rate_change_percent', ''))
            row[value_name] = peak['values'][0]
            row[error_name] = peak['values'][1]
            if selection is not None:
                full=peak.get('all_points_fit'); candidate=peak.get('candidate_fit')
                row['all_points_rate_s-1'] = full['values'][0] if full else ''
                row['all_points_error_s-1'] = full['values'][1] if full else ''
                row['without_shortest_rate_s-1'] = candidate['values'][0] if candidate else ''
                row['without_shortest_error_s-1'] = candidate['values'][1] if candidate else ''
            for index, (time, point, used) in enumerate(zip(report['times'], peak['points'], report['used']), 1):
                row.update({f'plane_{index:02d}_delay_s': time,
                            f'plane_{index:02d}_intensity': point['intensity'],
                            f'plane_{index:02d}_noise_sd': point['noise'],
                            f'plane_{index:02d}_used_in_fit': used})
            writer.writerow(row)


def write_support_archive(path, prefix, stage, native_source, native_names, internal_names):
    """Bundle technical records into one optional support file."""
    import zipfile
    archive_path=Path(path)
    readme=(f'Relaxation analysis support files for {prefix}.\n'
            'The adjacent *_numbers.csv is the one-row-per-peak analysis table; the PDF is the report.\n'
            'This archive keeps the saved fitter inputs/outputs, warnings, extraction and plotting records\n'
            'needed to audit or rebuild the report. Raw spectra, peak lists and CCPN projects\n'
            'are not included.\n')
    with zipfile.ZipFile(archive_path,'w',compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr('README.txt',readme)
        seen={'README.txt'}
        for name in list(internal_names)+list(native_names):
            name=Path(name).name
            if name in seen or name.endswith('_numbers.csv') or name.endswith('_intensity.pdf'):
                continue
            candidate=Path(stage)/name
            if not candidate.is_file(): candidate=Path(native_source)/name
            if candidate.is_file():
                archive.write(candidate,arcname=name)
                seen.add(name)


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
        input_fields = reader.fieldnames
        if input_fields not in (CSV_FIELDS, EARLY_CSV_FIELDS, LEGACY_CSV_FIELDS, LEGACY_EARLY_CSV_FIELDS):
            raise ValueError('Unsupported or incomplete analysis CSV schema')
        schema_fields = reader.fieldnames
        early = schema_fields in (EARLY_CSV_FIELDS, LEGACY_EARLY_CSV_FIELDS)
        records = list(reader)
    if not records or any(None in row or None in row.values() for row in records):
        raise ValueError('Empty or incomplete analysis CSV')
    legacy = input_fields in (LEGACY_CSV_FIELDS, LEGACY_EARLY_CSV_FIELDS)
    trace_fields = (LEGACY_TRACE_FIELDS if legacy else TRACE_FIELDS) + (['early_point_audit_json'] if early else [])
    if records[0]['schema_version'] != ('2' if legacy and early else '1' if legacy else REPORT_SCHEMA):
        raise ValueError('Unsupported analysis CSV schema version')
    context = {key: records[0][key] for key in (EARLY_TRACE_FIELDS if early else TRACE_FIELDS) if key in records[0]}
    if legacy:
        context.update(schema_version=REPORT_SCHEMA, analysis_id='', report_version=REPORT_VERSION)
    else:
        validate_analysis_id(context['analysis_id'])
    safe_output_name(context['output_prefix']); safe_output_name(context['pdf_name'])
    if context['fit_model'] != ('noe' if context['experiment_type'] == 'noe' else 'exp2'):
        raise ValueError('CSV model/experiment type mismatch')
    grouped = {}
    for record in records:
        if any(record[field] != records[0][field] for field in trace_fields):
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
    extra_fit_lines = {'all_points':[], 'candidate':[]}
    extra_warnings = {'all_points':{}, 'candidate':{}}
    varying = {'plane_index','delay_s','plane_code','intensity','noise_sd','prediction','residual',
               'residual_noise_sd','plane_median_sd','plane_q1_sd','plane_q3_sd',
               'plane_sign_fraction','plane_bias','delay_median_sd','delay_q1_sd','delay_q3_sd',
               'repeat_delay_chi2','repeat_delay_dof','repeat_delay_pair_z', 'used_in_fit'}
    for index, slots in sorted(grouped.items()):
        if sorted(slots) != list(range(1,len(slots)+1)):
            raise ValueError('Missing CSV measurement plane')
        first = slots[1]
        for record in slots.values():
            if any(record[key] != first[key] for key in input_fields if key not in varying):
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
        if early:
            for group in extra_fit_lines:
                parts=[first[group+'_'+key] for key in FIT_VALUE_FIELDS]
                if any(parts) and not all(parts):
                    raise ValueError('Incomplete saved early-point fitted parameters')
                if all(parts):
                    extra_fit_lines[group].append(tag+' '+' '.join(str(float(v)) for v in parts))
                    extra_warnings[group][tag]=first[group+'_warnings']
    counts = Counter(p['assignment'] for p in peaks)
    if any(p['duplicate_assignment'] != (counts[p['assignment']]>1) for p in peaks):
        raise ValueError('CSV duplicate-assignment flags disagree with identities')
    bundle = dict(peaks=peaks, extracted=extracted, fit_text='\n'.join(fit_lines)+'\n',
                  times=list(map(float,times)), noises=list(map(float,noises)), context=context, warning_map=warning_map)
    # Recompute diagnostics and verify every exported numerical value to reject
    # inconsistent derived fields.
    selection = None
    if early:
        audit=json.loads(context['early_point_audit_json'])
        audit.pop('archive_sha256', None)
        context['early_point_audit_json']=json.dumps(audit,sort_keys=True,separators=(',',':'))
        if len(extra_fit_lines['all_points']) != len(peaks):
            raise ValueError('Missing all-point reference fits')
        if bool(extra_fit_lines['candidate']) != audit['candidate_attempted']:
            raise ValueError('Saved candidate presence disagrees with audit')
        selection=dict(audit=audit,all_fit_text='\n'.join(extra_fit_lines['all_points'])+'\n',
                       candidate_fit_text=('\n'.join(extra_fit_lines['candidate'])+'\n' if extra_fit_lines['candidate'] else None),
                       all_warning_map=extra_warnings['all_points'],
                       candidate_warning_map=extra_warnings['candidate'] if extra_fit_lines['candidate'] else None)
        bundle['context']['_early_point_selection']=selection
    rows = fitted_rows(peaks,bundle['fit_text'],context['experiment_type'],warning_map)
    computed = build_assessment(rows,extracted,times,noises,context,selection)
    expected = {(str(r['peak_index']),str(r['plane_index'])): r for r in analysis_csv_rows(computed)}
    for record in records:
        wanted = expected[(record['peak_index'],record['plane_index'])]
        for key in (EARLY_CSV_FIELDS if early else CSV_FIELDS):
            if legacy and key in ('schema_version','analysis_id','report_version'): continue
            actual = record[key]; value = wanted.get(key)
            if legacy and key=='early_point_audit_json':
                previous=json.loads(actual);previous.pop('archive_sha256',None)
                if previous!=json.loads(value):
                    raise ValueError('CSV early-point audit disagrees with saved measurements/fit')
                continue
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
    if peak_list_name is None:
        outputs=[line for line in (source/'plot.gnu').read_text().splitlines() if line.startswith('set output ')]
        if len(outputs)!=1:raise ValueError('Cannot identify saved peak-list name; provide --peak-list-name')
        pdf_name=shlex.split(outputs[0])[2]
        suffix='_intensity.pdf'
        if not pdf_name.endswith(suffix):raise ValueError('Cannot identify saved peak-list name; provide --peak-list-name')
        peak_list_name=pdf_name[:-len(suffix)]
    safe_output_name(peak_list_name)
    inferred_mode = infer_experiment_type(peak_list_name)
    if mode and mode != inferred_mode:
        raise ValueError(f'--experiment-type {mode} conflicts with peak-list name; inferred {inferred_mode}')
    mode = inferred_mode
    diagnostic = source / (prefix+'_fit_diagnostics.tsv')
    with diagnostic.open(newline='') as handle:
        mappings = list(csv.DictReader(handle,delimiter='\t'))
    if not mappings: raise ValueError('Empty saved diagnostics mapping')
    if mappings[0].get('early_point_audit_json'):
        return load_selected_saved_analysis(source,prefix,mode,peak_list_name,mappings)
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
    context=report_context(mode,prefix,peak_list_name+'_intensity.pdf',source.parent.name+'/'+source.name,
                           peak_list=peak_list_name)
    if 'schema_version' in mappings[0]:
        schema=mappings[0]['schema_version']
        if schema not in ('1', REPORT_SCHEMA):
            raise ValueError('Unsupported saved diagnostics schema')
        trace_fields=LEGACY_TRACE_FIELDS if schema=='1' else TRACE_FIELDS
        if any(any(key not in row for key in trace_fields) for row in mappings):
            raise ValueError('Incomplete saved diagnostics provenance')
        previous={key:mappings[0][key] for key in trace_fields}
        if any(any(row[key]!=previous[key] for key in trace_fields) for row in mappings):
            raise ValueError('Inconsistent saved diagnostics provenance')
        if previous['fit_model']!=context['fit_model'] or previous['output_prefix']!=prefix:
            raise ValueError('Saved diagnostics model/prefix mismatch')
        context.update({key:previous[key] for key in ('source_collection','source_experiment','fit_origin')})
        if schema==REPORT_SCHEMA:
            context['analysis_id']=validate_analysis_id(previous['analysis_id'])
    return dict(peaks=peaks,extracted=extracted,fit_text=fit_text,times=times,noises=noises,context=context)


def early_archive_names():
    return [stem+suffix for stem in ('rlx_dfit_all_points','rlx_dfit_without_shortest')
            for suffix in ('.inp','.out','_noise.txt')] + ['rlx_dfit_noise.txt']


def validate_early_archives(source, bundle, mode):
    """Compare saved native values to the CSV, including both archived fits."""
    source=Path(source);peaks=bundle['peaks'];times=bundle['times'];noises=bundle['noises']
    selection=bundle['context']['_early_point_selection'];audit=selection['audit']
    def validate(stem, indices, fit_text, warning_map):
        lines=meaningful_lines(source/(stem+'.inp'))
        if not lines:
            raise ValueError('Missing archived native fitter input: '+stem)
        head=lines[0].split()
        if int(head[0])!=len(indices) or list(map(float,head[1:]))!=[times[i] for i in indices]:
            raise ValueError('Archived native delay vector disagrees with portable data: '+stem)
        measurements=parse_extracted('\n'.join(lines[1:]),peaks,len(indices))
        if any(measurements[tag]!=[bundle['extracted'][tag][i] for i in indices] for tag in measurements):
            raise ValueError('Archived native measurements disagree with portable data: '+stem)
        if counted_vector(source/(stem+'_noise.txt'),positive=True)!=[noises[i] for i in indices]:
            raise ValueError('Archived native noise disagrees with portable data: '+stem)
        native=fitted_rows(peaks,(source/(stem+'.out')).read_text(),mode)
        reference=fitted_rows(peaks,fit_text,mode,warning_map)
        for actual,wanted in zip(native,reference):
            if actual['warnings']!=wanted['warnings'] or any(not close_saved(a,b) for a,b in zip(actual['values'],wanted['values'])):
                raise ValueError('Archived native fit disagrees with portable data: '+stem)
    all_indices=list(range(len(times)))
    retained=[i for i in all_indices if i+1 not in audit['excluded_planes_1based']]
    validate('rlx_dfit',retained,bundle['fit_text'],bundle['warning_map'])
    validate('rlx_dfit_all_points',all_indices,selection['all_fit_text'],selection.get('all_warning_map'))
    if audit['candidate_attempted']:
        candidate=[i for i in all_indices if times[i]!=min(times)]
        validate('rlx_dfit_without_shortest',candidate,selection['candidate_fit_text'],selection.get('candidate_warning_map'))
    else:
        expected='# Candidate not attempted: '+audit['reason']+'\n'
        for suffix in ('.inp','.out','_noise.txt'):
            if (source/('rlx_dfit_without_shortest'+suffix)).read_text()!=expected:
                raise ValueError('Unexpected stale candidate archive')


def load_selected_saved_analysis(source,prefix,mode,peak_list_name,mappings):
    """Validate native selected artifacts against the portable, full-observation CSV."""
    bundle=load_analysis_csv(source/(prefix+'_analysis.csv'))
    context=bundle['context'];selection=context['_early_point_selection'];audit=selection['audit']
    if context['experiment_type'] == 'noe' or mode == 'noe':
        raise ValueError('Selected decay cannot be relabelled as NOE')
    with (source/(prefix+'_analysis.csv')).open(newline='') as handle:
        original=next(csv.DictReader(handle))
    traces=(LEGACY_TRACE_FIELDS if 'report_version' not in original else TRACE_FIELDS)+['early_point_audit_json']
    for row in mappings:
        if any(row.get(key) != original.get(key) for key in traces):
            raise ValueError('Saved diagnostics provenance disagrees with CSV')
    if [r['peak_id'] for r in mappings] != [p['peak_id'] for p in bundle['peaks']]:
        raise ValueError('Saved diagnostics peak identities disagree with CSV')
    records,_=parse_fit(bundle['fit_text'],mode)
    for row in mappings:
        values=records[row['peak_id']]
        if any(not close_saved(float(row[key]),value) for key,value in zip(
                ['value','error','i0','i0_error','chi2','chi2_threshold'],values)):
            raise ValueError('Saved diagnostics parameters disagree with CSV')
        peak=next(p for p in bundle['peaks'] if p['peak_id']==row['peak_id'])
        if (row['assignment'] != peak['assignment'] or int(row['residue']) != peak['residue'] or
            strict_bool(row['overlap']) != peak['overlap'] or
            strict_bool(row['duplicate_assignment']) != peak['duplicate_assignment'] or
            [row[f'w{i}'] for i in (1,2,3) if row.get(f'w{i}')] != peak['coordinates']):
            raise ValueError('Saved diagnostics peak mapping disagrees with CSV')
    validate_early_archives(source,bundle,mode)
    if counted_vector(source/'relax_header.txt') != bundle['times'] or counted_vector(
            source/'vEstNoise.txt',positive=True) != bundle['noises']:
        raise ValueError('Original processing vectors disagree with portable data')
    used=[i for i in range(len(bundle['times'])) if i+1 not in audit['excluded_planes_1based']]
    inp=meaningful_lines(source/'rlx_dfit.inp')
    head=inp[0].split()
    if int(head[0]) != len(used) or list(map(float,head[1:])) != [bundle['times'][i] for i in used]:
        raise ValueError('Selected native delay vector disagrees with mask')
    actual=parse_extracted('\n'.join(inp[1:]),bundle['peaks'],len(used))
    if any(actual[tag] != [bundle['extracted'][tag][i] for i in used] for tag in actual):
        raise ValueError('Selected native measurements disagree with mask')
    if counted_vector(source/'rlx_dfit_noise.txt',positive=True) != [bundle['noises'][i] for i in used]:
        raise ValueError('Selected native noise vector disagrees with mask')
    # Reconstruct the plotting table from saved values only, and compare bytes.
    # This detects mixed generations without refitting or touching the source.
    with tempfile.TemporaryDirectory(prefix='relax-verify-') as temp:
        stage=Path(temp)
        prepare_reports(stage,bundle['peaks'],bundle['extracted'],bundle['fit_text'],bundle['times'],
                        bundle['noises'],mode,prefix,context['pdf_name'],
                        context['source_collection']+'/'+context['source_experiment'],context=context,
                        warning_map=bundle['warning_map'])
        def numeric_blocks(path):
            return [list(map(float,l.split())) for l in Path(path).read_text().splitlines()
                    if l.strip() and not l.lstrip().startswith('#')]
        a,b=numeric_blocks(stage/'intensity.out'),numeric_blocks(source/'intensity.out')
        if len(a)!=len(b) or any(len(x)!=len(y) or any(not close_saved(u,v) for u,v in zip(x,y)) for x,y in zip(a,b)):
            raise ValueError('Saved plotting table disagrees with selected fit')
    if peak_list_name and peak_list_name != context['peak_list']:
        raise ValueError('Saved peak-list name disagrees with portable data')
    context['experiment_type']=mode
    return bundle


def extract_support_archive(archive_path, destination):
    import zipfile
    destination=Path(destination)
    with zipfile.ZipFile(archive_path) as archive:
        for info in archive.infolist():
            name=Path(info.filename)
            if name.is_absolute() or len(name.parts)!=1 or name.name in ('','.','..'):
                raise ValueError('Unsafe support archive member')
            if not info.is_dir():
                archive.extract(info,destination)


def refresh_report(args):
    extraction=None
    if args.csv:
        input_path=Path(args.csv).resolve()
        bundle=load_analysis_csv(input_path)
        source=input_path.parent
        native_source=None
        watched=[input_path]
        native_names=[]
        if args.experiment_type and args.experiment_type!=bundle['context']['experiment_type']:
            raise ValueError('CSV experiment type is fixed; use saved-result refresh to relabel legacy fits')
    else:
        source=Path.cwd().resolve()
        support=source/(args.saved_prefix+'_support.zip')
        if (source/'rlx_dfit.inp').is_file():
            native_source=source
            native_names=['rlx_dfit.inp','rlx_dfit.out','vEstNoise.txt','relax_header.txt',
                          'data_xtr.out']+early_archive_names()
            watched=[source/name for name in native_names if (source/name).is_file()]
            watched += [source/(args.saved_prefix+'_fit_diagnostics.tsv'),source/'intensity.out',source/'plot.gnu']
        elif support.is_file():
            extraction=tempfile.TemporaryDirectory(prefix='relax-support-')
            native_source=Path(extraction.name)
            extract_support_archive(support,native_source)
            native_names=[p.name for p in native_source.iterdir() if p.is_file()]
            watched=[support]
        else:
            raise ValueError('No saved fitter files or support archive found')
        bundle=load_saved_analysis(native_source,args.saved_prefix,args.experiment_type,args.peak_list_name)
    destination=Path(args.output_dir).resolve() if args.output_dir else source
    if not destination.is_dir():raise ValueError('Report output directory must already exist')
    try:
        with run_directory(destination,'report') as stage:
            before={p:p.read_bytes() for p in watched}
            if args.csv:
                bundle=load_analysis_csv(input_path)
            elif native_source == source:
                bundle=load_saved_analysis(source,args.saved_prefix,args.experiment_type,args.peak_list_name)
            else:
                bundle=load_saved_analysis(native_source,args.saved_prefix,args.experiment_type,args.peak_list_name)
            c=bundle['context'];mode=c['experiment_type']
            internal_names=prepare_reports(stage,bundle['peaks'],bundle['extracted'],bundle['fit_text'],
                              bundle['times'],bundle['noises'],mode,c['output_prefix'],c['pdf_name'],
                              c['source_collection']+'/'+c['source_experiment'],context=c,
                              warning_map=bundle.get('warning_map'))
            checked_run(['gnuplot','plot.gnu'],stage)
            nonempty(stage/c['pdf_name'])
            if any(p.read_bytes()!=contents for p,contents in before.items()):
                raise ValueError('Saved report inputs changed during rendering')
            numbers_name=c['output_prefix']+'_numbers.csv'
            support_name=c['output_prefix']+'_support.zip'
            write_support_archive(stage/support_name,c['output_prefix'],stage,native_source or source,
                                  native_names,internal_names)
            with (stage/(c['output_prefix']+'_analysis.csv')).open(newline='') as handle:
                report_id=next(csv.DictReader(handle))['analysis_id']
            outputs=[c['pdf_name'],numbers_name,support_name]
            existing=[name for name in outputs if (destination/name).is_file()]
            if existing:
                stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
                backup=destination/'.report-backups'/report_id/stamp
                backup.mkdir(parents=True)
                for name in existing:shutil.copy2(destination/name,backup/name)
            publish(stage,destination,outputs)
            print(f'Report {report_id}: {destination/numbers_name} (one row per peak)',flush=True)
    finally:
        if extraction is not None:
            extraction.cleanup()



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


def early_report_intro(selection):
    audit=selection['audit']
    commands=[]
    if audit['candidate_bins']:
        before={b['delay_s']:b for b in audit['full_bins']}
        commands+=list(AUTO_PLOT_MARGINS)+['set title "Shared residuals on the same retained delays"',
             'set xlabel "Delay (s)"','set ylabel "Median residual / fixed all-point I0 (%)"',
             'set autoscale','set key top right',
             'plot "-" using 1:2 with linespoints title "All-point fit", "-" using 1:2 with linespoints title "Shortest-delay candidate", "-" using 1:2 with lines dt 2 title "+limit", "-" using 1:2 with lines dt 2 title "-limit"']
        bins=audit['candidate_bins']
        for values in ([before[b['delay_s']]['median'] for b in bins],
                       [b['median'] for b in bins],[b['limit'] for b in bins],[-b['limit'] for b in bins]):
            commands += [f"{b['delay_s']:.17g} {v*100:.17g}" for b,v in zip(bins,values)]+['e']
    return commands


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
    footer = c['source_collection'] + '/' + c['source_experiment'] + ' | report v' + c['report_version']
    commands = ['set label 9000 ' + gnuplot_string(footer) + ' at screen .02,.022 font "Sans,5.5" front']
    lines = [CONCLUSIONS[report['conclusion']],
             'Dataset: ' + c['source_collection'] + '/' + c['source_experiment'],
             'Peak list: ' + c['peak_list'],
             f"{len(rows)} peaks; {len(report['times'])} acquisition planes; {len(set(report['times']))} distinct plane values"]
    if mode != 'noe':
        lines += [f"Assessment: {report['primary_count']} eligible; {len(rows)-report['primary_count']} excluded",
                  f"Holm-adjusted departures (alpha=0.05): {report['rejected']}",
                  f"Descriptive shared-plane flags: {report['coherent'] or 'none'}",
                  'Assessment uses supplied noise SDs and assumes independent Gaussian measurement noise.',
                  'Shared-plane screen: >=10 eligible peaks, >=80% same sign, |median residual/SD| > 2.',
                  'This is a descriptive screen. Departures do not establish a biexponential mechanism.']
    if report.get('selection'):
        audit=report['selection']['audit']
        lines.insert(4,f"Fit basis: {audit['fit_basis']}; observations retained {audit['fit_count']}/{audit['original_count']}")
        if not all(report['used']):
            lines[0]+=' (retained observations; conditional on selection)'
            lines.append('Open orange markers are excluded observations; all measurements remain in the CSV.')
    commands += text_page(name + ' | analysis summary',lines)
    if mode != 'noe':
        eligible = [r for r in rows if r['primary_eligible']]
        if eligible:
            commands += list(AUTO_PLOT_MARGINS) + ['set title "Original-noise goodness of fit (eligible peaks)"',
                        'set xlabel "Residue (eligible peaks)"', 'set ylabel "Original-noise reduced chi2"',
                        'set autoscale', 'set key off', 'plot "-" using 1:2 with points pt 7 lc rgb "#285879"']
            commands += [f"{r['residue']} {r['reduced_chi2_input']:.17g}" for r in eligible] + ['e']
            for kind in ('plane','delay'):
                summaries = ([(i+1,s) for i,s in enumerate(report['planes'])] if kind=='plane' else sorted(report['delays'].items()))
                summaries = [(x,s) for x,s in summaries if s['median'] is not None]
                commands += ['set title ' + gnuplot_string('Shared residual summaries by acquisition plane' if kind=='plane' else 'Shared residual summaries by delay'),
                             'set xlabel ' + gnuplot_string('Acquisition plane (original order)' if kind=='plane' else 'Delay (s); repeated delays pooled'),
                             'set ylabel "Residual / input SD (median and IQR)"', 'set autoscale', 'set key top right',
                             'plot "-" using 1:2:3:4 with yerrorbars pt 7 title "Median, IQR", 2 dt 2 title "+2 SD", -2 dt 2 title "-2 SD"']
                # Pad the axis without conflating original-order planes and delays.
                # Insert the range before this plot's command and inline data.
                span = summaries[-1][0] - summaries[0][0] or 1.
                commands.insert(len(commands)-1, f'set xrange [{summaries[0][0]-.04*span:.17g}:{summaries[-1][0]+.04*span:.17g}]')
                commands += [f"{x:.17g} {s['median']:.17g} {s['q1']:.17g} {s['q3']:.17g}" for x,s in summaries]+['e']
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


def prepare_reports(stage, peaks, extracted, fit_text, times, noises, mode, prefix, pdf_name, run_tag, context=None, warning_map=None, early_selection=None):
    rows = fitted_rows(peaks, fit_text, mode, warning_map)
    context = dict(context) if context is not None else report_context(mode, prefix, pdf_name, run_tag)
    early_selection = early_selection or context.pop('_early_point_selection',None)
    context = {key:value for key,value in context.items() if key in TRACE_FIELDS}
    context['schema_version'] = REPORT_SCHEMA
    context['analysis_id'] = validate_analysis_id(context.get('analysis_id') or str(uuid.uuid4()))
    context['report_version'] = REPORT_VERSION
    report = build_assessment(rows, extracted, times, noises, context, early_selection)
    context = report['context']
    traces=EARLY_TRACE_FIELDS if early_selection is not None else TRACE_FIELDS
    tag_header = f'# {run_tag}\n' + ''.join(f'# {key}: {context[key]}\n' for key in traces if key in context)
    if early_selection is not None:
        audit=early_selection['audit']
        tag_header+=''.join(f'# {key}: {audit[key]}\n' for key in ('policy','fit_basis','reason',
            'excluded_delay_s','excluded_planes_1based','original_count','fit_count'))
    csv_name = prefix + '_analysis.csv'
    write_analysis_csv(stage / csv_name, report)
    numbers_name = prefix + '_numbers.csv'
    write_numbers_csv(stage / numbers_name, report)
    names = [csv_name, numbers_name]
    if early_selection is not None:
        summary_name=prefix+'_early_point_summary.json'
        (stage/summary_name).write_text(json.dumps(dict(context=context,selection=early_selection['audit']),
                                                  sort_keys=True,indent=2,allow_nan=False)+'\n')
        names.append(summary_name)
    diagnostics = f'{prefix}_fit_diagnostics.tsv'
    fields = ['peak_id', 'assignment', 'w1', 'w2', 'w3', 'residue', 'value', 'error', 'i0', 'i0_error',
              'chi2', 'chi2_threshold', 'overlap', 'duplicate_assignment', 'valid_fit', 'flags', 'warnings'] + traces
    extra_fields = ['fit_basis','selection_reason','rate_change_percent']+[
        group+'_'+field for group in ('all_points','candidate') for field in FIT_VALUE_FIELDS]+[
        'all_points_warnings','candidate_warnings'] if early_selection is not None else []
    fields += extra_fields
    csv_extra={r['peak_id']:r for r in analysis_csv_rows(report) if r['plane_index']==1}
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
            info.update({key:context[key] for key in traces})
            info.update({key:csv_extra[row['peak_id']][key] for key in extra_fields})
            writer.writerow(info)
    names.append(diagnostics)
    plot = [tag_header.rstrip(), 'set terminal pdfcairo noenhanced color size 6,5 font "Sans,9"',
            'set output ' + gnuplot_string(pdf_name), 'set grid', 'set key off', 'set autoscale']
    plot += report_intro(report)
    if early_selection is not None: plot += early_report_intro(early_selection)
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
            for plane, (t, y, sd) in enumerate(zip(times, values, noises)):
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
                if early_selection is not None:
                    block[-1] += ' '+str(int(report['used'][plane]))
            label = row['assignment'] + ' [' + row['peak_id'] + ']'
            label += '\n' + fit_annotation(row, mode)
            if early_selection is not None:
                label += '\nFit basis: '+early_selection['audit']['fit_basis']
                if not all(report['used']):
                    label += f"; excluded {early_selection['audit']['excluded_delay_s']:.6g} s (open orange)"
            if row['flags']:
                label += '\n' + textwrap.fill(', '.join(row['flags']), width=76)
            span = max(times) - min(times) or max(abs(max(times)), .01)
            lo, hi = min(times) - .03 * span, max(times) + .03 * span
            plot += [f'set xrange [{lo:.17g}:{hi:.17g}]', 'set multiplot title ' + gnuplot_string(label) + ' font "Sans,8"']
            plot += ['set autoscale y'] + decay_panel_layout(top=.82, bottom=.11)
            plot.append('set ylabel ' + gnuplot_string('Intensity / fitted I0' if normalized else 'Raw intensity (invalid fit/I0)'))
            if early_selection is not None and not all(report['used']):
                command = (f'plot "intensity.out" index {index} using 1:($7==1?$2:1/0):3 with yerrorbars pt 7 lc rgb "#285879", '
                           f'"intensity.out" index {index} using 1:($7==0?$2:1/0):3 with yerrorbars pt 6 lc rgb "#c17b19"')
            else:
                command = f'plot "intensity.out" index {index} using 1:2:3 with yerrorbars pt 7'
            if normalized:
                command += f', exp(-({rate:.17g})*x) with lines'
            plot += [command] + decay_panel_layout(residual=True, top=.82, bottom=.11) + ['set ylabel "Residual / SD"']
            if normalized:
                if early_selection is not None and not all(report['used']):
                    plot.append(f'plot "intensity.out" index {index} using 1:($7==1?$6:1/0) with points pt 7 lc rgb "#285879", '
                                f'"intensity.out" index {index} using 1:($7==0?$6:1/0) with points pt 6 lc rgb "#c17b19", 0 with lines')
                else:
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


def infer_experiment_type(peak_list_name):
    """Infer the relaxation experiment from a peak-list filename token."""
    name = Path(peak_list_name).name.lower()
    # Match the specific T1rho spelling before the shorter T1 token. Separators
    # are filename punctuation; alphanumeric adjacency prevents accidental
    # matches inside sample names or unrelated numbers.
    patterns = (
        ('r1rho', r'(?<![a-z0-9])(?:t1rho|r1rho)(?![a-z0-9])'),
        ('noe', r'(?<![a-z0-9])(?:het[-_ ]?)?noe(?![a-z0-9])'),
        ('r1', r'(?<![a-z0-9])(?:t1|r1)(?![a-z0-9])'),
        ('r2', r'(?<![a-z0-9])(?:t2|r2)(?![a-z0-9])'),
    )
    found = [mode for mode, pattern in patterns if re.search(pattern, name)]
    if len(found) != 1:
        raise ValueError('Cannot infer one experiment type from peak-list name '
                         f'{Path(peak_list_name).name!r}; include one of T1, T1rho, T2 or NOE')
    return found[0]


def analyze(args):
    source = Path.cwd().resolve()
    matches = sorted(Path(p).resolve() for p in glob.glob(args.peaks) if Path(p).is_file())
    if len(matches) != 1:
        raise ValueError('Select exactly one peak list; set SPARKY_PEAKS to its filename')
    peak_path = matches[0]
    mode = infer_experiment_type(peak_path.name)
    if args.mode and args.mode != mode:
        raise ValueError(f'--mode {args.mode} conflicts with peak-list name; inferred {mode}')
    data_path, header_path, noise_path = (Path(p).resolve() for p in (args.data, args.header, args.noise))
    times = counted_vector(header_path, nonnegative=mode != 'noe')
    noises = counted_vector(noise_path, positive=True)
    if len(times) != len(noises):
        raise ValueError('Delay and noise counts differ')
    if args.format == 'nmrpipe' and plane_count(data_path) != len(times):
        raise ValueError('Spectrum plane, delay and noise counts differ')
    if mode != 'noe' and len(set(times)) < 2:
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
        early_selection=None
        artifact_names=[]
        if mode == 'noe':
            checked_run(['rlx_dfit','-mode','noe','-rms_noise','-data','rlx_dfit.inp',
                         '-d_esd','vEstNoise.txt','-out','rlx_dfit.out'],stage)
            nonempty(stage/'rlx_dfit.out')
        else:
            early_selection, artifact_names=fit_with_early_test(stage,peaks,extracted,times,noises,
                                                               mode,args.early_point_policy)
        if not any(row['valid'] for row in fitted_rows(peaks,(stage/'rlx_dfit.out').read_text(),mode)):
            raise ValueError('No valid fitted peaks; previous results retained')
        names = prepare_reports(stage, peaks, extracted, (stage / 'rlx_dfit.out').read_text(),
                                times, noises, mode, prefix, pdf_name, run_tag,
                                context=report_context(mode, prefix, pdf_name, run_tag,
                                    peak_list=peak_path.name, fit_origin='analysis'),
                                early_selection=early_selection)
        restored = ''.join(p['assignment'] + ' ' + ' '.join(f'{v:.17g}' for v in extracted[p['peak_id']]) + '\n' for p in peaks)
        (stage / 'data_xtr.out').write_text(restored)
        checked_run(['gnuplot', 'plot.gnu'], stage)
        nonempty(stage / pdf_name)
        numbers_name=prefix+'_numbers.csv'
        support_name=prefix+'_support.zip'
        native_names=['data_xtr.out','rlx_dfit.inp','rlx_dfit.out','vEstNoise.txt',
                      'relax_header.txt']+artifact_names
        write_support_archive(stage/support_name,prefix,stage,stage,native_names,names)
        publish(stage,source,[pdf_name,numbers_name,support_name])


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
    p.add_argument('--format', default='nmrpipe'); p.add_argument('--mode', choices=('r1', 'r1rho', 'r2', 'noe'),
                    help='deprecated consistency check; type is inferred from the peak-list filename')
    p.add_argument('--peak-format', choices=('sparky2d', 'sparky3d'), default='sparky3d')
    p.add_argument('--header', required=True); p.add_argument('--noise', required=True); p.add_argument('--prefix', default='')
    p.add_argument('--early-point-policy', choices=('auto','report','off'),default='auto')
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
