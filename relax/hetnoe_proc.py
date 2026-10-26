#!/usr/bin/env python3
"""Process reviewed interleaved HetNOE acquisitions without new Bruker expnos.

Two matched signed planes: reference, saturated. See hetnoe_proc.txt.
"""
import argparse
import array
import math
from pathlib import Path
import re
import subprocess
import sys

import relax_tools as rt


ROLES = ('reference', 'saturated')
OUTPUTS = ['planes', 'm.ft3', 'proc.ucsf', 'vEstNoise.txt', 'relax_header.txt', 'processing.txt']


def scalars(path):
    """Read scalar JCAMP parameters; unrelated arrays may contain quoted spaces."""
    result = {}
    for key, value in re.findall(r'^##\$(\w+)=(.*?)(?=^##|\Z)', Path(path).read_text(), re.M | re.S):
        if key in result:
            raise ValueError(f'{path}: duplicate parameter {key}')
        value = '\n'.join(line.split('$$', 1)[0] for line in value.splitlines()).strip()
        if not value.startswith('('):
            result[key] = value.strip('<>')
    return result


def number(parameters, key, positive=False, integer=False):
    try:
        value = float(parameters[key])
    except (KeyError, ValueError) as error:
        raise ValueError(f'Missing or invalid acquisition parameter {key}') from error
    if not math.isfinite(value) or (positive and value <= 0) or (integer and value != int(value)):
        raise ValueError(f'Invalid acquisition parameter {key}: {value}')
    return int(value) if integer else value


def signature(path):
    st = Path(path).stat()
    return st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns, st.st_ctime_ns


def profile(folder, raw_rows=None):
    folder = Path(folder).resolve()
    acq, indirect = scalars(folder/'acqus'), scalars(folder/'acqu2s')
    program = acq.get('PULPROG')
    pulse = (folder/'pulseprogram').read_text()
    active = re.sub(r'\s+', '', ''.join(line.split(';', 1)[0] for line in pulse.splitlines()
                                      if not line.lstrip().startswith('#')))
    if program == '15N_HETNOE':
        required = ['"l0=1"', 'if"l0%2==1"{d1', '}else{d2',
                    '(p0ph1)5mloto4timesl4', 'wr#0if#0zdiu0',
                    'lotoLBLF1I1times2', '"ph1loop+=1"', 'lotoLBLF1times2',
                    '"t1loop+=1"', 'lotoLBLST1timesST1CNT', 'trunc(td1/(2*2))']
        ordered = ['wr#0if#0zdiu0', 'lotoLBLF1I1times2', '"ph1loop+=1"',
                   'lotoLBLF1times2', '"t1loop+=1"', 'lotoLBLST1timesST1CNT']
        counter = 'l0'
    elif program == 'N15HetNoe_rk.nl':
        required = ['"l2=1"', 'if"l2==1"goto9', '8'+'11m(p1*2ph0):f1',
                    'loto8timesl8', '9'+'11md8*211m', 'loto9timesl8',
                    'wr#0if#0zd', 'iu2', 'loto3times2', 'ip8*2', 'igradEA',
                    'ru2', 'loto4times2', 'id0', 'loto5timesl3']
        ordered = ['wr#0if#0zd', 'iu2', 'loto3times2', 'ip8*2', 'igradEA',
                   'ru2', 'loto4times2', 'id0', 'loto5timesl3']
        counter = 'l2'
    else:
        raise ValueError(f'Unsupported pulse program {program!r}; plane order has not been reviewed')
    if any(part not in active for part in required) or not re.search('.*?'.join(map(re.escape, ordered)), active):
        raise ValueError('Saved pulse program differs from the reviewed reference/saturation and echo/antiecho loops')
    if len(re.findall(r'"'+counter+r'=1"', active)) != 1 or re.search(r'"'+counter+r'=(?!=)(?!1")', active):
        raise ValueError('Unrecognised initial saturation counter')
    audit = folder/'audita.txt'
    if audit.exists() and re.search(r'created by split\s+\d+:', audit.read_text()):
        raise ValueError('This is a TopSpin split copy; use the original interleaved experiment')
    td = number(acq, 'TD', positive=True, integer=True)
    metadata_rows = number(indirect, 'TD', positive=True, integer=True)
    if (number(acq, 'PARMODE', integer=True) != 1 or number(acq, 'AQ_mod', integer=True) != 3 or
            number(indirect, 'FnMODE', integer=True) != 6 or td % 2 or number(acq, 'FnTYPE', integer=True) != 0):
        raise ValueError('Expected uniformly sampled packed 2D DQD / echo-antiecho acquisition')
    datatype = number(acq, 'DTYPA', integer=True)
    byte_order = number(acq, 'BYTORDA', integer=True)
    if datatype not in (0, 2) or byte_order not in (0, 1):
        raise ValueError('Supported Bruker raw types: int32 or float64, known byte order')
    word_bytes = 4 if datatype == 0 else 8
    row_bytes = ((td*word_bytes+1023)//1024)*1024
    size = (folder/'ser').stat().st_size
    observed, remainder = divmod(size, row_bytes)
    rows = metadata_rows if raw_rows is None else raw_rows
    if rows <= 0 or rows % 4:
        raise ValueError('Raw row count must be a positive multiple of four (two states x echo/antiecho)')
    if remainder or observed != rows:
        raise ValueError(f'Raw ser has {observed} full rows and {remainder} trailing bytes; acqu2s TD={metadata_rows}, '
                         f'selected rows={rows}. Review acquisition status; --raw-rows N explicitly selects the complete file, never truncates it')
    for parameters in (acq, indirect):
        for key in ('SFO1', 'SW_h'):
            number(parameters, key, positive=True)
        number(parameters, 'O1')
    for key in ('DECIM', 'DSPFVS'):
        number(acq, key, positive=True)
    number(acq, 'GRPDLY')
    return dict(folder=folder, acq=acq, indirect=indirect, program=program, td=td,
                rows=rows, metadata_rows=metadata_rows, row_bytes=row_bytes,
                datatype=datatype, byte_order=byte_order, raw_rows_override=raw_rows is not None)


def split_raw(p, directory):
    """Alternate whole FIDs before echo/antiecho decoding; never mix states.

    Float64 raw values are converted to native float32 for NMRPipe without
    intensity normalization. Padding is recreated for the 4-byte word size.
    """
    paths = [directory/f'plane_{i:02d}.ser' for i in (1, 2)]
    handles = []
    try:
        handles = [path.open('wb') for path in paths]
        with (p['folder']/'ser').open('rb') as source:
            for index in range(p['rows']):
                record = source.read(p['row_bytes'])
                if len(record) != p['row_bytes']:
                    raise ValueError('Raw ser changed while reading')
                if p['datatype'] == 2:
                    values = array.array('d'); values.frombytes(record[:p['td']*8])
                    if (sys.byteorder == 'little') != (p['byte_order'] == 0):
                        values.byteswap()
                    if any(not math.isfinite(v) or abs(v) > 3.4028234663852886e38 for v in values):
                        raise ValueError('Nonfinite or out-of-range float64 raw intensity')
                    record = array.array('f', values).tobytes()
                    record += bytes((-len(record)) % 1024)
                handles[index % 2].write(record)
            if source.read(1):
                raise ValueError('Raw ser grew while reading')
    finally:
        for handle in handles:
            handle.close()
    return paths


def pipeline(commands, directory):
    running = []
    try:
        for index, command in enumerate(commands):
            print('+ '+' '.join(map(str, command)), flush=True)
            previous = running[-1].stdout if running else None
            process = subprocess.Popen(list(map(str, command)), cwd=directory, stdin=previous,
                                       stdout=subprocess.PIPE if index+1 < len(commands) else subprocess.DEVNULL)
            if previous is not None:
                previous.close()
            running.append(process)
        for process in reversed(running):
            process.wait()
        for command, process in zip(commands, running):
            if process.returncode:
                raise subprocess.CalledProcessError(process.returncode, command)
    finally:
        for process in running:
            if process.poll() is None:
                process.terminate()
        for process in running:
            process.wait()


def processing_commands(fid, out, args):
    """POLY runs along the active last dimension: N before TP, H after TP."""
    return [['nmrPipe', '-in', fid, '-fn', 'SOL'],
            ['nmrPipe', '-fn', 'SP', '-off', '.4', '-end', '.98', '-pow', '2', '-c', '.5'],
            ['nmrPipe', '-fn', 'ZF', '-size', args.h_zf], ['nmrPipe', '-fn', 'FT'],
            ['nmrPipe', '-fn', 'EXT', '-left', '-sw'],
            ['nmrPipe', '-fn', 'PS', '-p0', args.h_p0, '-p1', args.h_p1, '-di'],
            ['nmrPipe', '-fn', 'TP'],
            ['nmrPipe', '-fn', 'SP', '-off', '.4', '-end', '.98', '-pow', '2', '-c', '.5'],
            ['nmrPipe', '-fn', 'ZF', '-size', args.n_zf], ['nmrPipe', '-fn', 'FT'],
            ['nmrPipe', '-fn', 'PS', '-p0', args.n_p0, '-p1', args.n_p1, '-di'],
            ['nmrPipe', '-fn', 'POLY', '-auto', '-ord', args.n_baseline_order],
            ['nmrPipe', '-fn', 'TP'],
            ['nmrPipe', '-fn', 'POLY', '-auto', '-ord', args.h_baseline_order, '-out', out, '-ov']]


def process(args):
    p = profile(args.experiment, args.raw_rows)
    for key in ('h_p0', 'h_p1', 'n_p0', 'n_p1'):
        if not math.isfinite(getattr(args, key)):
            raise ValueError('Phase settings must be finite')
    if args.h_zf < max(2, p['td']//2) or args.n_zf < max(2, p['rows']//4):
        raise ValueError('Zero-fill sizes must accommodate acquired complex points')
    if not 0 <= args.n_baseline_order <= 6 or not 0 <= args.h_baseline_order <= 6:
        raise ValueError('Baseline polynomial orders must be between zero and six')
    destination = Path(args.output_dir).resolve() if args.output_dir else p['folder']/'hetnoe_processed'
    if destination == p['folder'] or destination in p['folder'].parents:
        raise ValueError('Use a dedicated processing directory, not an acquisition or its parent')
    destination.mkdir(parents=True, exist_ok=True)
    inputs = [p['folder']/name for name in ('ser', 'acqus', 'acqu2s', 'pulseprogram')]
    if (p['folder']/'audita.txt').exists():
        inputs.append(p['folder']/'audita.txt')
    initial = {path: signature(path) for path in inputs}
    with rt.run_directory(destination, 'hetnoe-process') as stage:
        raw = stage/'raw'; raw.mkdir()
        raw_paths = split_raw(p, raw)
        planes = stage/'planes'; planes.mkdir()
        swap = '-noswap' if p['datatype'] == 2 or ((p['byte_order'] == 0) == (sys.byteorder == 'little')) else '-swap'
        for index, path in enumerate(raw_paths, 1):
            fid = stage/f'plane_{index:02d}.fid'
            command = ['bruk2pipe', '-in', path, '-bad', '0.0', '-ext', swap, '-AMX']
            if p['datatype'] == 2:
                command.append('-noi2f')
            command += ['-decim', str(number(p['acq'], 'DECIM')), '-dspfvs', str(number(p['acq'], 'DSPFVS')),
                        '-grpdly', str(number(p['acq'], 'GRPDLY')),
                        '-xN', str(p['td']), '-xT', str(p['td']//2),
                        '-yN', str(p['rows']//2), '-yT', str(p['rows']//4),
                        '-xMODE', 'DQD', '-yMODE', 'Echo-AntiEcho']
            for axis, parameters in [('x', p['acq']), ('y', p['indirect'])]:
                obs = number(parameters, 'SFO1', positive=True)
                command += ['-'+axis+'SW', str(number(parameters, 'SW_h', positive=True)),
                            '-'+axis+'OBS', str(obs), '-'+axis+'CAR', str(number(parameters, 'O1')/obs)]
            command += ['-xLAB', 'HN', '-yLAB', '15N', '-ndim', '2', '-aq2D', 'States', '-out', fid, '-ov']
            # bruk2pipe applies the Rance-Kay macro itself for Echo-AntiEcho.
            # Do not run a second shuffling macro or use the etaXY axis corrections.
            rt.checked_run(command, stage); rt.nonempty(fid)
            pipeline(processing_commands(fid, planes/f'plane_{index:02d}.ft2', args), stage)
        paths = [planes/f'plane_{i:02d}.ft2' for i in (1, 2)]
        headers = [rt.pipe_header(path) for path in paths]
        if any(h[9] != 2 or (h[99], h[219]) != (headers[0][99], headers[0][219]) for h in headers):
            raise ValueError('Processed planes have inconsistent grids')
        rt.set_series_metadata(paths)
        with (stage/'m.ft3').open('wb') as handle:
            rt.checked_run(['xyz2pipe', '-in', str(planes/'plane_%02d.ft2')], stage, handle)
        if rt.plane_count(stage/'m.ft3') != 2:
            raise ValueError('Processed stack must have two planes')
        rt.checked_run(['pipe2ucsf', 'm.ft3', 'proc.ucsf'], stage); rt.nonempty(stage/'proc.ucsf')
        noises = []
        for path in paths:
            rt.checked_run(['sethdr', path, '-ndim', '2', '-zN', '1', '-zT', '1', '-zMODE', 'Real'], stage)
            rt.checked_run(['setfdata', path, '-fileCount', '1'], stage)
            with (stage/'noise_stat.txt').open('w') as handle:
                rt.checked_run(['specStat.com', '-in', path, '-stat', 'vEstNoise'], stage, handle)
            noises.append(rt.noise_value(stage/'noise_stat.txt'))
        # rlx_dfit NOE convention: 0 = reference, 1 = saturated (not delays).
        (stage/'relax_header.txt').write_text('2 0 1\n')
        (stage/'vEstNoise.txt').write_text('2 '+' '.join(map(str, noises))+'\n')
        (stage/'processing.txt').write_text(
            f"source_folder: {p['folder'].parent.name}\nexpno: {p['folder'].name}\npulse_program: {p['program']}\n"
            f"plane_roles: {', '.join(ROLES)}\nplane_codes: 0, 1\n"
            f"raw_type: {'int32' if p['datatype'] == 0 else 'float64 -> native float32, no normalization'}\n"
            f"metadata_TD1: {p['metadata_rows']}\nraw_rows: {p['rows']}\nraw_rows_override: {p['raw_rows_override']}\n"
            f"h_phase: {args.h_p0}, {args.h_p1}\nn_phase: {args.n_p0}, {args.n_p1}\n"
            f"zero_fill: {args.h_zf}, {args.n_zf}\n"
            f"15N_baseline: POLY -auto -ord {args.n_baseline_order} (after N PS -di, before TP)\n"
            f"H_baseline: POLY -auto -ord {args.h_baseline_order} (after final TP)\n"
            'Whole-record alternating demultiplex; Echo-AntiEcho conversion (one Rance-Kay decode).\n'
            'One common phase/window/baseline recipe; SOL; SP .4/.98/2/.5; H left-half extraction.\n'
            'Signed intensities retained; no magnitude calculation, peak clipping or plane normalization.\n')
        if any(signature(path) != old for path, old in initial.items()):
            raise ValueError('Acquisition inputs changed during processing')
        rt.publish(stage, destination, OUTPUTS)
    print(f'Processed reference (plane 1) and saturated (plane 2): {destination}')
    print('Export one assigned peak list from the well-phased reference plane; use fixed coordinates for both.')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('inspect', 'process'))
    parser.add_argument('--experiment', default='.')
    parser.add_argument('--raw-rows', type=int, help='Explicit complete-file row count when status TD disagrees; never truncates')
    parser.add_argument('--output-dir')
    for key in ('h-p0', 'h-p1', 'n-p0', 'n-p1'):
        parser.add_argument('--'+key, type=float, default=0.)
    parser.add_argument('--h-zf', type=int, default=4096)
    parser.add_argument('--n-zf', type=int, default=1024)
    parser.add_argument('--n-baseline-order', type=int, default=4)
    parser.add_argument('--h-baseline-order', type=int, default=4)
    args = parser.parse_args(argv)
    try:
        if args.command == 'inspect':
            p = profile(args.experiment, args.raw_rows)
            print(f"{p['folder']}\nPulse program: {p['program']}\n"
                  f"Raw rows: {p['rows']} (acqu2s TD={p['metadata_rows']}; explicit override={p['raw_rows_override']})\n"
                  f"{p['rows']//4} complex N points per plane; echo/antiecho\n"
                  'Plane 1: reference; plane 2: saturated.\n'
                  'Plan: common phasing/windows, explicit 15N baseline correction, then H baseline correction.')
        else:
            process(args)
    except (OSError, ValueError, OverflowError, subprocess.SubprocessError) as error:
        print(f'HetNOE processing: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
