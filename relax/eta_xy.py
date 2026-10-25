#!/usr/bin/env python3
"""Four-pathway etaXY analysis and processing for nico_nnh_tr_noSE.

Inspect is read-only. Process and analyze run only when explicitly requested.
Keep this file beside relax_tools.py. See eta_xy.txt for the acquisition profile.
"""
import argparse
import csv
from datetime import datetime, timezone
from decimal import Decimal
import glob
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import zipfile

import relax_tools as rt


ROLES = ('auto_A', 'cross_AB', 'cross_BA', 'auto_B')
METHOD_URL = 'https://doi.org/10.1016/S1090-7807(02)00190-8'


def jcamp(path):
    text = Path(path).read_text()
    result = {}
    for match in re.finditer(r'^##\$(\w+)=(.*?)(?=^##|\Z)', text, re.M | re.S):
        key, value = match.groups()
        if key in result:
            raise ValueError(f'{path}: duplicate parameter {key}')
        value = '\n'.join(line.split('$$', 1)[0] for line in value.splitlines()).strip()
        array = re.match(r'^\(0\.\.(\d+)\)\s*(.*)', value, re.S)
        if array:
            value = array[2].split()
            if len(value) != int(array[1])+1:
                raise ValueError(f'{path}: incomplete array {key}')
        else:
            value = value.strip('<>')
        result[key] = value
    return result


def numeric(parameters, key, index=None, positive=False, integer=False):
    try:
        if index is not None and not isinstance(parameters[key], list):
            raise TypeError('Expected array parameter')
        value = parameters[key] if index is None else parameters[key][index]
        number = float(value)
    except (KeyError, IndexError, TypeError, ValueError) as error:
        raise ValueError(f'Missing or invalid acquisition parameter {key}[{index}]') from error
    if not math.isfinite(number) or (positive and number <= 0) or (integer and number != int(number)):
        raise ValueError(f'Invalid acquisition parameter {key}: {number}')
    return int(number) if integer else number


def experiment_profile(folder):
    folder = Path(folder).resolve()
    acq, indirect = jcamp(folder/'acqus'), jcamp(folder/'acqu2s')
    if acq.get('PULPROG') != 'nico_nnh_tr_noSE':
        raise ValueError('Supported pulse program is nico_nnh_tr_noSE; refusing to assume another sequence has the same equation')
    pulse = (folder/'pulseprogram').read_text()
    active = ''.join(line.split(';', 1)[0] for line in pulse.splitlines()
                     if not line.lstrip().startswith('#'))
    active = re.sub(r'\s+', '', active)
    required = ['"l5=1"', 'if"l5%4==1"goto20', 'if"l5%4==2"goto20',
                'if"l5%4==1"goto40', 'if"l5%4==3"goto40',
                'lotoLBLF1I1times4', 'zdiu5', '"d15=5u"', '"d16=300u"',
                '"d24=d20*0.5-d15-p16-d16-p2-p2*0.637"',
                'ph10={0}*4{2}*4', 'ph11={0}*1{2}*1', 'ph14={0}*2{2}*2']
    if any(part not in active for part in required) or '"d20=' in active:
        raise ValueError('Saved pulse program differs from the reviewed four-plane nico profile')
    if active.count('d24') != 4:  # definition, d27 expression, and two relaxation delays
        raise ValueError('Unrecognised relaxation block in saved pulse program')
    audit = folder/'audita.txt'
    split = re.search(r'created by split 4: dataset (\d+) of 4 from EXPNO\s+(\d+)',
                      audit.read_text() if audit.exists() else '')
    if split:
        raise ValueError(f'This is split plane {split[1]}; select original experiment {split[2]}')
    td = numeric(acq, 'TD', positive=True, integer=True)
    total_rows = numeric(indirect, 'TD', positive=True, integer=True)
    if (numeric(acq, 'PARMODE', integer=True) != 1 or
            numeric(acq, 'AQ_mod', integer=True) != 3 or
            numeric(indirect, 'FnMODE', integer=True) != 5 or td % 2 or total_rows % 8):
        raise ValueError('Expected packed 2D DQD/States-TPPI acquisition with four interleaved pathways')
    if numeric(acq, 'DTYPA', integer=True) != 0:
        raise ValueError('Only the reviewed int32 Bruker raw-data format is supported')
    byte_order = numeric(acq, 'BYTORDA', integer=True)
    if byte_order not in (0, 1):
        raise ValueError('Unknown Bruker byte order')
    time = numeric(acq, 'D', 20, positive=True)
    d24 = time*.5-5e-6-numeric(acq, 'P', 16, positive=True)*1e-6-300e-6-1.637e-6*numeric(acq, 'P', 2, positive=True)
    if d24 <= 0:
        raise ValueError('D20 produces a nonpositive relaxation delay in this sequence')
    return dict(folder=folder, acq=acq, indirect=indirect, time=time, td=td,
                rows=total_rows, row_bytes=((td*4+1023)//1024)*1024,
                byte_order=byte_order)


def calculate(intensities, noises, time):
    """Magnitude estimator and first-order independent-noise propagation.

    q = sqrt(I_AB I_BA / (I_AA I_BB)); |eta| = atanh(q)/T.
    The product ratio must be positive; its sign is never repaired with abs().
    """
    if len(intensities) != 4 or len(noises) != 4 or not math.isfinite(time) or time <= 0:
        raise ValueError('Four intensities/noise SDs and positive total relaxation time are required')
    if any(not math.isfinite(v) for v in intensities) or any(not math.isfinite(s) or s <= 0 for s in noises):
        raise ValueError('Nonfinite intensity or invalid noise SD')
    flags = []
    if any(abs(intensities[i])/noises[i] < 3 for i in (0, 3)):
        flags.append('low_auto_snr')
    if any(abs(intensities[i])/noises[i] < 3 for i in (1, 2)):
        flags.append('low_cross_snr')
    result = dict(ratio=None, ratio_error=None, eta_xy_abs_s_1=None,
                  eta_xy_abs_error_s_1=None, valid=False, flags=flags)
    if any(v == 0 for v in intensities):
        flags.append('zero_intensity'); return result
    auto_product_positive = (intensities[0] > 0) == (intensities[3] > 0)
    cross_product_positive = (intensities[1] > 0) == (intensities[2] > 0)
    if auto_product_positive != cross_product_positive:
        flags.append('inconsistent_product_sign'); return result
    log_q = .5*(math.log(abs(intensities[1]))+math.log(abs(intensities[2]))-
               math.log(abs(intensities[0]))-math.log(abs(intensities[3])))
    try:
        q = math.exp(log_q)
        sq = .5*q*math.hypot(*(s/abs(v) for s,v in zip(noises,intensities)))
    except OverflowError:
        flags.append('nonfinite_ratio'); return result
    if not math.isfinite(q) or not math.isfinite(sq) or q == 0:
        flags.append('nonfinite_ratio'); return result
    result.update(ratio=q, ratio_error=sq)
    if q >= 1:
        flags.append('ratio_outside_model_domain'); return result
    eta = math.atanh(q)/time
    error = sq/(time*(1-q*q))
    if not math.isfinite(eta) or not math.isfinite(error):
        flags.append('nonfinite_rate'); return result
    if q-2*sq <= 0 or q+2*sq >= 1:
        flags.append('nonlinear_uncertainty')
    result.update(eta_xy_abs_s_1=eta, eta_xy_abs_error_s_1=error, valid=True)
    return result


def numeric_option(text):
    value = float(text)
    if not math.isfinite(value):
        raise argparse.ArgumentTypeError('Value must be finite')
    return value


def file_signature(path):
    stat = Path(path).stat()
    return stat.st_dev,stat.st_ino,stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns


def peak_list(path):
    matches = sorted(Path(p).resolve() for p in glob.glob(path) if Path(p).is_file())
    if len(matches) != 1:
        raise ValueError('Select exactly one peak list')
    path = matches[0]
    fields = rt.meaningful_lines(path)[0].lower().split()
    fmt = 'sparky3d' if 'w3' in fields else 'sparky2d'
    return path, fmt, rt.read_peaks(path, fmt)


def write_numbers(path, peaks, extracted, noises, profile, peak_name, spectrum_name, peak_notes=None):
    coordinate_fields = [f'w{i}' for i in range(1,len(peaks[0]['coordinates'])+1)]
    shared = {field: peaks[0]['coordinates'][i] for i,field in enumerate(coordinate_fields)
              if all(p['coordinates'][i] == peaks[0]['coordinates'][i] for p in peaks)}
    fields = ['assignment','residue','eta_xy_abs_s-1','eta_xy_abs_error_s-1','overlap',
              'valid','flags','ratio','ratio_error','duplicate_assignment','peak_id']
    fields += [f for f in coordinate_fields if f not in shared]
    fields += [f'plane_{i:02d}_intensity' for i in range(1,5)]
    results = []
    with Path(path).open('w',newline='') as handle:
        writer = csv.DictWriter(handle,fieldnames=fields,quoting=csv.QUOTE_NONNUMERIC)
        writer.writeheader()
        for peak in peaks:
            values = extracted[peak['peak_id']]
            fit = calculate(values,noises,profile['time'])
            if peak['overlap']: fit['flags'].append('overlap')
            if peak['duplicate_assignment']: fit['flags'].append('duplicate_assignment')
            results.append(dict(peak=peak, **fit))
            row = {key: peak[key] for key in ('assignment','residue','overlap','peak_id','duplicate_assignment')}
            row.update({key.replace('_s_1','_s-1'): value for key,value in fit.items() if key != 'flags'})
            row['flags'] = ';'.join(fit['flags'])
            row.update({field: peak['coordinates'][i] for i,field in enumerate(coordinate_fields) if field not in shared})
            row.update({f'plane_{i:02d}_intensity': v for i,v in enumerate(values,1)})
            writer.writerow(row)
        metadata = dict(source_folder=profile.get('source_folder',profile['folder'].parent.name),
                        expno=profile.get('expno',profile['folder'].name),
                        experiment_type='eta_xy', peak_list=peak_name, source_spectrum=spectrum_name,
                        pulse_program='nico_nnh_tr_noSE', evolution_time_s=profile['time'],
                        plane_roles=', '.join(ROLES), plane_noise_sd=', '.join(map(str,noises)),
                        rate_convention='magnitude', noise_propagation='first_order_independent_planes',
                        **shared)
        metadata.update(peak_notes or {})
        for key,value in metadata.items():
            text = str(value).replace('\r','\\r').replace('\n','\\n')
            handle.write(f'# {key}: {text}\n')
    return results


def plot_script(stage, prefix, results, profile, peak_name):
    q = rt.gnuplot_string
    residues = [r['peak']['residue'] for r in results]
    low, high = min(residues)-2, max(residues)+2
    title = profile.get('source_folder',profile['folder'].parent.name)+'/'+profile.get('expno',profile['folder'].name)+' | '+peak_name
    valid = [r for r in results if r['valid']]
    note = (f"T = {profile['time']:.6g} s (D20); {len(valid)}/{len(results)} valid; "
            f"{sum(r['peak']['overlap'] for r in results)} overlaps. |etaXY| = atanh(q)/T (magnitude).\n"
            'Errors: independent plane noise, first-order propagation. Weak-signal flags are in the CSV.')
    lines = ['set terminal pdfcairo size 11,8 noenhanced font "Sans,10"',
             'set output '+q(prefix+'.pdf'), 'set multiplot layout 2,1 title '+q(title),
             f'set xrange [{low}:{high}]', 'set grid ytics', 'set xlabel "Residue"',
             'set ylabel "|etaXY| (s^-1)"', 'set title '+q(note),
             'set key outside right']
    for what in ('rate','ratio'):
        if what == 'ratio':
            lines += ['set ylabel "Intensity-product ratio q"',
                      'set title "q = sqrt(I2 I3 / (I1 I4)); expected domain: 0 < q < 1. See CSV flags."',
                      'set arrow 1 from graph 0,first 1 to graph 1,first 1 nohead dt 2 lc rgb "gray"']
        chosen = valid if what == 'rate' else [r for r in results if r['ratio'] is not None]
        if not chosen:
            lines += ['set yrange [0:1]', 'plot 1/0 notitle', 'set yrange [*:*]']
            continue
        lines += ['plot \'-\' using 1:2:3 with yerrorbars pt 7 lc rgb "#23628b" title "No overlap", '
                  '\'-\' using 1:2:3 with yerrorbars pt 6 lc rgb "#b44731" title "Overlap"']
        value,error = ('eta_xy_abs_s_1','eta_xy_abs_error_s_1') if what == 'rate' else ('ratio','ratio_error')
        for overlap in (False,True):
            lines += [f"{r['peak']['residue']} {r[value]:.17g} {r[error]:.17g}" for r in chosen if r['peak']['overlap'] == overlap]
            lines.append('e')
    lines += ['unset multiplot']
    (stage/'plot.gnu').write_text('\n'.join(lines)+'\n')


def analyze(args):
    profile = experiment_profile(args.experiment)
    peak_path, fmt, peaks = peak_list(args.peaks)
    noises_path = Path(args.noise).resolve()
    noises = rt.counted_vector(noises_path, positive=True)
    if len(noises) != 4:
        raise ValueError('etaXY requires four per-plane noise SDs, not a single peak-list noise estimate')
    data_path = Path(args.extracted or args.data).resolve()
    if not args.extracted and rt.plane_count(data_path) != 4:
        raise ValueError('Expected exactly four processed planes, ordered auto_A, cross_AB, cross_BA, auto_B')
    destination = Path(args.output_dir).resolve()
    prefix = args.prefix or peak_path.name.split('.')[0]+'_etaXY'
    rt.safe_output_name(prefix)
    names = [prefix+'.pdf',prefix+'_numbers.csv',prefix+'_support.zip']
    if not destination.is_dir(): raise ValueError('Output directory must already exist')
    with rt.run_directory(destination,'eta-analyze') as stage:
        watched = [peak_path,noises_path,profile['folder']/'acqus',profile['folder']/'acqu2s',profile['folder']/'pulseprogram']
        before = {p:p.read_bytes() for p in watched}
        processing_record = data_path.parent/'processing.txt'
        if processing_record.is_file():
            before[processing_record] = processing_record.read_bytes()
            settings = dict(line.split(': ',1) for line in processing_record.read_text().splitlines() if ': ' in line)
            if settings.get('indirect_ft_negated_planes') != '2, 4':
                raise ValueError('Processing record lacks the required nitrogen-direction correction for planes 2/4; reprocess this series')
            if float(settings.get('evolution_time_s','nan')) != profile['time']:
                raise ValueError('Processing record and acquisition D20 disagree')
        data_stat = file_signature(data_path)
        clean = stage/'peaks.list'
        clean.write_text('assignment '+ ' '.join(f'w{i}' for i in range(1,len(peaks[0]['coordinates'])+1))+'\n'+
            ''.join(p['peak_id']+' '+' '.join(format(Decimal(c),'f') for c in p['coordinates'])+'\n' for p in peaks))
        raw = stage/'extracted_ids.out'
        peak_name = peak_path.name
        spectrum_name = str(Path(args.data))
        if args.extracted:
            mapping = data_path.with_name('peaks.list')
            watched.append(mapping);before[mapping]=mapping.read_bytes()
            if rt.meaningful_lines(mapping) != rt.meaningful_lines(clean):
                raise ValueError('Saved extraction mapping differs from peak-list identities/coordinates')
            shutil.copyfile(data_path,raw)
            origin = {}
            for line in raw.read_text().splitlines():
                if line.startswith('# '):
                    key,separator,value = line[2:].partition(': ')
                    if separator:
                        if key in origin:raise ValueError('Duplicate saved extraction metadata')
                        origin[key] = value
            if not {'source_folder','expno','peak_list','source_spectrum','evolution_time_s'}.issubset(origin):
                raise ValueError('Saved extraction lacks its source/timing metadata')
            if float(origin['evolution_time_s']) != profile['time']:
                raise ValueError('Saved extraction time disagrees with acquisition D20')
            if (profile['folder']/'ser').exists() and (origin['source_folder'],origin['expno']) != (
                    profile['folder'].parent.name,profile['folder'].name):
                raise ValueError('Saved extraction belongs to a different acquisition')
            profile.update(source_folder=origin['source_folder'],expno=origin['expno'])
            peak_name,spectrum_name = origin['peak_list'],origin['source_spectrum']
        else:
            with raw.open('w') as handle:
                rt.checked_run(['data_xtr','-data',data_path,'-dfmt','nmrpipe','-peak',clean,'-pfmt',fmt,'-afmt','full'],stage,handle)
            with raw.open('a') as handle:
                for key,value in dict(source_folder=profile['folder'].parent.name,expno=profile['folder'].name,
                                      peak_list=peak_name,source_spectrum=spectrum_name,
                                      evolution_time_s=profile['time']).items():
                    text=str(value).replace('\r','\\r').replace('\n','\\n')
                    handle.write(f'# {key}: {text}\n')
        extracted = rt.parse_extracted(raw.read_text(),peaks,4)
        peak_notes = {}
        for line in before[peak_path].decode().splitlines():
            if line.startswith('# '):
                key,separator,value = line[2:].partition(': ')
                if separator and key in ('coordinates_source','common_coordinate_shift_ppm','purpose'):
                    if key in peak_notes:raise ValueError('Duplicate peak-list coordinate metadata')
                    peak_notes[key] = value
        results = write_numbers(stage/names[1],peaks,extracted,noises,profile,peak_name,spectrum_name,peak_notes)
        plot_script(stage,prefix,results,profile,peak_name)
        rt.checked_run(['gnuplot','plot.gnu'],stage)
        rt.nonempty(stage/names[0])
        with zipfile.ZipFile(stage/names[2],'w',compression=zipfile.ZIP_DEFLATED) as archive:
            for path in (clean,raw,stage/'plot.gnu'):
                archive.write(path,path.name)
            archive.writestr('original_peak_list.list',before[peak_path])
            archive.writestr('vEstNoise.txt',noises_path.read_bytes())
            if processing_record in before:
                archive.writestr('processing.txt',before[processing_record])
            for filename in ('acqus','acqu2s','pulseprogram'):
                archive.writestr(filename,before[profile['folder']/filename])
            archive.writestr('README.txt','Four-plane nico_nnh_tr_noSE analysis. '+METHOD_URL+'\n'
                'Plane order: auto_A, cross_AB, cross_BA, auto_B. T is D20, in seconds.\n'
                'The square-root product determines |etaXY|, not its sign.\n'
                'Noise SDs must match the extracted intensity units. Errors assume independent Gaussian plane noise.\n'
                'First-order errors can be unreliable for weak cross peaks or ratios near 0 or 1; inspect CSV flags.\n'
                'No intensity magnitudes, clipping, plane-wise normalization or empirical uncertainty rescaling are applied.\n'
                'The CSV contains peak identities/coordinates; peaks.list maps them to extraction tags.\n'
                'Rebuild using the original list, saved extracted_ids.out, vEstNoise.txt and acquisition files.\n')
        if any(p.read_bytes() != value for p,value in before.items()) or file_signature(data_path) != data_stat:
            raise ValueError('Analysis inputs changed during the run')
        backup_reports(destination,names)
        rt.publish(stage,destination,names)
    print(f'{len(results)} peaks; {sum(r["valid"] for r in results)} valid. {destination/names[1]}')


def backup_reports(destination,names):
    if any((destination/name).is_symlink() for name in names):
        raise ValueError('Refusing to replace report symlinks')
    existing = [name for name in names if (destination/name).exists()]
    if existing:
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        backup = destination/'.report-backups'/('eta-xy-'+stamp)
        backup.mkdir(parents=True)
        for name in existing: shutil.copy2(destination/name,backup/name)


def split_raw(profile, output):
    """Deinterleave whole Bruker FID records; retain the States-TPPI pairs."""
    raw = profile['folder']/'ser'
    expected = profile['rows']*profile['row_bytes']
    if raw.stat().st_size != expected:
        raise ValueError('Raw ser size differs from the reviewed packed four-plane layout')
    paths = [output/f'plane_{i:02d}.ser' for i in range(1,5)]
    handles = []
    try:
        for path in paths:handles.append(path.open('wb'))
        with raw.open('rb') as handle:
            for index in range(profile['rows']):
                record = handle.read(profile['row_bytes'])
                if len(record) != profile['row_bytes']:raise ValueError('Raw ser changed while reading')
                handles[index%4].write(record)
            if handle.read(1):raise ValueError('Unexpected trailing raw data')
    finally:
        for handle in handles:handle.close()
    return paths


def pipeline(commands,stage):
    """Check every subprocess, including intermediate failures."""
    running = []
    try:
        for index,command in enumerate(commands):
            print('+ '+' '.join(map(str,command)),flush=True)
            previous = running[-1].stdout if running else None
            process = subprocess.Popen(list(map(str,command)),cwd=stage,stdin=previous,
                                       stdout=subprocess.PIPE if index+1 < len(commands) else subprocess.DEVNULL)
            if previous is not None:previous.close()
            running.append(process)
        for process in reversed(running):process.wait()
        for command,process in zip(commands,running):
            if process.returncode:raise subprocess.CalledProcessError(process.returncode,command)
    finally:
        for process in running:
            if process.poll() is None:process.terminate()
        for process in running:process.wait()


def process(args):
    profile = experiment_profile(args.experiment)
    destination = Path(args.output_dir).resolve() if args.output_dir else profile['folder']/'eta_xy_processed'
    if destination == profile['folder']:raise ValueError('Use a dedicated processing subdirectory')
    if args.h_zf < profile['td']//2 or args.n_zf < profile['rows']//8:
        raise ValueError('Zero-fill sizes must accommodate the acquired complex points')
    destination.mkdir(parents=True,exist_ok=True)
    with rt.run_directory(destination,'eta-process') as stage:
        initial = {p:file_signature(p) for p in [profile['folder']/'ser',profile['folder']/'acqus',
                                       profile['folder']/'acqu2s',profile['folder']/'pulseprogram']}
        raw_dir=stage/'raw';raw_dir.mkdir()
        raw_paths=split_raw(profile,raw_dir)
        planes=stage/'planes';planes.mkdir()
        acq,indirect=profile['acq'],profile['indirect']
        swap='-noswap' if (profile['byte_order']==0) == (sys.byteorder=='little') else '-swap'
        for index,raw in enumerate(raw_paths,1):
            fid=stage/f'plane_{index:02d}.fid'
            converted=['bruk2pipe','-in',raw,'-bad','0.0','-ext',swap,'-AMX',
                       '-decim',str(numeric(acq,'DECIM',positive=True)),
                       '-dspfvs',str(numeric(acq,'DSPFVS',positive=True)),
                       '-grpdly',str(numeric(acq,'GRPDLY')),
                       '-xN',str(profile['td']),'-xT',str(profile['td']//2),
                       '-yN',str(profile['rows']//4),'-yT',str(profile['rows']//8),
                       '-xMODE','DQD','-yMODE','States-TPPI']
            for axis,parameters in [('x',acq),('y',indirect)]:
                obs=numeric(parameters,'SFO1',positive=True)
                converted += ['-'+axis+'SW',str(numeric(parameters,'SW_h',positive=True)),
                              '-'+axis+'OBS',str(obs),'-'+axis+'CAR',str(numeric(parameters,'O1')/obs)]
            converted += ['-xLAB','HN','-yLAB','15N','-ndim','2','-aq2D','States','-out',fid,'-ov']
            rt.checked_run(converted,stage);rt.nonempty(fid)
            out=planes/f'plane_{index:02d}.ft2'
            # Decode the opposite 15N frequency sense in readout branches 2/4
            # before applying the same phase/window recipe to every plane.
            commands=[['nmrPipe','-in',fid,'-fn','SOL'],
                      ['nmrPipe','-fn','SP','-off','.4','-end','.98','-pow','2','-c','.5'],
                      ['nmrPipe','-fn','ZF','-size',args.h_zf],['nmrPipe','-fn','FT'],
                      ['nmrPipe','-fn','EXT','-left','-sw'],
                      ['nmrPipe','-fn','PS','-p0',args.h_p0,'-p1',args.h_p1,'-di'],
                      ['nmrPipe','-fn','TP'],
                      ['nmrPipe','-fn','SP','-off','.4','-end','.98','-pow','2','-c','.5'],
                      ['nmrPipe','-fn','ZF','-size',args.n_zf],
                      ['nmrPipe','-fn','FT','-alt']+(['-neg'] if index in (2,4) else []),
                      ['nmrPipe','-fn','PS','-p0',args.n_p0,'-p1',args.n_p1,'-di'],
                      ['nmrPipe','-fn','TP'],['nmrPipe','-fn','POLY','-auto','-ord','4','-out',out,'-ov']]
            pipeline(commands,stage);rt.nonempty(out)
        paths=[planes/f'plane_{i:02d}.ft2' for i in range(1,5)]
        headers = [rt.pipe_header(path) for path in paths]
        if any(h[9] != 2 or (h[99],h[219]) != (headers[0][99],headers[0][219]) for h in headers):
            raise ValueError('Processed planes have inconsistent spectral grids')
        rt.set_series_metadata(paths)
        with (stage/'m.ft3').open('wb') as handle:
            rt.checked_run(['xyz2pipe','-in',str(planes/'plane_%02d.ft2')],stage,handle)
        if rt.plane_count(stage/'m.ft3') != 4:raise ValueError('Processed stack does not contain four planes')
        rt.checked_run(['pipe2ucsf','m.ft3','proc.ucsf'],stage);rt.nonempty(stage/'proc.ucsf')
        noises=[]
        for path in paths:
            rt.checked_run(['sethdr',path,'-ndim','2','-zN','1','-zT','1','-zMODE','Real'],stage)
            rt.checked_run(['setfdata',path,'-fileCount','1'],stage)
            with (stage/'noise_stat.txt').open('w') as handle:
                rt.checked_run(['specStat.com','-in',path,'-stat','vEstNoise'],stage,handle)
            noises.append(rt.noise_value(stage/'noise_stat.txt'))
        (stage/'vEstNoise.txt').write_text('4 '+' '.join(map(str,noises))+'\n')
        settings = (f"source_folder: {profile['folder'].parent.name}\nexpno: {profile['folder'].name}\n"
                    f"pulse_program: nico_nnh_tr_noSE\nevolution_time_s: {profile['time']}\n"
                    f"plane_roles: {', '.join(ROLES)}\nh_phase: {args.h_p0}, {args.h_p1}\n"
                    f"n_phase: {args.n_p0}, {args.n_p1}\nzero_fill: {args.h_zf}, {args.n_zf}\n"
                    'indirect_ft_negated_planes: 2, 4\n'
                    'States-TPPI; one common phase/window/baseline recipe for all four planes.\n'
                    'SOL; SP .4/.98/2/.5; H left-half extraction; POLY -auto -ord 4.\n')
        (stage/'processing.txt').write_text(settings)
        if any(file_signature(p)!=value for p,value in initial.items()):raise ValueError('Acquisition changed during processing')
        rt.publish(stage,destination,['planes','m.ft3','proc.ucsf','vEstNoise.txt','processing.txt'])
    print(f'Processed four pathways in {destination}. Pick/export one list on planes/plane_01.ft2.')


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('inspect',help='Validate saved acquisition metadata; read-only')
    p.add_argument('--experiment',default='.')
    p=sub.add_parser('process',help='Process the packed raw acquisition as four matched planes')
    p.add_argument('--experiment',default='.');p.add_argument('--output-dir')
    for key in ('h-p0','h-p1','n-p0','n-p1'):p.add_argument('--'+key,type=numeric_option,default=0.)
    p.add_argument('--h-zf',type=int,default=4096);p.add_argument('--n-zf',type=int,default=1024)
    p=sub.add_parser('analyze',help='Extract four intensities per peak and calculate |etaXY|')
    p.add_argument('--experiment',default='.');p.add_argument('--peaks',default='*.list')
    p.add_argument('--data',default='eta_xy_processed/m.ft3')
    p.add_argument('--extracted',help='Saved data_xtr intensities instead of extracting spectra')
    p.add_argument('--noise',default='eta_xy_processed/vEstNoise.txt')
    p.add_argument('--prefix',default='');p.add_argument('--output-dir',default='.')
    args=parser.parse_args(argv)
    try:
        if args.command=='inspect':
            profile=experiment_profile(args.experiment)
            print(f"{profile['folder']}\nPulse program: nico_nnh_tr_noSE\n"
                  f"D20 total relaxation time: {profile['time']} s\n"
                  f"Four pathways: {', '.join(ROLES)}\n"
                  f"{profile['rows']//8} complex 15N points per plane; States-TPPI\n"
                  'Plan: retain one acquisition; extract one common peak list across four processed planes.')
        elif args.command=='process':process(args)
        else:analyze(args)
    except (OSError,ValueError,subprocess.SubprocessError) as error:
        print(f'etaXY: {error}',file=sys.stderr);return 1
    return 0


if __name__=='__main__':sys.exit(main())
