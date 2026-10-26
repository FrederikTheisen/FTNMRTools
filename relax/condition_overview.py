#!/usr/bin/env python3
"""Build condition TSVs and refresh pulse-program reports from saved results only.

No spectrum processing, intensity extraction, peak matching or fitting occurs.
See condition_overview.txt. Python standard library plus gnuplot for PDFs.
"""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import io
import json
import math
import os
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import sys
import zipfile

import relax_tools as rt


MODES = {'15N T1': 'r1', '15N T1rho': 'r1rho', '1H-15N hetNOE': 'noe',
         'r1': 'r1', 'r1rho': 'r1rho', 'noe': 'noe'}
SAMPLE_NAMES = {'ANP32A': 'ANP32aIDR', 'SeV Ntail': 'SeVNtail',
                'HeTDR1 (folder/title evidence)': 'HeTDR1', 'heTDR1': 'HeTDR1'}


def load_json(path):
    return json.loads(Path(path).read_text())


def read_csv(path):
    with Path(path).open(newline='') as handle:
        return list(csv.DictReader(handle))


def read_numbers_bytes(data):
    """Separate metadata after parsing CSV records, preserving multiline warnings."""
    reader = csv.reader(io.StringIO(data.decode('utf-8'), newline=''))
    fields = next(reader)
    rows, metadata, footer = [], {}, False
    for values in reader:
        if not values:
            continue
        if values[0].startswith('# '):
            footer = True
            key, sep, value = ','.join(values)[2:].partition(': ')
            if sep:
                if key in metadata:
                    raise ValueError('Duplicate numbers metadata: '+key)
                metadata[key] = value
        else:
            if footer or len(values) != len(fields):
                raise ValueError('Invalid numbers table/footer')
            rows.append(dict(zip(fields, values)))
    if not rows:
        raise ValueError('Empty numbers table')
    return fields, rows, metadata


def add_pulse_footer(data, pulse):
    _, _, metadata = read_numbers_bytes(data)
    if 'pulse_program' in metadata:
        if metadata['pulse_program'] == pulse:
            return data
        if metadata['pulse_program'] not in ('unknown', ''):
            raise ValueError('Existing numbers pulse program disagrees with acquisition')
        lines = data.splitlines(keepends=True)
        reader = csv.reader(io.StringIO(data.decode('utf-8'), newline=''))
        start = 0
        for record in reader:
            end = reader.line_num
            if record and record[0].startswith('# pulse_program:'):
                return b''.join(lines[:start])+f'# pulse_program: {pulse}\n'.encode()+b''.join(lines[end:])
            start = end
        raise ValueError('Cannot locate existing pulse footer')
    # Retain every existing byte, including original table precision and quoting.
    return data+(b'' if data.endswith(b'\n') else b'\n')+f'# pulse_program: {pulse}\n'.encode()


def numeric(value):
    try:
        value = float(value)
        return value if math.isfinite(value) else None
    except (TypeError, ValueError):
        return None


def truth(value):
    return str(value).lower() in ('true', '1', 'yes')


def fmt(value):
    value = numeric(value)
    return format(value, '.12g') if value is not None else ''


def safe_source(name):
    path = Path(name)
    if path.is_absolute() or len(path.parts) != 2 or '..' in path.parts or '211@950' in path.parts:
        raise ValueError('Invalid/protected acquisition path: '+name)
    return path


def confirmed_condition(row):
    c = json.loads(row.get('conditions') or '{}')
    if isinstance(c, str):
        raise ValueError('Condition record must be an object')
    salt = numeric(row.get('NaCl_mM'))
    temperature = numeric(row.get('temperature_C'))
    if temperature is None:
        temperature = numeric(c.get('nominal_temperature_C'))
    field = numeric(row.get('field_MHz'))
    issues = []
    if 'provisional' in row.get('temperature_basis', '').lower():
        issues.append('nominal temperature provisional: '+row['temperature_basis'])
    if salt is None: issues.append('salt unconfirmed')
    if temperature is None: issues.append('nominal temperature unconfirmed')
    if field is None: issues.append('field unconfirmed')
    if c.get('condition_status', '').startswith('proposal') or c.get('review_unresolved_metadata'):
        issues.append('condition review pending')
    # Never replace unconfirmed reviewed values with controller TE or folder tokens.
    if 'review_nominal_temperature' in c and numeric(c['review_nominal_temperature']) is None:
        issues.append('reviewed nominal temperature unresolved')
    if 'review_salt_concentration_mM' in c and numeric(c['review_salt_concentration_mM']) is None:
        issues.append('reviewed salt unresolved')
    protein = SAMPLE_NAMES.get(row['protein'], row['protein'])
    construct = row.get('construct', '')
    additives = c.get('additives') or c.get('additive_condition') or c.get('DNA_condition') or ''
    if isinstance(additives, dict): additives = json.dumps(additives, sort_keys=True)
    if not construct: issues.append('construct unconfirmed')
    if issues:
        return None, '; '.join(dict.fromkeys(issues))
    return (protein, construct, MODES[row['experiment_type']], fmt(salt), fmt(temperature), fmt(field), str(additives)), ''


def qc_records(fid_inventory=None, r1_qc=None, decisions=None):
    records, protocols, exclusions = {}, {}, {}
    if fid_inventory:
        for row in read_csv(fid_inventory):
            name = row['relative_path']; safe_source(name)
            status = row.get('FID_screen_status', '')
            label = 'screened_no_major_flags' if status == 'no_screen_flag' else 'review'
            records[name] = (label, '; '.join(filter(None, [status, row.get('notes', '')])))
            protocols[name] = row
    if r1_qc:
        for row in read_csv(r1_qc):
            name = row['relative_path']; safe_source(name)
            if row.get('error'):
                records[name] = ('review', row['error']); continue
            summary = json.loads(row.get('event_summary') or '{}')
            flags = summary.get('flag_counts_native_processed_range', {})
            count = summary.get('flagged_trace_count_native_processed_range')
            if count is None:
                records[name] = ('review', 'screened; retained-range flags unavailable')
            else:
                records[name] = ('review' if count else 'screened_no_major_flags',
                                 'retained-range trace flags: '+(str(flags) if flags else 'none'))
    if decisions:
        data = load_json(decisions)
        # Explicit reviewed findings, not a new threshold applied to raw metrics.
        for name in data.get('known_strong_growth_experiments', []):
            safe_source(name)
            exclusions[name] = 'reviewed strong late-FID growth: '+Path(decisions).name
        for row in data.get('exclusions', []):
            if row.get('exclude', True):
                safe_source(row['relative_path']); exclusions[row['relative_path']] = row['reason']
    return records, protocols, exclusions


def discover(catalog, project, data_root, fid_inventory=None, r1_qc=None, decisions=None):
    project, data_root = Path(project).resolve(), Path(data_root).resolve()
    qcs, protocols, excluded = qc_records(fid_inventory, r1_qc, decisions)
    result, seen, directories = [], set(), set()
    for source in read_csv(catalog):
        if source.get('status') != 'analysed' or source.get('experiment_type') not in MODES:
            continue
        relative = str(safe_source(source['relative_path']))
        if relative in seen: raise ValueError('Duplicate current acquisition: '+relative)
        seen.add(relative)
        destination = Path(source['output_directory']).resolve()
        try:
            relative_directory = destination.relative_to(project)
        except ValueError as error:
            raise ValueError('Report directory is outside project: '+str(destination)) from error
        if any(part.startswith('.') or part in ('pilots', 'workers', 'previous_driver_versions')
               for part in relative_directory.parts):
            raise ValueError('Catalog points to archived/diagnostic output: '+relative)
        if destination in directories: raise ValueError('Duplicate current report directory')
        directories.add(destination)
        numbers = sorted(destination.glob('*_numbers.csv'))
        if len(numbers) != 1: raise ValueError('Expected one numbers file in '+str(destination))
        fields, rows, metadata = read_numbers_bytes(numbers[0].read_bytes())
        mode = MODES[source['experiment_type']]
        custom_noe = mode == 'noe' and 'hetNOE' in fields
        if metadata.get('source_folder') and (metadata['source_folder'], metadata.get('expno')) != tuple(Path(relative).parts):
            raise ValueError('Numbers metadata/source acquisition mismatch: '+relative)
        condition, pending = confirmed_condition(source)
        investigator = json.loads(source.get('investigator_exclusion') or '{}')
        exclusion = excluded.get(relative, '')
        if investigator.get('exclude_from_primary_analysis'):
            exclusion = investigator.get('assessment') or investigator.get('observation') or 'investigator exclusion'
        pulse = rt.acquisition_pulse_program(data_root/relative)
        recorded = metadata.get('pulse_program', 'unknown')
        if recorded not in ('unknown', '', pulse) and pulse != 'unknown':
            raise ValueError('Pulse-program conflict: '+relative)
        pulse = pulse if pulse != 'unknown' else recorded
        job = load_json(destination/'job.json') if (destination/'job.json').is_file() else {}
        if job and job.get('spec', {}).get('relative_path', relative) != relative:
            raise ValueError('Job/source acquisition mismatch: '+relative)
        row = dict(source, relative_path=relative, directory=destination, numbers=numbers[0], rows=rows,
                   metadata=metadata, mode=mode, custom_noe=custom_noe, condition=condition, pending=pending,
                   exclusion=exclusion, pulse=pulse, qc=qcs.get(relative, ('not_assessed', 'no saved acquisition-level FID screen')),
                   protocol=protocols.get(relative, {}), job=job, acqus_path=data_root/relative/'acqus')
        result.append(row)
    return sorted(result, key=lambda r: (r.get('acquisition_date', ''), r['relative_path']))


def read_user_reviews(path, records, project):
    """Read editable input; never rewrite it during overview generation."""
    if not Path(path).is_file():
        return
    source_targets = {row['relative_path'] for row in records}
    condition_targets = {str(Path(*row['directory'].relative_to(project).parts[:3])) for row in records}
    reviews = {}
    for review in read_tsv(path):
        target = review['target'].strip()
        if target not in source_targets | condition_targets:
            raise ValueError('Unknown review target: '+target)
        if target in reviews:
            raise ValueError('Duplicate review target: '+target)
        assessment = review['assessment'].strip().lower()
        if assessment not in ('', 'usable', 'review', 'exclude'):
            raise ValueError('Review assessment must be blank, usable, review or exclude: '+target)
        reviews[target] = (assessment, review['comment'].strip())
    for row in records:
        condition_target = str(Path(*row['directory'].relative_to(project).parts[:3]))
        applicable = [reviews[target] for target in (condition_target, row['relative_path']) if target in reviews]
        assessment = next((a for a,c in reversed(applicable) if a), '')
        comments = [c for a,c in applicable if c]
        row['user_assessment'] = assessment
        row['user_comment'] = ' | '.join(comments)
        if assessment == 'exclude':
            row['exclusion'] = 'user excluded'+(': '+row['user_comment'] if row['user_comment'] else '')
        elif assessment == 'usable':
            row['exclusion'] = ''  # Explicit user acceptance overrides acquisition exclusions, not condition metadata.
        if assessment or comments:
            automatic = row['qc'][0]+'; '+row['qc'][1]
            row['qc'] = ('user_'+assessment if assessment else row['qc'][0],
                         ('user: '+row['user_comment']+'; ' if comments else '')+'automated: '+automatic)


def read_tsv(path):
    with Path(path).open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames != ['target', 'assessment', 'comment']:
            raise ValueError('Review TSV columns must be target, assessment, comment')
        rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError('Incomplete review TSV row')
    return rows


def initialise_reviews(path, records):
    path = Path(path)
    if path.exists():
        return
    # Exclusive creation: even a file created since inventory is never overwritten.
    try:
        with path.open('x', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(['target', 'assessment', 'comment'])
            for row in records:
                writer.writerow([row['relative_path'], '', ''])
    except FileExistsError:
        return


def native_members(archive):
    """Bytes of saved native fitter records; reporting must never alter these."""
    return {n: archive.read(n) for n in archive.namelist()
            if Path(n).name.startswith('rlx_dfit') or Path(n).name in ('data_xtr.out', 'vEstNoise.txt', 'relax_header.txt')}


def render_noe(stage, row, pdf_name):
    q = rt.gnuplot_string
    title = row['relative_path']+' | Pulse sequence: '+row['pulse']
    commands = ['set terminal pdfcairo size 11,5 noenhanced font "Sans,10"',
                'set output '+q(pdf_name), 'set title '+q(title), 'set xlabel "Residue"',
                'set ylabel "Signed Isat / Iref"', 'set grid ytics', 'set key outside right',
                'set label 1 "Bars: saved independent-noise SD; placement/overlap/processing systematics excluded." at screen .07,.025 font ",8"',
                'set bmargin 4', 'set arrow from graph 0,first 0 to graph 1,first 0 nohead lc rgb "gray"',
                'plot "-" using 1:2:3 with yerrorbars pt 7 lc rgb "#285879" title "Eligible annotation screen", "-" using 1:2:3 with yerrorbars pt 6 lc rgb "#c17b19" title "Flagged"']
    for eligible in (True, False):
        for peak in row['rows']:
            if truth(peak.get('comparison_eligible')) == eligible:
                value, error = numeric(peak['hetNOE']), numeric(peak['hetNOE_sd'])
                if value is not None and error is not None:
                    commands.append(f"{int(peak['residue'])} {peak['hetNOE']} {peak['hetNOE_sd']}")
        commands.append('e')
    (stage/'plot.gnu').write_text('\n'.join(commands)+'\n')
    rt.checked_run(['gnuplot', 'plot.gnu'], stage)
    rt.nonempty(stage/pdf_name)


def refresh(row):
    """One locked, backed-up transaction per acquisition; rerender saved results."""
    destination = row['directory']; prefix = row['numbers'].name[:-len('_numbers.csv')]
    support = destination/(prefix+'_support.zip')
    pulse = row['pulse'] or 'unknown'
    job_path = destination/'job.json'
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    with rt.run_directory(destination, 'overview-refresh') as stage:
        watched = {path: path.read_bytes() for path in (row['numbers'], support)}
        for path in (job_path, destination/'provenance.json'):
            if path.is_file(): watched[path] = path.read_bytes()
        acq_path = row.get('acqus_path')
        if acq_path and acq_path.is_file():
            watched[acq_path] = acq_path.read_bytes()
            if rt.pulse_program_from_text(watched[acq_path].decode(), acq_path) != pulse:
                raise ValueError('Acquisition pulse program changed since inventory')
        _, row['rows'], row['metadata'] = read_numbers_bytes(watched[row['numbers']])
        current_job = json.loads(watched[job_path]) if job_path in watched else {}
        if row['metadata'].get('pulse_program') == pulse and current_job.get('report_refresh', {}).get('pulse_program') == pulse:
            return 'already_current'
        record = dict(kind='pulse-program metadata refresh from saved results', pulse_program=pulse,
                      refreshed_utc=datetime.now(timezone.utc).isoformat(), refitted=False,
                      numerical_table_unchanged=True, original_fit_provenance_preserved=True)
        with zipfile.ZipFile(io.BytesIO(watched[support])) as original:
            old_native = native_members(original)
            replacements = {}
            if row['custom_noe']:
                pdfs = sorted(destination.glob('*_intensity.pdf'))
                if len(pdfs) != 1: raise ValueError('Cannot identify saved NOE PDF')
                pdf_name = pdfs[0].name
                render_noe(stage, row, pdf_name)
                provenance_path = destination/'provenance.json'
                provenance = json.loads(watched[provenance_path])
                provenance.update(pulse_program=pulse, report_refresh=record)
                replacements['provenance.json'] = (json.dumps(provenance, indent=2, ensure_ascii=False)+'\n').encode()
                (stage/'provenance.json').write_bytes(replacements['provenance.json'])
            else:
                csv_names = [n for n in original.namelist() if n.endswith('_analysis.csv') and '/' not in n]
                if len(csv_names) != 1: raise ValueError('Missing portable saved analysis: '+row['relative_path'])
                saved_csv = stage/'saved_analysis.csv'; saved_csv.write_bytes(original.read(csv_names[0]))
                bundle = rt.load_analysis_csv(saved_csv); context = bundle['context']
                if (context['source_collection'], context['source_experiment']) != tuple(Path(row['relative_path']).parts):
                    raise ValueError('Portable analysis/source mismatch')
                if context.get('pulse_program', 'unknown') not in ('unknown', '', pulse):
                    raise ValueError('Portable pulse-program conflict')
                if context['experiment_type'] != row['mode'] or context['output_prefix'] != prefix:
                    raise ValueError('Portable observable/prefix mismatch')
                context['pulse_program'] = pulse
                pdf_name = rt.safe_output_name(context['pdf_name'])
                generated = rt.prepare_reports(stage, bundle['peaks'], bundle['extracted'], bundle['fit_text'],
                    bundle['times'], bundle['noises'], row['mode'], prefix, pdf_name, row['relative_path'],
                    context=context, warning_map=bundle.get('warning_map'))
                # New code must reproduce the saved numbers before any publication.
                old_table = read_numbers_bytes(watched[row['numbers']])[:2]
                new_table = read_numbers_bytes((stage/row['numbers'].name).read_bytes())[:2]
                if old_table != new_table:
                    raise ValueError('Report rebuild would change saved numerical cells: '+row['relative_path'])
                replacements = {name: (stage/name).read_bytes() for name in generated
                                if name not in (row['numbers'].name, pdf_name)}
                rt.checked_run(['gnuplot', 'plot.gnu'], stage); rt.nonempty(stage/pdf_name)
            watched[destination/pdf_name] = (destination/pdf_name).read_bytes()
            (stage/row['numbers'].name).write_bytes(add_pulse_footer(watched[row['numbers']], pulse))
            replacements['report_refresh.json'] = (json.dumps(record, indent=2)+'\n').encode()
            with zipfile.ZipFile(stage/support.name, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
                for info in original.infolist():
                    archive.writestr(info, replacements.pop(info.filename, original.read(info.filename)))
                for name, data in replacements.items(): archive.writestr(name, data)
        with zipfile.ZipFile(stage/support.name) as archive:
            if native_members(archive) != old_native:
                raise ValueError('Native fitter/extraction bytes changed')
        names = [row['numbers'].name, pdf_name, support.name]
        if row['custom_noe']: names.append('provenance.json')
        if job_path in watched:
            job = json.loads(watched[job_path]); job['report_refresh'] = record
            # Existing frozen-driver resume checks require CURRENT output checksums.
            # Update only fields already present; expose no new hashes in reports.
            import hashlib
            for path in list(job.get('output_hashes', {})):
                staged = stage/Path(path).name
                if staged.name in names:
                    job['output_hashes'][path] = hashlib.sha256(staged.read_bytes()).hexdigest()
            (stage/'job.json').write_text(json.dumps(job, indent=2, ensure_ascii=False)+'\n'); names.append('job.json')
        if any(path.read_bytes() != data for path, data in watched.items()):
            raise ValueError('Saved report input changed during rendering')
        backup = destination/'.report-backups'/('condition-overview-'+stamp)
        backup.mkdir(parents=True)
        for name in names:
            if (destination/name).is_file(): shutil.copy2(destination/name, backup/name)
        rt.publish(stage, destination, names)
    _, row['rows'], row['metadata'] = read_numbers_bytes(row['numbers'].read_bytes())
    return 'refreshed'


def row_flags(peak):
    flags = [s for field in ('flags', 'review_flags') for s in peak.get(field, '').split(';') if s]
    if truth(peak.get('uncertainty_rescaled')): flags.append('uncertainty_rescaled')
    if peak.get('warnings'): flags.append('fitter_warning')
    if numeric(peak.get('p_holm')) is not None and float(peak['p_holm']) < .05: flags.append('lack_of_fit')
    if truth(peak.get('duplicate_assignment')): flags.append('duplicate_assignment')
    return sorted({'error_adjusted_by_fitter' if flag == 'uncertainty_rescaled' else
                   'lack_of_fit' if flag == 'original_noise_lack_of_fit' else flag for flag in flags})


def peak_cells(peaks, mode):
    width = 4 if mode == 'noe' else 5
    if not peaks: return ['']*(width-1)+['missing']
    if len(peaks) > 1 or any(truth(p.get('duplicate_assignment')) for p in peaks): return ['']*(width-1)+['multiple_peaks']
    peak = peaks[0]; flags = row_flags(peak)
    value = peak.get('hetNOE', peak.get('noe_ratio')) if mode == 'noe' else peak.get('rate_s-1')
    error = peak.get('hetNOE_sd', peak.get('noe_ratio_error')) if mode == 'noe' else peak.get('rate_error_s-1')
    valid = numeric(value) is not None and numeric(error) is not None and float(error) >= 0
    if 'valid_fit' in peak: valid = valid and truth(peak['valid_fit'])
    if not valid: flags.append('invalid_fit')
    overlap = any(truth(peak.get(key)) for key in ('overlap', 'ref_overlap', 'NOE_overlap'))
    if overlap: flags.append('overlap')
    cells = [value if valid else '', error if valid else '', str(overlap)]
    if mode != 'noe': cells.append(peak.get('reduced_chi2_input', ''))
    return cells+[';'.join(sorted(set(flags)))]


def first_height(peak):
    return numeric(peak.get('plane_01_intensity', peak.get('I_ref')))


def relative_noise_header(row):
    counts = Counter(int(p['residue']) for p in row['rows'])
    heights = []
    for peak in row['rows']:
        height = first_height(peak)
        overlap = any(truth(peak.get(k)) for k in ('overlap', 'ref_overlap', 'NOE_overlap'))
        if (height is not None and height != 0 and not overlap and counts[int(peak['residue'])] == 1 and
                not truth(peak.get('duplicate_assignment')) and truth(peak.get('valid_fit', True))):
            heights.append(abs(height))
    if not heights:
        return 'unavailable: no usable first-plane peak heights'
    reference = statistics.median(heights)
    raw = row['metadata'].get('plane_noise_sd', '')
    if raw:
        noises = [numeric(value.strip()) for value in raw.split(',')]
        if any(v is None or v < 0 for v in noises):
            raise ValueError('Invalid saved noise vector')
        percentages = ', '.join(format(100*v/reference, '.3g')+'%' for v in noises)
    elif row['mode'] == 'noe':
        parts = []
        for label, key in [('ref', 'I_ref_sd'), ('sat', 'I_NOE_sd')]:
            noises = sorted({numeric(p.get(key)) for p in row['rows'] if numeric(p.get(key)) is not None})
            if not noises:
                parts.append(label+'=unknown')
            else:
                value = format(100*noises[0]/reference, '.3g')+'%'
                if len(noises) > 1: value += '..'+format(100*noises[-1]/reference, '.3g')+'%'
                parts.append(label+'='+value)
        percentages = ', '.join(parts)
    else:
        return 'unavailable: no saved noise SD'
    return percentages+'; median |first-plane height|='+format(reference, '.6g')+f' (n={len(heights)})'


def header_cells(row, project):
    meta, peaks, protocol = row['metadata'], row['rows'], row['protocol']
    source = row['relative_path']
    pulse = row['pulse']
    if row['mode'] == 'r1rho':
        pulse += '; spinlock_Hz='+str(protocol.get('spinlock_nominal_Hz') or 'unknown')
        if protocol.get('spinlock_offset_from_carrier_Hz'):
            pulse += '; offset_Hz='+protocol['spinlock_offset_from_carrier_Hz']
    basis = meta.get('fit_basis', 'unknown')
    if row['mode'] == 'noe':
        selection = 'not applicable (ratio)'
    elif basis == 'all_points':
        selection = 'all delays used'
    elif basis == 'exclude_shortest':
        delay = numeric(meta.get('excluded_shortest_delay_s'))
        label = format(delay*1000, '.6g')+' ms' if delay is not None else 'shortest'
        selection = label+' delay omitted'
    else:
        selection = 'fit selection unknown'
    noise = relative_noise_header(row)
    qc = row['qc'][0]+'; '+row['qc'][1]
    flagged = sum(bool(row_flags(p)) or any(truth(p.get(k)) for k in ('overlap','ref_overlap','NOE_overlap')) for p in peaks)
    valid = sum('invalid_fit' not in peak_cells([p],row['mode'])[-1] for p in peaks)
    fit = f'{valid}/{len(peaks)} valid peaks; {flagged} flagged'
    if row['mode'] == 'noe':
        fit += '; ratio analysis; exponential goodness of fit not applicable'
    else:
        fit += '; '+meta.get('dataset_conclusion', 'assessment unavailable')
        fit += '; Holm='+meta.get('holm_rejection_count', '?')+'/'+meta.get('primary_peak_count', '?')
    return [source, pulse, selection, noise, qc, fit]


def write_overview(path, group, excluded, project, review_path=None):
    condition = group[0]['condition']; protein, construct, mode, salt, temp, field, additives = condition
    group = sorted(group, key=lambda r: (r.get('acquisition_date', ''),r['relative_path']))
    per = 4 if mode == 'noe' else 5
    maps = []
    for acquisition in group:
        residues = defaultdict(list)
        for row in acquisition['rows']:
            residues[int(row['residue'])].append(row)
        maps.append(residues)
    with Path(path).open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['# generated', 'Automatically generated by condition_overview.py. DO NOT EDIT; use the persistent review TSV for assessments.'])
        writer.writerow(['# condition',f'{protein}; construct={construct}; {mode}; NaCl={salt}mM; nominal_T={temp}C; field={field}MHz'+(' ; '+additives if additives else '')])
        writer.writerow(['# units', 'ratio' if mode == 'noe' else 's^-1'])
        headers = [header_cells(row, project) for row in group]
        for index,label in enumerate(['source', 'pulse / protocol', 'delays used', 'noise / median first-plane height', 'acquisition QC', 'fit assessment']):
            writer.writerow([label]+[cell for info in headers for cell in [info[index]]+['']*(per-1)])
        parameter = {'r1': 'R1', 'r1rho': 'R1rho', 'noe': 'hetNOE'}[mode]
        labels = [parameter, 'error', 'overlap']+([] if mode == 'noe' else ['reduced_chi2'])+['flags']
        writer.writerow(['residue']+labels*len(group))
        observed = set().union(*(set(m) for m in maps))
        for residue in range(min(observed), max(observed)+1):
            writer.writerow([residue]+[cell for m in maps for cell in peak_cells(m.get(residue, []), mode)])
        writer.writerow(['# flags','missing = no saved peak; multiple_peaks = ambiguous residue; blank value/error for invalid fits; review/not_assessed QC is not a clean-data certificate'])
        writer.writerow(['# error', 'Reported fitter SD. error_adjusted_by_fitter means the error was adjusted using residual scatter; the fitted value is unchanged. This does not fix a poor fit.'])
        writer.writerow(['# noise', 'Percentages = 100 x plane noise SD / median absolute first-plane height, using unique valid non-overlapping peaks; first acquisition plane even if omitted from fitting.'])
        for row in group:
            writer.writerow(['# numbers_file', row['relative_path'], os.path.relpath(row['numbers'], project), 'date='+str(row.get('acquisition_date') or 'unknown')])
        if review_path:
            writer.writerow(['# user_reviews', os.path.relpath(review_path, project), 'Editable input; never overwritten by this generator.'])
        for row in excluded:
            writer.writerow(['# excluded',row['relative_path'],row['exclusion']])


def destination_for(group, project):
    # The directory hierarchy is used only to PLACE already-reviewed groups.
    # Conditions themselves never come from these historical directory names.
    parents = set()
    for row in group:
        parts = row['directory'].relative_to(project).parts
        if len(parts) != 6: raise ValueError('Expected protein/temperature/condition/run/collection/exp output layout')
        parents.add(Path(project).joinpath(*parts[:3]))
    if len(parents) != 1:
        raise ValueError('Reviewed condition spans different output directories; explicit placement needed')
    return next(iter(parents))/'overview.tsv'


def run(args):
    project = Path(args.project_root).resolve()
    records = discover(args.catalog,project,args.data_root,args.fid_inventory,args.r1_qc,args.qc_decisions)
    review_path = Path(getattr(args, 'user_reviews', None) or project/'condition-reviews.tsv').resolve()
    read_user_reviews(review_path, records, project)
    groups, excluded = defaultdict(list), defaultdict(list)
    for row in records:
        if row['exclusion']:
            if row['condition']: excluded[row['condition']].append(row)
        elif row['condition']:
            groups[row['condition']].append(row)
    destinations = {key:destination_for(group,project) for key,group in groups.items()}
    if len(set(destinations.values())) != len(destinations):
        raise ValueError('Distinct reviewed conditions would overwrite one overview')
    counts = dict(acquisitions=len(records), included=sum(map(len,groups.values())),
                  excluded=sum(bool(r['exclusion']) for r in records),
                  condition_review=sum(not r['condition'] and not r['exclusion'] for r in records), overviews=len(groups))
    print(json.dumps(counts,sort_keys=True),flush=True)
    if args.dry_run:
        for row in records:
            print(row['relative_path'], 'excluded: '+row['exclusion'] if row['exclusion'] else 'condition_review: '+row['pending'] if row['pending'] else 'included',flush=True)
        return counts
    with rt.run_directory(project,'overviews') as stage:
        initialise_reviews(review_path, records)
        review_before = review_path.read_bytes()
        for row in records:
            status = refresh(row) if args.refresh_reports else (
                'already_current' if row['metadata'].get('pulse_program') == row['pulse'] and
                row['job'].get('report_refresh', {}).get('pulse_program') == row['pulse'] else 'not_requested')
            row['refresh_status'] = status
            print(row['relative_path']+': '+status,flush=True)
        for index,(key,group) in enumerate(sorted(groups.items()),1):
            path = destinations[key]; name = f'overview_{index}.tsv'
            write_overview(stage/name, group, excluded[key], project, review_path)
            with rt.run_directory(path.parent,'overview-publish') as local:
                shutil.copyfile(stage/name,local/'overview.tsv')
                if path.is_file():
                    backup = path.parent/'.report-backups'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
                    backup.mkdir(parents=True);shutil.copy2(path,backup/'overview.tsv')
                rt.publish(local,path.parent,['overview.tsv'])
        index_name = 'condition-overviews.tsv'
        with (stage/index_name).open('w',newline='') as handle:
            fields=['acquisition','protein','observable','disposition','reason','acquisition_QC','pulse_program','overview','numbers_file','report_refresh','user_assessment','user_comment']
            writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
            for row in records:
                writer.writerow(dict(acquisition=row['relative_path'],protein=row['protein'],observable=row['mode'],
                    disposition='excluded' if row['exclusion'] else 'condition_review' if row['pending'] else 'included',
                    reason=row['exclusion'] or row['pending'],acquisition_QC=row['qc'][0],pulse_program=row['pulse'],
                    overview=os.path.relpath(destinations[row['condition']],project) if row['condition'] in destinations else '',
                    numbers_file=os.path.relpath(row['numbers'],project),report_refresh=row['refresh_status'],user_assessment=row.get('user_assessment',''),user_comment=row.get('user_comment','')))
        if (project/index_name).is_file():
            backup=project/'.report-backups'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
            backup.mkdir(parents=True);shutil.copy2(project/index_name,backup/index_name)
        if review_path.read_bytes() != review_before:
            raise ValueError('User review input changed while generating; rerun to incorporate changes')
        rt.publish(stage,project,[index_name])
    print('Delivery index: '+str(project/'condition-overviews.tsv'),flush=True)
    return counts


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',required=True)
    parser.add_argument('--data-root',required=True)
    parser.add_argument('--catalog',required=True,help='Reviewed master acquisition CSV; current analysed reports only')
    parser.add_argument('--fid-inventory',help='Reviewed R1rho FID inventory CSV')
    parser.add_argument('--r1-qc',help='Saved retained-range R1 FID screening CSV')
    parser.add_argument('--qc-decisions',help='Reviewed findings JSON with known_strong_growth_experiments / explicit exclusions')
    parser.add_argument('--user-reviews',help='Editable target/assessment/comment TSV; defaults to PROJECT/condition-reviews.tsv; never overwritten')
    parser.add_argument('--refresh-reports',action='store_true')
    parser.add_argument('--dry-run',action='store_true',help='Read-only inventory and disposition summary')
    args=parser.parse_args(argv)
    try:
        run(args)
    except (OSError,ValueError,KeyError,zipfile.BadZipFile,subprocess.SubprocessError) as error:
        print('Condition overview: '+str(error),file=sys.stderr)
        return 1
    return 0


if __name__=='__main__':sys.exit(main())
