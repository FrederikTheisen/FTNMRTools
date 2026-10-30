#!/usr/bin/env python3
"""Build condition TSVs and refresh pulse-program reports from saved results only.

No spectrum processing, intensity extraction, peak matching or fitting occurs.
See condition_overview.txt. Python standard library plus gnuplot for PDFs.
"""
import argparse
import ast
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import io
import json
import math
import os
from pathlib import Path
try:
    from zoneinfo import ZoneInfo
except ImportError:  # Python 3.8 can still inventory records; dates remain unknown.
    ZoneInfo = None
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import zipfile

import relax_tools as rt


MODES = {'15N T1': 'r1', '15N T1rho': 'r1rho', '1H-15N hetNOE': 'noe',
         'r1': 'r1', 'r1rho': 'r1rho', 'noe': 'noe'}
REVIEW_EXTRA_MODES = {'etaXY': 'eta_xy', 'water 1H T1': 'water_t1sat',
                      '15N CPMG': 'cpmg', '15N T2': 'r2',
                      'unresolved relaxation candidate': 'unresolved'}
REVIEW_MODES = dict(MODES, **REVIEW_EXTRA_MODES)
REVIEW_LABELS = {'r1': 'R1', 'r1rho': 'R1rho', 'noe': 'hetNOE',
                 'eta_xy': 'etaXY', 'water_t1sat': 'water T1sat',
                 'cpmg': '15N CPMG', 'r2': 'R2', 'unresolved': 'unresolved relaxation'}
SAMPLE_NAMES = {'ANP32A': 'ANP32aIDR', 'SeV Ntail': 'SeVNtail',
                'HeTDR1 (folder/title evidence)': 'HeTDR1', 'heTDR1': 'HeTDR1'}
REVIEW_FIELDS = ['target', 'protein', 'experiment_type', 'pulse_program',
                 'acquisition_date', 'NaCl_mM', 'nominal_temperature_C', 'field',
                 'machine_disposition', 'machine_assessment', 'machine_reason',
                 'assessment', 'comment']

RETIRED_REVIEW_FIELDS = {'field_MHz', 'delays_used', 'fit_assessment'}

REVIEW_LOOKUP = [
    ('edit', 'Edit assessment and comment. Machine columns are saved evidence, not your decisions; ordinary regeneration leaves this file unchanged.'),
    ('assessment', 'Blank keeps saved decisions; usable accepts an acquisition; review requests inspection; exclude removes it from condition overviews.'),
    ('included', 'Condition is assigned and no reviewed exclusion applies. Inclusion does not mean QC passed.'),
    ('condition_review', 'Condition assignment is unresolved or provisional; not included in a confirmed-condition overview.'),
    ('processed_only_review', 'Retained processed data have no original raw acquisition available; provenance needs review. This is not an additional confirmed raw acquisition.'),
    ('excluded', 'A reviewed exclusion applies. A blank machine_reason for an unassessed acquisition does not cancel that decision.'),
    ('not_assessed', 'No saved acquisition-level QC screen is available. machine_reason is intentionally blank.'),
    ('screened_no_major_flags', 'The saved acquisition-level screen found no major flags; this is not proof of clean intensities.'),
    ('screened_no_retained_flags', 'No remaining flags in the saved acquisition screen after omitting the initial-FID check.'),
    ('omitted check', 'Initial-FID startup screening is omitted from this assessment. Original screening records are preserved.'),
    ('screen scope', 'Automatic review of stored time-domain FIDs. Flags describe signal features; their cause and effect on amide intensities are not established.'),
    ('trace counts', 'A trace is one stored direct-dimension FID at a relaxation delay, 15N evolution increment and echo/antiecho channel. Counts are flagged traces, not residues or independent repeats.'),
    ('later FID discontinuity', 'An unusually large adjacent-sample change is flagged after 0.5 ms of the screened direct-time window; cause is unassigned.'),
    ('water-region frequency outlier', 'The strongest Fourier feature near water has an unusual frequency across 15N increments at the same delay and echo/antiecho channel. This is a feature outlier, not a confirmed physical frequency jump.'),
    ('water-region amplitude outlier', 'The water-region signal amplitude is unusual relative to other 15N increments at the same delay and echo/antiecho channel.'),
    ('trace review', 'Isolated or modest trace features need review; they are not treated as confirmed corruption.'),
    ('late-FID growth', 'Reviewed strong water-dominated growth late in the FID; an existing exclusion applies. Any stated milliseconds identify affected sampled delays.'),
    ('reviewed exclusion', 'An explicit saved exclusion applies; consult the reviewed acquisition records for its original reason.'),
    ('water growth', 'Flat or growing water FIDs were noted at some delays; the effect on short-delay intensities is unresolved.'),
    ('ADC overflow', 'ADC overflow was logged; its buffer number has not been assigned to a delay plane.'),
    ('water increment', 'A shared water-amplitude increment is a candidate observation; cause and physical time are unresolved.'),
    ('condition pending', 'Salt, nominal temperature, field, construct or another reviewed condition remains unresolved/provisional. Unknown values are not inferred from controller settings.'),
    ('QC review', 'Other saved screening notes require review; consult the original FID-screen records for the detailed evidence.'),
    ('field', 'Field is rounded to the nearest MHz for display only; condition grouping retains the reviewed numerical value.'),
]



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
    return (protein, construct, REVIEW_MODES[row['experiment_type']], fmt(salt), fmt(temperature), fmt(field), str(additives)), ''


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


def discover_reviews(catalog, project, data_root, fid_inventory=None, r1_qc=None,
                     decisions=None, analysed_records=None):
    """Inventory every primary relaxation record, independent of analysis status."""
    project, data_root = Path(project).resolve(), Path(data_root).resolve()
    qcs, _, excluded = qc_records(fid_inventory, r1_qc, decisions)
    result = list(analysed_records or [])
    existing = {row['relative_path'] for row in result}
    seen = set()
    sources = read_csv(catalog)
    dispositions = Path(catalog).with_name('source_record_dispositions.csv')
    if dispositions.is_file():
        known = {row['relative_path'] for row in sources}
        sources += [row for row in read_csv(dispositions) if row['relative_path'] not in known]
    for source in sources:
        if source.get('experiment_type') not in REVIEW_MODES:
            continue
        kind = source.get('primary_record_status', '')
        processed_only = kind.startswith('processed only')
        legacy_analysed = not kind and source.get('status') == 'analysed'
        if ((kind != 'primary raw acquisition' and not processed_only and not legacy_analysed) or
                source.get('related_primary_path') or
                source.get('status') in ('excluded_duplicate_raw', 'excluded_nonprimary')):
            continue  # Split components, setup records and copies are not extra acquisitions.
        relative = str(safe_source(source['relative_path']))
        if relative in seen:
            raise ValueError('Duplicate assessment acquisition: '+relative)
        seen.add(relative)
        if relative in existing:
            continue  # Reuse report records so user reviews reach the numerical overviews.
        row = dict(source)
        conditions = json.loads(source.get('conditions') or '{}')
        # The inventory-temperature column can represent controller grouping.
        # Prefer the saved nominal annotation, retaining provisional conditions.
        nominal = conditions.get('review_nominal_temperature', conditions.get('nominal_temperature_C'))
        if nominal is not None:
            row['temperature_C'] = nominal or 'unknown'
        condition, pending = confirmed_condition(row)
        pulse = rt.acquisition_pulse_program(data_root/relative)
        recorded = source.get('pulse_program') or 'unknown'
        if pulse != 'unknown' and recorded not in ('unknown', pulse):
            raise ValueError('Pulse-program conflict: '+relative)
        exclusion = excluded.get(relative, '')
        investigator = json.loads(source.get('investigator_exclusion') or '{}')
        if investigator.get('exclude_from_primary_analysis'):
            exclusion = investigator.get('assessment') or investigator.get('observation') or 'investigator exclusion'
        directory = Path(source['output_directory']).resolve() if source.get('output_directory') else None
        if directory is not None:
            directory.relative_to(project)
        mode = REVIEW_MODES[source['experiment_type']]
        row.update(relative_path=relative, mode=mode, directory=directory, condition=condition,
                   pending=pending, processed_only=processed_only, exclusion=exclusion,
                   pulse=pulse if pulse != 'unknown' else recorded,
                   qc=qcs.get(relative, ('not_assessed', 'no saved acquisition-level FID screen')))
        result.append(row)
    result.extend(discover_new_relaxation(data_root, sources))
    return sorted(result, key=lambda r: (r.get('acquisition_date', ''), r['relative_path']))


def discover_new_relaxation(data_root, sources):
    """Find acquired relaxation datasets added since the reviewed catalog."""
    known = {row['relative_path'] for row in sources}
    by_collection = {}
    for row in sources:
        by_collection.setdefault(Path(row['relative_path']).parts[0], row)
    result = []
    for collection in sorted(Path(data_root).iterdir()):
        if ('211@950' in collection.name or collection.name.startswith('.') or
                collection.name in ('inventory', 'scripts') or collection.is_symlink() or not collection.is_dir()):
            continue
        for folder in sorted(collection.iterdir()):
            if ('211@950' in folder.name or not folder.name.isdigit() or folder.is_symlink() or not folder.is_dir()):
                continue
            relative = str(folder.relative_to(data_root))
            if relative in known or not (folder/'acqus').is_file():
                continue
            pulse = rt.acquisition_pulse_program(folder)
            lower = pulse.lower()
            mode = ('noe' if 'hetnoe' in lower else 'eta_xy' if 'nnh_tr' in lower
                    else 'water_t1sat' if 't1sat' in lower else 'cpmg' if 'cpmg' in lower or 'rexet' in lower
                    else 'r1rho' if 't1rho' in lower or 'hsqctret' in lower
                    else 'r2' if 't2' in lower else 'r1' if 't1' in lower else None)
            if mode is None or not any((folder/name).is_file() and (folder/name).stat().st_size > 0 for name in ('ser', 'fid')):
                continue
            audit = (folder/'audita.txt').read_text(errors='replace') if (folder/'audita.txt').is_file() else ''
            if re.search(r'created by split\s+\d+:.*?from EXPNO\s+\d+', audit, re.I | re.S):
                continue
            text = (folder/'acqus').read_text(errors='replace')
            def parameter(key):
                match = re.search(r'^##\$'+key+r'=\s*([^\n]*)', text, re.M)
                return match[1].strip() if match else ''
            date = ''
            timestamp = numeric(parameter('DATE_START'))
            if timestamp is not None and timestamp > 0 and ZoneInfo is not None:
                date = datetime.fromtimestamp(timestamp, ZoneInfo('Europe/Paris')).date().isoformat()
            sibling = by_collection.get(collection.name, {})
            protein = sibling.get('protein') or next((name for token, name in
                [('anp32', 'ANP32A'), ('sev', 'SeV Ntail'), ('tdr1', 'heTDR1'), ('trd1', 'heTDR1')]
                if token in collection.name.lower()), 'unknown')
            result.append(dict(relative_path=relative, protein=protein, construct=sibling.get('construct', ''),
                               mode=mode, pulse=pulse, acquisition_date=date, field_MHz=parameter('BF1'),
                               NaCl_mM='unknown', temperature_C='unknown', directory=None, condition=None,
                               pending='condition metadata not yet reviewed', processed_only=False,
                               exclusion='', qc=('not_assessed', 'no saved acquisition-level FID screen')))
    return result


def review_condition_target(row, project):
    directory = row.get('directory')
    return str(Path(*directory.relative_to(project).parts[:3])) if directory is not None else None


def read_user_reviews(path, records, project):
    """Read editable input; never rewrite it during overview generation."""
    for row in records:
        row.setdefault('machine_review_context', machine_review_context(row))
    if not Path(path).is_file():
        return
    source_targets = {row['relative_path'] for row in records}
    condition_targets = {target for row in records if (target := review_condition_target(row, project)) is not None}
    reviews = {}
    for review in read_tsv(path):
        target = review['target'].strip()
        if target not in source_targets | condition_targets:
            raise ValueError('Unknown review target: '+target)
        if target in reviews:
            raise ValueError('Duplicate review target: '+target)
        assessment = review['assessment'].strip()
        reviews[target] = (assessment, review['comment'].strip())
    for row in records:
        condition_target = review_condition_target(row, project)
        applicable = [reviews[target] for target in (condition_target, row['relative_path']) if target in reviews]
        assessment = next((a for a,c in reversed(applicable) if a), '')
        comments = [c for a,c in applicable if c]
        row['user_assessment'] = assessment
        row['user_comment'] = ' | '.join(comments)
        if assessment.lower() == 'exclude':
            row['exclusion'] = 'user excluded'+(': '+row['user_comment'] if row['user_comment'] else '')
        elif assessment.lower() == 'usable':
            row['exclusion'] = ''  # Explicit user acceptance overrides acquisition exclusions, not condition metadata.
        if assessment or comments:
            automatic = row['qc'][0]+'; '+row['qc'][1]
            row['qc'] = ('user_'+assessment if assessment else row['qc'][0],
                         ('user: '+row['user_comment']+'; ' if comments else '')+'automated: '+automatic)


def read_review_table(path):
    with Path(path).open(newline='') as handle:
        records = csv.reader(handle, delimiter='\t')
        # Read the preamble as TSV records, so a quoted multiline user comment
        # beginning with "#" inside the table is never mistaken for a lookup row.
        for fields in records:
            if any(value.strip() for value in fields) and not fields[0].lstrip('\ufeff').startswith('# '):
                break
        else:
            fields = []
        if fields:
            fields[0] = fields[0].lstrip('\ufeff')
        if (not fields or len(set(fields)) != len(fields) or
                not {'target', 'assessment', 'comment'}.issubset(fields)):
            raise ValueError('Review TSV must contain target, assessment, comment; context columns are optional')
        rows = []
        for values in records:
            if not any(value.strip() for value in values):
                continue
            if len(values) != len(fields):
                raise ValueError('Incomplete review TSV row')
            rows.append(dict(zip(fields, values)))
    return fields, rows


def read_tsv(path):
    return read_review_table(path)[1]


def review_machine_assessment(row):
    """Omit startup-only warnings from the sheet without changing saved QC."""
    status, note = row['qc']
    if status != 'review' or row['exclusion']:
        return status
    if note.split(';', 1)[0].strip() == 'startup_pattern_only':
        return 'screened_no_retained_flags'
    prefix = 'retained-range trace flags: '
    if note.startswith(prefix):
        try:
            flags = ast.literal_eval(note[len(prefix):])
            if (isinstance(flags, dict) and flags.get('direct_startup_discontinuity_review', 0) > 0 and
                    all(isinstance(n, (int, float)) and n >= 0 for n in flags.values()) and
                    not any(n > 0 for flag, n in flags.items() if flag != 'direct_startup_discontinuity_review')):
                return 'screened_no_retained_flags'
        except (ValueError, SyntaxError, TypeError):
            pass  # Unrecognised notes retain their existing review status.
    return status


def brief_machine_reason(row):
    """Display recurring screen findings using the lookup terms, not long notes."""
    if row['qc'][0] == 'not_assessed':
        return ''
    reasons = []
    note = row['qc'][1]
    if row['exclusion']:
        if 'late-FID growth' in row['exclusion']:
            growth = re.search(r'Strong water-dominated late-FID growth at ([\d, ]+) ms', note)
            reasons.append('late-FID growth'+(' ('+growth[1]+' ms)' if growth else ''))
        else:
            reasons.append('reviewed exclusion')
    if row['pending']:
        reasons.append('condition pending')
    for flag, label in [('isolated_direct_discontinuity_review', 'later FID discontinuity'),
                        ('water_frequency_change_review', 'water-region frequency outlier'),
                        ('isolated_water_amplitude_change_review', 'water-region amplitude outlier')]:
        count = re.search(r"['\"]"+flag+r"['\"]\s*:\s*(\d+)", note)
        if count:
            reasons.append(label+' ('+count[1]+(' trace)' if int(count[1]) == 1 else ' traces)'))
    for fragment, label in [('isolated_or_modest_signature_review', 'trace review'),
                             ('Median late/early >1', 'water growth'),
                             ('Persistent near-flat/slightly growing water FID', 'water growth'),
                             ('ADC overflow', 'ADC overflow'),
                             ('Shared water-amplitude increment', 'water increment')]:
        if fragment in note:
            reasons.append(label)
    if not reasons and review_machine_assessment(row) == 'review':
        reasons.append('QC review')
    return '; '.join(dict.fromkeys(reasons))


def machine_review_context(row):
    """Snapshot of automated evidence, kept separate from user decisions."""
    disposition = ('excluded' if row['exclusion'] else 'processed_only_review' if row.get('processed_only')
                   else 'condition_review' if row['pending'] else 'included')
    return dict(target=row['relative_path'], protein=row['protein'],
                experiment_type=REVIEW_LABELS[row['mode']],
                pulse_program=row['pulse'], acquisition_date=row.get('acquisition_date') or 'unknown',
                NaCl_mM=row.get('NaCl_mM') or 'unknown',
                nominal_temperature_C=row.get('temperature_C') or 'unknown',
                field=(f"{numeric(row.get('field_MHz')):.0f} MHz" if numeric(row.get('field_MHz')) is not None else 'unknown'),
                machine_disposition=disposition,
                machine_assessment=review_machine_assessment(row), machine_reason=brief_machine_reason(row))


def append_missing_reviews(path, records):
    """Preserve the entire edited document byte-for-byte; append new targets."""
    path = Path(path)
    before = path.read_bytes()
    fields, existing = read_review_table(path)
    targets = [row['target'].strip() for row in existing]
    contexts = []
    seen = set(targets)
    for row in records:
        target = row['relative_path']
        if target not in seen:
            context = row.get('machine_review_context') or machine_review_context(row)
            contexts.append(dict(context, assessment='', comment=''))
            seen.add(target)
    if not contexts:
        return
    newline = '\r\n' if b'\r\n' in before else '\n'
    added = io.StringIO(newline='')
    writer = csv.DictWriter(added, fieldnames=fields, delimiter='\t', lineterminator=newline, extrasaction='ignore')
    writer.writerows(contexts)
    after = before + (b'' if before.endswith((b'\n', b'\r')) else newline.encode()) + added.getvalue().encode('utf-8')
    with tempfile.TemporaryDirectory(prefix='.review-append-', dir=path.parent) as folder:
        stage = Path(folder)
        (stage/path.name).write_bytes(after)
        if path.read_bytes() != before:
            raise ValueError('User review input changed during append; rerun to include latest edits')
        backup = path.parent/'.report-backups'/('review-append-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
        backup.mkdir(parents=True)
        (backup/path.name).write_bytes(before)
        rt.publish(stage, path.parent, [path.name])


def initialise_reviews(path, records, update_context=False, append_missing=False):
    path = Path(path)
    if path.exists() and not update_context:
        if append_missing:
            append_missing_reviews(path, records)
        return
    before = path.read_bytes() if path.exists() else None
    existing = read_tsv(path) if before is not None else []
    contexts = {row['relative_path']: row.get('machine_review_context') or machine_review_context(row)
                for row in records}
    fields = list(REVIEW_FIELDS)
    for row in existing:
        fields += [key for key in row if key not in fields and key not in RETIRED_REVIEW_FIELDS]
    output, seen = [], set()
    for row in existing:
        target = row['target'].strip()
        if target in seen:
            raise ValueError('Duplicate review target: '+target)
        seen.add(target)
        # Preserve user assessment/comment and any custom columns verbatim.
        enriched = dict(row)
        old_field = enriched.pop('field_MHz', '')  # Previous generated display column.
        enriched.pop('delays_used', None)
        enriched.pop('fit_assessment', None)
        if old_field and not enriched.get('field'):
            value = numeric(old_field)
            enriched['field'] = f'{value:.0f} MHz' if value is not None else old_field
        if target in contexts:
            enriched.update(contexts[target])
        output.append(enriched)
    for row in records:
        target = row['relative_path']
        if target not in seen:
            output.append(dict(contexts[target], assessment='', comment=''))
    def write_sheet(handle):
        lookup = csv.writer(handle, delimiter='\t', lineterminator='\n')
        lookup.writerows([['# '+term, explanation] for term, explanation in REVIEW_LOOKUP])
        lookup.writerow([])
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(output)
    if before is None:
        # Exclusive creation never clobbers a concurrently created user input.
        try:
            with path.open('x', newline='') as handle:
                write_sheet(handle)
        except FileExistsError:
            return
    else:
        # Only an explicit context update rewrites an existing review sheet.
        with tempfile.TemporaryDirectory(prefix='.review-context-', dir=path.parent) as folder:
            stage = Path(folder)
            with (stage/path.name).open('w', newline='') as handle:
                write_sheet(handle)
            if path.read_bytes() != before:
                raise ValueError('User review input changed during context update')
            backup = path.parent/'.report-backups'/('review-context-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
            backup.mkdir(parents=True)
            shutil.copy2(path, backup/path.name)
            rt.publish(stage, path.parent, [path.name])


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


def overview_filename(condition):
    """Name a table using reviewed conditions, independent of its directory."""
    protein, construct, mode, salt, temperature, field, additives = condition
    def token(value):
        text = re.sub(r'[^A-Za-z0-9_.+-]+', '-', str(value)).strip('.')
        if not text or text in ('.', '..'):
            raise ValueError('Condition cannot be represented in an overview filename')
        return text
    sample = token(protein)
    if construct and not re.sub(r'\W', '', protein).lower().endswith(re.sub(r'\W', '', construct).lower()):
        sample += '_'+token(construct)
    observable = {'r1': 'R1', 'r1rho': 'R1rho', 'noe': 'hetNOE'}[mode]
    parts = [sample, observable, token(salt)+'mMNaCl', token(temperature)+'C', token(field)+'MHz']
    if additives:
        parts.append(token(additives))
    name = '_'.join(parts)+'_overview.tsv'
    if len(name.encode('utf-8')) > 255:
        raise ValueError('Condition overview filename exceeds filesystem limit')
    return name


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
    conditions = {row['condition'] for row in group}
    if len(conditions) != 1 or None in conditions:
        raise ValueError('Overview filename requires one confirmed condition')
    return next(iter(parents))/overview_filename(next(iter(conditions)))


def run(args):
    project = Path(args.project_root).resolve()
    reviews_only = getattr(args, 'reviews_only', False)
    if reviews_only and args.refresh_reports:
        raise ValueError('--reviews-only cannot refresh numerical reports')
    records = [] if reviews_only else discover(args.catalog,project,args.data_root,args.fid_inventory,args.r1_qc,args.qc_decisions)
    review_records = discover_reviews(args.catalog, project, args.data_root, args.fid_inventory,
                                      args.r1_qc, args.qc_decisions, None if reviews_only else records)
    grist_snapshot = getattr(args, 'grist_review_snapshot', None)
    if grist_snapshot and (getattr(args, 'user_reviews', None) or getattr(args, 'update_review_context', False)):
        raise ValueError('--grist-review-snapshot cannot be combined with --user-reviews or --update-review-context')
    review_path = Path(grist_snapshot or getattr(args, 'user_reviews', None) or project/'condition-reviews.tsv').resolve()
    if grist_snapshot and not review_path.is_file():
        raise ValueError('Grist review snapshot does not exist: '+str(review_path))
    if reviews_only:
        counts = dict(review_records=len(review_records), **dict(Counter(
            machine_review_context(row)['machine_disposition'] for row in review_records)))
        print(json.dumps(counts, sort_keys=True), flush=True)
        if not args.dry_run and not grist_snapshot:
            with rt.run_directory(project, 'reviews'):
                initialise_reviews(review_path, review_records,
                                   update_context=getattr(args, 'update_review_context', False), append_missing=True)
            print('Assessment file: '+str(review_path), flush=True)
        return counts
    read_user_reviews(review_path, review_records, project)
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
        if not grist_snapshot:
            initialise_reviews(review_path, review_records, update_context=getattr(args, 'update_review_context', False), append_missing=True)
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
                shutil.copyfile(stage/name,local/path.name)
                if path.is_file():
                    backup = path.parent/'.report-backups'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
                    backup.mkdir(parents=True);shutil.copy2(path,backup/path.name)
                rt.publish(local,path.parent,[path.name])
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
    parser.add_argument('--catalog',required=True,help='Reviewed acquisition CSV; all primary relaxation records enter the assessment regardless of analysis status')
    parser.add_argument('--fid-inventory',help='Reviewed R1rho FID inventory CSV')
    parser.add_argument('--r1-qc',help='Saved retained-range R1 FID screening CSV')
    parser.add_argument('--qc-decisions',help='Reviewed findings JSON with known_strong_growth_experiments / explicit exclusions')
    parser.add_argument('--grist-review-snapshot',help='Read-only Grist review snapshot; disables legacy review writes')
    parser.add_argument('--user-reviews',help='Editable review TSV; defaults to PROJECT/condition-reviews.tsv; normal overview regeneration leaves it unchanged')
    parser.add_argument('--update-review-context',action='store_true',help='Explicitly refresh informational review columns, preserving all user assessments/comments and backing up the sheet')
    parser.add_argument('--reviews-only',action='store_true',help='Append missing relaxation acquisitions to the assessment, preserving the edited document without requiring analysis results')
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
