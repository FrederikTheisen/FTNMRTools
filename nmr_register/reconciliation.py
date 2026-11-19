"""Explicit, conservative duplicate reconciliation; no discovery or fuzzy matching.

Snapshots map table IDs to REST ``records`` lists. Columns are normalized metadata
objects with table, colId, type, isFormula. Caller must include every table and
column (including hidden ones), download a full document backup, and exclude
other writers throughout preflight/apply. Grist has no compare-and-swap here.
"""
from copy import deepcopy
from dataclasses import dataclass
from .schema import CANONICAL_ACQUISITION_FIELDS, assessment_from_values


class ReconciliationError(ValueError):
    pass


def empty(value):
    # False and zero are data, not empty cells.
    return value is None or value == '' or value == ['L']


def semantic_type(value):
    key = str(value).lower().replace(' ', '').replace('_', '').replace('-', '')
    return {'15nhsqc': 'hsqc', '1h15nhsqc': 'hsqc', 't1rho': 'r1rho',
            't1ρ': 'r1rho', 'r1': 't1', 'r2': 't2'}.get(key, key)


def _resolved_migration_issues(value, merged):
    """Remove only complete migration-generated clauses proven obsolete here."""
    import re
    exact = {
        'Data not found in scanned roots; recorded from workbook',
        'Inherited exact workbook path; data not found in scanned roots',
        'Incomplete acquisition identity; source details retained',
        'Planned or incomplete record; identity and sample unresolved',
        'Review-only record; sample and registration unresolved',
    }
    # Match complete clauses with boundaries so genuine warning substrings survive.
    for phrase in exact:
        value = re.sub(r'(?:(?<=^)|(?<=[;\n]))[ \t]*' + re.escape(phrase) +
                       r'[ \t]*(?=;|\n|$)', '', value)
    parts = re.split(r'[;\n]', value)
    kept = []
    skip_condition = False
    for part in parts:
        clause = part.strip()
        if not clause:
            continue
        if re.fullmatch(r'Potential matching row [0-9]+ has unresolved acquisition identity', clause):
            continue
        if merged.get('Sample') and re.fullmatch(r'Sample link unresolved: .+', clause):
            skip_condition = True
            continue
        if skip_condition and re.fullmatch(r'(?:workbook condition|condition): .+', clause, re.IGNORECASE):
            skip_condition = False
            continue
        skip_condition = False
        if clause == 'Assessment reason missing' and str(merged.get('AssessmentReason') or '').strip():
            continue
        kept.append(clause)
    return '; '.join(kept)


@dataclass
class MergePlan:
    old_id: int
    target_id: int
    key: str
    before: dict
    columns: list
    actions: list
    after: dict


def _rows(snapshot, table):
    return {r['id']: r['fields'] for r in snapshot[table]}


def _ref(value, old, target, multiple):
    if multiple:
        if not isinstance(value, list) or not value or value[0] != 'L':
            if empty(value):
                return value
            raise ReconciliationError('Malformed reference list')
        return ['L'] + list(dict.fromkeys(target if x == old else x for x in value[1:]))
    return target if value == old else value


def build_plan(snapshot, columns, old_id, target_id, target_key):
    """Build a reviewable plan for ONE explicitly reviewed old/target pair.

    Every nonformula field is preserved, combined if text evidence, or conflicts.
    Observed acquisition facts and acquisition identity stay with linked target.
    No migration warning is automatically discarded.
    """
    if old_id == target_id:
        raise ReconciliationError('Pair must contain distinct records')
    before = deepcopy(snapshot)
    rows = _rows(before, 'Experiments')
    if old_id not in rows or target_id not in rows:
        raise ReconciliationError('Both records must exist when planning')
    old, target = rows[old_id], rows[target_id]
    if not target_key or target.get('AcquisitionKey') != target_key:
        raise ReconciliationError('Reviewed target key does not match')
    if str(old.get('Instrument')) == '950' and str(old.get('Expno')) == '211':
        raise ReconciliationError('Protected acquisition 211@950')
    if str(target.get('Instrument')) == '950' and str(target.get('Expno')) == '211':
        raise ReconciliationError('Protected acquisition 211@950')
    if any(rid != target_id and f.get('AcquisitionKey') == target_key
           for rid, f in rows.items()):
        raise ReconciliationError('Target acquisition key is not unique')
    after = deepcopy(before)
    surviving = _rows(after, 'Experiments')[target_id]
    actions = []
    patch = {}
    identities = {'AcquisitionKey', 'SourcePath', 'SourceID', 'Folder', 'Expno',
                  'Instrument', 'RawPresent', 'RawBytes', 'TransferStable',
                  'ScanError', 'LastChangedAt'}
    combined = {'Notes', 'SourceEvidence', 'MigrationKey', 'MigrationIssue',
                'AcquisitionDateText', 'AssessmentReason'}
    def unset(name, value, typ=''):
        return empty(value) or (name.startswith('QC') and value == 'Not assessed') or (
            typ in ('Ref:Samples', 'Ref:Series') and value == 0)
    # Preserve both the displayed decision and a fresh calculation from inputs.
    prior_decisions = set()
    for record in (old, target):
        for decision in (record.get('Assessment'), assessment_from_values(record)):
            if not empty(decision):
                prior_decisions.add(decision)
    if len(prior_decisions) > 1:
        raise ReconciliationError('Conflicting effective Assessment decisions')
    # Decisions must agree before explanations may be combined.
    for name in ('ImportedAssessment', 'AnalysisStatus'):
        a, b = old.get(name), target.get(name)
        if not empty(a) and not empty(b) and a != b:
            raise ReconciliationError('Conflicting preserved field: ' + name)
    if not empty(old.get('ResultCSV')) and not empty(target.get('ResultCSV')) and old['ResultCSV'] != target['ResultCSV']:
        raise ReconciliationError('Conflicting preserved field: ResultCSV')
    evidence = []
    provenance = {k: old[k] for k in ('AcquisitionKey', 'SourcePath', 'SourceID',
                  'Folder', 'Expno', 'Instrument') if not empty(old.get(k))}
    import json
    evidence.append('Reviewed alias provenance: legacy record ' + str(old_id) +
                    ' ' + json.dumps(provenance, sort_keys=True) +
                    ' -> linked record ' + str(target_id) + ' ' + target_key)
    for c in columns:
        if c['table'] != 'Experiments' or c.get('isFormula') or c['type'] == 'ManualSortPos':
            continue
        name = c['colId']
        a, b = old.get(name), target.get(name)
        if unset(name, a, c['type']) or a == b:
            continue
        if name in (CANONICAL_ACQUISITION_FIELDS - {'ExperimentType'}) or name in ('Title', 'TitleSource'):
            # Reviewed exact identity establishes source ownership. A missing raw
            # value is kept missing; retain workbook facts as explicit evidence.
            evidence.append('Preserved legacy ' + name + ': ' + str(a) +
                            '; linked acquisition ' + name + ': ' + str(b))
            continue
        if name.startswith('Observed') or name in identities:
            if empty(b) and name.startswith('Observed'):
                patch[name] = deepcopy(a)
            elif empty(b):
                evidence.append('Preserved legacy ' + name + ': ' + str(a))
            continue
        if unset(name, b, c['type']):
            if name == 'SourceEvidence' and name in patch:
                patch[name] += '\n' + str(a)
            else:
                patch[name] = deepcopy(a)
        elif name == 'AcquiredAt':
            evidence.append('Preserved legacy ' + name + ': ' + str(a) +
                            '; linked record ' + name + ': ' + str(b))
        elif name == 'ExperimentType' and semantic_type(a) == semantic_type(b):
            continue
        elif name in combined and isinstance(a, str) and isinstance(b, str):
            base = patch.get(name, b)
            patch[name] = base + '\n' + a if a not in base.split('\n') else base
        elif name != 'ResultCSV' and c['type'] == 'Attachments' and isinstance(a, list) and isinstance(b, list) and a[:1] == b[:1] == ['L']:
            patch[name] = ['L'] + list(dict.fromkeys(b[1:] + a[1:]))
        else:
            raise ReconciliationError('Conflicting preserved field: ' + name)
    evidence_base = patch.get('SourceEvidence', target.get('SourceEvidence') or '')
    patch['SourceEvidence'] = '\n'.join([v for v in [evidence_base] + evidence if v])
    merged_values = dict(target, **patch)
    issues = merged_values.get('MigrationIssue')
    if isinstance(issues, str):
        cleaned = _resolved_migration_issues(issues, merged_values)
        if cleaned != issues:
            patch['MigrationIssue'] = cleaned
            merged_values['MigrationIssue'] = cleaned
    if prior_decisions and assessment_from_values(merged_values) not in prior_decisions:
        raise ReconciliationError('Merge would mask or change effective Assessment')
    if patch:
        actions.append(['UpdateRecord', 'Experiments', target_id, patch])
        surviving.update(deepcopy(patch))
    for c in columns:
        typ = c['type']
        if typ not in ('Ref:Experiments', 'RefList:Experiments'):
            continue
        table, name = c['table'], c['colId']
        if table not in after:
            raise ReconciliationError('Missing reference table: ' + table)
        if c.get('isFormula'):
            # Formula references recompute; they are never manually overwritten.
            continue
        for row in after[table]:
            if table == 'Experiments' and row['id'] == old_id:
                continue
            value = row['fields'].get(name)
            new = _ref(value, old_id, target_id, typ.startswith('RefList:'))
            if new != value:
                actions.append(['UpdateRecord', table, row['id'], {name: new}])
                row['fields'][name] = new
    after['Experiments'] = [r for r in after['Experiments'] if r['id'] != old_id]
    actions.append(['RemoveRecord', 'Experiments', old_id])
    return MergePlan(old_id, target_id, target_key, before, deepcopy(columns), actions, after)


def apply_plan(client, plan, *, read_snapshot, read_columns, backup_path,
               exclusive_write=False):
    """Apply once; fresh full-state equality guards edits and uncertain retries.

    Snapshot readers must fetch fresh server state. A timeout is propagated; retry
    this same plan to recognize a committed bundle without repeating deletion.
    """
    from pathlib import Path
    if not exclusive_write:
        raise ReconciliationError('An exclusive editing window is required')
    backup = Path(backup_path)
    if not backup.is_file() or backup.stat().st_size == 0:
        raise ReconciliationError('Full document backup required')
    def comparable(snapshot):
        result = deepcopy(snapshot)
        formulas = {(c['table'], c['colId']) for c in plan.columns if c.get('isFormula')}
        for table, records in result.items():
            for record in records:
                record['fields'] = {k: v for k, v in record['fields'].items()
                                    if (table, k) not in formulas}
            records.sort(key=lambda r: r['id'])
        return result
    current = comparable(read_snapshot())
    if read_columns() != plan.columns:
        raise ReconciliationError('Schema changed after review')
    if current == comparable(plan.after):
        return 'already_applied'
    if current != comparable(plan.before):
        raise ReconciliationError('Document changed after review; rebuild plan')
    client.apply(deepcopy(plan.actions))
    if comparable(read_snapshot()) != comparable(plan.after):
        raise ReconciliationError('Unexpected post-apply state; inspect before retry')
    return 'applied'
