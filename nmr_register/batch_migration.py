"""Reviewable, insert-only acquisition batch migration.

Caller supplies an approved manifest hash and a fresh full observation inventory.
A single writer must be used: Grist does not provide an atomic uniqueness constraint.
"""
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

from .split import is_derived_split
from .schema import CONDITIONS
from .discovery import safe_read_path
from .sync import _canonical_payload


class BatchUncertainError(RuntimeError):
    """A POST failed and the full outcome is not established by a fresh read."""
    def __init__(self, message, committed, remaining):
        super().__init__(message)
        self.committed = committed
        self.remaining = remaining


def _identity(fields, observed=False):
    prefix = 'Observed' if observed else ''
    source, folder, expno = (fields.get(prefix + key) for key in ('SourceID','Folder','Expno'))
    if (not isinstance(source,str) or not source or any(c.isspace() or c in ':/' + chr(92) for c in source) or not isinstance(folder,str) or expno is None):
        return None
    expno=str(expno)
    if (not expno.isascii() or not expno.isdigit() or folder.startswith('/') or chr(92) in folder or (folder and any(not part or part in ('.','..') or '211@950' in part for part in folder.split('/')))):
        return None
    return source,folder,expno


def _partial_matches(fields, identity):
    source,folder,expno=identity
    actual=[fields.get(k) for k in ('SourceID','Folder','Expno')]
    # A row with a populated folder or experiment number can denote this
    # acquisition even if the other components are still missing.
    if _identity(fields) or not any(v not in (None,'') for v in actual[1:]):return False
    return all(v in (None,'') or str(v)==wanted for v,wanted in zip(actual,(source,folder,expno)))


def _key(identity):
    source,folder,expno=identity
    return source+':'+(folder+'/' if folder else '')+expno


def _checked_bytes(path):
    # Map only the user's home-directory alias; never resolve input descendants.
    value=Path(path)
    if not value.is_absolute() or any('211@950' in part or part in ('.','..') for part in value.parts):
        raise ValueError('Path is relative, traversing, or protected')
    home=Path.home()
    try: relative=value.relative_to(home)
    except ValueError: pass
    else: value=home.resolve()/relative
    return safe_read_path(value).read_bytes()


def _manifest(path, expected_sha256):
    data=_checked_bytes(path); actual=hashlib.sha256(data).hexdigest()
    if not expected_sha256 or actual != expected_sha256:
        raise ValueError('Manifest hash does not match the reviewed batch')
    try: manifest=json.loads(data)
    except (ValueError,UnicodeError): raise ValueError('Invalid batch manifest') from None
    rows=manifest.get('records') if isinstance(manifest,dict) else None
    if not isinstance(rows,list) or not rows:raise ValueError('Batch manifest is empty or invalid')
    return actual,rows,manifest


def _provenance(manifest, path):
    path=str(Path(path));expected=manifest.get('source_sha256',{}).get(path)
    data=_checked_bytes(path)
    if not expected or hashlib.sha256(data).hexdigest()!=expected:raise ValueError('Source reconciliation hash changed')
    try:source=json.loads(data)
    except (ValueError,UnicodeError):raise ValueError('Invalid source reconciliation') from None
    files=source.get('source_files');entries=source.get('planned_targets')
    if not isinstance(files,list) or not files or not isinstance(entries,list):raise ValueError('Incomplete source reconciliation')
    for item in files:
        try:actual=hashlib.sha256(_checked_bytes(item['path'])).hexdigest()
        except (KeyError,OSError,TypeError):raise ValueError('Original source input unavailable') from None
        if actual!=item.get('sha256'):raise ValueError('Original source input hash changed')
    by_key=defaultdict(list);identities=set()
    for entry in entries:
        if isinstance(entry,dict):
            fields=entry.get('fields') or {};by_key[fields.get('MigrationKey')].append(entry)
            for identity in (_identity(fields),_identity(fields,True)):
                if identity:identities.add(identity)
    return by_key,identities


def _legacy(row,by_key):
    matches=by_key.get(row.get('MigrationKey'),[])
    if len(matches)!=1:raise ValueError('Legacy source coordinate is missing or ambiguous')
    fields=matches[0]['fields'];context=row.get('legacy_source_context') or {}
    if (_identity(fields)!=_identity(row) or fields.get('SourceEvidence','')!=(row.get('SourceEvidence') or '')
        or any((fields.get(k) or '')!=(context.get(k) or '') for k in ('Notes','MigrationIssue'))):
        raise ValueError('Legacy source evidence differs from reviewed manifest')
    if fields.get('Assessment') or 'Conflict:' in (fields.get('MigrationIssue') or '') or 'duplicate workbook target' in (fields.get('MigrationIssue') or ''):
        raise ValueError('Legacy decision or conflict requires separate review')
    return fields


def _payload(row, obs, legacy, registration_timestamp):
    # The source adapter validates the canonical allowlist and keeps observed
    # facts out of the destination's retired columns.
    if any(k in obs for k in ('ResultCSV','ResultSHA256','ResultSourcePath')):
        raise ValueError('Observation may not set result attachments')
    machine,facts=_canonical_payload({k:v for k,v in obs.items() if not k.startswith('_')})
    numeric=('AcquiredAt','Field_MHz','Temperature_C','D1_s','SpinLock_Hz','SpinLockOffset_Hz','RawBytes')
    values={name:None for name in numeric}
    values.update({'Snapshot'+name:None for name,typ in CONDITIONS.items() if typ in ('Numeric','Int')})
    values.update(machine);values.update(facts)
    # Use an explicit recorded workbook fact only where the collector has no
    # value. Status is excluded: observed Unknown is an actual current state.
    for name in ('Instrument','ExperimentType','AcquiredAt','Field_MHz','Temperature_C','PulseProgram','D1_s','DelayList_s','SpinLock_Hz','SpinLockOffset_Hz'):
        if name in facts:continue
        value=legacy.get(name)
        if value not in (None,''):values[name]=value
    values['RegisteredAt']=registration_timestamp
    for name in ('SourceID','Folder','Expno'):values[name]=row[name]
    context=row.get('legacy_source_context') or {}
    if not isinstance(context,dict):raise ValueError('Invalid legacy source context')
    values['SourceEvidence']=row.get('SourceEvidence') or ''
    values['MigrationKey']=row.get('MigrationKey') or ''
    values['Notes']=context.get('Notes') or ''
    values['MigrationIssue']=context.get('MigrationIssue') or ''
    for name in (*numeric,*('Snapshot'+name for name,typ in CONDITIONS.items() if typ in ('Numeric','Int'))):
        value=values.get(name)
        if value is not None and (type(value) not in (int,float) or not math.isfinite(value)):
            raise ValueError('Non-finite or nonnumeric acquisition value: '+name)
    return values


def plan_batch(manifest_path, observations, client, *, expected_manifest_sha256, source_reconciliation_path, registration_timestamp):
    """Read current Grist records and return an exact, guarded insert plan.

    observations must be a fresh full inventory, not only selected rows: this
    detects fingerprints shared with acquisitions outside the proposed batch.
    No writes are made. Previously inserted identical-identity rows are skipped.
    """
    if type(registration_timestamp) not in (int,float) or not math.isfinite(registration_timestamp) or registration_timestamp<=0:raise ValueError('A fixed registration timestamp is required')
    digest,rows,manifest=_manifest(manifest_path,expected_manifest_sha256)
    by_source,legacy_identities=_provenance(manifest,source_reconciliation_path)
    if not isinstance(observations,list):raise ValueError('A fresh full observation inventory is required')
    by_key=defaultdict(list);fingerprints=Counter()
    for obs in observations:
        if not isinstance(obs,dict):raise ValueError('Invalid observation')
        if is_derived_split(obs):continue
        k=obs.get('AcquisitionKey')
        if k:by_key[k].append(obs)
        if obs.get('ObservedFingerprint'):fingerprints[obs['ObservedFingerprint']]+=1
    existing=client.records('Experiments');keys=defaultdict(list);tuples=defaultdict(list);migration=defaultdict(list);remote_fp=defaultdict(list)
    for record in existing:
        f=record['fields'];k=f.get('AcquisitionKey')
        if k:keys[k].append(record)
        for identity in set(filter(None,(_identity(f),_identity(f,True)))):tuples[identity].append(record)
        if f.get('MigrationKey'):migration[f['MigrationKey']].append(record)
        if f.get('ObservedFingerprint'):remote_fp[f['ObservedFingerprint']].append(record)
    seen=set(); planned=[];skipped=[]
    for row in rows:
        if not isinstance(row,dict) or row.get('exclusions') != []:raise ValueError('Manifest contains unresolved candidate')
        ident=_identity(row);key=row.get('AcquisitionKey');finger=row.get('ObservedFingerprint');mkey=row.get('MigrationKey')
        if not ident or key!=_key(ident) or not isinstance(finger,str) or len(finger)!=64 or not mkey:
            raise ValueError('Manifest has invalid exact identity or provenance')
        origin=row.get('origin','legacy')
        if origin=='legacy':legacy=_legacy(row,by_source)
        elif origin=='inventory':
            if (mkey!='discovery:'+key or row.get('SourceEvidence')!='Observed acquisition: '+key
                or row.get('legacy_source_context') not in (None,{})):
                raise ValueError('Inventory provenance differs from deterministic source identity')
            if ident in legacy_identities:raise ValueError('Inventory target has legacy evidence requiring review')
            legacy={}
        else:raise ValueError('Unknown batch source origin')
        if key in seen:raise ValueError('Manifest contains duplicate acquisition key')
        seen.add(key)
        match=by_key.get(key,[])
        if len(match)!=1:raise ValueError('Fresh observation missing or ambiguous: '+key)
        obs=match[0]
        if _identity(obs,True)!=ident or obs.get('ObservedFingerprint')!=finger or fingerprints[finger]!=1:
            raise ValueError('Fresh acquisition identity/fingerprint requires review: '+key)
        if obs.get('ScanError') or obs.get('ReconciliationReason') or obs.get('ObservedSampleTagIssue') or obs.get('ObservedSampleID'):
            raise ValueError('Fresh acquisition needs reconciliation: '+key)
        for name in ('SourcePath','ObservedPulseProgram','ObservedType','ObservedAcquisitionStatus'):
            if obs.get(name)!=row.get(name):raise ValueError('Fresh observation differs from reviewed manifest: '+key)
        if any(_partial_matches(r['fields'],ident) for r in existing):raise ValueError('Partial manual identity requires review: '+key)
        conflicts={r['id']:r for r in keys.get(key,[])+tuples.get(ident,[])+migration.get(mkey,[])+remote_fp.get(finger,[])}
        if conflicts:
            if len(conflicts)!=1:raise ValueError('Registered identity is ambiguous: '+key)
            f=next(iter(conflicts.values()))['fields']
            if f.get('AcquisitionKey')!=key or _identity(f,True)!=ident or _identity(f)!=ident or f.get('ObservedFingerprint')!=finger or f.get('MigrationKey')!=mkey or f.get('ReconciliationReason'):
                raise ValueError('Registered identity or provenance conflicts: '+key)
            skipped.append(key);continue
        planned.append({'key':key,'fields':_payload(row,obs,legacy,registration_timestamp)})
    return {'manifest_sha256':digest,'insert':planned,'skipped':skipped}


def apply_batch(manifest_path, observations, client, *, expected_manifest_sha256, source_reconciliation_path, registration_timestamp, batch_size=50, dry_run=True):
    """Preview or insert; after an uncertain POST, re-read before a later retry.

    Re-run with the same reviewed manifest to resume. Existing matched records
    are never updated, including records whose human fields have since changed.
    """
    if not isinstance(batch_size,int) or isinstance(batch_size,bool) or not 1<=batch_size<=100:raise ValueError('Invalid batch size')
    plan=plan_batch(manifest_path,observations,client,expected_manifest_sha256=expected_manifest_sha256,source_reconciliation_path=source_reconciliation_path,registration_timestamp=registration_timestamp)
    if dry_run:return plan
    committed=[]
    for start in range(0,len(plan['insert']),batch_size):
        # A fresh preflight catches edits since the previous batch. A single
        # writer remains necessary between this read and Grist's POST.
        current=plan_batch(manifest_path,observations,client,expected_manifest_sha256=expected_manifest_sha256,source_reconciliation_path=source_reconciliation_path,registration_timestamp=registration_timestamp)
        pending={r['key']:r for r in current['insert']}
        batch=[pending[r['key']] for r in plan['insert'][start:start+batch_size] if r['key'] in pending]
        if not batch:continue
        try:client.add_records('Experiments',[r['fields'] for r in batch])
        except Exception as error:
            try:
                checked=plan_batch(manifest_path,observations,client,expected_manifest_sha256=expected_manifest_sha256,source_reconciliation_path=source_reconciliation_path,registration_timestamp=registration_timestamp)
                remaining={r['key'] for r in checked['insert']}
            except Exception:
                raise BatchUncertainError('POST outcome could not be reconciled; stop and inspect before retry',committed, [r['key'] for r in batch]) from error
            done=[r['key'] for r in batch if r['key'] not in remaining]
            committed.extend(done)
            unresolved=[r['key'] for r in batch if r['key'] in remaining]
            if unresolved:raise BatchUncertainError('POST outcome partially reconciled; re-run only after review',committed,unresolved) from error
            continue
        committed.extend(r['key'] for r in batch)
    return {'manifest_sha256':plan['manifest_sha256'],'inserted':committed,'skipped':plan['skipped']}
