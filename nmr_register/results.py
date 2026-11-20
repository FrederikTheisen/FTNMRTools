"""Attach explicitly selected original numbers CSVs, verified by hash and source."""
import csv
from dataclasses import dataclass
import hashlib
import io
import json
from pathlib import Path
import re
import time
from .discovery import acquisition_key,safe_read_path
from .grist import GristError


@dataclass
class ResultSelection:
    acquisition_key: str
    path: str
    sha256: str
    content: bytes


def read_numbers_bytes(data):
    """Mirror the existing CSV-record-aware reader, including metadata footers."""
    reader = csv.reader(io.StringIO(data.decode('utf-8'),newline=''))
    try:
        fields = next(reader)
    except StopIteration:
        raise ValueError('Empty numbers CSV') from None
    rows,metadata,footer = [],{},False
    for values in reader:
        if not values:
            continue
        if values[0].startswith('# '):
            footer = True
            key,sep,value = ','.join(values)[2:].partition(': ')
            if sep:
                if key in metadata:
                    raise ValueError('Duplicate numbers metadata: '+key)
                metadata[key] = value
        else:
            if footer or len(values) != len(fields):
                raise ValueError('Invalid numbers table/footer')
            rows.append(dict(zip(fields,values)))
    if not rows:
        raise ValueError('Empty numbers table')
    return fields,rows,metadata


def _read_selected(path,sha256,key,max_bytes):
    try:
        path = safe_read_path(path)
    except (OSError,ValueError):
        raise ValueError('Selected results must be an explicit unprotected regular absolute path without symlink components') from None
    if not path.is_file() or path.stat().st_size > max_bytes:
        raise ValueError('Selected results are missing or exceed the configured size bound')
    with path.open('rb') as handle:
        content = handle.read(max_bytes+1)
    if len(content)>max_bytes:
        raise ValueError('Selected results exceed the configured size bound')
    if not re.fullmatch(r'[0-9a-f]{64}',sha256) or hashlib.sha256(content).hexdigest()!=sha256:
        raise ValueError('Selected results SHA256 mismatch')
    _,_,metadata = read_numbers_bytes(content)
    relative = key.split(':',1)[1]
    folder,_,expno = relative.rpartition('/')
    if metadata.get('source_folder') != folder or metadata.get('expno') != expno:
        raise ValueError('Selected results source identity differs from manifest acquisition')
    return ResultSelection(key,str(path),sha256,content)


def select_results(manifest,source_id,max_bytes=2_000_000,targets=None):
    manifest = Path(manifest)
    with manifest.open(newline='') as handle:
        reader = csv.DictReader(handle,delimiter='\t')
        required = {'acquisition','numbers_file','numbers_sha256'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError('Manifest requires acquisition, numbers_file and numbers_sha256')
        selections,keys,hashes = [],set(),{}
        for row in reader:
            relative = row['acquisition'].strip()
            if targets is not None and relative not in set(targets):
                continue
            parts = Path(relative).parts
            if not parts or Path(relative).is_absolute() or '211@950' in parts or any(p in ('.','..') for p in parts):
                raise ValueError('Invalid or protected manifest acquisition')
            folder,_,expno = relative.rpartition('/')
            key = acquisition_key(source_id,folder,expno)
            if key in keys:
                raise ValueError('Multiple CURRENT result selections for one acquisition')
            selection = _read_selected(row['numbers_file'],row['numbers_sha256'].strip().lower(),key,max_bytes)
            if selection.sha256 in hashes and hashes[selection.sha256] != key:
                raise ValueError('Result hash shared by different acquisitions requires reconciliation')
            hashes[selection.sha256] = key
            keys.add(key)
            selections.append(selection)
    return selections


def queue_results(state,selections):
    with state.db:
        state.db.execute('''CREATE TABLE IF NOT EXISTS result_queue
          (key TEXT PRIMARY KEY,path TEXT NOT NULL,sha256 TEXT NOT NULL)''')
        for selection in selections:
            sent = state.get_meta('result_sent:'+selection.acquisition_key)
            signature = [selection.path,selection.sha256]
            if sent != signature:
                state.db.execute('INSERT OR REPLACE INTO result_queue VALUES (?,?,?)',
                                 (selection.acquisition_key,selection.path,selection.sha256))


def _attachment_records(client):
    records = client.attachments()
    return [{'id':r.get('id'),**r.get('fields',r)} for r in records]


def _attachment_id(state,client,selection):
    filename = 'nmr-result-'+selection.sha256+'.csv'
    cached = state.db.execute('SELECT * FROM uploads WHERE hash=?',(selection.sha256,)).fetchone()
    if cached and cached['attachment']:
        # Numeric attachment IDs are scoped to one document. Recheck bytes on a
        # recovery path so stale or restored local state cannot link another file.
        content = client.download_attachment(cached['attachment'])
        if hashlib.sha256(content).hexdigest() != selection.sha256:
            raise GristError('Cached result attachment hash mismatch; reconcile before linking')
        return cached['attachment']
    # Reconcile even the very first upload by deterministic filename so recovery
    # from local state loss does not create another copy.
    candidates = [r for r in _attachment_records(client) if r.get('fileName') == filename]
    if len(candidates)>1:
        raise GristError('Multiple result attachments require reconciliation')
    if candidates:
        attachment = candidates[0]['id']
        content = client.download_attachment(attachment)
        if hashlib.sha256(content).hexdigest() != selection.sha256:
            raise GristError('Existing result attachment hash mismatch; reconcile manually')
    else:
        attempts = state.get_meta('upload_attempts:'+selection.sha256,0)
        if attempts >= 3:
            raise GristError('Result upload retry bound reached; reconcile attachments before resetting')
        state.set_meta('upload_attempts:'+selection.sha256,attempts+1)
        with state.db:
            state.db.execute('INSERT OR REPLACE INTO uploads VALUES (?,NULL,?,?,?)',
                             (selection.sha256,'uploading',filename,time.time()))
        ids = client.upload(filename,selection.content)
        if len(ids)!=1:
            raise GristError('Unexpected attachment upload response')
        attachment = ids[0]
    with state.db:
        state.db.execute('INSERT OR REPLACE INTO uploads VALUES (?,?,?,?,?)',
                         (selection.sha256,attachment,'uploaded',filename,time.time()))
    return attachment


def attach_results(state,client,selections=None,max_bytes=2_000_000,discovery_issues=(),eligible_keys=None):
    if getattr(client,'base',None):
        state.bind_document(client.base)
    if selections is not None:
        queue_results(state,selections)
    state.db.execute('''CREATE TABLE IF NOT EXISTS result_queue
        (key TEXT PRIMARY KEY,path TEXT NOT NULL,sha256 TEXT NOT NULL)''')
    queued = state.db.execute('SELECT * FROM result_queue ORDER BY key').fetchall()
    if not queued:
        state.set_meta('result_issues',[])
        return {'attached':0,'pending_results':0,'blocked_results':0,'result_issues':[]}
    # A conflicted local version must also suppress an older queued attachment.
    suppressed = {issue['key']: issue['reason'] for issue in discovery_issues}
    blocked, selected = [], []
    # Re-verify immutable input before upload, even after restart. A local file
    # problem blocks that acquisition without poisoning the remaining queue.
    for row in queued:
        reason = suppressed.get(row['key'])
        if eligible_keys is not None and row['key'] not in eligible_keys:
            reason = 'Acquisition not in the current scan scope; result selection retained'
        if reason:
            blocked.append({'key':row['key'], 'reason':reason})
            continue
        try:
            selected.append(_read_selected(row['path'],row['sha256'],row['key'],max_bytes))
        except (OSError, ValueError) as error:
            blocked.append({'key':row['key'], 'reason':str(error)})
    rows = client.records('Experiments') if selected else []
    by_key = {}
    for row in rows:
        by_key.setdefault(row['fields'].get('AcquisitionKey'),[]).append(row)
    attached = 0
    def block(selection,reason):
        blocked.append({'key':selection.acquisition_key,'reason':reason})
    for selection in selected:
        matching = by_key.get(selection.acquisition_key,[])
        if len(matching)!=1:
            block(selection,'Results require one exact registered acquisition; selection retained')
            continue
        row = matching[0]
        fields = row['fields']
        if fields.get('ReconciliationReason'):
            block(selection,'Acquisition requires reconciliation before attaching results')
            continue
        old_hash = fields.get('ResultSHA256')
        acknowledged = state.get_meta('result_sent:'+selection.acquisition_key)
        if old_hash and old_hash != selection.sha256 and (not acknowledged or acknowledged[1] != old_hash):
            block(selection,'A different result hash is already attached; reconcile explicitly before replacement')
            continue
        if old_hash != selection.sha256 or not fields.get('ResultCSV'):
            attachment = _attachment_id(state,client,selection)
            # An explicitly changed CURRENT selection replaces only the result we
            # previously acknowledged. Recheck after upload before moving pointer.
            fresh = [r for r in client.records('Experiments') if r['fields'].get('AcquisitionKey') == selection.acquisition_key]
            if len(fresh) != 1 or fresh[0]['id'] != row['id'] or fresh[0]['fields'].get('ResultSHA256') != old_hash:
                block(selection,'Result pointer changed during upload; reconcile before replacement')
                continue
            if fresh[0]['fields'].get('ReconciliationReason'):
                block(selection,'Acquisition became conflicted during upload; reconcile before linking results')
                continue
            client.update_records('Experiments',[{'id':row['id'],'fields':{
                'ResultCSV':['L',attachment],'ResultSHA256':selection.sha256,'ResultSourcePath':selection.path}}])
            attached += 1
        elif fields.get('ResultSourcePath') != selection.path:
            block(selection,'Existing result source pointer differs; reconcile explicitly')
            continue
        # Timeout after PATCH leaves queue; GET on next run recognizes committed hash.
        with state.db:
            state.db.execute('DELETE FROM result_queue WHERE key=?',(selection.acquisition_key,))
        state.set_meta('result_sent:'+selection.acquisition_key,[selection.path,selection.sha256])
    state.set_meta('result_issues',blocked)
    return {'attached':attached,'pending_results':state.db.execute('SELECT count(*) FROM result_queue').fetchone()[0],
            'blocked_results':len(blocked),'result_issues':blocked}
