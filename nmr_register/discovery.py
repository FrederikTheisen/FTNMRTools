"""Directory traversal with protected paths pruned before any inspection."""
import os
import math
import stat
from pathlib import Path
from .bruker import observe
from .split import is_derived_split


def safe_read_path(path):
    """Reject symlinks in every component before inspecting any child of them."""
    path = Path(path)
    if not path.is_absolute() or any(part in ('.','..') or '211@950' in part for part in path.parts):
        raise ValueError('Path is relative, traversing, or protected')
    current = Path(path.anchor)
    for part in path.parts[1:]:
        current = current/part
        try:
            info = current.lstat()
        except FileNotFoundError:
            # No descendant can exist without its parent. Caller handles missing
            # roots/files; importantly no symlink was followed to get this far.
            break
        if stat.S_ISLNK(info.st_mode):
            raise ValueError('Path contains a symlink; configure its canonical physical path')
    return path


def acquisition_key(source_id, folder, expno):
    if not source_id or ':' in source_id or '/' in source_id:
        raise ValueError('Invalid source ID')
    folder = str(folder).replace('\\','/').strip('/')
    if any(part in ('..','.') for part in folder.split('/')) or not str(expno).isdigit():
        raise ValueError('Invalid exact acquisition path')
    return source_id+':'+(folder+'/' if folder else '')+str(expno)


def _block_key(record):
    source, folder, expno = (record.get(key) for key in
                           ('ObservedSourceID', 'ObservedFolder', 'ObservedExpno'))
    if (not isinstance(source, str) or not source or not isinstance(folder, str)
            or not isinstance(expno, (str, int)) or isinstance(expno, bool)
            or not str(expno).isascii() or not str(expno).isdigit()):
        return None
    return source, folder, int(expno)//100


def resolve_block_sample_defaults(observations):
    """Resolve hundred-block opt-ins using only the acquisitions already observed."""
    groups = {}
    for record in observations:
        # Inherited values are never direct evidence on a subsequent resolution,
        # including when records from separately scanned scopes are combined.
        evidence = record.get('ObservedSampleTagEvidence', '')
        if evidence.startswith(('SAMPLE_ID_FOLDER default from ', 'SAMPLE_ID_BLOCK default from ')):
            record.update(ObservedSampleID='', ObservedSampleTagEvidence='', ObservedSampleTagIssue='')
        key = _block_key(record)
        if key is not None:
            groups.setdefault(key, []).append(record)
    for records in groups.values():
        declared = [r for r in records if r.get('_BlockSampleEvidence')]
        if not declared:
            continue
        defaults = {r['_BlockSampleID'] for r in declared if r.get('_BlockSampleID')}
        problems = []
        if len(defaults) > 1:
            problems.append('Conflicting SAMPLE_ID_BLOCK defaults in monitored hundred block')
        for record in records:
            issue = record.get('_BlockSampleIssue', '') or record.get('ObservedSampleTagIssue', '')
            if issue:
                problems.append('EXPNO '+str(record.get('ObservedExpno', ''))+': '+issue)
            direct = record.get('ObservedSampleID')
            if direct and len(defaults) == 1 and direct not in defaults:
                problems.append('EXPNO '+str(record.get('ObservedExpno', ''))+': conflicting direct SAMPLE_ID in SAMPLE_ID_BLOCK')
        origin = '; '.join('EXPNO '+str(r.get('ObservedExpno', ''))+': '+r['_BlockSampleEvidence'] for r in declared[:8])
        if len(declared) > 8:
            origin += '; ... ('+str(len(declared))+' declaring acquisitions)'
        evidence = 'SAMPLE_ID_BLOCK default from '+origin
        for record in records:
            if record.get('ObservedSampleID') or record.get('ObservedSampleTagIssue'):
                continue  # Preserve direct valid IDs and direct errors.
            record['ObservedSampleTagEvidence'] = evidence
            if problems or len(defaults) != 1:
                record['ObservedSampleID'] = ''
                record['ObservedSampleTagIssue'] = '; '.join(problems[:8]) or 'No coherent SAMPLE_ID_BLOCK default'
            else:
                record['ObservedSampleID'] = next(iter(defaults))
    return observations


def resolve_folder_sample_defaults(observations):
    """Resolve explicit defaults only across acquisitions already in this scan.

    No filesystem inspection, cached inheritance, parent-folder matching, sample
    creation, or human-link writes happen here. A direct valid tag always wins.
    """
    resolve_block_sample_defaults(observations)
    blocked = {_block_key(r) for r in observations if r.get('_BlockSampleEvidence')}
    groups = {}
    for record in observations:
        key = (record.get('ObservedSourceID'),record.get('ObservedFolder'))
        if not isinstance(key[0],str) or not key[0] or not isinstance(key[1],str):
            continue
        groups.setdefault(key,[]).append(record)
    for records in groups.values():
        declared = [r for r in records if r.get('_FolderSampleEvidence')]
        if not declared:
            continue  # Preserve legacy observations/signatures byte-for-byte.
        defaults = {r['_FolderSampleID'] for r in declared if r.get('_FolderSampleID')}
        problems = []
        for record in records:
            issue = record.get('_FolderSampleIssue','')
            # Unreadable titles could hide another explicit declaration. Missing
            # titles alone remain eligible for a readable sibling's default.
            if (not issue and 'Title' not in record and record.get('ObservedSampleTagIssue')
                    and not record.get('ObservedSampleTagEvidence', '').startswith('SAMPLE_ID_BLOCK default from ')):
                issue = record['ObservedSampleTagIssue']
            if issue:
                problems.append('EXPNO '+str(record.get('ObservedExpno',''))+': '+issue)
        if len(defaults)>1:
            problems.append('Conflicting SAMPLE_ID_FOLDER defaults in monitored folder')
        origin = '; '.join('EXPNO '+str(r.get('ObservedExpno',''))+': '+r['_FolderSampleEvidence'] for r in declared[:8])
        if len(declared)>8:
            origin += '; ... ('+str(len(declared))+' declaring acquisitions)'
        evidence = 'SAMPLE_ID_FOLDER default from '+origin
        for record in records:
            if _block_key(record) in blocked:
                continue  # Explicit blocks supersede folder defaults, even when broken.
            if record.get('ObservedSampleID') and not record.get('ObservedSampleTagIssue'):
                continue  # Acquisition-local SAMPLE_ID has explicit priority.
            if record.get('ObservedSampleTagIssue'):
                continue  # Invalid direct evidence is never replaced by a guess.
            if problems or len(defaults)!=1:
                record['ObservedSampleID'] = ''
                record['ObservedSampleTagIssue'] = '; '.join(problems[:8]) or 'No coherent SAMPLE_ID_FOLDER default'
                record['ObservedSampleTagEvidence'] = evidence
            else:
                record['ObservedSampleID'] = next(iter(defaults))
                record['ObservedSampleTagEvidence'] = evidence
    return observations


def match_instrument(frequency, profiles=()):
    """Return the unique raw proton-frequency match; ambiguity stays unresolved."""
    if isinstance(frequency, bool) or not isinstance(frequency, (int, float)) or not math.isfinite(frequency) or frequency <= 0:
        return None
    matches = [profile for profile in profiles
               if abs(frequency - profile.proton_mhz) <= profile.tolerance_mhz]
    return matches[0] if len(matches) == 1 else None


def scan_source(source, instrument_profiles=()):
    """Return observations and errors. Missing roots never mean deletion.

    scandir names are filtered BEFORE is_dir/is_symlink/stat or descent. Every
    symlink is excluded, including roots and metadata files.
    """
    try:
        root = safe_read_path(source.path)
    except (OSError,ValueError):
        return [], ['source root is protected, unavailable, or contains a symlink: '+source.id]
    if not root.is_dir():
        return [], ['source unavailable: '+source.id]
    observations, errors = [], []
    targets = set(source.targets) if source.targets is not None else None
    excluded = set(source.exclude) | {'211@950','.git','.svn','inventory','scripts'}
    def walk(directory):
        try:
            with os.scandir(directory) as entries:
                children = sorted(entries,key=lambda item:item.name)
        except OSError as error:
            errors.append(f'{source.id}: cannot read {directory.relative_to(root)} ({error.__class__.__name__})')
            return
        for entry in children:
            # Hard stop even if an exclusion is nested. Do not inspect the entry.
            if entry.name in excluded or '211@950' in entry.name or entry.name.startswith('.'):
                continue
            try:
                if entry.is_symlink() or not entry.is_dir(follow_symlinks=False):
                    continue
                path = Path(entry.path)
                relative_path = path.relative_to(root).as_posix()
                if targets is not None and not any(target == relative_path or target.startswith(relative_path+'/') for target in targets):
                    continue
                if entry.name.isdigit():
                    # Bruker experiment directories may be partial: record if any
                    # acquisition marker exists, without assuming finished transfer.
                    markers = (path/'acqus',path/'acqu',path/'fid',path/'ser')
                    if any(not p.is_symlink() and p.is_file() for p in markers):
                        relative = path.relative_to(root)
                        folder = relative.parent.as_posix()
                        if folder == '.': folder = ''
                        fields = observe(path)
                        if is_derived_split(fields):
                            continue
                        fields.update(AcquisitionKey=acquisition_key(source.id,folder,entry.name),
                                      SourcePath=str(path),ObservedSourceID=source.id,
                                      ObservedFolder=folder,ObservedExpno=entry.name)
                        if source.instrument:
                            fields['ObservedInstrument'] = source.instrument
                        else:
                            profile = match_instrument(fields.get('ObservedField_MHz'), instrument_profiles)
                            if profile is not None:
                                fields['ObservedInstrument'] = profile.id
                                evidence = (f'instrument={profile.id} from proton frequency '
                                            f'{fields["ObservedField_MHz"]:g} MHz; configured match '
                                            f'±{profile.tolerance_mhz:g} MHz')
                                previous = fields.get('ObservedProtocolEvidence', '')
                                fields['ObservedProtocolEvidence'] = previous+'; '+evidence if previous else evidence
                        observations.append(fields)
                        continue
                if entry.name not in ('pdata',):
                    walk(path)
            except OSError as error:
                errors.append(f'{source.id}: cannot observe {entry.name} ({error.__class__.__name__})')
    walk(root)
    resolve_folder_sample_defaults(observations)
    return observations, errors
