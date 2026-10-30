"""Directory traversal with protected paths pruned before any inspection."""
import os
import stat
from pathlib import Path
from .bruker import observe


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


def scan_source(source):
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
                        fields.update(AcquisitionKey=acquisition_key(source.id,folder,entry.name),
                                      SourcePath=str(path),ObservedSourceID=source.id,
                                      ObservedFolder=folder,ObservedExpno=entry.name)
                        if source.instrument:
                            fields['ObservedInstrument'] = source.instrument
                        observations.append(fields)
                        continue
                if entry.name not in ('pdata',):
                    walk(path)
            except OSError as error:
                errors.append(f'{source.id}: cannot observe {entry.name} ({error.__class__.__name__})')
    walk(root)
    return observations, errors
