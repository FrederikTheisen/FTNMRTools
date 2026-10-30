"""Validated, atomic, seven-day Grist document backups with internal attachments."""
from datetime import datetime,timezone
import os
from pathlib import Path
import sqlite3
import tempfile

REQUIRED_TABLES = ('_grist_Tables','_grist_Tables_column','Experiments','Samples','Series')


def validate_backup(path,required_tables=REQUIRED_TABLES):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError('Backup must be a regular file')
    with path.open('rb') as handle:
        if handle.read(16) != b'SQLite format 3\x00':
            raise ValueError('Backup is not a SQLite Grist document')
    connection = sqlite3.connect(path.resolve().as_uri()+'?mode=ro&immutable=1',uri=True)
    try:
        if connection.execute('PRAGMA integrity_check').fetchone()[0] != 'ok':
            raise ValueError('Backup SQLite integrity check failed')
        tables = {row[0] for row in connection.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        missing = set(required_tables)-tables
        if missing:
            raise ValueError('Backup missing required tables: '+', '.join(sorted(missing)))
        # External attachment references cannot constitute a self-contained backup.
        count = 0
        if '_grist_Attachments' in tables:
            cols = {row[1] for row in connection.execute('PRAGMA table_info(_grist_Attachments)')}
            if 'storageId' in cols and connection.execute("SELECT count(*) FROM _grist_Attachments WHERE storageId IS NOT NULL AND storageId != ''").fetchone()[0]:
                raise ValueError('Backup contains external attachments; switch document storage to internal first')
            count = connection.execute('SELECT count(*) FROM _grist_Attachments').fetchone()[0]
            if count:
                if '_gristsys_Files' not in tables:
                    raise ValueError('Backup attachment content table is missing')
                file_cols = {row[1] for row in connection.execute('PRAGMA table_info(_gristsys_Files)')}
                if {'ident','data'}.issubset(file_cols) and 'fileIdent' in cols:
                    absent = connection.execute('''SELECT count(*) FROM _grist_Attachments a
                      LEFT JOIN _gristsys_Files f ON a.fileIdent=f.ident
                      WHERE f.ident IS NULL OR f.data IS NULL''').fetchone()[0]
                    if absent:
                        raise ValueError('Backup is missing internal attachment content')
        return {'valid':True,'tables':sorted(tables),'attachments':count,'bytes':path.stat().st_size}
    finally:
        connection.close()


def _atomic_bytes(content,destination,required_tables=REQUIRED_TABLES):
    destination = Path(destination)
    destination.parent.mkdir(parents=True,exist_ok=True,mode=0o700)
    fd,name = tempfile.mkstemp(prefix='.nmr-backup-',suffix='.grist',dir=destination.parent)
    try:
        with os.fdopen(fd,'wb') as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        validate_backup(name,required_tables)
        os.replace(name,destination)
        directory_fd = os.open(destination.parent,os.O_RDONLY)
        try: os.fsync(directory_fd)
        finally: os.close(directory_fd)
    finally:
        if os.path.exists(name): os.unlink(name)
    return destination


def create_backup(client,directory,premigration=False,now=None,required_tables=REQUIRED_TABLES):
    now = now or datetime.now(timezone.utc)
    directory = Path(directory)
    directory.mkdir(parents=True,exist_ok=True,mode=0o700)
    # Daily slots retain seven distinct successful days; premigration snapshots
    # are immutable and never participate in rotation.
    name = ('premigration-'+now.strftime('%Y%m%dT%H%M%S%fZ') if premigration else 'daily-'+now.strftime('%Y-%m-%d'))+'.grist'
    path = directory/name
    _atomic_bytes(client.download(),path,required_tables)
    if not premigration:
        daily = sorted(directory.glob('daily-????-??-??.grist'),reverse=True)
        for old in daily[7:]:
            if not old.is_symlink(): old.unlink()
    return {'path':str(path),**validate_backup(path,required_tables)}


def restore_copy(path,destination,required_tables=REQUIRED_TABLES):
    """Verify and restore to a new local file, suitable for a Grist import test.

    Never changes a remote document and never overwrites a destination.
    """
    validate_backup(path,required_tables)
    destination = Path(destination)
    if destination.exists() or destination.is_symlink():
        raise ValueError('Restore destination must not exist')
    return _atomic_bytes(Path(path).read_bytes(),destination,required_tables)
