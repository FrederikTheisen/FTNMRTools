"""Local credential enrollment. Secrets are never accepted as CLI arguments."""
import getpass
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys


def credential_path(home):
    return Path(home)/'credential'


def _private_file(path):
    info = path.lstat()
    if stat.S_ISLNK(info.st_mode) or not stat.S_ISREG(info.st_mode) or info.st_mode & 0o777 != 0o600:
        raise ValueError('Credential file must be a regular file with mode 0600')
    if info.st_uid != os.getuid():
        raise ValueError('Credential must be owned by the current user')
    return path


def read_token(home,doc_url=None,credential_file=None):
    path = credential_path(home)
    try:
        _private_file(path)
    except FileNotFoundError:
        # Reuse the user's existing enrollment, without copying or displaying it.
        from .config import Config
        path = Path(credential_file) if credential_file else Path.home()/'.config/nmr-register/grist.json'
        _private_file(path)
        try:
            value = json.loads(path.read_text())
            token = value['api_key']
            enrolled_url = value['document_url']
        except (ValueError,KeyError,TypeError):
            raise ValueError('Invalid enrolled Grist credential JSON') from None
        if not doc_url or Config(doc_url=enrolled_url).api_base != Config(doc_url=doc_url).api_base:
            raise ValueError('Enrolled credential document host and ID must match collector configuration')
    else:
        token = path.read_text().strip()
    if not isinstance(token,str) or not token.strip() or '\n' in token.strip() or '\r' in token.strip():
        raise ValueError('Invalid credential file')
    return token.strip()


def prompt_token():
    if sys.stdin.isatty():
        return getpass.getpass('Grist API key (stored locally, never echoed): ').strip()
    if shutil.which('zenity') and (os.environ.get('DISPLAY') or os.environ.get('WAYLAND_DISPLAY')):
        result = subprocess.run(['zenity','--password','--title=Grist API key'],capture_output=True,text=True)
        if result.returncode == 0:
            return result.stdout.strip()
        raise ValueError('Credential enrollment cancelled')
    raise ValueError('Credential enrollment requires a terminal or a graphical session with zenity')


def enroll(home, client_factory=None):
    home = Path(home).expanduser().resolve()
    if any((parent/'.git').exists() for parent in [home,*home.parents]):
        raise ValueError('Credential storage must be outside a source repository')
    token = prompt_token()
    if not token or '\n' in token or '\r' in token:
        raise ValueError('Empty or invalid API key')
    # Test access read-only before replacing an existing credential.
    if client_factory:
        client_factory(token).test_connection()
    home.mkdir(parents=True,exist_ok=True,mode=0o700)
    path = credential_path(home)
    tmp = path.with_name('credential.tmp')
    fd = os.open(tmp,os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o600)
    try:
        with os.fdopen(fd,'w') as handle:
            handle.write(token+'\n')
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp,path)
    finally:
        if tmp.exists(): tmp.unlink()
    return path
