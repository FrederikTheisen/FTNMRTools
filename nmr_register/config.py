"""Non-secret collector configuration."""
from dataclasses import dataclass, field, asdict
import json
import os
from pathlib import Path
from urllib.parse import urlparse

DEFAULT_DOC = 'https://docs.getgrist.com/e259CYFK6mGK/registernmrtest'

def default_home():
    return Path(os.environ.get('NMR_REGISTER_HOME', Path.home()/'.local/share/nmr-register'))

@dataclass
class Source:
    id: str
    path: str
    instrument: str = ''
    exclude: list = field(default_factory=lambda: ['211@950'])
    targets: list | None = None

@dataclass
class Config:
    doc_url: str = DEFAULT_DOC
    sources: list = field(default_factory=list)
    state_dir: str = ''
    results_manifest: str = ''
    result_source_id: str = ''
    max_result_bytes: int = 2_000_000
    credential_file: str = ''
    timeout: float = 30
    batch_size: int = 100

    @property
    def home(self):
        return Path(self.state_dir) if self.state_dir else default_home()

    @property
    def api_base(self):
        url = urlparse(self.doc_url)
        parts = [part for part in url.path.split('/') if part]
        if url.scheme != 'https' and not (url.scheme == 'http' and url.hostname in ('localhost','127.0.0.1')):
            raise ValueError('Grist URL must use HTTPS')
        if not parts or url.username or url.password or url.query or url.fragment:
            raise ValueError('Expected a Grist document URL without credentials or query')
        # Grist UI may prefix the document with /o/<organization>.
        doc = parts[2] if parts[0] == 'o' and len(parts) > 2 else parts[0]
        return f'{url.scheme}://{url.netloc}/api/docs/{doc}'

    def validate(self):
        self.api_base
        ids = [s.id for s in self.sources]
        if len(ids) != len(set(ids)) or any(not i or ':' in i or '/' in i for i in ids):
            raise ValueError('Source IDs must be unique, nonempty and contain no colon or slash')
        for source in self.sources:
            if not Path(source.path).is_absolute():
                raise ValueError('Source paths must be absolute')
            if '211@950' not in source.exclude:
                source.exclude.append('211@950')
        if not 1 <= self.batch_size <= 500 or self.max_result_bytes < 1 or self.timeout <= 0:
            raise ValueError('Invalid collector bounds')
        return self


def load_config(path=None):
    path = Path(path) if path else default_home()/'config.json'
    value = json.loads(path.read_text())
    value['sources'] = [Source(**source) for source in value.get('sources', [])]
    return Config(**value).validate()


def save_config(config, path=None):
    config.validate()
    path = Path(path) if path else config.home/'config.json'
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    temporary = path.with_name(path.name+'.tmp')
    fd = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with os.fdopen(fd,'w') as handle:
        json.dump(asdict(config),handle,indent=2)
        handle.write('\n')
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary,path)
    return path
