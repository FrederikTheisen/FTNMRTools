"""Select current local numbers tables for sync without walking result archives."""
import csv
from dataclasses import dataclass
import hashlib
from pathlib import Path
from .discovery import acquisition_key, safe_read_path
from .results import _read_selected
from .split import is_derived_split


@dataclass
class ResultPlan:
    selections: list
    issues: list
    discovered: int
    eligible_keys: set


def prepare_results(config, observations):
    """Prefer one source-verified local result; retain valid manifest fallbacks.

    Only acquisitions found by this scan are eligible. Targets, exclusions,
    protected paths and derived splits retain the acquisition scanner's scope.
    Differing local versions block the key instead of guessing from timestamps.
    A bad result cannot make another acquisition's valid results disappear.
    """
    sources = {source.id: source for source in config.sources}
    eligible = {}
    for record in observations:
        if is_derived_split(record):
            continue
        source_id = record.get('ObservedSourceID')
        if source_id not in sources:
            continue
        source = sources[source_id]
        folder, expno = record.get('ObservedFolder', ''), str(record.get('ObservedExpno', ''))
        key = acquisition_key(source_id, folder, expno)
        if record.get('AcquisitionKey') != key:
            raise ValueError('Result discovery acquisition key disagrees with source identity')
        relative = folder+'/'+expno if folder else expno
        if source.targets is not None and relative not in source.targets:
            continue
        parts = Path(relative).parts
        if any(part in set(source.exclude) | {'inventory', 'scripts'} or
               '211@950' in part or part.startswith('.') for part in parts):
            continue
        path = safe_read_path(Path(source.path)/relative)
        eligible[key] = path

    baselines = {}
    if config.results_manifest:
        if config.result_source_id not in sources:
            raise ValueError('results_manifest requires a configured result_source_id')
        with Path(config.results_manifest).open(newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            if not {'acquisition', 'numbers_file', 'numbers_sha256'}.issubset(reader.fieldnames or []):
                raise ValueError('Manifest requires acquisition, numbers_file and numbers_sha256')
            for row in reader:
                relative = row['acquisition'].strip()
                # Filter by current scan keys before inspecting any selected file.
                key = config.result_source_id+':'+relative
                if key in eligible:
                    baselines.setdefault(key, []).append(row)

    selections, issues, discovered = [], [], 0
    for key, folder in sorted(eligible.items()):
        rows = baselines.get(key, [])
        local = []
        try:
            candidates = sorted(p for p in folder.iterdir()
                                if p.name.endswith('_numbers.csv') and not p.name.startswith('.'))
            for path in candidates:
                safe_read_path(path)
                if not path.is_file() or path.stat().st_size > config.max_result_bytes:
                    raise ValueError('Local result is missing, not regular, or exceeds configured size: '+str(path))
                with path.open('rb') as handle:
                    content = handle.read(config.max_result_bytes+1)
                if len(content) > config.max_result_bytes:
                    raise ValueError('Local result exceeds configured size: '+str(path))
                digest = hashlib.sha256(content).hexdigest()
                local.append(_read_selected(str(path), digest, key, config.max_result_bytes))
            discovered += len(local)
            if len({s.sha256 for s in local}) > 1:
                raise ValueError('Multiple different local result files require an explicit current version: '+
                                 ', '.join(s.path for s in local))
            if local:
                selected = local[0]  # Multiple byte-identical files are equivalent.
                # Preserve an already selected archived path for unchanged bytes.
                # A stale manifest hash must not veto a newly analysed local file.
                if len(rows) == 1 and rows[0]['numbers_sha256'].strip().lower() == selected.sha256:
                    try:
                        selected = _read_selected(rows[0]['numbers_file'], selected.sha256,
                                                  key, config.max_result_bytes)
                    except (OSError, ValueError):
                        pass  # The verified local copy remains available.
                selections.append(selected)
            elif len(rows) > 1:
                raise ValueError('Multiple manifest result selections require reconciliation')
            elif rows:
                selections.append(_read_selected(rows[0]['numbers_file'],
                                  rows[0]['numbers_sha256'].strip().lower(), key,
                                  config.max_result_bytes))
        except (OSError, ValueError) as error:
            issues.append({'key': key, 'reason': str(error)})

    # Keep the existing rule that content must not be shared by different keys.
    by_hash = {}
    for selection in selections:
        by_hash.setdefault(selection.sha256, []).append(selection.acquisition_key)
    shared = {key for keys in by_hash.values() if len(keys) > 1 for key in keys}
    for key in sorted(shared):
        issues.append({'key': key, 'reason': 'Result hash shared by different acquisitions requires reconciliation'})
    return ResultPlan([s for s in selections if s.acquisition_key not in shared], issues, discovered, set(eligible))
