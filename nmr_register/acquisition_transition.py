"""Local evidence for generation-preserving collector metadata transitions."""
CORE_KEYS = ('DATE_START','PULPROG','BF1','SFO1','NUC1','TD')


def raw_ids(context):
    return {v[0]:tuple(v[3:5]) for v in context.get('_RawFileStats', []) if len(v) >= 5 and v[1] > 0}


def clean_setup(record, context):
    """Positive raw absence with no retained evidence of a prior acquisition."""
    identity = context.get('_AcquisitionIdentity') or {}
    stats = context.get('_RawFileStats')
    if (context.get('_RawGenerationFiles') or record.get('RawPresent') is not False
            or record.get('RawBytes') != 0 or not isinstance(stats, (list, tuple))):
        return False
    if any(not isinstance(item, (list, tuple)) or len(item) != 5
           or item[0] not in ('fid', 'ser') or item[1] != 0 for item in stats):
        return False
    if identity:
        return bool(identity.get('version') == 1 and identity.get('raw_clean'))
    # Older collectors supplied explicit file stats and Setup status, but no
    # acquisition identity. A present, unclean identity never falls back here.
    return record.get('ObservedAcquisitionStatus') == 'Setup'


def permitted_transition(previous, old_context, current, context):
    """Return the explicit transition kind, never authorize a new raw generation."""
    old = old_context.get('_AcquisitionIdentity') or {}
    new = context.get('_AcquisitionIdentity') or {}
    if new.get('version') != 1 or not new.get('raw_clean'):
        return None
    old_setup = clean_setup(previous, old_context)
    if old_setup and clean_setup(current, context):
        return 'setup-to-setup'
    if old_setup and current.get('RawPresent') and raw_ids(context):
        return 'setup-to-raw'
    if not (old.get('version') == 1 and old.get('raw_clean') and previous.get('RawPresent')
            and current.get('RawPresent') and raw_ids(old_context) == raw_ids(context)
            and raw_ids(context)):
        return None
    a, b = old.get('parameters', {}), new.get('parameters', {})
    try:
        positive_start = float(a.get('DATE_START') or 0) > 0
    except (TypeError, ValueError):
        positive_start = False
    if not positive_start or a.get('DATE_START') != b.get('DATE_START'):
        return None
    if new.get('latest_in_progress') and previous.get('ObservedAcquisitionStatus') in ('Completed', 'Aborted', 'Finished', 'Stopped'):
        return None
    if old.get('deletion_seen') != new.get('deletion_seen') or new.get('latest_deleted'):
        return None
    # Every raw transition needs positive same-generation audit evidence. A
    # stable inode and integer DATE_START alone cannot distinguish quick restarts.
    terminal = current.get('ObservedAcquisitionStatus') in ('Completed', 'Aborted')
    old_anchors = old.get('progress_anchors') or []
    new_anchors = new.get('progress_anchors') or []
    old_anchor = old_anchors[-1] if old_anchors else None
    finalized_progress = (not old.get('uuid') and old.get('progress_verified')
        and old.get('latest_in_progress') and old_anchor is not None
        and terminal and bool(new.get('uuid'))
        and bool(new_anchors) and old_anchor == new_anchors[-1]
        and _anchor_time(old_anchor) == new.get('terminal_start'))
    same_terminal = (bool(old.get('uuid')) and terminal
        and old['uuid'] == new.get('uuid'))
    same_progress = (not old.get('uuid') and not new.get('uuid')
        and old.get('progress_verified') and new.get('progress_verified')
        and old.get('latest_in_progress') and new.get('latest_in_progress')
        and old_anchor is not None and bool(new_anchors)
        and old_anchor == new_anchors[-1])
    if not (same_terminal or finalized_progress or same_progress):
        return None
    if any(a.get(k) != b.get(k) for k in CORE_KEYS):
        if not finalized_progress:
            return None
        return 'verified-progress-finalization'
    return 'mutable-end-date'


def _anchor_time(anchor):
    from .bruker import _audit_timestamp
    try:
        return _audit_timestamp(anchor[1].strip())
    except (ValueError, TypeError, IndexError, OverflowError):
        return None


def reviewed(record):
    return (bool(record.get('Assessment') or record.get('ImportedAssessment') or record.get('AssessmentReason'))
            or any(record.get(k) not in (None, '', 'Not assessed') for k in
                   ('QCAcquisitionCoverage','QCSettings','QCFID','QCSNR','QCSample','QCAnalysis'))
            or any(record.get(k) for k in ('ResultCSV','ResultSHA256','ResultSourcePath')))
