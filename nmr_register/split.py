"""Positive evidence for TopSpin derived datasets; never infer from paths."""
import re

SPLIT_SCAN_ERROR = 'derived split experiment; original acquisition requires review'
SPLIT_EVIDENCE = 'audita.txt: derived split; original acquisition outcome is not this dataset'


def _audit_is_single_fid_extraction(audit):
    """Require structured rser provenance and no recorded actual acquisition."""
    # Import at call time: bruker uses this module while defining its parser.
    from .bruker import _AUDIT_RECORD
    if len(re.findall(r'^##TITLE=\s*Audit trail, TopSpin [^\r\n]+$', audit, re.M)) != 1:
        return False
    if len(re.findall(r'^##AUDIT TRAIL=', audit, re.M)) != 1 or len(re.findall(r'^##END=', audit, re.M)) != 1:
        return False
    section = re.search(r'^##AUDIT TRAIL=[^\n]*\n(.*?)^##END=[^\n]*(?:\n|\Z)', audit, re.M | re.S)
    if not section:
        return False
    def comments_only(text):
        return all(not line.strip() or line.lstrip().startswith('$$') for line in text.splitlines())
    if not comments_only(audit[section.end():]):
        return False
    records = []; end = 0; previous = 0
    for record in _AUDIT_RECORD.finditer(section[1]):
        if not comments_only(section[1][end:record.start()]):
            return False
        if int(record[1]) <= previous or re.search(r'^\s*\(\s*\d+\s*,\s*<', record[7], re.M):
            return False
        records.append(record); end = record.end(); previous = int(record[1])
    if not comments_only(section[1][end:]):
        return False
    # Even an earlier/later acquisition event makes exclusion unsafe: the copy
    # could have been reused as an acquisition dataset.
    if any(re.fullmatch(r'go[0-9]+', r[5].strip(), re.I) for r in records):
        return False
    provenance = re.compile(
        r'rser [1-9][0-9]* from\s*\n[ \t]*[1-9][0-9]* "[^"\r\n]+" "[^"\r\n]+"\s*\n'
        r'[ \t]*data hash MD5: [^\r\n]+\s*\n[ \t]*(?:[0-9A-Fa-f]{2}[ \t]*){16}')
    return any(r[5].strip() == 'mvrow' and provenance.fullmatch(r[7].strip()) for r in records)


def audit_is_derived_split(audit):
    """Recognize the existing explicit TopSpin split provenance statement."""
    return _audit_is_single_fid_extraction(audit) or bool(re.search(r'created by split\s+\d+:.*?from EXPNO\s+\d+', audit, re.I | re.S))


def is_derived_split(fields):
    """Use only generated machine evidence, never titles or reconciliation notes."""
    return (fields.get('ObservedAcquisitionEvidence') == SPLIT_EVIDENCE
            or bool(re.search(r'(?:^|; )'+re.escape(SPLIT_SCAN_ERROR)+r'(?:; |$)', fields.get('ScanError', ''))))
