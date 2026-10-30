"""Bounded Bruker metadata observations; raw bytes never imply completeness."""
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

MAX_METADATA_BYTES = 256_000


def text_file(path):
    if path.is_symlink() or not path.is_file():
        return ''
    with path.open('rb') as handle:
        return handle.read(MAX_METADATA_BYTES).decode('utf-8', errors='replace')


def parameters(text):
    return {m[1]: m[2].strip() for m in re.finditer(r'^##\$([A-Za-z0-9_]+)=\s*([^\n]*)',text,re.M)}


def numeric(value):
    try:
        result = float(value)
        return result if result == result and abs(result) != float('inf') else None
    except (ValueError, TypeError):
        return None


def pulse_type(pulse):
    lower = pulse.lower()
    for token, kind in [('hetnoe','noe'),('nnh_tr','eta_xy'),('t1sat','water_t1sat'),('cpmg','cpmg'),('rexet','cpmg'),('t1rho','r1rho'),('hsqctret','r1rho'),('t2','r2'),('t1','r1'),('hsqc','hsqc'),('hmbc','hmbc'),('cosy','cosy'),('noesy','noesy')]:
        if token in lower:
            return kind
    return 'unknown'


def observe(path):
    path = Path(path)
    acqus = text_file(path/'acqus')
    params = parameters(acqus)
    audit = text_file(path/'audita.txt')
    pulse = params.get('PULPROG','').strip('<> ')
    if not pulse:
        pulse_text = text_file(path/'pulseprogram')
        match = re.search(r'^;\s*([^\s;]+)',pulse_text,re.M)
        pulse = match[1] if match else ''
    raw_bytes = 0
    raw_files = []
    for name in ('fid','ser'):
        raw = path/name
        if not raw.is_symlink() and raw.is_file():
            info = raw.stat()
            size = info.st_size
            raw_bytes += size
            raw_files.append((name,size,info.st_mtime_ns,info.st_ino,info.st_dev))
    errors = []
    if not acqus:
        errors.append('acqus missing or empty')
    if re.search(r'created by split\s+\d+:.*?from EXPNO\s+\d+',audit,re.I|re.S):
        errors.append('derived split experiment; original acquisition requires review')
    # Stable acquisition identity excludes file mtime, raw bytes and location. It is
    # a clue for reconciliation, never authorization to copy human assessments.
    identity = {k:params.get(k,'') for k in ('DATE_START','DATE','PULPROG','BF1','SFO1','NUC1','TD')}
    # audita.txt continues changing during processing and is deliberately excluded.
    acquired = numeric(params.get('DATE_START')) or numeric(params.get('DATE'))
    fingerprint = hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest() if acquired and acquired > 0 and pulse else ''
    fields = {'ObservedPulseProgram':pulse, 'ObservedType':pulse_type(pulse),
              'RawPresent':raw_bytes > 0, 'RawBytes':raw_bytes,
              'ObservedFingerprint':fingerprint, 'ScanError':'; '.join(errors), '_RawFileStats':raw_files}
    if acquired and acquired > 0:
        fields['ObservedAcquiredAt'] = acquired
    field = numeric(params.get('BF1'))
    if field is not None:
        fields['ObservedField_MHz'] = field
    temperature = numeric(params.get('TE'))
    if temperature is not None and 0 < temperature <= 1000:
        fields['ObservedTemperature_C'] = round(temperature-273.15,6)
    return fields
