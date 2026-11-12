"""Bounded Bruker metadata observations; raw bytes never imply completeness."""
from nmr_register.split import audit_is_derived_split
import hashlib
import json
import re
import stat
import math
import os
from datetime import datetime
from decimal import Decimal, InvalidOperation
from pathlib import Path

MAX_METADATA_BYTES = 256_000
MAX_LIST_VALUES = 4096


def _safe_path(path):
    """Check ancestors before inspecting children; never resolve a symlink."""
    path = Path(path)
    if not path.is_absolute() or any(part in ('.', '..') or '211@950' in part for part in path.parts):
        raise ValueError('metadata path is relative, traversing, or protected')
    current = Path(path.anchor)
    for part in path.parts[1:]:
        current = current/part
        try:
            info = current.lstat()
        except FileNotFoundError:
            break
        if stat.S_ISLNK(info.st_mode):
            raise ValueError('metadata path contains a symlink')
    return path


def _read_metadata(path):
    path = _safe_path(path)
    try:
        info = path.lstat()
    except FileNotFoundError:
        return None
    if not stat.S_ISREG(info.st_mode):
        raise ValueError('metadata is not a regular file')
    if info.st_size > MAX_METADATA_BYTES:
        raise ValueError('metadata exceeds size limit')
    with path.open('rb') as handle:
        data = handle.read(MAX_METADATA_BYTES+1)
    if len(data) > MAX_METADATA_BYTES:
        raise ValueError('metadata exceeds size limit')
    return data


def text_file(path):
    try:
        data = _read_metadata(path)
    except (OSError, ValueError):
        return ''
    return data.decode('utf-8', errors='replace') if data is not None else ''


def parameters(text):
    # Keep this legacy scalar representation for acquisition identity. Protocol
    # additions must not change the established reconciliation fingerprint.
    return {m[1]: m[2].strip() for m in re.finditer(r'^##\$([A-Za-z0-9_]+)=\s*([^\n]*)', text, re.M)}


def protocol_parameters(text):
    """Read complete JCAMP blocks without scientific processing dependencies."""
    result = {}
    for match in re.finditer(r'^##\$(\w+)=(.*?)(?=^##|\Z)', text, re.M | re.S):
        key, value = match.groups()
        if key in result:
            raise ValueError('duplicate acquisition parameter '+key)
        value = '\n'.join(line.split('$$', 1)[0] for line in value.splitlines()).strip()
        if value.startswith('('):
            array = re.fullmatch(r'\(0\.\.(\d+)\)\s*(.*)', value, re.S)
            if not array or int(array[1])+1 > MAX_LIST_VALUES:
                raise ValueError('invalid acquisition array '+key)
            values = re.findall(r'<[^>]*>|[^\s]+', array[2])
            if len(values) != int(array[1])+1:
                raise ValueError('incomplete acquisition array '+key)
            value = [item[1:-1] if item.startswith('<') and item.endswith('>') else item for item in values]
        else:
            value = value[1:-1] if value.startswith('<') and value.endswith('>') else value
        result[key] = value
    return result


def numeric(value):
    try:
        result = float(value)
        return result if result == result and abs(result) != float('inf') else None
    except (ValueError, TypeError, OverflowError):
        return None


def _indexed(params, key, index):
    values = params.get(key)
    return numeric(values[index]) if isinstance(values, list) and index < len(values) else None


def pulse_type(pulse):
    lower = pulse.strip().lower()
    if lower == 'ibs_1d_sculpting':
        return '1d'
    if lower in ('zg','zg30','15n_calib','13c_calib'):
        return 'calibration'
    for token, kind in [('hetnoe','noe'),('hsqcnoe','noe'),('nnh_tr','eta_xy'),('t1sat','water_t1sat'),('cpmg','cpmg'),('rexet','cpmg'),('t1rho','r1rho'),('hsqctret','r1rho'),('t2','r2'),('t1','r1'),('hsqc','hsqc'),('hmbc','hmbc'),('cosy','cosy'),('noesy','noesy')]:
        if token in lower:
            return kind
    return 'unknown'


def _active_program(text):
    return re.sub(r'\s+', '', ''.join(line.split(';', 1)[0] for line in text.splitlines()
                                    if not line.lstrip().startswith('#'))).lower()


def _profile(pulse, program):
    """Names plus saved active-code proof; comments alone never establish units."""
    active = _active_program(program)
    profiles = {
        'hsqct1etf3gpsitc3d.3.ts3': ('vd', ('definelist<delay>vd_list=<$vdlist>',
            '"d31=vd_list[l11]"', '"counter1=(d31/40m)"', 'loto8timescounter1')),
        'hsqctretf3gpsitc3d.3.ts3': ('vd', ('definelist<delay>vd_list=<$vdlist>',
            '"d31=vd_list[l11]"', '"p25=d31"',
            '"plw25=plw3*pow((p21*cnst27*4/1000000),2)"', 'fq=cnst30(sfohz):f3', 'pl25:f3', '(p25ph1):f3')),
        '15n_t1rho': ('vd', ('definelist<delay>vd_list=<$vdlist>',
            '"p15=vd_list[t1loop]"', '"plw8=plw3*(pow((p21*4*cnst8/1000000),2))"', 'pl8p15ph0')),
        't1rhon15_rk.nl': ('vp', ('pl8vpph0', 'ivp', '"p17=p18+1m-vp"',
            '"d17=d1-p17"', '"d8=p8+vp*0.25-p1*2-2u"', '"d9=vp*0.5-p1*4-24u"')),
    }
    name = pulse.lower()
    profile = profiles.get(name)
    if profile and all(token in active for token in profile[1]):
        # Reject competing assignments to a nominal field calibration. Initial
        # delay/pulse zeros are intentional in the reviewed programs.
        output = {'hsqctretf3gpsitc3d.3.ts3':'plw25', '15n_t1rho':'plw8'}.get(name)
        if output and len(re.findall(r'"'+output+r'=', active)) != 1:
            return None
        calibration = ('cnst27','cnst30','p21','plw3') if name == 'hsqctretf3gpsitc3d.3.ts3' else (('cnst8','p21','plw3') if name == '15n_t1rho' else ())
        if any(re.search(r'"'+key+r'=', active) for key in calibration):
            return None  # Saved code overrides acqus: its scalar is not authoritative.
        return profile[0]
    return None


def _reference_name(reference):
    # Archived absolute instrument paths supply only a filename, never a path to
    # traverse. Reject traversal even though only the basename would be used.
    reference = str(reference).strip().replace('\\', '/')
    parts = reference.split('/')
    if not reference or any(part in ('.', '..') or '211@950' in part for part in parts) or '\x00' in reference:
        raise ValueError('unsafe or empty list reference')
    name = parts[-1]
    if not name or name.startswith('.') or not re.fullmatch(r'[\w .+@=-]+', name):
        raise ValueError('unsafe list reference')
    return name


def _selected_list(path, kind, reference):
    """Select acquisition-local declared lists; conflicting copies fail closed."""
    name = _reference_name(reference)
    candidates = [path/(kind+'list'), path/name, path/'lists'/kind/name]
    found = []
    for candidate in dict.fromkeys(candidates):
        data = _read_metadata(candidate)
        if data is not None:
            found.append((candidate, data))
    if not found:
        raise ValueError('referenced '+kind+' list is missing')
    # Permit whitespace/comment differences, but never choose one different list
    # merely because it is first, newest, or the sole unrelated archive file.
    def tokens(data):
        text = data.decode('utf-8')
        return [item for line in text.splitlines()
                for item in re.split(r'\s+', re.split(r'\$\$|#|;', line, 1)[0].strip()) if item]
    selected = tokens(found[0][1])
    if any(tokens(data) != selected for _, data in found[1:]):
        raise ValueError('conflicting '+kind+' list candidates')
    if not selected or len(selected) > MAX_LIST_VALUES:
        raise ValueError('empty or oversized '+kind+' list')
    chosen, data = found[0]
    return selected, chosen.relative_to(path).as_posix(), hashlib.sha256(data).hexdigest()


_DELAY = re.compile(r'([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)(s|m|ms|u|us|µs|μs)?')


def _delay_seconds(token):
    match = _DELAY.fullmatch(token)
    if not match:
        raise ValueError('invalid programmed delay')
    scale = {None:1, 's':1, 'm':.001, 'ms':.001, 'u':1e-6, 'us':1e-6, 'µs':1e-6, 'μs':1e-6}[match[2]]
    value = numeric(match[1])
    if value is not None:
        try:
            value = float(Decimal(match[1])*Decimal(str(scale)))
        except (InvalidOperation,OverflowError):
            value = None
    if value is None or numeric(value) is None or value < 0:
        raise ValueError('invalid programmed delay')
    return value


def _protocol(path, acqus, pulse):
    fields = {'ObservedD1_s':None, 'ObservedDelayList_s':'',
              'ObservedSpinLock_Hz':None, 'ObservedSpinLockOffset_Hz':None, 'ObservedField_MHz':None}
    notes = []
    try:
        params = protocol_parameters(acqus)
    except ValueError as error:
        fields['ObservedProtocolEvidence'] = str(error)+'; protocol settings unknown'
        return fields
    if acqus:
        notes.append('acqus#'+hashlib.sha256(acqus.encode()).hexdigest()[:12])
    d1 = _indexed(params, 'D', 1)
    if d1 is not None and d1 >= 0:
        fields['ObservedD1_s'] = d1
        notes.append('D1=D[1] s')
    elif 'D' in params:
        notes.append('D[1] invalid; D1 unknown')
    try:
        data = _read_metadata(path/'pulseprogram')
        program = data.decode('utf-8', errors='replace') if data is not None else ''
    except (OSError, ValueError) as error:
        program = ''
        notes.append('pulseprogram unavailable: '+str(error))
    kind = _profile(pulse, program)
    if program:
        notes.append('pulseprogram#'+hashlib.sha256(program.encode()).hexdigest()[:12])
    if kind:
        notes.append('reviewed saved-code profile='+pulse)
    elif pulse_type(pulse) in ('r1','r1rho','r2','cpmg'):
        notes.append('no reviewed saved-code profile; relaxation list seconds unverified')
    for source_kind in ('vd', 'vp', 'vc'):
        reference = params.get(source_kind.upper()+'LIST', '')
        if not reference:
            if kind == source_kind:
                notes.append(source_kind.upper()+'LIST not declared; delay list unknown')
            continue
        try:
            values, selected, digest = _selected_list(path, source_kind, reference)
            notes.append(source_kind.upper()+'LIST '+str(reference)+' => '+selected+'#'+digest[:12])
            if kind == source_kind:
                delays = [_delay_seconds(value) for value in values]
                fields['ObservedDelayList_s'] = json.dumps(delays, separators=(',', ':'), allow_nan=False)
                notes.append('programmed seconds, ordered with duplicates')
                if pulse.lower() == 'hsqct1etf3gpsitc3d.3.ts3':
                    notes.append('40 ms quantization: programmed list, not measured duration')
            else:
                shown = ','.join(values[:24])+(',...' if len(values) > 24 else '')
                notes.append(('cycle counts' if source_kind == 'vc' else 'uninterpreted list')+'=['+shown+']; not converted to seconds')
        except (OSError, ValueError, UnicodeError) as error:
            notes.append(source_kind.upper()+'LIST unresolved: '+str(error))
    spin_index = {'hsqctretf3gpsitc3d.3.ts3':27, '15n_t1rho':8}.get(pulse.lower()) if kind else None
    if spin_index is not None:
        frequency = _indexed(params, 'CNST', spin_index)
        calibration = (_indexed(params, 'P', 21), _indexed(params, 'PLW', 3))
        isotope = str(params.get('NUC3','')).upper().replace(' ','')
        if isotope == '15N' and frequency is not None and 0 <= frequency <= 1e6 and all(value is not None and value > 0 for value in calibration):
            fields['ObservedSpinLock_Hz'] = frequency
            notes.append('nominal spin-lock=CNST['+str(spin_index)+'] Hz, P[21]/PLW[3] calibration')
        else:
            notes.append('nominal spin-lock unknown: invalid CNST/calibration')
        if pulse.lower() == 'hsqctretf3gpsitc3d.3.ts3':
            offset = _indexed(params, 'CNST', 30)
            if isotope == '15N' and offset is not None and abs(offset) <= 1e9:
                fields['ObservedSpinLockOffset_Hz'] = offset
                notes.append('spin-lock offset=CNST[30] Hz relative to carrier')
        else:
            notes.append('spin-lock offset units unverified; unknown')
    elif pulse_type(pulse) == 'r1rho':
        settings = ['PLW['+str(i)+']='+str(value)+' W' for i in (3,8,25)
                    if (value := _indexed(params, 'PLW', i)) is not None]
        notes.append('; '.join(settings)+'; RF power is not spin-lock field Hz')
    proton = []
    for channel in range(1, 9):
        nucleus = str(params.get('NUC'+str(channel), '')).upper().replace(' ', '')
        field = numeric(params.get('BF'+str(channel)))
        if nucleus == '1H' and field is not None and 0 < field <= 10000:
            proton.append((channel,field))
    if proton and all(abs(value-proton[0][1]) < .01 for _,value in proton):
        fields['ObservedField_MHz'] = proton[0][1]
        notes.append('proton MHz=BF'+str(proton[0][0])+' (NUC=1H)')
    else:
        notes.append('proton MHz unknown: no unique valid 1H channel')
    fields['ObservedProtocolEvidence'] = '; '.join(notes)
    return fields



# DATE/DATE_START are integer Unix seconds. Audit timestamps carry fractions and
# explicit offsets. Saved TopSpin 3.7/4.4.1 examples show DATE written up to 1.238s
# before the terminal audit entry: allow at most 2s, never an arbitrary clock gap.
AUDIT_TERMINAL_TOLERANCE_S = 2
STATUS_CONTRADICTION_PREFIX = "audita.txt: latest attempt unresolved: "
_AUDIT_TIME = r'\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}:\d{2}(?:\.\d{1,6})? ?[+-]\d{4}'
_AUDIT_RECORD = re.compile(
    r'^[ \t]*\(\s*(\d+),\s*<([^<>\r\n]*)>,\s*<([^<>\r\n]*)>,\s*<([^<>\r\n]*)>,'
    r'\s*<([^<>\r\n]*)>,\s*<([^<>\r\n]*)>,\s*<(.*?)>\)[ \t]*(?=\r?\n|\Z)', re.S | re.M)


def _audit_timestamp(value):
    if not re.fullmatch(_AUDIT_TIME, value):
        raise ValueError('invalid audit timestamp or missing timezone')
    # Both "2025-11-07 18:41:32.195 +0100" and TopSpin4's
    # "2025-05-28T15:31:41.371+0000" occur in saved audit trails.
    offset = value[-5:]
    hours, minutes = int(offset[1:3]), int(offset[3:5])
    if minutes >= 60 or hours > 14 or (hours == 14 and minutes):
        raise ValueError('invalid audit timezone offset')
    compact = value[:10]+'T'+value[11:].replace(' ', '')
    result = datetime.fromisoformat(compact)
    if abs(result.utcoffset().total_seconds()) > 14*3600:
        raise ValueError('invalid audit timezone offset')
    return result.timestamp()


def _supported_audit_version(value):
    return bool(re.fullmatch(r'TopSpin (?:3\.6\.5|3\.7(?:\.\d+)?|4\.4\.1|5\.0\.0)', value.strip()))


def acquisition_status(acqus, audit, raw_present=None):
    """Instrument-run outcome from the latest matching acquisition audit record.

    This never establishes live running, transfer completeness, or assessment.
    Only the reviewed TopSpin 3.6.5/3.7/4.4.1/5.0.0 formats are accepted; missing timezone,
    metadata, terminal events or uncertain generations produce Unknown. Raw file
    presence is deliberately independent from the instrument's logged outcome.
    With clean metadata, raw_present=False permits Setup: this describes only the
    last observed local copy and never proves transfer completeness.
    """
    def result(status, evidence):
        return {'ObservedAcquisitionStatus':status,
                'ObservedAcquisitionEvidence':'audita.txt: '+evidence}
    setup = False
    if raw_present is False:
        try:
            parsed = protocol_parameters(acqus)
            pulse = parsed.get('PULPROG')
            setup = isinstance(pulse,str) and bool(pulse.strip())
        except (ValueError,TypeError):
            pass
    def unknown(reason, contradictory=False, allow_setup=True):
        if contradictory:
            reason = 'latest attempt unresolved: '+reason
        if setup and allow_setup:
            return result('Setup', reason+'; no nonempty fid/ser in local copy; transfer not established')
        return result('Unknown', reason)
    if not audit.strip():
        return unknown('missing, unreadable, empty, or oversized audit')
    if audit_is_derived_split(audit):
        return unknown('derived split; original acquisition outcome is not this dataset', allow_setup=False)
    try:
        params = protocol_parameters(acqus)
        start, end = numeric(params.get('DATE_START')), numeric(params.get('DATE'))
        if start is None or end is None or start <= 0 or end < start or not start.is_integer() or not end.is_integer():
            return unknown('valid integer DATE_START and DATE required')
        titles = re.findall(r'^##TITLE=\s*Audit trail,\s*([^\r\n]+)', audit, re.M)
        if len(titles) != 1 or not _supported_audit_version(titles[0]):
            return unknown('unsupported or missing audit format/version')
        if len(re.findall(r'^##AUDIT TRAIL=', audit, re.M)) != 1 or len(re.findall(r'^##END=', audit, re.M)) != 1:
            return unknown('missing or conflicting audit boundaries')
        section = re.search(r'^##AUDIT TRAIL=[^\n]*\n(.*?)^##END=[^\n]*(?:\n|\Z)', audit, re.M | re.S)
        if not section:
            return unknown('incomplete audit records')
        if any(line.strip() and not line.lstrip().startswith('$$') for line in audit[section.end():].splitlines()):
            return unknown('unparsed content after audit end')
        records, previous_end, previous_number = [], 0, 0
        for record in _AUDIT_RECORD.finditer(section[1]):
            if any(line.strip() and not line.lstrip().startswith('$$')
                   for line in section[1][previous_end:record.start()].splitlines()):
                return unknown('malformed audit record')
            if re.search(r'^\s*\(\s*\d+\s*,\s*<', record[7], re.M):
                return unknown('malformed nested audit record; previous record is not closed')
            number = int(record[1])
            if number <= previous_number:
                return unknown('conflicting audit record numbers')
            previous_number, previous_end = number, record.end()
            if record[5].strip().lower() == 'go4':
                if not _supported_audit_version(record[6]):
                    return unknown('unsupported acquisition record version')
                when = _audit_timestamp(record[2].strip())
                records.append((number, when, record[7].strip()))
        if any(line.strip() and not line.lstrip().startswith('$$')
               for line in section[1][previous_end:].splitlines()):
            return unknown('malformed or truncated audit record')
        if not records:
            return unknown('no acquisition-process audit records')
        number, when, body = records[-1]
        if body == 'acquisition in progress':
            return unknown('entry '+str(number)+' records in progress; copied audit cannot establish live running', contradictory=True)
        if body == 'raw data deleted by zg':
            return unknown('entry '+str(number)+' deletes raw data; no later matching terminal event', contradictory=True)
        if not re.match(r'^created by zg on\s', body):
            return unknown('latest acquisition event is unrecognized or aborted')
        identifiers = re.findall(r"\bon data set with UUID '([0-9a-f]{8}(?:-[0-9a-f]{4}){3}-[0-9a-f]{12})'", body, re.I)
        started = re.findall(r'^\s*started at ('+_AUDIT_TIME+r'),', body, re.M)
        finished = re.findall(r'^\s*completed at ('+_AUDIT_TIME+r')\s*$', body, re.M)
        interrupted = re.findall(r"^\s*terminated by command '(stop|halt)' at ("+_AUDIT_TIME+r')\s*$', body, re.M)
        stopped = [item[1] for item in interrupted]
        if len(identifiers) != 1 or len(started) != 1 or len(finished)+len(stopped) != 1:
            return unknown('latest acquisition lacks one unambiguous start and recognized terminal event')
        # Extra start/terminal wording is a conflict, even if one line happens
        # to fit our known format. Never find a completion inside a comment.
        if len(re.findall(r'\bstarted at\b', body)) != 1 or len(re.findall(r'\bcompleted at\b|\bterminated by command\b|\baborted\b', body)) != 1:
            return unknown('conflicting acquisition terminal events')
        start_time = _audit_timestamp(started[0])
        terminal_text = (finished or stopped)[0]
        terminal_time = _audit_timestamp(terminal_text)
        if math.floor(start_time) != start:
            return unknown('latest acquisition start does not match DATE_START; stale or restarted generation', contradictory=True)
        if terminal_time < start_time or abs(terminal_time-end) > AUDIT_TERMINAL_TOLERANCE_S:
            return unknown('terminal event does not match DATE within 2 s')
        if abs(when-terminal_time) > AUDIT_TERMINAL_TOLERANCE_S:
            return unknown('terminal record time conflicts with its event')
        # Two terminal acquisitions beginning within the same stored second
        # cannot be distinguished by DATE_START alone.
        for _, _, earlier in records[:-1]:
            earlier_start = re.findall(r'^\s*started at ('+_AUDIT_TIME+r'),', earlier, re.M)
            if earlier.startswith('created by zg on ') and earlier_start and math.floor(_audit_timestamp(earlier_start[0])) == start:
                return unknown('multiple acquisition generations share DATE_START second')
        status = 'Completed' if finished else 'Aborted'
        verb = 'completed' if finished else "terminated by command '"+interrupted[0][0]+"'"
        return result(status, 'entry '+str(number)+' zg '+verb+' at '+terminal_text+
                      '; DATE_START matched, DATE within 2 s')
    except (ValueError, TypeError, OverflowError, OSError):
        return unknown('invalid acquisition metadata or audit timestamps')


MAX_SAMPLE_ID_CHARS = 256


def _title_sample_ids(text, marker='SAMPLE_ID'):
    """Only whole, explicitly marked lines declare a sample; other prose is inert."""
    identifiers, malformed = [], False
    for line in text.splitlines():
        if line.startswith(marker+'='):
            value = line[len(marker)+1:]
            if (not value or value != value.strip() or len(value) > MAX_SAMPLE_ID_CHARS
                    or not value.isprintable()):
                malformed = True
            else:
                identifiers.append(value)
        elif re.match(r'^\s*'+re.escape(marker)+r'\b', line, re.I):
            malformed = True
    return identifiers, malformed


def title_observation(path):
    """Read only expno/pdata/1/title; processing titles never supply acquisition facts.

    Both sample markers are optional. Explicit malformed/conflicting declarations
    fail closed, while ordinary comments and absent markers require no action.
    """
    fields = {'ObservedSampleID':'', 'ObservedSampleTagEvidence':'', 'ObservedSampleTagIssue':''}
    try:
        path = _safe_path(path)
        data = _read_metadata(path/'pdata'/'1'/'title')
        if data is None:
            return fields
        text = data.decode('utf-8')
        fields['Title'] = ''.join(line for line in text.splitlines(keepends=True)
            if not any(_title_sample_ids(line, marker)[0] and not _title_sample_ids(line, marker)[1]
                       for marker in ('SAMPLE_ID', 'SAMPLE_ID_FOLDER')))
        fields['TitleSource'] = 'pdata/1/title'
        values, malformed = _title_sample_ids(text)
        issues = []
        if malformed:
            issues.append('pdata/1/title: malformed SAMPLE_ID line')
        if len(set(values)) > 1:
            issues.append('Conflicting SAMPLE_ID declarations in acquisition title')
        if values or malformed:
            fields['ObservedSampleTagEvidence'] = 'Explicit SAMPLE_ID in pdata/1/title'
        if len(set(values)) == 1 and not issues:
            fields['ObservedSampleID'] = values[0]
        fields['ObservedSampleTagIssue'] = '; '.join(issues)
        folder_values, folder_malformed = _title_sample_ids(text, 'SAMPLE_ID_FOLDER')
        if folder_values or folder_malformed:
            folder_issues = []
            if folder_malformed:
                folder_issues.append('pdata/1/title: malformed SAMPLE_ID_FOLDER line')
            if len(set(folder_values)) > 1:
                folder_issues.append('Conflicting SAMPLE_ID_FOLDER declarations in acquisition title')
            fields['_FolderSampleID'] = folder_values[0] if folder_values and not folder_issues else ''
            fields['_FolderSampleEvidence'] = 'pdata/1/title'
            fields['_FolderSampleIssue'] = '; '.join(folder_issues)
    except (OSError, ValueError, UnicodeError) as error:
        issue = 'Acquisition title unreadable or unsafe ('+error.__class__.__name__+')'
        fields['ObservedSampleTagIssue'] = issue
        fields['_TitleIssue'] = issue+'; comments unavailable'
    return fields


def sample_tag_observation(path):
    # Compatibility surface for callers interested only in explicit sample tags.
    return {key:value for key,value in title_observation(path).items() if key.startswith('ObservedSample')}


def observe(path):
    path = _safe_path(path)
    acqus = text_file(path/'acqus')
    params = parameters(acqus)
    try:
        audit_bytes = _read_metadata(path/'audita.txt')
        audit = audit_bytes.decode('utf-8') if audit_bytes is not None else ''
        audit_inspection_clean = True
    except (OSError,ValueError,UnicodeError):
        audit, audit_inspection_clean = '', False
    pulse = params.get('PULPROG','').strip('<> ')
    if not pulse:
        pulse_text = text_file(path/'pulseprogram')
        match = re.search(r'^;\s*([^\s;]+)', pulse_text, re.M)
        pulse = match[1] if match else ''
    raw_bytes = 0
    raw_files = []
    raw_inspection_clean = True
    for name in ('fid', 'ser'):
        raw = path/name
        try:
            info = raw.lstat()
        except FileNotFoundError:
            continue
        if not stat.S_ISREG(info.st_mode):
            raw_inspection_clean = False
            continue
        size = info.st_size
        raw_bytes += size
        raw_files.append((name,size,info.st_mtime_ns,info.st_ino,info.st_dev))
    errors = []
    if not acqus:
        errors.append('acqus missing, empty, unreadable, or oversized')
    if audit_is_derived_split(audit):
        errors.append('derived split experiment; original acquisition requires review')
    # Stable acquisition identity excludes file mtime, raw bytes and location. It is
    # a clue for reconciliation, never authorization to copy human assessments.
    identity = {k:params.get(k,'') for k in ('DATE_START','DATE','PULPROG','BF1','SFO1','NUC1','TD')}
    # audita.txt continues changing during processing and is deliberately excluded.
    acquired = numeric(params.get('DATE_START')) or numeric(params.get('DATE'))
    fingerprint = hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest() if acquired and acquired > 0 and pulse else ''
    observed_type = pulse_type(pulse)
    if observed_type == 'hsqc':
        # The evolved nucleus distinguishes carbon HSQC from nitrogen HSQC;
        # a carbon decoupling channel or sample labelling alone does not.
        indirect = parameters(text_file(path/'acqu2s')).get('NUC1','').strip('<> ').upper()
        if indirect == '13C':
            observed_type = '13c_hsqc'
    fields = {'ObservedPulseProgram':pulse, 'ObservedType':observed_type,
              'RawPresent':raw_bytes > 0, 'RawBytes':raw_bytes,
              'ObservedFingerprint':fingerprint, 'ScanError':'; '.join(errors), '_RawFileStats':raw_files,
              'ObservedField_MHz':None, 'ObservedTemperature_C':None, 'ObservedAcquiredAt':None}
    if acquired and acquired > 0:
        fields['ObservedAcquiredAt'] = acquired
    protocol = _protocol(path, acqus, pulse)
    temperature = numeric(params.get('TE'))
    if temperature is not None and 0 < temperature <= 1000:
        fields['ObservedTemperature_C'] = round(temperature-273.15, 6)
    fields.update(protocol)
    # A missing audit can accompany a setup copy; an unsafe/unreadable audit
    # cannot establish a clean setup observation. No raw bytes are opened.
    fields.update(acquisition_status(acqus, audit,
        raw_present=(raw_bytes > 0) if acqus and not errors and raw_inspection_clean and audit_inspection_clean else None))
    titles = title_observation(path)
    title_issue = titles.pop('_TitleIssue', '')
    if title_issue:
        fields['ScanError'] = '; '.join(filter(None, [fields['ScanError'], title_issue]))
    fields.update(titles)
    return fields
