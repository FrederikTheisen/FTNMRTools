"""Changed-only registration, with stable-key reconciliation after ambiguity."""
from collections import defaultdict
import time
from datetime import datetime, timezone
from zoneinfo import ZoneInfo,ZoneInfoNotFoundError
from .discovery import scan_source
from .schema import MACHINE_FIELDS,CONDITIONS
from .grist import GristError


def _tuple(fields, observed=False):
    prefix = 'Observed' if observed else ''
    values = tuple(str(fields.get(prefix+k,'') if fields.get(prefix+k,'') is not None else '').strip() for k in ('SourceID','Folder','Expno'))
    # Grist numeric Int values may arrive as floats.
    if values[2].endswith('.0'):
        values = (*values[:2],values[2][:-2])
    return values if values[0] and values[2] else None


def _differences(machine, remote):
    result = {key:value for key,value in machine.items() if remote.get(key) != value}
    return result


def _type(value):
    normalized = ''.join(str(value).lower().split()).replace('-', '').replace('_','')
    return {'t1':'r1','t2':'r2','t1rho':'r1rho','hetnoe':'noe','heteronuclearnoe':'noe','15nhsqc':'hsqc'}.get(normalized,normalized)


def _acquisition_zone():
    try:
        return ZoneInfo('Europe/Paris')
    except ZoneInfoNotFoundError:
        # Conda may omit system zoneinfo directories from its search path.
        with open('/usr/share/zoneinfo/Europe/Paris','rb') as handle:
            return ZoneInfo.from_file(handle,key='Europe/Paris')


def _consistent(human,confirmed,value):
    if human == 'ExperimentType':
        return _type(confirmed) == _type(value)
    if str(confirmed) == str(value):
        return True
    try:
        number,observed = float(confirmed),float(value)
        if human == 'AcquiredAt' and number % 86400 == 0:
            # Workbook dates have day precision, represented as UTC midnight.
            # Bruker timestamps have second precision in the acquisition zone.
            return datetime.fromtimestamp(number,timezone.utc).date() == datetime.fromtimestamp(observed,_acquisition_zone()).date()
        if human == 'Field_MHz' and number.is_integer():
            return round(observed) == int(number)
        return abs(number-observed) < 0.01
    except (TypeError,ValueError,OverflowError):
        return False


METADATA_PAIRS = (('Instrument','ObservedInstrument'),('ExperimentType','ObservedType'),
                  ('AcquiredAt','ObservedAcquiredAt'),('Field_MHz','ObservedField_MHz'),
                  ('PulseProgram','ObservedPulseProgram'),('Temperature_C','ObservedTemperature_C'))


def _conflicts(machine, remote):
    differences = []
    for human, observed in METADATA_PAIRS:
        confirmed = remote.get(human)
        value = machine.get(observed)
        if confirmed not in ('',None) and value not in ('',None,'unknown') and not _consistent(human,confirmed,value):
            differences.append('observed '+human+' differs from confirmed value')
    return differences


def _identity_conflicts(machine,remote):
    expected,confirmed = _tuple(machine,observed=True),_tuple(remote)
    return ['observed source/folder/expno differs from confirmed identity'] if expected and confirmed and expected != confirmed else []


def _legacy_metadata_only(reason):
    allowed = {'observed '+human+' differs from confirmed value' for human,_ in METADATA_PAIRS}
    return bool(reason) and all(piece.strip() in allowed for piece in reason.split(';'))


def _metadata_fields(observed,remote):
    result = {'MetadataDifference':'; '.join(_conflicts(observed,remote))}
    if _legacy_metadata_only(remote.get('ReconciliationReason','')):
        result['ReconciliationReason'] = ''
    return result


def flush(state, client, batch_size=100):
    if getattr(client,'base',None):
        state.bind_document(client.base)
    pending = state.pending()
    # One Experiments read per run also notices human corrections when raw files
    # are unchanged. Samples and Series are never polled by the collector.
    rows = client.records('Experiments')
    keys, manual = defaultdict(list), defaultdict(list)
    metadata_updates = {}
    metadata_warnings = set()
    for row in rows:
        fields = row['fields']
        if fields.get('AcquisitionKey'):
            keys[fields['AcquisitionKey']].append(row)
        identity = _tuple(fields)
        if identity:
            manual[identity].append(row)
        computed = _metadata_fields(fields,fields)
        if computed['MetadataDifference']:
            metadata_warnings.add(row['id'])
        # Grist defaults blank Text to ''. Avoid creating hourly empty writes.
        changed = {k:v for k,v in computed.items() if (fields.get(k) or '') != v}
        if changed:
            metadata_updates[row['id']] = {'id':row['id'],'fields':changed}
    create, update, acknowledged = [], [], []
    conflict_updates = {}
    conflicts = 0
    for key,payload in pending:
        forbidden = set(payload)-set(MACHINE_FIELDS['Experiments'])
        if forbidden:
            raise ValueError('Refusing non-machine collector fields: '+', '.join(sorted(forbidden)))
        matches = keys.get(key,[])
        identity = _tuple(payload,observed=True)
        exact = manual.get(identity,[]) if identity else []
        if len(matches)>1 or len(exact)>1 or (matches and exact and matches[0]['id'] != exact[0]['id']):
            conflicts += 1
            for item in matches+exact:
                conflict_updates[item['id']] = {'id':item['id'],'fields':{'ReconciliationReason':'duplicate exact acquisition identity; reconcile rows before binding collector key'}}
            continue
        row = matches[0] if matches else None
        if not row and len(exact)==1:
            if exact[0]['fields'].get('AcquisitionKey') not in ('',None,key):
                conflicts += 1
                conflict_updates[exact[0]['id']] = {'id':exact[0]['id'],'fields':{'ReconciliationReason':'exact manual identity already has a different collector key; reconcile before binding'}}
                continue
            row = exact[0]
        machine = dict(payload)
        if row:
            metadata_updates.pop(row['id'],None)
            computed = _metadata_fields(machine,row['fields'])
            machine['MetadataDifference'] = computed['MetadataDifference']
            if machine['MetadataDifference']: metadata_warnings.add(row['id'])
            else: metadata_warnings.discard(row['id'])
            reason = machine.get('ReconciliationReason','')
            if _legacy_metadata_only(reason): reason = ''
            remote_reason = row['fields'].get('ReconciliationReason','')
            if not reason and remote_reason and not _legacy_metadata_only(remote_reason):
                reason = remote_reason  # Never erase a genuine identity ambiguity.
            identity_reasons = _identity_conflicts(machine,row['fields'])
            if identity_reasons:
                reason = '; '.join(filter(None,[reason,*identity_reasons]))
                conflicts += 1
            machine['ReconciliationReason'] = reason
            changed = _differences(machine,row['fields'])
            if changed:
                update.append((key,{'id':row['id'],'fields':changed}))
            else:
                acknowledged.append(key)
        else:
            initialize = {'Field_MHz':None,'Temperature_C':None,
                'ObservedField_MHz':None,'ObservedTemperature_C':None,'RawBytes':None,'MetadataDifference':''}
            initialize.update({'Snapshot'+name:None for name,typ in CONDITIONS.items() if typ in ('Numeric','Int')})
            create.append((key,{**initialize,**machine}))
    state.acknowledge(acknowledged)
    for identifier,record in conflict_updates.items():
        if identifier in metadata_updates:
            metadata_updates[identifier]['fields'].update(record['fields'])
        else:
            metadata_updates[identifier] = record
    remote_only = list(metadata_updates.values())
    for start in range(0,len(remote_only),batch_size):
        client.update_records('Experiments',remote_only[start:start+batch_size])
    for items,is_create in ((create,True),(update,False)):
        for start in range(0,len(items),batch_size):
            batch = items[start:start+batch_size]
            # Every subgroup stays pending after an ambiguous response. The next
            # GET reconciles committed keys before another create is attempted.
            if is_create:
                client.add_records('Experiments',[fields for _,fields in batch])
            else:
                client.update_records('Experiments',[record for _,record in batch])
            state.acknowledge([key for key,_ in batch])
    return {'created':len(create),'updated':len(update)+len(remote_only),'conflicts':conflicts,
            'metadata_warnings':len(metadata_warnings)}


def sync_scan_status(state,client,config,scan_error=''):
    """One status upsert per hourly run, two calls maximum."""
    scanner = 'nmr-register'
    fields = {'Scanner':scanner,'LastScanAt':state.get_meta('last_scan'),
              'LastSyncAt':state.get_meta('last_sync'),'LastError':scan_error,
              'ScannedRoot':'; '.join(s.path for s in config.sources)}
    rows = client.records('ScanStatus')
    matches = [row for row in rows if row['fields'].get('Scanner')==scanner]
    if len(matches)>1:
        raise GristError('Duplicate ScanStatus collector rows require reconciliation')
    if matches:
        client.update_records('ScanStatus',[{'id':matches[0]['id'],'fields':fields}])
    else:
        client.add_records('ScanStatus',[fields])


def run(config,state,client=None,dry_run=False):
    observations,errors = [],[]
    for source in config.sources:
        found,failures = scan_source(source)
        observations.extend(found)
        errors.extend(failures)
    if dry_run:
        # Dry runs do not alter local observations/outbox or contact Grist.
        return {'observed':len(observations),'errors':errors,'dry_run':True}
    state.bind_document(config.api_base)
    state.observe(observations)
    state.set_meta('last_scan',time.time())
    state.set_meta('last_error','; '.join(errors))
    result = {'observed':len(observations),'errors':errors,'pending':len(state.pending())}
    if client is None:
        return result
    try:
        result.update(flush(state,client,config.batch_size))
        if result['conflicts']:
            errors.append(str(result['conflicts'])+' acquisition identities require reconciliation')
            state.set_meta('last_error','; '.join(errors))
        state.set_meta('last_sync',time.time())
        sync_scan_status(state,client,config,'; '.join(errors))
    except (GristError,OSError,ValueError) as error:
        state.set_meta('last_error',str(error))
        raise
    result['pending'] = len(state.pending())
    return result
