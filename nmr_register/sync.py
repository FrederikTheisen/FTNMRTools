"""Changed-only registration, with stable-key reconciliation after ambiguity."""
from collections import defaultdict
import json
import math
import time
from datetime import datetime, timezone
from zoneinfo import ZoneInfo,ZoneInfoNotFoundError
from .discovery import scan_source
from .schema import MACHINE_FIELDS,CONDITIONS,canonical_experiment_type,FILL_ONLY_IDENTITY_FIELDS
from .grist import GristError
from .state import observation_warning,sufficient_fingerprint,MISSING_FINGERPRINT_WARNING
from .acquisition_transition import reviewed
from .experiment_taxonomy import classify_pulse


def _tuple(fields, observed=False):
    if not observed and 'incomplete acquisition identity' in str(fields.get('MigrationIssue') or '').lower():
        return None  # Blank workbook folder is unresolved, not an explicit root acquisition.
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
    if human == 'DelayList_s':
        try:
            left, right = json.loads(confirmed), json.loads(value)
            return (isinstance(left,list) and isinstance(right,list) and len(left) == len(right)
                    and all(type(a) in (int,float) and type(b) in (int,float)
                        and math.isfinite(a) and math.isfinite(b) and a >= 0 and b >= 0
                        and math.isclose(a,b,rel_tol=1e-8,abs_tol=1e-9) for a,b in zip(left,right)))
        except (TypeError,ValueError,OverflowError):
            return False
    if str(confirmed) == str(value):
        return True
    try:
        number,observed = float(confirmed),float(value)
        if human == 'AcquiredAt' and number % 86400 == 0:
            # Workbook dates have day precision, represented as UTC midnight.
            # Bruker timestamps have second precision in the acquisition zone.
            return datetime.fromtimestamp(number,timezone.utc).date() == datetime.fromtimestamp(observed,_acquisition_zone()).date()
        if human == 'D1_s':
            return math.isclose(number,observed,rel_tol=1e-8,abs_tol=1e-6)
        if human in ('SpinLock_Hz','SpinLockOffset_Hz'):
            return math.isclose(number,observed,rel_tol=1e-6,abs_tol=.1)
        if human == 'Field_MHz' and number.is_integer():
            return round(observed) == int(number)
        return abs(number-observed) < 0.01
    except (TypeError,ValueError,OverflowError):
        return False


METADATA_PAIRS = (('Instrument','ObservedInstrument'),('ExperimentType','ObservedType'),
                  ('AcquiredAt','ObservedAcquiredAt'),('Field_MHz','ObservedField_MHz'),
                  ('PulseProgram','ObservedPulseProgram'),('Temperature_C','ObservedTemperature_C'),
                  ('D1_s','ObservedD1_s'),('DelayList_s','ObservedDelayList_s'),
                  ('SpinLock_Hz','ObservedSpinLock_Hz'),('SpinLockOffset_Hz','ObservedSpinLockOffset_Hz'),
                  ('SpectralWidths_ppm','ObservedSpectralWidths_ppm'),('SpectralOffsets_ppm','ObservedSpectralOffsets_ppm'))


SOURCE_FACTS = dict(METADATA_PAIRS + (('AcquisitionStatus','ObservedAcquisitionStatus'),))


def _status_facts_for_remote(facts, machine, remote):
    # Missing copied metadata does not undo an already documented terminal run.
    # An explicit newer attempt/deletion does invalidate that prior outcome.
    if remote.get("AcquisitionStatus") == "Dummy":
        return {**facts, "AcquisitionStatus": "Dummy"}
    from .bruker import STATUS_CONTRADICTION_PREFIX
    terminal = {'Completed','Aborted','Finished','Stopped','Done','Acquired'}
    if (facts.get('AcquisitionStatus') in ('Unknown','Setup')
            and remote.get('AcquisitionStatus') in terminal
            and not machine.get('ObservedAcquisitionEvidence','').startswith(STATUS_CONTRADICTION_PREFIX)):
        retained = {'Finished':'Completed','Done':'Completed','Acquired':'Completed','Stopped':'Aborted'}.get(
            remote['AcquisitionStatus'],remote['AcquisitionStatus'])
        return {**facts,'AcquisitionStatus':retained}
    return facts


def _canonical_payload(payload):
    """Translate local observations; unknown facts are absent, never blank overrides."""
    facts = {}
    numeric = {'AcquiredAt','Field_MHz','Temperature_C','D1_s','SpinLock_Hz','SpinLockOffset_Hz'}
    for canonical, observed in SOURCE_FACTS.items():
        value = payload.get(observed, payload.get(canonical))
        if canonical == 'AcquisitionStatus':
            value = {'Finished':'Completed','Stopped':'Aborted'}.get(value,value)
            if value in ('Completed','Aborted','Setup','Unknown'):
                facts[canonical] = value
        elif canonical == 'ExperimentType':
            # Cached observations from older collectors may contain only the
            # broad Calibration token or an unknown assignment classification.
            # Reclassify those from the unchanged raw pulse name so retries
            # cannot restore a coarse label after the taxonomy upgrade.
            if value in (None, '', 'unknown', 'Unknown', 'calibration', 'Calibration'):
                pulse = payload.get('ObservedPulseProgram', payload.get('PulseProgram', ''))
                value = classify_pulse(pulse) or value
            value = canonical_experiment_type(value)
            if value:
                facts[canonical] = value
        elif canonical in numeric:
            if type(value) in (int,float) and math.isfinite(value):
                if canonical in ('AcquiredAt','Field_MHz') and value <= 0:
                    continue
                facts[canonical] = math.floor(value+.5) if canonical == 'Field_MHz' else value
        elif isinstance(value,str) and value.strip():
            if canonical != 'Instrument' or value.strip().lower() != 'unknown':
                facts[canonical] = value
    for name in ('Title','TitleSource'):
        if name in payload and isinstance(payload[name],str):
            facts[name] = payload[name]  # A readable empty source explicitly clears comments.
    machine = {key:value for key,value in payload.items()
               if key not in SOURCE_FACTS and key not in SOURCE_FACTS.values() and key not in ('Title','TitleSource')}
    forbidden = (set(machine)|set(facts))-set(MACHINE_FIELDS['Experiments'])
    if forbidden:
        raise ValueError('Refusing non-machine collector fields: '+', '.join(sorted(forbidden)))
    return machine, facts


def _conflicts(machine, remote):
    differences = []
    for human, observed in METADATA_PAIRS:
        confirmed = remote.get(human)
        value = machine.get(observed)
        if confirmed not in ('',None) and value not in ('',None,'unknown') and not _consistent(human,confirmed,value):
            differences.append('observed '+human+' differs from confirmed value')
    return differences


def _identity_conflicts(machine,remote,accepted=(),transition_kind=None):
    source,folder,expno = (machine.get('Observed'+name) for name in ('SourceID','Folder','Expno'))
    if not isinstance(source,str) or not source or not isinstance(folder,str) or expno is None or not str(expno).isdigit():
        return ['observed acquisition identity is missing or invalid; reconcile before binding']
    expected = (source,folder,str(expno))
    expected_key = source+':'+(folder+'/' if folder else '')+str(expno)
    if machine.get('AcquisitionKey') != expected_key:
        return ['collector key differs from observed acquisition identity; reconcile before binding']
    reasons = []
    if any(remote.get(name) not in (None,'') and str(remote[name]) != expected[index]
           for index,name in enumerate(('SourceID','Folder','Expno'))):
        reasons.append('observed source/folder/expno differs from confirmed identity')
    fingerprint, previous = machine.get('ObservedFingerprint'), remote.get('ObservedFingerprint')
    if fingerprint and previous and fingerprint != previous and previous not in accepted:
        reasons.append('acquisition fingerprint differs from registered identity; reconcile before reusing assessment')
    if transition_kind in ('setup-to-setup','setup-to-raw') and previous != fingerprint:
        if reviewed(remote):
            reasons.append(transition_kind+' transition has existing review or result; reconcile before reusing assessment')
        if remote.get('RawPresent') or (remote.get('RawBytes') or 0) > 0:
            reasons.append(transition_kind+' transition conflicts with previously registered raw data; reconcile before reusing assessment')
    return reasons


def _identity_defaults(machine, remote,accepted=(),transition_kind=None):
    """Fill blanks only after exact observed-key and existing-identity validation."""
    if _identity_conflicts(machine, remote,accepted,transition_kind):
        return {}
    return {name:machine['Observed'+name] for name in FILL_ONLY_IDENTITY_FIELDS
            if remote.get(name) in (None,'') and remote.get(name) != machine['Observed'+name]}


def _legacy_metadata_only(reason):
    allowed = {'observed '+human+' differs from confirmed value' for human,_ in METADATA_PAIRS}
    return bool(reason) and all(piece.strip() in allowed for piece in reason.split(';'))


def _source_review_reason(local,remote):
    # A local copy warning cannot erase a serious warning already in Grist.
    if remote and not _legacy_metadata_only(remote) and not observation_warning(remote):
        if local and not observation_warning(local) and local not in remote:
            return local if remote in local else remote+'; '+local
        return remote
    return local or (remote if not _legacy_metadata_only(remote) else '')


def _metadata_fields(observed,remote):
    # Former override-versus-observation differences are obsolete. Preserve any
    # unrelated diagnostic and every genuine identity reconciliation warning.
    difference = remote.get('MetadataDifference','') or ''
    result = {'MetadataDifference':'' if _legacy_metadata_only(difference) else difference}
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
    # The same hourly GET catches edits to source-owned canonical facts without
    # rescanning files or issuing any additional remote reads. Queue corrections
    # before transport so ambiguous responses keep the existing durable retry.
    already_pending = {key for key,_ in pending}
    requeue = []
    for key,cached in state.cached_observations():
        if key in already_pending:
            continue
        machine,facts = _canonical_payload(cached)
        accepted,transition_kind = state.transition_lineage(key,machine.get('ObservedFingerprint'))
        matches = keys.get(key,[])
        exact = manual.get(_tuple(machine,observed=True),[])
        if len(matches)>1 or len(exact)>1 or (matches and exact and matches[0]['id']!=exact[0]['id']):
            requeue.append(key)
            continue
        row = matches[0] if matches else exact[0] if len(exact)==1 else None
        if row and (bool(_differences(_status_facts_for_remote(facts,machine,row['fields']),row['fields'])) or _identity_conflicts(machine,row['fields'],accepted,transition_kind) or _identity_defaults(machine,row['fields'],accepted,transition_kind)
                    or (row['fields'].get('ReconciliationReason')==MISSING_FINGERPRINT_WARNING and sufficient_fingerprint(machine.get('ObservedFingerprint')))
                    or state.resolved_copy_warning(key,machine.get('ObservedFingerprint'),row['fields'].get('ReconciliationReason',''))
                    or any(reason and not _legacy_metadata_only(reason) and not observation_warning(reason) for reason in
                        (machine.get('ReconciliationReason',''),row['fields'].get('ReconciliationReason','')))):
            requeue.append(key)
    state.requeue_observations(requeue)
    pending = state.pending()
    create, update, acknowledged = [], [], []
    blocked = set()
    conflict_updates = {}
    conflicts = 0
    for key,payload in pending:
        machine, facts = _canonical_payload(payload)
        accepted,transition_kind = state.transition_lineage(key,machine.get('ObservedFingerprint'))
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
        if row:
            if not machine.get('ObservedFingerprint') and row['fields'].get('ObservedFingerprint'):
                machine.pop('ObservedFingerprint',None)  # A partial read must not erase the last identity guard.
            metadata_updates.pop(row['id'],None)
            computed = _metadata_fields(machine,row['fields'])
            machine['MetadataDifference'] = computed['MetadataDifference']
            if machine['MetadataDifference']: metadata_warnings.add(row['id'])
            else: metadata_warnings.discard(row['id'])
            reason = machine.get('ReconciliationReason','')
            if _legacy_metadata_only(reason): reason = ''
            remote_reason = row['fields'].get('ReconciliationReason','')
            if not reason and state.resolved_copy_warning(key,machine.get('ObservedFingerprint'),remote_reason):
                remote_reason = ''  # State recomputed copy groups from current observations.
            if sufficient_fingerprint(machine.get('ObservedFingerprint')):
                if reason==MISSING_FINGERPRINT_WARNING:reason=''
                if remote_reason==MISSING_FINGERPRINT_WARNING:remote_reason=''
            reason = _source_review_reason(reason,remote_reason)
            identity_reasons = _identity_conflicts(machine,row['fields'],accepted,transition_kind)
            if identity_reasons:
                reason = '; '.join(filter(None,[reason,*(r for r in identity_reasons if r not in reason)]))
            machine['ReconciliationReason'] = reason
            if reason and not observation_warning(reason):
                blocked.add(key)
                conflicts += 1
            else:
                machine.update(_identity_defaults(machine,row['fields'],accepted,transition_kind))
                machine.update(_status_facts_for_remote(facts,machine,row['fields']))
            if key not in blocked:
                state.mark_transition_attempt(key,machine.get('ObservedFingerprint'),row['fields'].get('ObservedFingerprint'))
            changed = _differences(machine,row['fields'])
            if changed:
                update.append((key,{'id':row['id'],'fields':changed}))
            elif key not in blocked:
                acknowledged.append(key)
        else:
            reason = machine.get('ReconciliationReason','')
            if _legacy_metadata_only(reason):
                reason = ''
            reason = '; '.join(filter(None,[reason,*_identity_conflicts(machine,{})]))
            machine['ReconciliationReason'] = reason
            if reason and not observation_warning(reason):
                blocked.add(key)
                conflicts += 1
            else:
                machine.update(_identity_defaults(machine,{}))
                machine.update(facts)
            initialize = {'Field_MHz':None,'Temperature_C':None,'AcquiredAt':None,
                'RawBytes':None,'MetadataDifference':''}
            initialize.update({name:None for name in ('D1_s','SpinLock_Hz','SpinLockOffset_Hz')})
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
            state.acknowledge([key for key,_ in batch if key not in blocked])
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
        found,failures = scan_source(source, config.instrument_profiles)
        observations.extend(found)
        errors.extend(failures)
    from .result_discovery import prepare_results
    plan = prepare_results(config, observations)
    result_summary = {
        'results_discovered': plan.discovered,
        'results_selected': len(plan.selections),
        'results_changed': sum(state.get_meta('result_sent:'+selection.acquisition_key) !=
                               [selection.path, selection.sha256] for selection in plan.selections),
        'result_discovery_issues': plan.issues,
    }
    if dry_run:
        # Preview never queues results, changes observations or contacts Grist.
        return {'observed':len(observations),'errors':errors,'dry_run':True, **result_summary}
    state.bind_document(config.api_base)
    state.observe(observations)
    state.set_meta('last_scan',time.time())
    state.set_meta('last_error','; '.join(errors))
    result = {'observed':len(observations),'errors':errors,'pending':len(state.pending()), **result_summary}
    if client is None:
        return result
    try:
        result.update(flush(state,client,config.batch_size))
        if result['conflicts']:
            errors.append(str(result['conflicts'])+' acquisition identities require reconciliation')
            state.set_meta('last_error','; '.join(errors))
        from .results import attach_results
        result.update(attach_results(state, client, plan.selections, config.max_result_bytes,
                                     discovery_issues=plan.issues, eligible_keys=plan.eligible_keys))
        state.set_meta('last_sync',time.time())
        sync_scan_status(state,client,config,'; '.join(errors))
    except (GristError,OSError,ValueError) as error:
        state.set_meta('last_error',str(error))
        raise
    result['pending'] = len(state.pending())
    return result
