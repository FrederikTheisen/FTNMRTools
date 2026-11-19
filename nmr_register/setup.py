"""Produce setup actions; caller must provide authenticated Grist transport."""
import json
from .schema import schema_actions, trigger_actions, DATE_ONLY_OPTIONS, DATE_ONLY_EXPERIMENT_COLUMNS, LABELLING_CHOICES, acquisition_columns, spectral_columns, PRE_ASSIGNED_SAMPLE_ATTENTION_FORMULA, PRE_LIVE_SAMPLE_ATTENTION_FORMULA, PRE_STUDY_ASSESSMENT_ATTENTION_FORMULA, PRE_CALIBRATION_ATTENTION_FORMULA, PRE_DUMMY_ATTENTION_FORMULA, PRE_DUMMY_EFFECTIVE_FORMULA, PRE_CATEGORY_ATTENTION_FORMULA, PRE_REGISTRATION_REMOVAL_ATTENTION_FORMULA

# Reviewed one-time conversion from split observed/confirmed acquisition facts.
LEGACY_OBSERVED_FACTS = {
    'ObservedInstrument':'Instrument','ObservedType':'ExperimentType','ObservedAcquiredAt':'AcquiredAt',
    'ObservedAcquisitionStatus':'AcquisitionStatus','ObservedField_MHz':'Field_MHz',
    'ObservedTemperature_C':'Temperature_C','ObservedPulseProgram':'PulseProgram',
    'ObservedD1_s':'D1_s','ObservedDelayList_s':'DelayList_s',
    'ObservedSpinLock_Hz':'SpinLock_Hz','ObservedSpinLockOffset_Hz':'SpinLockOffset_Hz',
}
LEGACY_EFFECTIVE_FORMULAS = {
    'EffectiveInstrument':'$Instrument or $ObservedInstrument',
    'EffectiveType':'$ExperimentType or $ObservedType',
    'EffectiveAcquiredAt':'$AcquiredAt or $ObservedAcquiredAt',
    'EffectiveAcquisitionStatus':'"Finished" if $AcquisitionStatus in ("Done", "Acquired") else ($AcquisitionStatus or $ObservedAcquisitionStatus or "Unknown")',
    'EffectiveField_MHz':'$Field_MHz if $Field_MHz is not None else (round($ObservedField_MHz) if $ObservedField_MHz is not None else None)',
    'EffectiveTemperature_C':'$Temperature_C if $Temperature_C is not None else $ObservedTemperature_C',
    'EffectivePulseProgram':'$PulseProgram or $ObservedPulseProgram',
    'EffectiveD1_s':'$D1_s if $D1_s not in (None, "") else $ObservedD1_s',
    'EffectiveDelayList_s':'$DelayList_s if $DelayList_s not in (None, "") else $ObservedDelayList_s',
    'EffectiveSpinLock_Hz':'$SpinLock_Hz if $SpinLock_Hz not in (None, "") else $ObservedSpinLock_Hz',
    'EffectiveSpinLockOffset_Hz':'$SpinLockOffset_Hz if $SpinLockOffset_Hz not in (None, "") else $ObservedSpinLockOffset_Hz',
}

def canonical_experiment_migration(client, approved_digest=None, dry_run=True):
    """Produce an evidence audit; apply only the exact reviewed, conflict-free plan.

    Source metadata may replace recorded facts only for a unique exact-bound
    acquisition. Never writes sample, assessments, links, titles or human notes.
    Existing source facts remain visible in the returned audit and caller backup.
    """
    import hashlib, re
    from collections import Counter, defaultdict
    from .schema import columns, canonical_experiment_type
    tables=client.records('_grist_Tables');meta=client.records('_grist_Tables_column');rows=client.records('Experiments')
    found=[r for r in tables if r['fields'].get('tableId')=='Experiments']
    if len(found)!=1:raise ValueError('Experiments table missing or ambiguous')
    table=found[0]['id'];records=[r for r in meta if r['fields'].get('parentId')==table]
    byname={r['fields'].get('colId'):r for r in records}
    if len(byname)!=len(records):raise ValueError('Duplicate Experiments column')
    present=[k for k in LEGACY_OBSERVED_FACTS if k in byname]
    if not present:
        expected={c['id']:c for c in columns()['Experiments']}
        required=set(LEGACY_OBSERVED_FACTS.values())|set(LEGACY_EFFECTIVE_FORMULAS)|{'Title','TitleSource'}
        if not required <= set(byname):raise ValueError('Incomplete canonical experiment schema')
        if any(not byname[k]['fields'].get('isFormula') or byname[k]['fields'].get('formula')!=expected[k]['formula'] for k in LEGACY_EFFECTIVE_FORMULAS):raise ValueError('Canonical experiment formula requires review')
        return {'digest':None,'audit':[],'actions':[],'already_migrated':True}
    if len(present)!=len(LEGACY_OBSERVED_FACTS):raise ValueError('Partial legacy observation schema requires review')
    expected={c['id']:c for c in columns()['Experiments']};actions=[]
    for name in ('Title','TitleSource'):
        if name not in byname:actions.append(['AddColumn','Experiments',name,{k:v for k,v in expected[name].items() if k!='id'}])
        elif byname[name]['fields'].get('type')!='Text' or byname[name]['fields'].get('isFormula'):raise ValueError('Unexpected '+name+' schema')
    for name,old in LEGACY_EFFECTIVE_FORMULAS.items():
        entry=byname.get(name);wanted=expected[name]['formula']
        if not entry or not entry['fields'].get('isFormula') or entry['fields'].get('formula') not in (old,wanted):
            raise ValueError('Unexpected formula for '+name)
    # Permit only the known auto-generated summary sum as an outside dependency.
    token=re.compile(r'\b('+'|'.join(map(re.escape,LEGACY_OBSERVED_FACTS))+r')\b')
    summary_ids={r['id'] for r in tables if r['fields'].get('summarySourceTable')==table}
    remove_summary=[]
    for entry in meta:
        f=entry['fields'];formula=f.get('formula') or ''
        if not token.search(formula):continue
        if f.get('parentId')==table and f.get('colId') in LEGACY_EFFECTIVE_FORMULAS:continue
        if (f.get('parentId') in summary_ids and f.get('isFormula') and f.get('colId')=='ObservedField_MHz'
            and formula=='SUM($group.ObservedField_MHz)'):
            table_name=next(r['fields']['tableId'] for r in tables if r['id']==f['parentId'])
            remove_summary.append(['RemoveColumn',table_name,'ObservedField_MHz']);continue
        raise ValueError('Unfamiliar observed-field formula dependency: '+str(f.get('colId')))
    def tup(f,prefix=''):
        source=f.get(prefix+'SourceID');folder=f.get(prefix+'Folder');exp=f.get(prefix+'Expno')
        if not isinstance(source,str) or not source or folder is None or not isinstance(folder,str) or exp is None:return None
        exp=str(exp)
        if not exp.isdigit():return None
        return source,folder,exp
    keys=Counter(r['fields'].get('AcquisitionKey') for r in rows if r['fields'].get('AcquisitionKey'))
    identities=Counter(tup(r['fields'],'Observed') for r in rows if tup(r['fields'],'Observed'))
    manual=defaultdict(list)
    for r in rows:
        identity=tup(r['fields'])
        if identity:manual[identity].append(r['id'])
    names=list(LEGACY_OBSERVED_FACTS.values()); audit=[];updates=[]
    if len({r['id'] for r in rows})!=len(rows):raise ValueError('Duplicate experiment row ID')
    recognized={'observed '+name+' differs from confirmed value' for name in names}
    for row in sorted(rows,key=lambda x:x['id']):
        f=row['fields'];errors=row.get('errors',{})
        if any(k in errors for k in (*LEGACY_OBSERVED_FACTS,*names,'AcquisitionKey','SourceID','Folder','Expno','ObservedSourceID','ObservedFolder','ObservedExpno','ReconciliationReason','MetadataDifference')):
            raise ValueError('Cannot read acquisition facts for row '+str(row['id']))
        observed=tup(f,'Observed');key=f.get('AcquisitionKey') or ''
        expected_key=(observed[0]+':'+(observed[1]+'/' if observed[1] else '')+observed[2]) if observed else None
        disagreements=any(f.get(k) not in (None,'') and str(f[k])!=observed[i] for i,k in enumerate(('SourceID','Folder','Expno'))) if observed else False
        coherent=bool(observed and key and key==expected_key and keys[key]==1 and identities[observed]==1
                      and all(rid==row['id'] for rid in manual[observed]) and not disagreements and not f.get('ReconciliationReason'))
        facts={k:f.get(k) for k in LEGACY_OBSERVED_FACTS};old={k:f.get(k) for k in names};change={}
        if not coherent and any(v not in (None,'') for v in facts.values()):
            raise ValueError('Observed facts require identity reconciliation before removal: row '+str(row['id']))
        if coherent:
            applied_sources=set()
            for src,dst in LEGACY_OBSERVED_FACTS.items():
                v=facts[src]
                if v is None or v=='':continue
                if isinstance(v,str) and v.strip().lower()=='unknown' and src!='ObservedAcquisitionStatus':continue
                if src=='ObservedType':
                    v=canonical_experiment_type(v)
                    if not v:continue
                if src=='ObservedAcquisitionStatus' and v not in ('Finished','Stopped','Unknown'):raise ValueError('Unexpected observed acquisition status at row '+str(row['id']))
                if src=='ObservedField_MHz':
                    if not isinstance(v,(int,float)) or isinstance(v,bool) or v<=0:raise ValueError('Invalid observed field at row '+str(row['id']))
                    v=round(v)
                if src=='ObservedAcquiredAt' and v==0:continue
                applied_sources.add(dst)
                if old[dst]!=v:change[dst]=v
            warning=f.get('MetadataDifference') or ''
            if warning and all(piece.strip() in recognized and piece.strip()[9:-29] in applied_sources for piece in warning.split(';')):change['MetadataDifference']=''
        audit.append({'id':row['id'],'identity':{'AcquisitionKey':key,'SourceID':f.get('SourceID'),'Folder':f.get('Folder'),'Expno':f.get('Expno'),
             'ObservedSourceID':f.get('ObservedSourceID'),'ObservedFolder':f.get('ObservedFolder'),'ObservedExpno':f.get('ObservedExpno'),
             'ReconciliationReason':f.get('ReconciliationReason'),'MetadataDifference':f.get('MetadataDifference')},
             'coherent':coherent,'recorded':old,'observed':facts,'changes':change})
        if change:updates.append(['UpdateRecord','Experiments',row['id'],change])
    actions.extend(updates)
    for name in LEGACY_EFFECTIVE_FORMULAS:
        entry=byname[name];wanted=expected[name]['formula']
        if entry['fields'].get('formula')!=wanted:actions.append(['UpdateRecord','_grist_Tables_column',entry['id'],{'formula':wanted}])
    for name in ('AcquisitionStatus','PulseProgram','D1_s','DelayList_s','SpinLock_Hz','SpinLockOffset_Hz'):
        entry=byname.get(name)
        if not entry:raise ValueError('Missing canonical field: '+name)
        target=expected[name];changes={}
        if entry['fields'].get('label')!=target.get('label',name):changes.update(label=target.get('label',name),untieColIdFromLabel=True)
        if name=='AcquisitionStatus' and entry['fields'].get('description'):changes['description']=''
        if changes:actions.append(['UpdateRecord','_grist_Tables_column',entry['id'],changes])
    actions.extend(remove_summary)
    actions.extend(['RemoveColumn','Experiments',name] for name in LEGACY_OBSERVED_FACTS)
    payload={'audit':audit,'actions':actions}
    digest=hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(',',':'),ensure_ascii=False,allow_nan=False).encode()).hexdigest()
    if not dry_run:
        if not approved_digest or approved_digest!=digest:raise ValueError('Canonical migration audit changed or was not approved')
        client.apply(actions)
    return {'digest':digest,'audit':audit,'actions':actions,'already_migrated':False}


LEGACY_TEMPERATURE_FORMULA = '$Temperature_C if $Temperature_C is None or abs($Temperature_C - round($Temperature_C)) > 0.1 else round($Temperature_C)'


def upgrade_temperature_grouping(client, *, dry_run=True, tables=None, metadata=None):
    """Plan or apply the one known display-formula change; never rewrite temperatures."""
    from .schema import columns
    if tables is None:
        found=[r for r in client.records('_grist_Tables') if r['fields'].get('tableId')=='Experiments']
        if len(found)!=1:raise ValueError('Experiments table missing or ambiguous')
        table=found[0]['id']
    else:
        table=tables.get('Experiments')
        if not table:raise ValueError('Experiments table missing or ambiguous')
    if metadata is None:metadata=client.records('_grist_Tables_column')
    matches=[r for r in metadata if r['fields'].get('parentId')==table and r['fields'].get('colId')=='EffectiveTemperature_C']
    if len(matches)!=1:raise ValueError('Effective temperature column missing or ambiguous')
    row=matches[0];fields=row['fields'];wanted=next(c['formula'] for c in columns()['Experiments'] if c['id']=='EffectiveTemperature_C')
    if fields.get('type')!='Numeric' or not fields.get('isFormula') or fields.get('formula') not in (LEGACY_TEMPERATURE_FORMULA,wanted):
        raise ValueError('Effective temperature formula requires review')
    changes={}
    if fields['formula']!=wanted:changes['formula']=wanted
    if fields.get('label') in ('Temperature (°C)','Temperature group (°C)'):
        if fields.get('label')!='Temperature group (°C)':changes['label']='Temperature group (°C)'
        if not fields.get('untieColIdFromLabel'):changes['untieColIdFromLabel']=True
    actions=[['UpdateRecord','_grist_Tables_column',row['id'],changes]] if changes else []
    for summary in metadata:
        f=summary['fields']
        if f.get('summarySourceCol')==row['id'] and f.get('label') in ('Temperature (°C)','Temperature group (°C)'):
            update={}
            if f.get('label')!='Temperature group (°C)':update['label']='Temperature group (°C)'
            if not f.get('untieColIdFromLabel'):update['untieColIdFromLabel']=True
            if update:actions.append(['UpdateRecord','_grist_Tables_column',summary['id'],update])
    if actions and not dry_run:client.apply(actions)
    return actions


LEGACY_SOURCE_LABEL_FORMULA = '($Folder + "/" + $Expno) if $Folder and $Expno else ($Folder or $Expno or "")'


def upgrade_source_label(client, *, dry_run=True, tables=None, metadata=None):
    """Update only the shipped source label formula and its default label."""
    from .schema import SOURCE_LABEL_FORMULA
    if tables is None:
        found=[r for r in client.records('_grist_Tables') if r['fields'].get('tableId')=='Experiments']
        if len(found)!=1:raise ValueError('Experiments table missing or ambiguous')
        table=found[0]['id']
    else:
        table=tables.get('Experiments')
        if not table:raise ValueError('Experiments table missing or ambiguous')
    if metadata is None:metadata=client.records('_grist_Tables_column')
    found=[r for r in metadata if r['fields'].get('parentId')==table and r['fields'].get('colId')=='SourceLabel']
    if len(found)!=1:raise ValueError('Source label missing or ambiguous')
    row=found[0];f=row['fields']
    if f.get('type')!='Text' or not f.get('isFormula') or f.get('formula') not in (LEGACY_SOURCE_LABEL_FORMULA,SOURCE_LABEL_FORMULA):
        raise ValueError('Source label formula requires review')
    changes={};actions=[]
    if f['formula']!=SOURCE_LABEL_FORMULA:changes['formula']=SOURCE_LABEL_FORMULA
    if f.get('label') in ('Acquisition source','EXPNO · Folder'):
        # Grist may rename the column ID if label and untie are changed in
        # the same action; untie first in a separate ordered action.
        if not f.get('untieColIdFromLabel'):
            actions.append(['UpdateRecord','_grist_Tables_column',row['id'],{'untieColIdFromLabel':True}])
        if f.get('label')!='EXPNO · Folder':changes['label']='EXPNO · Folder'
    if changes:actions.append(['UpdateRecord','_grist_Tables_column',row['id'],changes])
    if actions and not dry_run:client.apply(actions)
    return actions


def upgrade_attention_identity(client, *, dry_run=True, tables=None, metadata=None):
    """Add only the two derived identity displays used by Needs attention."""
    from .schema import columns
    if tables is None:
        found=[r for r in client.records('_grist_Tables') if r['fields'].get('tableId')=='Experiments']
        if len(found)!=1:raise ValueError('Experiments table missing or ambiguous')
        table=found[0]['id']
    else:
        table=tables.get('Experiments')
        if not table:raise ValueError('Experiments table missing or ambiguous')
    if metadata is None:metadata=client.records('_grist_Tables_column')
    expected={c['id']:c for c in columns()['Experiments']};actions=[]
    for name in ('AttentionFolder','AttentionExpno'):
        found=[r for r in metadata if r['fields'].get('parentId')==table and r['fields'].get('colId')==name]
        if len(found)>1:raise ValueError('Attention identity column ambiguous: '+name)
        want=expected[name]
        if not found:actions.append(['AddColumn','Experiments',name,{k:v for k,v in want.items() if k!='id'}])
        else:
            f=found[0]['fields']
            if f.get('type')!='Text' or not f.get('isFormula') or f.get('formula')!=want['formula']:
                raise ValueError('Attention identity formula requires review: '+name)
    if actions and not dry_run:client.apply(actions)
    return actions


def create_schema(client):
    """Create tables in an empty document, then configure snapshot triggers.

    client.apply(actions) and client.records(table) are required; intentionally no
    credentials or network handling. Fails if tables already exist.
    """
    existing={r['fields']['tableId']:r['id'] for r in client.records('_grist_Tables')}
    actions=schema_actions(); wanted={a[1] for a in actions}
    metadata=client.records('_grist_Tables_column')
    if 'Experiments' in existing and any(r['fields'].get('parentId')==existing['Experiments'] and r['fields'].get('colId') in LEGACY_OBSERVED_FACTS for r in metadata):
        raise ValueError('Review and apply the canonical experiment migration first')
    if any(name in existing for name in ('Series','Samples')):
        upgrade_labelling_choices(client,existing,metadata)
        metadata=client.records('_grist_Tables_column')
        upgrade_protein_choices(client,existing,metadata)
        metadata=client.records('_grist_Tables_column')
    if 'Samples' in existing:
        upgrade_sample_status(client,existing,metadata)
        metadata=client.records('_grist_Tables_column')
    if 'Samples' in existing and 'StorageGrid' in existing:
        upgrade_single_rack(client,existing,metadata)
        metadata=client.records('_grist_Tables_column')
    if 'Experiments' in existing:
        upgrade_acquisition_workflow_status(client,dry_run=False)
        upgrade_completion_schema(client)
        upgrade_acquisition_schema(client)
        metadata=client.records('_grist_Tables_column')
        upgrade_temperature_grouping(client,dry_run=False,tables=existing,metadata=metadata)
        if any(r['fields'].get('parentId')==existing['Experiments'] and r['fields'].get('colId')=='SourceLabel' for r in metadata):
            upgrade_source_label(client,dry_run=False,tables=existing,metadata=metadata)
        upgrade_electrostatics_schema(client)
        metadata=client.records('_grist_Tables_column')
        upgrade_attention_identity(client,dry_run=False,tables=existing,metadata=metadata)
        upgrade_sample_tag_schema(client)
        upgrade_assessment_schema(client)
        metadata=client.records('_grist_Tables_column')
    for action in actions:
        name=action[1]
        if name not in existing: continue
        actual={r['fields']['colId']:r['fields'] for r in metadata if r['fields']['parentId']==existing[name]}
        for expected in action[2]:
            found=actual.get(expected['id'])
            if not found or found.get('type') != expected['type'] or found.get('formula','') != expected['formula'] or bool(found.get('isFormula')) != expected['isFormula']:
                raise ValueError('Existing register schema differs at '+name+'.'+expected['id'])
    missing=[a for a in actions if a[1] not in existing]
    if missing:
        client.apply(missing)
        if any(a[1]=='Experiments' for a in missing): upgrade_acquisition_schema(client)
    tables={r['fields']['tableId']:r['id'] for r in client.records('_grist_Tables')}
    inverse={v:k for k,v in tables.items()}
    cols={(inverse[r['fields']['parentId']],r['fields']['colId']):r['id']
          for r in client.records('_grist_Tables_column') if r['fields']['parentId'] in inverse}
    client.apply(trigger_actions(tables,cols))
    apply_experiment_date_format(client)
    # AddTable creates a primary page for every table. Hide only the pages just
    # created here; the underlying raw views and tables remain available.
    created=[r for r in client.records('_grist_Tables') if r['fields']['tableId'] in wanted]
    primary={r['fields'].get('primaryViewId') for r in created} - {0,None}
    cleanup=[['UpdateRecord','_grist_Views',r['fields']['primaryViewId'],
              {'name':'Register raw '+r['fields']['tableId']}]
             for r in created if r['fields'].get('primaryViewId')]
    for meta in ('_grist_Pages','_grist_TabBar'):
        cleanup.extend(['RemoveRecord',meta,r['id']] for r in client.records(meta)
                       if r['fields'].get('viewRef') in primary)
    if cleanup: client.apply(cleanup)

def upgrade_single_rack(client,tables=None,metadata=None):
    """Migrate only known rack display formulas; retain all sample/box values."""
    from .samples_ui import STORAGE_FORMULA,LEGACY_STORAGE_FORMULA
    if tables is None: tables={r['fields']['tableId']:r['id'] for r in client.records('_grist_Tables')}
    if 'Samples' not in tables or 'StorageGrid' not in tables: raise ValueError('Sample storage tables missing')
    if metadata is None: metadata=client.records('_grist_Tables_column')
    bycol={(r['fields']['parentId'],r['fields']['colId']):r for r in metadata}
    rows=client.records('StorageGrid'); letters=[r['fields'].get('Row') for r in rows]
    if any(not isinstance(row,str) or len(row)!=1 or row not in 'ABCDEFGH' for row in letters) or len(set(letters))!=len(letters):
        raise ValueError('Storage grid rows are invalid or duplicated; review before continuing')
    actions=[]
    display=bycol.get((tables['Samples'],'StorageDisplay'))
    if display:
        f=display['fields']
        if not f.get('isFormula') or f.get('formula') not in (STORAGE_FORMULA,LEGACY_STORAGE_FORMULA):
            raise ValueError('Existing storage display differs; review before continuing')
        if f['formula']!=STORAGE_FORMULA:
            actions.append(['UpdateRecord','_grist_Tables_column',display['id'],{'formula':STORAGE_FORMULA}])
    for number in range(1,13):
        column=bycol.get((tables['StorageGrid'],'Col'+str(number)))
        if not column: raise ValueError('Missing rack grid column '+str(number))
        old='", ".join(s.SampleID for s in Samples.lookupRecords(StorageBox=$StorageBox, StorageRow=$Row, StorageColumn='+str(number)+'))'
        new='", ".join(s.SampleID for s in Samples.lookupRecords(StorageRow=$Row, StorageColumn='+str(number)+'))'
        f=column['fields']
        if not f.get('isFormula') or f.get('formula') not in (old,new):
            raise ValueError('Existing rack grid formula differs: '+str(number))
        if f['formula']!=new: actions.append(['UpdateRecord','_grist_Tables_column',column['id'],{'formula':new}])
    actions.extend(['AddRecord','StorageGrid',None,{'Row':row}] for row in 'ABCDEFGH' if row not in letters)
    if actions: client.apply(actions)
    return actions


LEGACY_WORKFLOW_EFFECTIVE = '"Finished" if $AcquisitionStatus in ("Done", "Acquired") else ($AcquisitionStatus or "Unknown")'
LEGACY_WORKFLOW_RECORDS = '"%s / %s" % (len($group), len([r for r in $group if r.EffectiveAcquisitionStatus == "Finished"]))'

def upgrade_acquisition_workflow_status(client, *, dry_run=True):
    """Update known status presentation without rewriting stored human values."""
    import json
    from .schema import columns, ATTENTION_FORMULA, ACQUISITION_STATUSES
    tables=client.records('_grist_Tables'); metadata=client.records('_grist_Tables_column')
    ids=[r['id'] for r in tables if r['fields'].get('tableId')=='Experiments']
    if len(ids)!=1:raise ValueError('Experiments table missing or ambiguous')
    table=ids[0]; cols={r['fields'].get('colId'):r for r in metadata if r['fields'].get('parentId')==table}
    if not {'AcquisitionStatus','EffectiveAcquisitionStatus','AttentionReason'} <= cols.keys():raise ValueError('Acquisition status columns missing')
    wanted={c['id']:c for c in columns()['Experiments']};actions=[]
    effective=cols['EffectiveAcquisitionStatus']; f=effective['fields']
    if not f.get('isFormula') or f.get('formula') not in (PRE_DUMMY_EFFECTIVE_FORMULA,LEGACY_WORKFLOW_EFFECTIVE,wanted['EffectiveAcquisitionStatus']['formula']):raise ValueError('Status formula differs; review')
    if f.get('formula')!=wanted['EffectiveAcquisitionStatus']['formula']:actions.append(['UpdateRecord','_grist_Tables_column',effective['id'],{'formula':wanted['EffectiveAcquisitionStatus']['formula']}])
    old_attention=PRE_DUMMY_ATTENTION_FORMULA.replace('if $EffectiveAcquisitionStatus == "Aborted": issues.append("Acquisition aborted")\nif $EffectiveAcquisitionStatus == "Unknown": issues.append("Acquisition status unknown")', 'if $EffectiveAcquisitionStatus == "Stopped": issues.append("Acquisition stopped")').replace('$EffectiveAcquisitionStatus == "Completed"','$EffectiveAcquisitionStatus == "Finished"')
    attention=cols['AttentionReason']; f=attention['fields']
    if not f.get('isFormula') or f.get('formula') not in (old_attention,PRE_ASSIGNED_SAMPLE_ATTENTION_FORMULA, PRE_LIVE_SAMPLE_ATTENTION_FORMULA,PRE_STUDY_ASSESSMENT_ATTENTION_FORMULA,PRE_DUMMY_ATTENTION_FORMULA,PRE_CALIBRATION_ATTENTION_FORMULA,PRE_CATEGORY_ATTENTION_FORMULA,PRE_REGISTRATION_REMOVAL_ATTENTION_FORMULA,ATTENTION_FORMULA):raise ValueError('Attention formula differs; review')
    if f.get('formula')!=ATTENTION_FORMULA:actions.append(['UpdateRecord','_grist_Tables_column',attention['id'],{'formula':ATTENTION_FORMULA}])
    status=cols['AcquisitionStatus']; f=status['fields']
    if f.get('isFormula') or f.get('formula','') or f.get('type') not in ('Text','Choice'):raise ValueError('Stored status column differs; review')
    options=json.loads(f.get('widgetOptions') or '{}'); oldchoices=options.get('choices',[])
    if not isinstance(oldchoices,list) or any(not isinstance(c,str) for c in oldchoices):raise ValueError('Status choices require review')
    choices=list(ACQUISITION_STATUSES)+[c for c in oldchoices if c not in ACQUISITION_STATUSES]
    if f.get('type')!='Choice' or oldchoices!=choices:
        options['choices']=choices
        actions.append(['ModifyColumn','Experiments','AcquisitionStatus',{'type':'Choice','widgetOptions':json.dumps(options)}])
    table_names={r['id']:r['fields'].get('tableId','') for r in tables}
    formulas=electrostatics_formulas()
    old_verified=formulas['VerifiedAcquisitions'].replace('r.AcquisitionStatus in ("Done", "Acquired", "Finished", "Completed")','r.AcquisitionStatus in ("Done", "Acquired")')
    for row in metadata:
        f=row['fields']; name=table_names.get(f.get('parentId'),'');column=f.get('colId')
        if not name.startswith('Experiments_summary_') or column not in ('RecordsDone','VerifiedAcquisitions'):continue
        target=formulas[column]; old=LEGACY_WORKFLOW_RECORDS if column=='RecordsDone' else old_verified
        if not f.get('isFormula') or f.get('formula') not in (old,target):raise ValueError('Summary status formula differs; review')
        if f.get('formula')!=target:actions.append(['UpdateRecord','_grist_Tables_column',row['id'],{'formula':target}])
        if column=='RecordsDone' and f.get('label')=='Records / finished':
            if not f.get('untieColIdFromLabel'):actions.append(['UpdateRecord','_grist_Tables_column',row['id'],{'untieColIdFromLabel':True}])
            actions.append(['UpdateRecord','_grist_Tables_column',row['id'],{'label':'Records / completed'}])
    views={r['id']:r['fields'].get('name') for r in client.records('_grist_Views')}
    grids={(views.get(r['fields'].get('parentId')),r['fields'].get('title')):r['id'] for r in client.records('_grist_Views_section') if r['fields'].get('parentKey')=='record'}
    targets=[grids[x] for x in (('Experiments','Experiments'),('Samples','Sample experiments'),('Needs attention','Needs attention'),('Electrostatics','Individual acquisitions')) if x in grids]
    section_fields=client.records('_grist_Views_section_field'); raw=cols['AcquisitionStatus']['id']; display=effective['id']
    for section in targets:
        rows=[r for r in section_fields if r['fields'].get('parentId')==section]
        old=[r for r in rows if r['fields'].get('colRef')==raw]; new=[r for r in rows if r['fields'].get('colRef')==display]
        if len(old)>1 or len(new)>1:raise ValueError('Ambiguous status field in grid')
        if old and new and new[0]['fields'].get('width')==0:
            actions.append(['RemoveRecord','_grist_Views_section_field',new[0]['id']])
            actions.append(['UpdateRecord','_grist_Views_section_field',old[0]['id'],{'colRef':display}])
        elif old and new:actions.append(['RemoveRecord','_grist_Views_section_field',old[0]['id']])
        elif old:actions.append(['UpdateRecord','_grist_Views_section_field',old[0]['id'],{'colRef':display}])
    if actions and not dry_run:client.apply(actions)
    return actions


def normalize_acquisition_status(client, *, dry_run=True, batch_size=200):
    """Convert only exact synonymous legacy labels, without touching other fields.

    A single-writer window is required between the final read and transaction.
    On a transport error inspect current records; never blindly retry.
    """
    if not isinstance(batch_size,int) or isinstance(batch_size,bool) or not 1<=batch_size<=500:raise ValueError('Invalid status batch size')
    mapping={'Finished':'Completed','Done':'Completed','Acquired':'Completed','Stopped':'Aborted'}
    rows=client.records('Experiments');changes=[];seen=set()
    for row in rows:
        rid=row['id']
        if rid in seen:raise ValueError('Duplicate experiment ID')
        seen.add(rid);old=row['fields'].get('AcquisitionStatus')
        if old in mapping:changes.append((rid,old,mapping[old]))
    actions=[['BulkUpdateRecord','Experiments',[rid for rid,_,_ in changes[i:i+batch_size]],{'AcquisitionStatus':[new for _,_,new in changes[i:i+batch_size]]}] for i in range(0,len(changes),batch_size)]
    if actions and not dry_run:
        current={r['id']:r['fields'].get('AcquisitionStatus') for r in client.records('Experiments')}
        if any(current.get(rid)!=old for rid,old,_ in changes):raise ValueError('Status changed since preflight')
        client.apply(actions)
    return {'actions':actions,'changed':len(changes),'counts':{old:sum(1 for _,value,_ in changes if value==old) for old in mapping}}


LEGACY_COMPLETION_FORMULAS = {'AttentionReason': 'issues = []\nif $MigrationIssue: issues.append($MigrationIssue)\nif $ReconciliationReason: issues.append($ReconciliationReason)\nif $MetadataDifference: issues.append($MetadataDifference)\nif $ScanError: issues.append($ScanError)\nif not $RegisteredAt: issues.append("Registration incomplete")\nif not $Sample: issues.append("Sample not assigned")\nif $Assessment and $Assessment not in ("Use", "Limited use", "Do not use"): issues.append("Invalid assessment")\nif $Assessment in ("Limited use", "Do not use") and not ($AssessmentReason or "").strip(): issues.append("Assessment reason missing")\nif not $Assessment and ($RawPresent or $AcquisitionStatus == "Acquired") and $EffectiveAcquiredAt and (NOW() - $EffectiveAcquiredAt).total_seconds() > 86400: issues.append("Acquired over 24 hours ago; unassessed")\nreturn "; ".join(issues)', 'RecordsDone': '"%s / %s" % (len($group), len([r for r in $group if r.AcquisitionStatus in ("Done", "Acquired")]))'}

def upgrade_completion_schema(client):
    """Add status observations and update only the known shipped formulas."""
    from .schema import columns
    tables=client.records('_grist_Tables')
    table=next((r['id'] for r in tables if r['fields']['tableId']=='Experiments'),None)
    if table is None: raise ValueError('Experiments table is missing')
    metadata=client.records('_grist_Tables_column')
    existing={r['fields']['colId']:r for r in metadata if r['fields'].get('parentId')==table}
    expected={c['id']:c for c in columns()['Experiments']};actions=[]
    if 'ObservedAcquisitionStatus' not in existing:
        status=existing.get('EffectiveAcquisitionStatus')
        if not status or status['fields'].get('formula') not in (PRE_DUMMY_EFFECTIVE_FORMULA,expected['EffectiveAcquisitionStatus']['formula']):
            raise ValueError('Canonical completion schema differs')
        return []
    raise ValueError('Legacy observed acquisition status requires reviewed canonical migration')


def upgrade_sample_tag_schema(client):
    """Add tag facts and an atomic fill-only Sample trigger, preserving data."""
    from .schema import columns, ATTENTION_FORMULA, LEGACY_QC_ATTENTION_FORMULA, SAMPLE_TAG_TRIGGER
    tables=client.records('_grist_Tables')
    ids=[r['id'] for r in tables if r['fields'].get('tableId')=='Experiments']
    if len(ids)!=1:raise ValueError('Experiments table missing or ambiguous')
    table=ids[0]
    def metadata():
        rows=[r for r in client.records('_grist_Tables_column') if r['fields'].get('parentId')==table]
        found={}
        for row in rows:
            name=row['fields'].get('colId')
            if name in found:raise ValueError('Duplicate Experiments column: '+str(name))
            found[name]=row
        return found
    cols=metadata();wanted={c['id']:c for c in columns()['Experiments']};actions=[]
    for name in ('ObservedSampleID','ObservedSampleTagEvidence','ObservedSampleTagIssue','SampleTagIssue'):
        spec=wanted[name];existing=cols.get(name)
        if not existing:actions.append(['AddColumn','Experiments',name,{k:v for k,v in spec.items() if k!='id'}])
        elif (existing['fields'].get('type')!=spec['type'] or bool(existing['fields'].get('isFormula'))!=spec['isFormula'] or existing['fields'].get('formula','')!=spec['formula']):
            raise ValueError('Existing sample tag column differs at '+name)
    sample=cols.get('Sample');attention=cols.get('AttentionReason')
    if not sample or not attention:raise ValueError('Missing Sample or AttentionReason column')
    samplefields=sample['fields'];formula=samplefields.get('formula','')
    if samplefields.get('type')!='Ref:Samples' or samplefields.get('isFormula') or formula not in ('',SAMPLE_TAG_TRIGGER):
        raise ValueError('Existing Sample formula requires review')
    if not formula and samplefields.get('recalcDeps') not in (None,[],['L']):
        raise ValueError('Existing Sample trigger dependencies require review')
    current_target=ATTENTION_FORMULA if 'AssessmentExplanation' in cols else LEGACY_QC_ATTENTION_FORMULA.replace('if not $RegisteredAt: issues.append("Registration incomplete")\n', '')
    legacy_without_tag=LEGACY_QC_ATTENTION_FORMULA.replace('if $SampleTagIssue: issues.append($SampleTagIssue)\n','')
    old=current_target.replace('if $SampleTagIssue: issues.append($SampleTagIssue)\n','')
    current=attention['fields'].get('formula','')
    if current not in (old,current_target,LEGACY_QC_ATTENTION_FORMULA,legacy_without_tag,PRE_ASSIGNED_SAMPLE_ATTENTION_FORMULA, PRE_LIVE_SAMPLE_ATTENTION_FORMULA,PRE_STUDY_ASSESSMENT_ATTENTION_FORMULA,PRE_CALIBRATION_ATTENTION_FORMULA,PRE_DUMMY_ATTENTION_FORMULA,PRE_CATEGORY_ATTENTION_FORMULA,PRE_REGISTRATION_REMOVAL_ATTENTION_FORMULA,PRE_DUMMY_ATTENTION_FORMULA.replace('if $SampleTagIssue: issues.append($SampleTagIssue)\n','')):raise ValueError('AttentionReason formula requires review')
    if current!=current_target:actions.append(['UpdateRecord','_grist_Tables_column',attention['id'],{'formula':current_target}])
    if actions:client.apply(actions)
    cols=metadata()
    deps=('ObservedSampleID','ObservedSampleTagIssue','AcquisitionKey','ObservedSourceID','ObservedFolder','ObservedExpno','SourceID','Folder','Expno','ReconciliationReason')
    expected=['L',*[cols[name]['id'] for name in deps]]
    fields=cols['Sample']['fields'];changes={}
    if fields.get('formula','')!=SAMPLE_TAG_TRIGGER:changes['formula']=SAMPLE_TAG_TRIGGER
    if fields.get('recalcWhen')!=0:changes['recalcWhen']=0
    if fields.get('recalcDeps')!=expected:changes['recalcDeps']=expected
    if changes:
        action=['UpdateRecord','_grist_Tables_column',cols['Sample']['id'],changes]
        client.apply([action]);actions.append(action)
    return actions


def upgrade_assessment_schema(client, reviewed=None):
    """Atomically install QC criteria, preserving every reviewed legacy decision.

    reviewed maps row IDs to {expected_assessment, expected_reason, checks}; every
    nonblank legacy decision requires an explicit evidence-reviewed checklist.
    Refuse unknown formulas/choices and unexpected existing checklist values.
    """
    from .schema import (columns, QC_CRITERIA, QC_DEFAULT, assessment_from_values,
                         ASSESSMENT_FORMULA, ASSESSMENT_EXPLANATION_FORMULA, PRE_COMPACT_ASSESSMENT_DISPLAY_FORMULA, LEGACY_ASSESSMENT_FORMULA,
                         PRE_ANALYSIS_ASSESSMENT_FORMULA, PRE_ANALYSIS_ASSESSMENT_EXPLANATION_FORMULA, PRE_ANALYSIS_ASSESSMENT_DISPLAY_FORMULA,
                         LEGACY_ASSESSMENT_EXPLANATION_FORMULA, PRE_TOLERANCE_ASSESSMENT_FORMULA,
                         PRE_TOLERANCE_ASSESSMENT_EXPLANATION_FORMULA, PRE_PARTIAL_ASSESSMENT_EXPLANATION_FORMULA, PRE_PARTIAL_ATTENTION_FORMULA, ATTENTION_FORMULA, LEGACY_QC_ATTENTION_FORMULA)
    reviewed=reviewed or {}
    tables=[r['id'] for r in client.records('_grist_Tables') if r['fields'].get('tableId')=='Experiments']
    if len(tables)!=1:raise ValueError('Experiments table missing or ambiguous')
    metadata=[r for r in client.records('_grist_Tables_column') if r['fields'].get('parentId')==tables[0]]
    cols={r['fields']['colId']:r for r in metadata}
    if len(cols)!=len(metadata):raise ValueError('Duplicate Experiments column')
    expected={c['id']:c for c in columns()['Experiments']}
    assessment=cols.get('Assessment');attention=cols.get('AttentionReason')
    if not assessment or not attention:raise ValueError('Missing Assessment or AttentionReason')
    af=assessment['fields'];old_mode=not af.get('isFormula') and af.get('formula','')==''
    if af.get('type')!='Choice' or (not old_mode and (not af.get('isFormula') or af.get('formula') not in (ASSESSMENT_FORMULA, PRE_ANALYSIS_ASSESSMENT_FORMULA, PRE_TOLERANCE_ASSESSMENT_FORMULA, LEGACY_ASSESSMENT_FORMULA))):
        raise ValueError('Assessment schema requires review')
    current_attention=attention['fields'].get('formula','')
    if current_attention not in (LEGACY_QC_ATTENTION_FORMULA,PRE_ASSIGNED_SAMPLE_ATTENTION_FORMULA, PRE_LIVE_SAMPLE_ATTENTION_FORMULA,PRE_STUDY_ASSESSMENT_ATTENTION_FORMULA,PRE_DUMMY_ATTENTION_FORMULA,PRE_CALIBRATION_ATTENTION_FORMULA,PRE_PARTIAL_ATTENTION_FORMULA,PRE_CATEGORY_ATTENTION_FORMULA,PRE_REGISTRATION_REMOVAL_ATTENTION_FORMULA,ATTENTION_FORMULA):raise ValueError('AttentionReason formula requires review')
    actions=[]
    for name in (*QC_CRITERIA,'AssessmentExplanation','ImportedAssessment','AssessmentDisplay'):
        spec=expected[name];existing=cols.get(name)
        if existing:
            f=existing['fields']
            allowed=(spec['formula'], PRE_ANALYSIS_ASSESSMENT_EXPLANATION_FORMULA, PRE_PARTIAL_ASSESSMENT_EXPLANATION_FORMULA, PRE_TOLERANCE_ASSESSMENT_EXPLANATION_FORMULA, LEGACY_ASSESSMENT_EXPLANATION_FORMULA) if name=='AssessmentExplanation' else (spec['formula'],)
            if name == 'AssessmentDisplay':allowed = (spec['formula'], PRE_ANALYSIS_ASSESSMENT_DISPLAY_FORMULA, PRE_COMPACT_ASSESSMENT_DISPLAY_FORMULA)
            if f.get('type')!=spec['type'] or bool(f.get('isFormula'))!=spec['isFormula'] or f.get('formula','') not in allowed:
                raise ValueError('Existing QC schema differs at '+name)
            if name == 'QCAnalysis':
                import json
                try:
                    options = json.loads(f.get('widgetOptions') or '{}')
                except (ValueError, TypeError):
                    raise ValueError('QCAnalysis choices require review')
                target_choices = [QC_DEFAULT, *QC_CRITERIA[name][1]]
                if not isinstance(options, dict) or options.get('choices') not in ([QC_DEFAULT, 'Reliable', 'Limited', 'Unreliable'], target_choices):
                    raise ValueError('QCAnalysis choices require review')
                if options['choices'] != target_choices:
                    options['choices'] = target_choices
                    actions.append(['UpdateRecord','_grist_Tables_column',existing['id'],{'widgetOptions':json.dumps(options)}])
            if name in QC_CRITERIA or name == 'AssessmentDisplay':
                presentation = {}
                if not f.get('untieColIdFromLabel'):
                    presentation['untieColIdFromLabel'] = True
                if f.get('label') != spec['label']:
                    presentation.update(label=spec['label'], untieColIdFromLabel=True)
                if presentation:
                    actions.append(['UpdateRecord','_grist_Tables_column',existing['id'],presentation])
            if name == 'AssessmentDisplay' and f.get('formula') in (PRE_ANALYSIS_ASSESSMENT_DISPLAY_FORMULA, PRE_COMPACT_ASSESSMENT_DISPLAY_FORMULA):
                actions.append(['UpdateRecord','_grist_Tables_column',existing['id'],{'formula':spec['formula']}])
            if name=='AssessmentExplanation' and f.get('formula') in (PRE_ANALYSIS_ASSESSMENT_EXPLANATION_FORMULA, PRE_PARTIAL_ASSESSMENT_EXPLANATION_FORMULA, PRE_TOLERANCE_ASSESSMENT_EXPLANATION_FORMULA, LEGACY_ASSESSMENT_EXPLANATION_FORMULA):
                actions.append(['UpdateRecord','_grist_Tables_column',existing['id'],{'formula':ASSESSMENT_EXPLANATION_FORMULA}])
        else:
            info={k:v for k,v in spec.items() if k!='id'}
            if info.get('recalcDeps')==['L']:info['recalcDeps']=[['L']]
            actions.append(['AddColumn','Experiments',name,info])
    if old_mode:
        records=client.records('Experiments'); seen=set(); updates={name:[] for name in QC_CRITERIA}; rowids=[]
        for row in records:
            rid=row['id'];fields=row['fields'];old=fields.get('Assessment') or '';note=fields.get('AssessmentReason') or ''
            if rid in seen:raise ValueError('Duplicate experiment row ID')
            seen.add(rid);mapping=reviewed.get(rid)
            if old and not mapping:raise ValueError('Existing assessment needs reviewed mapping: '+str(rid))
            if mapping and (mapping.get('expected_assessment')!=old or mapping.get('expected_reason')!=note):
                raise ValueError('Assessment changed since review: '+str(rid))
            checks=mapping.get('checks',{}) if mapping else {}
            if set(checks)-set(QC_CRITERIA):raise ValueError('Unknown checklist field: '+str(rid))
            values={}
            for name,(_,choices) in QC_CRITERIA.items():
                existing=fields.get(name) or QC_DEFAULT;choice=checks.get(name,existing)
                if choice not in (QC_DEFAULT,*choices):raise ValueError('Invalid checklist choice at '+str(rid)+': '+name)
                if existing not in (QC_DEFAULT,choice):raise ValueError('Existing checklist value requires review: '+str(rid)+': '+name)
                values[name]=choice;updates[name].append(choice)
            if assessment_from_values(values)!=old:raise ValueError('Checklist changes existing assessment: '+str(rid))
            rowids.append(rid)
        if set(reviewed)-seen:raise ValueError('Reviewed mapping includes missing row')
        if rowids:actions.append(['BulkUpdateRecord','Experiments',rowids,updates])
        actions.append(['ModifyColumn','Experiments','Assessment',{'isFormula':True,'formula':ASSESSMENT_FORMULA}])
    else:
        if any(name not in cols for name in (*QC_CRITERIA,'AssessmentExplanation')):
            raise ValueError('Computed Assessment exists with incomplete QC schema; review before continuing')
        if af.get('formula') in (PRE_ANALYSIS_ASSESSMENT_FORMULA, PRE_TOLERANCE_ASSESSMENT_FORMULA, LEGACY_ASSESSMENT_FORMULA):
            actions.append(['UpdateRecord','_grist_Tables_column',assessment['id'],{'formula':ASSESSMENT_FORMULA}])
    if current_attention!=ATTENTION_FORMULA:
        actions.append(['UpdateRecord','_grist_Tables_column',attention['id'],{'formula':ATTENTION_FORMULA}])
    if actions:client.apply(actions)
    return actions


def apply_sample_tag_ui(client):
    """Expose concise tag and source on existing experiment detail cards."""
    tables=[r['id'] for r in client.records('_grist_Tables') if r['fields'].get('tableId')=='Experiments']
    if len(tables)!=1:raise ValueError('Experiments table missing or ambiguous')
    cols={r['fields']['colId']:r['id'] for r in client.records('_grist_Tables_column') if r['fields'].get('parentId')==tables[0]}
    for name in ('ObservedSampleID','ObservedSampleTagEvidence'):
        if name not in cols:raise ValueError('Missing '+name)
    views={r['id']:r['fields'].get('name') for r in client.records('_grist_Views')}
    targets={('Experiments','Experiment details'),('Needs attention','Attention details'),('Electrostatics','Acquisition details')}
    sections=client.records('_grist_Views_section');fields=client.records('_grist_Views_section_field');actions=[]
    for section in sections:
        f=section['fields']
        if (views.get(f.get('parentId')),f.get('title')) not in targets:continue
        if f.get('tableRef')!=tables[0]:raise ValueError('Unexpected tag detail card table')
        current=[r for r in fields if r['fields'].get('parentId')==section['id']]
        refs={r['fields'].get('colRef') for r in current};position=max([r['fields'].get('parentPos') or 0 for r in current]+[0])
        for name in ('ObservedSampleID','ObservedSampleTagEvidence'):
            if cols[name] not in refs:
                position+=1;actions.append(['AddRecord','_grist_Views_section_field',None,{'parentId':section['id'],'colRef':cols[name],'parentPos':position,'width':150}])
    if actions:client.apply(actions)
    return actions


def upgrade_labelling_choices(client, tables=None, metadata=None):
    """Offer standard labels while preserving each stored text and custom choice.

    Only existing Series and Samples Isotope columns are changed. The Text to
    Choice conversion keeps Grist's string-valued cells; experiment snapshots
    and all formula/trigger settings are left alone.
    """
    if tables is None:
        tables={r['fields']['tableId']:r['id'] for r in client.records('_grist_Tables')}
    if metadata is None: metadata=client.records('_grist_Tables_column')
    descriptions=[]; choices=list(LABELLING_CHOICES)
    for table in ('Series','Samples'):
        if table not in tables: continue
        matches=[r for r in metadata if r['fields'].get('parentId')==tables[table] and r['fields'].get('colId')=='Isotope']
        if len(matches)!=1: raise ValueError('Missing or ambiguous '+table+'.Isotope')
        fields=matches[0]['fields']; typ=fields.get('type')
        if typ not in ('Text','Choice'): raise ValueError('Unexpected type for '+table+'.Isotope')
        options=json.loads(fields.get('widgetOptions') or '{}')
        if not isinstance(options,dict): raise ValueError('Invalid labelling options for '+table)
        existing=options.get('choices',[])
        if not isinstance(existing,list) or any(not isinstance(v,str) for v in existing):
            raise ValueError('Invalid labelling choices for '+table)
        values=[]
        for row in client.records(table):
            if 'Isotope' in row.get('errors',{}): raise ValueError('Cannot read labelling value for '+table)
            value=row['fields'].get('Isotope')
            if value is not None and not isinstance(value,str): raise ValueError('Non-text labelling value in '+table)
            if value: values.append(value)
        for value in existing+values:
            if value and value not in choices: choices.append(value)
        descriptions.append((table,typ,options,existing))
    actions=[]
    for table,typ,options,existing in descriptions:
        if typ!='Choice' or existing!=choices:
            options['choices']=choices
            changes={'widgetOptions':json.dumps(options)}
            if typ!='Choice': changes['type']='Choice'
            actions.append(['ModifyColumn',table,'Isotope',changes])
    if actions: client.apply(actions)
    return actions


def upgrade_sample_status(client, tables=None, metadata=None, reviewed=None):
    """Migrate only known tube statuses; refuse unknown meanings before writing."""
    from .schema import SAMPLE_STATUSES
    if tables is None:tables={r['fields']['tableId']:r['id'] for r in client.records('_grist_Tables')}
    if 'Samples' not in tables:raise ValueError('Samples table missing')
    if metadata is None:metadata=client.records('_grist_Tables_column')
    matches=[r for r in metadata if r['fields'].get('parentId')==tables['Samples'] and r['fields'].get('colId')=='Available']
    if len(matches)!=1:raise ValueError('Samples.Available column missing or ambiguous')
    f=matches[0]['fields'];formula=f.get('formula','')
    if f.get('type')!='Choice' or f.get('isFormula') or formula not in ('','"Usable"'):
        raise ValueError('Sample status schema requires review')
    if f.get('recalcDeps') not in (None,[],['L']):raise ValueError('Sample status trigger dependencies require review')
    try:options=json.loads(f.get('widgetOptions') or '{}')
    except (TypeError,ValueError):raise ValueError('Invalid sample status options') from None
    if not isinstance(options,dict):raise ValueError('Invalid sample status options')
    choices=options.get('choices',[])
    if not isinstance(choices,list) or any(not isinstance(x,str) for x in choices):raise ValueError('Invalid sample status choices')
    allowed=set(SAMPLE_STATUSES)|{'Available','Unavailable','Retired','Unknown',''}
    if any(choice not in allowed for choice in choices):raise ValueError('Custom sample status choice requires review')
    mapping={None:'Usable','':'Usable','Available':'Usable','Unknown':'Usable','Unavailable':'Dead','Retired':'Dead'}
    ids=[];values=[];reviewed=reviewed or {};seen=set()
    for row in client.records('Samples'):
        seen.add(row['id'])
        if 'Available' in row.get('errors',{}):raise ValueError('Cannot read sample status')
        value=row['fields'].get('Available')
        if value is not None and not isinstance(value,str):raise ValueError('Invalid sample status value')
        if value not in allowed and value is not None:raise ValueError('Custom sample status value requires review: '+str(value))
        updated=mapping.get(value,value)
        override=reviewed.get(row['id'])
        if override is not None:
            updated=override.get('status')
            if updated not in SAMPLE_STATUSES:raise ValueError('Invalid reviewed sample status: '+str(row['id']))
            if value not in (override.get('expected_status'),updated) or override.get('expected_notes')!=(row['fields'].get('Notes') or ''):
                raise ValueError('Sample status evidence changed since review: '+str(row['id']))
        if updated!=value:ids.append(row['id']);values.append(updated)
    if set(reviewed)-seen:raise ValueError('Reviewed sample status row missing')
    actions=[];opts=dict(options);opts['choices']=list(SAMPLE_STATUSES)
    changes={}
    if choices!=list(SAMPLE_STATUSES):changes['widgetOptions']=json.dumps(opts)
    old_description='Current tube condition. Usable is the working default, not an experiment assessment. Retired means permanently unavailable: used up, discarded or irreversibly unusable. Clear or move its recorded rack position when the tube is removed.'
    if f.get('description')==old_description:changes['description']=old_description.replace('Retired means','Dead means')
    if formula!='"Usable"':changes['formula']='"Usable"'
    if f.get('recalcWhen')!=0:changes['recalcWhen']=0
    if f.get('recalcDeps') not in (None,['L']):changes['recalcDeps']=[['L']]
    if changes:actions.append(['ModifyColumn','Samples','Available',changes])
    if ids:actions.append(['BulkUpdateRecord','Samples',ids,{'Available':values}])
    if actions:client.apply(actions)
    return actions


def upgrade_protein_choices(client, tables=None, metadata=None):
    """Seed open Choice suggestions without changing values or trigger formulas."""
    if tables is None:
        tables={r['fields']['tableId']:r['id'] for r in client.records('_grist_Tables')}
    if metadata is None:metadata=client.records('_grist_Tables_column')
    choices=[];descriptions=[]
    def add(value):
        if value is not None and not isinstance(value,str):raise ValueError('Non-text protein value')
        if value and value not in choices:choices.append(value)
    for table in ('Series','Samples'):
        if table not in tables:continue
        matches=[r for r in metadata if r['fields'].get('parentId')==tables[table] and r['fields'].get('colId')=='Protein']
        if len(matches)!=1:raise ValueError('Missing or ambiguous '+table+'.Protein')
        fields=matches[0]['fields'];typ=fields.get('type')
        if typ not in ('Text','Choice'):raise ValueError('Unexpected type for '+table+'.Protein')
        try:options=json.loads(fields.get('widgetOptions') or '{}')
        except (ValueError,TypeError):raise ValueError('Invalid protein options for '+table) from None
        if not isinstance(options,dict):raise ValueError('Invalid protein options for '+table)
        existing=options.get('choices',[])
        if not isinstance(existing,list) or any(not isinstance(v,str) for v in existing):raise ValueError('Invalid protein choices for '+table)
        for value in existing:add(value)
        for row in client.records(table):
            if 'Protein' in row.get('errors',{}):raise ValueError('Cannot read protein value for '+table)
            add(row['fields'].get('Protein'))
        descriptions.append((table,typ,options,existing))
    if 'Experiments' in tables and any(r['fields'].get('parentId')==tables['Experiments'] and r['fields'].get('colId')=='SnapshotProtein' for r in metadata):
        for row in client.records('Experiments'):
            if 'SnapshotProtein' in row.get('errors',{}):raise ValueError('Cannot read historical protein snapshot')
            add(row['fields'].get('SnapshotProtein'))
    actions=[]
    for table,typ,options,existing in descriptions:
        if typ!='Choice' or existing!=choices:
            options['choices']=list(choices);changes={'widgetOptions':json.dumps(options)}
            if typ!='Choice':changes['type']='Choice'
            actions.append(['ModifyColumn',table,'Protein',changes])
    if actions:client.apply(actions)
    return actions


LEGACY_ACQUISITION_FORMULAS = {'CoverageConditions': 'parts = []\nfor label, val, unit in [("Construct", $SnapshotConstruct, ""), ("Label", $SnapshotIsotope, ""), ("Buffer", $SnapshotBuffer, ""), ("pH", $SnapshotpH, ""), ("Protein", $SnapshotProtein_uM, "µM"), ("Notes", $SnapshotOtherConditions, ""), ("Instrument", $EffectiveInstrument, "")]:\n  if val is not None and val != "": parts.append(label + ": " + str(val) + unit)\nif $CoverageUncertaintyKey: parts.append("Conditions incomplete")\nreturn "; ".join(parts)', 'EffectiveField_MHz': '$Field_MHz if $Field_MHz is not None else $ObservedField_MHz'}

def upgrade_acquisition_schema(client):
    """Add protocol fields and initialize new stored numerics to null atomically."""
    tables={r['fields']['tableId']:r['id'] for r in client.records('_grist_Tables')}
    if 'Experiments' not in tables: raise ValueError('Experiments table is missing')
    rows=client.records('Experiments')
    metadata=client.records('_grist_Tables_column')
    existing={r['fields']['colId']:r for r in metadata if r['fields'].get('parentId')==tables['Experiments']}
    desired=acquisition_columns()
    actions=[]; fresh_numeric=[]
    for column in desired:
        found=existing.get(column['id'])
        if found:
            f=found['fields']
            if f.get('type') != column['type'] or bool(f.get('isFormula')) != column['isFormula']:
                raise ValueError('Existing acquisition schema differs at Experiments.'+column['id'])
            if column['type']=='Numeric' and not column['isFormula']:
                if f.get('formula','') not in ('', 'None'):
                    raise ValueError('Existing acquisition default differs at Experiments.'+column['id'])
                if f.get('formula','') != 'None' or f.get('recalcWhen') != 0 or f.get('recalcDeps') not in (None,['L']):
                    actions.append(['ModifyColumn','Experiments',column['id'],{'formula':'None','recalcWhen':0,'recalcDeps':[['L']]}])
            elif f.get('formula','') != column['formula']:
                raise ValueError('Existing acquisition schema differs at Experiments.'+column['id'])
        else:
            info={k:v for k,v in column.items() if k!='id'}
            # Grist's column-info metadata expects one nested encoded list;
            # direct _grist_Tables_column UpdateRecord uses the flat encoding.
            if info.get('recalcDeps') == ['L']: info['recalcDeps']=[['L']]
            actions.append(['AddColumn','Experiments',column['id'],info])
            if column['type']=='Numeric' and not column['isFormula']: fresh_numeric.append(column['id'])
        if column.get('description') and (not found or found['fields'].get('description') != column['description']):
            actions.append(['ModifyColumn','Experiments',column['id'],{'description':column['description']}])
    # Grist creates Numeric columns with zero defaults. Ensure old experiments
    # instead have unknown values; both operations are a single transaction.
    if fresh_numeric and rows:
        ids=[r['id'] for r in rows]
        actions.append(['BulkUpdateRecord','Experiments',ids,{name:[None]*len(ids) for name in fresh_numeric}])
    pulse=existing.get('PulseProgram')
    if not pulse: raise ValueError('Missing existing field: PulseProgram')
    if pulse['fields'].get('label') != 'Pulse program':
        actions.append(['UpdateRecord','_grist_Tables_column',pulse['id'],{'label':'Pulse program','untieColIdFromLabel':True}])
    # Only replace the specific formula versions shipped by this register.
    from .schema import columns
    wanted={c['id']:c for c in columns()['Experiments']}
    for name,old in LEGACY_ACQUISITION_FORMULAS.items():
        found=existing.get(name)
        if not found:
            if name=='CoverageConditions': continue  # Added by the subsequent Electrostatics upgrade.
            raise ValueError('Missing existing formula: '+name)
        actual=found['fields'].get('formula',''); new=wanted[name]['formula']
        if actual==old:
            actions.append(['UpdateRecord','_grist_Tables_column',found['id'],{'formula':new}])
        elif actual!=new: raise ValueError('Existing formula differs at Experiments.'+name)
    if actions: client.apply(actions)
    return actions


def upgrade_electrostatics_schema(client):
    """Add only this feature's missing columns; preserve all existing values."""
    tables={r['fields']['tableId']:r['id'] for r in client.records('_grist_Tables')}
    if 'Experiments' not in tables: raise ValueError('Experiments table is missing')
    existing={r['fields']['colId']:r['fields'] for r in client.records('_grist_Tables_column') if r['fields'].get('parentId')==tables['Experiments']}
    definitions=next(a[2] for a in schema_actions() if a[1]=='Experiments')
    names=('Study','AnalysisStatus','AnalysisDisplay','SourceLabel','CoverageConditions','CoverageUncertaintyKey')
    new=[]
    for column in definitions:
        if column['id'] not in names: continue
        found=existing.get(column['id'])
        if found:
            if found.get('type') != column['type'] or bool(found.get('isFormula')) != column['isFormula'] or found.get('formula','') != column['formula']:
                raise ValueError('Existing register schema differs at Experiments.'+column['id'])
        else:
            new.append(['AddVisibleColumn','Experiments',column['id'],{k:v for k,v in column.items() if k!='id'}])
    if new: client.apply(new)


def apply_experiment_date_format(client):
    """Show experiment dates without changing their stored DateTime values."""
    table = next(r['id'] for r in client.records('_grist_Tables')
                 if r['fields']['tableId'] == 'Experiments')
    actions = []
    for row in client.records('_grist_Tables_column'):
        fields = row['fields']
        if fields.get('parentId') != table or fields.get('colId') not in DATE_ONLY_EXPERIMENT_COLUMNS:
            continue
        options = json.loads(fields.get('widgetOptions') or '{}')
        if any(options.get(k) != v for k, v in DATE_ONLY_OPTIONS.items()):
            options.update(DATE_ONLY_OPTIONS)
            actions.append(['UpdateRecord', '_grist_Tables_column', row['id'],
                            {'widgetOptions': json.dumps(options)}])
    if actions:
        client.apply(actions)


def electrostatics_formulas():
    """Counts are descriptive and never infer biological replication."""
    base='r.SourceID and r.Folder and r.Expno and r.AcquisitionKey and r.SourceID == r.ObservedSourceID and r.Folder == r.ObservedFolder and r.Expno == r.ObservedExpno and (r.RawPresent or r.AcquisitionStatus in ("Done", "Acquired", "Finished", "Completed")) and not r.ReconciliationReason and len(Experiments.lookupRecords(SourceID=r.SourceID, Folder=r.Folder, Expno=r.Expno)) == 1'
    return {
      'VerifiedAcquisitions':'len([r for r in $group if '+base+'])',
      'KnownTubes':'len(set(r.Sample.id for r in $group if r.Sample))',
      'UnknownSamples':'len([r for r in $group if not r.Sample])',
      'SelectedCSVs':'len([r for r in $group if r.ResultCSV])',
      'RecordsDone': '''"%s / %s" % (len($group), len([r for r in $group if r.EffectiveAcquisitionStatus == "Completed"]))''',
      'AnalysisSummary': '''"%s / %s" % (len([r for r in $group if not r.ReconciliationReason and (not r.AcquisitionKey or (r.SourceID == r.ObservedSourceID and r.Folder == r.ObservedFolder and r.Expno == r.ObservedExpno)) and r.AnalysisStatus == "Analyzed"]), len([r for r in $group if not r.ReconciliationReason and (not r.AcquisitionKey or (r.SourceID == r.ObservedSourceID and r.Folder == r.ObservedFolder and r.Expno == r.ObservedExpno)) and r.ResultCSV]))''',
      'AssessmentSummary': '''"%s / %s / %s / %s" % (len([r for r in $group if r.Assessment == "Use"]), len([r for r in $group if r.Assessment == "Limited use"]), len([r for r in $group if r.Assessment == "Do not use"]), len([r for r in $group if not r.Assessment]))''',
    }

def _result(client, action):
    """Return a Grist apply endpoint's first action result."""
    result=client.apply([action])
    if isinstance(result, dict): result=result.get('retValues', result.get('results'))
    if not isinstance(result, list) or len(result) != 1 or not isinstance(result[0], dict):
        raise ValueError('Grist apply transport must return retValues for each user action')
    return result[0]

def _compact_experiment_rows(sections, grid):
    matches=[r for r in sections if r['id']==grid]
    if len(matches)!=1:raise ValueError('Experiment grid missing or ambiguous')
    options=json.loads(matches[0]['fields'].get('options') or '{}')
    if not isinstance(options,dict):raise ValueError('Invalid experiment grid options')
    if options.get('rowHeight')==1 and options.get('rowHeightUniform') is True:return None
    options.update(rowHeight=1,rowHeightUniform=True)
    return ['UpdateRecord','_grist_Views_section',grid,{'options':json.dumps(options)}]


def build_layout(client, electrostatics_only=False):
    """Create four pages and native linked widgets in a freshly-created register.

    Requires records(table), add(table, fields), and apply(actions), where apply
    returns Grist's retValues. Refuses existing named pages to avoid duplication.
    The caller is responsible for authenticating the transport and reviewing the
    resulting view in a disposable document before use.
    """
    from .schema import CONDITIONS
    import json
    tables={r['fields']['tableId']:r['id'] for r in client.records('_grist_Tables')}
    byid={v:k for k,v in tables.items()}
    cols={(byid[r['fields']['parentId']],r['fields']['colId']):r['id']
          for r in client.records('_grist_Tables_column') if r['fields']['parentId'] in byid}
    def page(name, table):
        matches=[r['id'] for r in client.records('_grist_Views') if r['fields'].get('name') == name]
        if len(matches)>1: raise ValueError('Ambiguous page: '+name)
        if matches: return matches[0]
        view=_result(client,['AddView',table,'empty',name])['id']
        return view
    def section(view, table, kind='record', grouping=None, title=''):
        sections=client.records('_grist_Views_section')
        matches=[r['id'] for r in sections if r['fields'].get('parentId') == view and r['fields'].get('title') == title]
        if len(matches)>1: raise ValueError('Ambiguous widget: '+title)
        if matches:
            if grouping is not None:
                current=next(r for r in sections if r['id']==matches[0]); ref=current['fields'].get('tableRef')
                existing={r['fields'].get('summarySourceCol') for r in client.records('_grist_Tables_column')
                          if r['fields'].get('parentId')==ref and r['fields'].get('summarySourceCol')}
                wanted={cols[(table,name)] for name in grouping}
                if existing != wanted:
                    client.apply([['UpdateSummaryViewSection',matches[0],[cols[(table,name)] for name in grouping]]])
            return matches[0]
        # Recover a section created immediately before an interrupted title update.
        candidates=[r for r in sections if r['fields'].get('parentId') == view
                    and not r['fields'].get('title') and r['fields'].get('parentKey') == kind]
        if grouping is None:
            unmatched=[r['id'] for r in candidates if r['fields'].get('tableRef') == tables[table]]
        else:
            source_tables={r['id']:r['fields'].get('summarySourceTable') for r in client.records('_grist_Tables')}
            all_columns=client.records('_grist_Tables_column')
            wanted={cols[(table,c)] for c in grouping}
            unmatched=[]
            for r in candidates:
                ref=r['fields'].get('tableRef')
                refs={c['fields'].get('summarySourceCol') for c in all_columns
                      if c['fields'].get('parentId') == ref and c['fields'].get('summarySourceCol')}
                if source_tables.get(ref) == tables[table] and refs == wanted: unmatched.append(r['id'])
        if len(unmatched)>1: raise ValueError('Ambiguous unnamed widget: '+title)
        if unmatched:
            update('_grist_Views_section',unmatched[0],{'title':title})
            return unmatched[0]
        refs=None if grouping is None else [cols[(table,c)] for c in grouping]
        result=_result(client,['CreateViewSection',tables[table],view,kind,refs,None])['sectionRef']
        update('_grist_Views_section',result,{'title':title})
        return result
    def update(table,row,fields):
        client.apply([['UpdateRecord',table,row,fields]])
    def link(target,source,targetcol=None):
        update('_grist_Views_section',target,{'linkSrcSectionRef':source,'linkTargetColRef':targetcol or 0})
    def arrange(view, sections):
        # Keep already installed registration forms when setup is rerun.
        forms=[r['id'] for r in client.records('_grist_Views_section')
               if r['fields'].get('parentId')==view and r['fields'].get('title') in ('Register samples','Register experiment')]
        all_sections=list(sections)+[ref for ref in forms if ref not in sections]
        update('_grist_Views',view,{'layoutSpec':json.dumps({'children':[{'leaf':s} for s in all_sections]})})
    def fields(section_id, table, names, widths=None):
        wanted=[cols[(table,name)] for name in names]
        current=[r for r in client.records('_grist_Views_section_field') if r['fields'].get('parentId')==section_id]
        bycol={r['fields'].get('colRef'):r['id'] for r in current}
        actions=[['RemoveRecord','_grist_Views_section_field',r['id']] for r in current if r['fields'].get('colRef') not in wanted]
        for i,colref in enumerate(wanted,1):
            values={'parentPos':i,'width':(widths or {}).get(names[i-1],150)}
            if colref in bycol: actions.append(['UpdateRecord','_grist_Views_section_field',bycol[colref],values])
            else: actions.append(['AddRecord','_grist_Views_section_field',None,dict(values,parentId=section_id,colRef=colref)])
        if actions: client.apply(actions)
    def pins(section_id, table, names):
        current={(r['fields'].get('viewSectionRef'),r['fields'].get('colRef')) for r in client.records('_grist_Filters')}
        actions=[]
        for name in names:
            ref=cols[(table,name)]
            if (section_id,ref) not in current:
                actions.append(['AddRecord','_grist_Filters',None,{'viewSectionRef':section_id,'colRef':ref,'filter':'{}','pinned':True}])
        if actions: client.apply(actions)
    experiment_columns=['AcquiredAt','Folder','Instrument','AcquisitionStatus','ExperimentCategory','ExperimentType','PulseProgram','Expno','Sample','Title','Assessment']
    experiment_detail=['Study','AnalysisStatus','SourceID','Folder','Expno','Series','Sample','Instrument','ExperimentCategory','ExperimentType','AcquiredAt','AcquisitionStatus','ObservedAcquisitionEvidence','Field_MHz','Temperature_C','PulseProgram','D1_s','DelayList_s','SpinLock_Hz','SpinLockOffset_Hz','SpectralWidths_ppm','SpectralOffsets_ppm','ObservedProtocolEvidence','ObservedSampleID','ObservedSampleTagEvidence','Title','SnapshotProtein','SnapshotConstruct','SnapshotBuffer','SnapshotpH','SnapshotNaCl_mM','SnapshotDNA_uM','SnapshotTFE_pct','SnapshotProtein_uM','SnapshotOtherConditions','Assessment','AssessmentReason','Notes','AttentionReason','ResultCSV']
    sample_columns=['SampleID','Protein','Series','NaCl_mM','DNA_uM','TFE_pct','StorageRow','StorageColumn','Available']
    sample_detail=['SampleID','PreparedAt','Isotope','Protein','Construct','DistinguishingCondition','SampleNumber','Series','Buffer','pH','NaCl_mM','DNA_uM','TFE_pct','Protein_uM','OtherConditions','StorageRow','StorageColumn','Available','CurrentLocation','Notes']
    if not electrostatics_only:
        v=page('Experiments','Experiments'); e=section(v,'Experiments',title='Experiments'); card=section(v,'Experiments','single',title='Experiment details'); link(card,e); arrange(v,[e,card])
        compact=_compact_experiment_rows(client.records('_grist_Views_section'),e)
        if compact:client.apply([compact])
        fields(e,'Experiments',experiment_columns); fields(card,'Experiments',experiment_detail)
        pins(e,'Experiments',['Study','Sample','SnapshotProtein','ExperimentType','Instrument','Field_MHz','Temperature_C','SnapshotNaCl_mM','SnapshotDNA_uM','SnapshotTFE_pct','Assessment'])
        stale=[r['id'] for r in client.records('_grist_Filters') if r['fields'].get('viewSectionRef')==e and r['fields'].get('colRef')==cols[('Experiments','Temperature_C')]]
        if stale: client.apply([['RemoveRecord','_grist_Filters',i] for i in stale])
        v=page('Samples','Samples'); series=section(v,'Series',title='Series'); sample=section(v,'Samples',title='Samples'); detail=section(v,'Samples','single',title='Sample details'); exp=section(v,'Experiments',title='Sample experiments'); grid=section(v,'StorageGrid',title='Storage grid')
        link(sample,series,cols[('Samples','Series')]); link(detail,sample); link(exp,sample,cols[('Experiments','Sample')]); arrange(v,[series,sample,detail,exp,grid])
        fields(series,'Series',['SeriesID','Protein','Construct','Isotope','Buffer','pH','NaCl_mM','DNA_uM','TFE_pct','Protein_uM','OtherConditions','Notes'])
        fields(sample,'Samples',sample_columns); fields(detail,'Samples',sample_detail); fields(exp,'Experiments',experiment_columns)
        rack=next((r['fields'] for r in client.records('_grist_Views_section') if r['id']==grid),None)
        if not rack:raise ValueError('Storage widget is missing')
        if rack.get('parentKey')=='record' and rack.get('tableRef')==tables['StorageGrid']:
            fields(grid,'StorageGrid',['Row']+['Col'+str(i) for i in range(1,13)])
        elif rack.get('parentKey') not in ('creator','custom') or rack.get('tableRef')!=tables['Samples']:
            raise ValueError('Storage widget type or table requires review')
        pins(sample,'Samples',['Protein','Series','NaCl_mM','DNA_uM','TFE_pct','Available'])
        v=page('Needs attention','Experiments'); attention=section(v,'Experiments',title='Needs attention'); card=section(v,'Experiments','single',title='Attention details'); scanner=section(v,'ScanStatus','single',title='Scanner status'); link(card,attention)
        fields(attention,'Experiments',['AttentionFolder','AttentionExpno']+[name for name in experiment_columns if name not in ('SourceLabel','Folder','Expno')]+['AttentionReason'], {'AttentionFolder':300,'AttentionExpno':100}); fields(card,'Experiments',experiment_detail+['MigrationIssue'])
        fields(scanner,'ScanStatus',['Scanner','LastScanAt','LastSyncAt','LastError','ScannedRoot'])
        filters=[r for r in client.records('_grist_Filters') if r['fields'].get('viewSectionRef')==attention and r['fields'].get('colRef')==cols[('Experiments','NeedsAttention')]]
        if not filters: client.apply([['AddRecord','_grist_Filters',None,{'viewSectionRef':attention,'colRef':cols[('Experiments','NeedsAttention')], 'filter':json.dumps({'included':[True]}),'pinned':True}]])
        arrange(v,[attention,card,scanner])
        v=page('Overview','Experiments'); summary=section(v,'Experiments','record',['SnapshotProtein','SnapshotConstruct','SnapshotBuffer','SnapshotpH','SnapshotNaCl_mM','SnapshotDNA_uM','SnapshotTFE_pct','SnapshotProtein_uM','EffectiveTemperature_C','EffectiveField_MHz'],title='Conditions and counts'); arrange(v,[summary])
        summary_record=next(r for r in client.records('_grist_Views_section') if r['id']==summary)
        summary_ref=summary_record['fields']['tableRef']
        summary_table=next(r['fields']['tableId'] for r in client.records('_grist_Tables') if r['id']==summary_ref)
        aggregate_formulas={
            'UseCount':'len([r for r in $group if r.Assessment == "Use"])',
            'LimitedCount':'len([r for r in $group if r.Assessment == "Limited use"])',
            'ExcludedCount':'len([r for r in $group if r.Assessment == "Do not use"])',
            'UnassessedCount':'len([r for r in $group if not r.Assessment])',
            'ResultsCount':'len([r for r in $group if r.ResultCSV])',
        }
        present={r['fields']['colId'] for r in client.records('_grist_Tables_column') if r['fields'].get('parentId')==summary_ref}
        for name,formula in aggregate_formulas.items():
            if name not in present:
                client.apply([['AddVisibleColumn',summary_table,name,{'type':'Text' if name in ('RecordsDone','AnalysisSummary','AssessmentSummary') else 'Int','isFormula':True,'formula':formula}]])
        for r in client.records('_grist_Tables_column'):
            if r['fields'].get('parentId')==summary_ref: cols[(summary_table,r['fields']['colId'])]=r['id']
        fields(summary,summary_table,['SnapshotProtein','SnapshotConstruct','SnapshotBuffer','SnapshotpH','SnapshotNaCl_mM','SnapshotDNA_uM','SnapshotTFE_pct','SnapshotProtein_uM','EffectiveTemperature_C','EffectiveField_MHz','count',*aggregate_formulas])
    # Summary keys include all frozen conditions even when not all are visible.
    # Each unique identity is counted only if its machine-observed tuple agrees,
    # there is no unresolved identity/migration issue, and no other row has that tuple.
    grouping=['Study','CoverageUncertaintyKey','CoverageConditions']+['Snapshot'+k for k in ('Protein','Construct','Isotope','Buffer','pH','NaCl_mM','DNA_uM','TFE_pct','Protein_uM','OtherConditions')]+['EffectiveType','EffectiveInstrument','EffectiveField_MHz','EffectiveTemperature_C','EffectivePulseProgram','EffectiveD1_s','EffectiveDelayList_s','EffectiveSpinLock_Hz','EffectiveSpinLockOffset_Hz']
    v=page('Electrostatics','Experiments')
    matrices=[r for r in client.records('_grist_Views_section') if r['fields'].get('parentId')==v and r['fields'].get('title')=='Electrostatics matrix']
    if len(matrices)>1:raise ValueError('Ambiguous Electrostatics matrix')
    if matrices:
        if matrices[0]['fields'].get('tableRef') != tables['Experiments']:
            raise ValueError('Unexpected Electrostatics matrix table')
        coverage=matrices[0]['id']
    else:
        coverage=section(v,'Experiments','record',grouping,title='Condition coverage')
    acquisitions=section(v,'Experiments','record',title='Individual acquisitions')
    link(acquisitions,coverage)
    details=section(v,'Experiments','single',title='Acquisition details')
    link(details,acquisitions)
    arrange(v,[coverage,acquisitions,details])
    if not matrices:
        update('_grist_Views_section',coverage,{'description':'Records / completed includes all records and those with completed acquisition status. Tubes are distinct linked physical samples. Analyzed / CSV means marked analyzed / selected CSV. Assessment shows Use / Limited use / Do not use / blank. Unknown conditions stay separate; counts are not biological replicates.'})
        sr=next(r for r in client.records('_grist_Views_section') if r['id']==coverage)['fields']['tableRef']
        st=next(r['fields']['tableId'] for r in client.records('_grist_Tables') if r['id']==sr)
        formulas=electrostatics_formulas()
        present={r['fields']['colId']:r for r in client.records('_grist_Tables_column') if r['fields'].get('parentId')==sr}
        labels={'VerifiedAcquisitions':'Matched','KnownTubes':'Tubes','UnknownSamples':'No tube',
                'RecordsDone': '''"%s / %s" % (len($group), len([r for r in $group if r.EffectiveAcquisitionStatus == "Completed"]))''',
          'RecordsDone':'Records / completed','AnalysisSummary':'Analyzed / CSV','AssessmentSummary':'Use / limited / no / blank'}
        for name,formula in formulas.items():
            opts={'type':'Text' if name in ('RecordsDone','AnalysisSummary','AssessmentSummary') else 'Int',
                  'isFormula':True,'formula':formula}
            if name in labels: opts['label']=labels[name]
            if name not in present:
                client.apply([['AddVisibleColumn',st,name,opts]])
            elif present[name]['fields'].get('formula') != formula:
                client.apply([['UpdateRecord','_grist_Tables_column',present[name]['id'],{'formula':formula}]])
        for r in client.records('_grist_Tables_column'):
            if r['fields'].get('parentId')==sr: cols[(st,r['fields']['colId'])]=r['id']
        fields(coverage,st,['SnapshotProtein','SnapshotNaCl_mM','EffectiveTemperature_C','EffectiveField_MHz','EffectiveType','CoverageConditions','RecordsDone','KnownTubes','AnalysisSummary','AssessmentSummary'], {'SnapshotProtein':100,'SnapshotNaCl_mM':80,'EffectiveTemperature_C':85,'EffectiveField_MHz':75,'EffectiveType':85,'CoverageConditions':205,'RecordsDone':100,'KnownTubes':65,'AnalysisSummary':110,'AssessmentSummary':145})
        filters=[r for r in client.records('_grist_Filters') if r['fields'].get('viewSectionRef')==coverage and r['fields'].get('colRef')==cols[(st,'Study')]]
        if not filters: client.apply([['AddRecord','_grist_Filters',None,{'viewSectionRef':coverage,'colRef':cols[(st,'Study')],'filter':json.dumps({'included':['Electrostatics']}),'pinned':False}]])
    fields(acquisitions,'Experiments',['EffectiveAcquiredAt','SnapshotSampleID','SourceLabel','EffectiveAcquisitionStatus','AnalysisDisplay','Assessment','ResultCSV'], {'EffectiveAcquiredAt':100,'SnapshotSampleID':210,'SourceLabel':230,'EffectiveAcquisitionStatus':100,'AnalysisDisplay':100,'Assessment':100,'ResultCSV':90})
    fields(details,'Experiments',['Study','AcquiredAt','Sample','SnapshotSampleID','ExperimentCategory','ExperimentType','SourceID','Folder','Expno','AcquisitionStatus','ObservedAcquisitionEvidence','AnalysisStatus','Assessment','AssessmentReason','ResultCSV','AttentionReason','SnapshotConstruct','SnapshotIsotope','SnapshotBuffer','SnapshotpH','SnapshotNaCl_mM','SnapshotDNA_uM','SnapshotTFE_pct','SnapshotProtein_uM','SnapshotOtherConditions','Instrument','Field_MHz','Temperature_C','PulseProgram','D1_s','DelayList_s','SpinLock_Hz','SpinLockOffset_Hz','SpectralWidths_ppm','SpectralOffsets_ppm','ObservedProtocolEvidence','ObservedSampleID','ObservedSampleTagEvidence','Title','Notes'])
    pins(acquisitions,'Experiments',['EffectiveType','Assessment'])
    if electrostatics_only:
        # Add the two editable human fields to the existing experiment detail
        # without resetting that page's arrangement or other visible fields.
        experiment_page=next(r['id'] for r in client.records('_grist_Views') if r['fields'].get('name')=='Experiments')
        exp_sections=client.records('_grist_Views_section')
        exp_card=next(r['id'] for r in exp_sections if r['fields'].get('parentId')==experiment_page and r['fields'].get('title')=='Experiment details')
        existing_fields=[r for r in client.records('_grist_Views_section_field') if r['fields'].get('parentId')==exp_card]
        existing_refs={r['fields'].get('colRef') for r in existing_fields}
        position=max((r['fields'].get('parentPos') or 0 for r in existing_fields),default=0)
        for name in ('Study','AnalysisStatus','EffectivePulseProgram','PulseProgram','EffectiveD1_s','D1_s','EffectiveDelayList_s','DelayList_s','EffectiveSpinLock_Hz','SpinLock_Hz','EffectiveSpinLockOffset_Hz','SpinLockOffset_Hz','SpectralWidths_ppm','SpectralOffsets_ppm','ObservedProtocolEvidence'):
            if cols[('Experiments',name)] not in existing_refs:
                position+=1
                client.apply([['AddRecord','_grist_Views_section_field',None,{'parentId':exp_card,'colRef':cols[('Experiments',name)],'parentPos':position,'width':150}]])
    if not electrostatics_only:
        existing=[r['fields'].get('Row') for r in client.records('StorageGrid')]
        if any(not isinstance(row,str) or len(row)!=1 or row not in 'ABCDEFGH' for row in existing) or len(set(existing)) != len(existing):
            raise ValueError('Storage grid rows are invalid or duplicated; review before continuing')
        for row in 'ABCDEFGH':
            if row not in existing: client.add('StorageGrid',{'Row':row})
    # The acquisition-field reset above removes nonstandard fields; restore the
    # native row link here, including on electrostatics-only setup runs.
    from .navigation import apply_experiment_navigation
    apply_experiment_navigation(client,dry_run=False)


def _mcp(client, name, arguments):
    """Call a Grist MCP tool and decode its JSON text result, failing closed."""
    import json
    result=client.mcp_call(name,arguments)
    if result.get('isError'):
        message=' '.join(b.get('text','') for b in result.get('content',[]) if b.get('type')=='text')
        if '429' in message or 'monthly API limit' in message:
            from .grist import GristError
            raise GristError('Grist HTTP 429: monthly API limit reached')
        raise ValueError('Grist MCP tool failed: '+name)
    blocks=[block.get('text','') for block in result.get('content',[]) if block.get('type')=='text']
    if len(blocks)!=1: raise ValueError('Unexpected Grist MCP response: '+name)
    try: return json.loads(blocks[0])
    except (TypeError,ValueError): raise ValueError('Unexpected Grist MCP result: '+name) from None

def build_registration_widgets(client,doc_id):
    """Create/update the two document-owned forms; safe to rerun by title.

    These forms need full document access to atomically create linked records.
    """
    from .registration_widget import registration_html
    pages=_mcp(client,'grist_get_pages',{'doc_id':doc_id})
    result={}
    for page_name,table,mode,title in (
        ('Samples','Samples','sample','Register samples'),
        ('Experiments','Experiments','experiment','Register experiment')):
        matches=[p for p in pages if p.get('name')==page_name]
        if len(matches)!=1: raise ValueError('Page is missing or ambiguous: '+page_name)
        page_id=matches[0]['id']
        widgets=_mcp(client,'grist_get_page_widgets',{'doc_id':doc_id,'page_id':page_id})
        matches=[w for w in widgets if w.get('title')==title]
        if len(matches)>1: raise ValueError('Registration widget is ambiguous: '+title)
        if matches:
            w=matches[0]
            if w.get('type')!='creator' or w.get('table_id')!=table:
                raise ValueError('Existing registration widget has unexpected type or table: '+title)
        else:
            w=_mcp(client,'grist_add_page_widget',{'doc_id':doc_id,'page_id':page_id,
                'widget_options':{'table_id':table,'type':'creator','title':title}})
        wid=w['id'] if 'id' in w else w['widget_id']
        html=registration_html(mode)
        current=_mcp(client,'grist_get_creator_html',{'doc_id':doc_id,'widget_id':wid})
        if current.get('html') != html:
            _mcp(client,'grist_set_creator_html',{'doc_id':doc_id,'widget_id':wid,'html':html})
        _mcp(client,'grist_set_custom_widget_settings',{'doc_id':doc_id,'widget_id':wid,'access':'full'})
        result[mode]=wid
    # build_layout also serves fresh installations; restore the compact Samples
    # presentation after both registration widgets have been installed or reused.
    from .samples_ui import apply_samples_ui
    result['assessment_layout']=apply_assessment_ui(client,doc_id)
    result['experiments_layout']=apply_experiments_ui(client,doc_id)
    result['samples_layout']=apply_samples_ui(client,doc_id)
    result['electrostatics_layout']=apply_electrostatics_ui(client,doc_id)
    return result


def apply_assessment_ui(client,doc_id):
    """Move the sole assessment editor into the Experiments widget tray."""
    import copy,re
    from .assessment_widget import assessment_html
    call=lambda name,**kwargs:_mcp(client,name,{'doc_id':doc_id,**kwargs})
    views=client.records('_grist_Views');sections=client.records('_grist_Views_section')
    tables=client.records('_grist_Tables');meta=client.records('_grist_Tables_column')
    exp=[r for r in views if r['fields'].get('name')=='Experiments']
    old=[r for r in views if r['fields'].get('name')=='Assessment']
    table=[r for r in tables if r['fields'].get('tableId')=='Experiments']
    if len(exp)!=1 or len(old)>1 or len(table)!=1:raise ValueError('Assessment pages or table are ambiguous')
    page=exp[0]['id'];table_id=table[0]['id'];old_id=old[0]['id'] if old else None
    on_page=[r for r in sections if r['fields'].get('parentId')==page]
    grids=[r for r in on_page if r['fields'].get('title')=='Experiments']
    cards=[r for r in sections if r['fields'].get('title') in ('Assessment checklist','Assess experiment') and r['fields'].get('parentId') in (page,old_id)]
    if len(grids)!=1 or grids[0]['fields'].get('tableRef')!=table_id or len(cards)>1:raise ValueError('Assessment selector or editor is ambiguous')
    grid=grids[0]['id']
    created=not cards
    if created:
        w=call('grist_add_page_widget',page_id=page,widget_options={'table_id':'Experiments','type':'creator','title':'Assess experiment'})
        card_id=w.get('id',w.get('widget_id'))
        if not isinstance(card_id,int):raise ValueError('Assessment editor creation failed')
        sections=client.records('_grist_Views_section');cards=[r for r in sections if r['id']==card_id]
    if len(cards)!=1:raise ValueError('Assessment editor missing')
    card=cards[0];card_id=card['id'];cf=card['fields'];source=cf.get('parentId')
    if cf.get('tableRef')!=table_id or cf.get('parentKey') not in ('custom','creator'):
        raise ValueError('Unexpected Assessment editor')
    old_sections=[r for r in sections if r['fields'].get('parentId')==old_id] if old_id else []
    old_grid=[r for r in old_sections if r['fields'].get('title')=='Experiments to assess']
    if len(old_grid)>1:raise ValueError('Assessment selector is ambiguous')
    if source==old_id and (len(old_grid)!=1 or cf.get('linkSrcSectionRef')!=old_grid[0]['id']):
        raise ValueError('Old assessment link requires review')
    if source==page and cf.get('linkSrcSectionRef') not in (0,None,grid):raise ValueError('Assessment editor has a custom link')
    if any(r['fields'].get('linkSrcSectionRef')==card_id for r in sections):raise ValueError('Assessment editor controls another widget')
    old_grid_id=old_grid[0]['id'] if old_grid else None
    clean_old=bool(old_id and old_grid_id and {r['id'] for r in old_sections} <= {old_grid_id,card_id})
    if clean_old:
        nav=client.records('_grist_Pages');entries=[r for r in nav if r['fields'].get('viewRef')==old_id]
        if len(entries)!=1:raise ValueError('Old assessment page navigation requires review')
        ordered=sorted(nav,key=lambda r:r['fields'].get('pagePos',0));idx=ordered.index(entries[0]);navf=entries[0]['fields']
        if navf.get('shareRef') or navf.get('options') or (idx+1<len(ordered) and ordered[idx+1]['fields'].get('indentation',0)>navf.get('indentation',0)):
            clean_old=False  # Preserve shared or parent pages.
    if clean_old:
        f=old_grid[0]['fields']
        if f.get('tableRef')!=table_id or f.get('parentKey')!='record':raise ValueError('Old assessment grid requires review')
        if any(r['id']!=card_id and r['fields'].get('linkSrcSectionRef')==old_grid_id for r in sections):raise ValueError('Old assessment grid controls another widget')
        for r in meta:
            f=r['fields'];formula=f.get('formula') or ''
            if f.get('colId')!='AssessmentLink' and (re.search(r'#a[12]\.s'+str(old_grid_id)+r'\b',formula) or re.search(r'page\s*=\s*'+str(old_id)+r'\b',formula)):raise ValueError('Old assessment grid has a formula dependency')
    def load_layout(row):
        try:layout=json.loads(row['fields'].get('layoutSpec') or '{}')
        except (TypeError,ValueError):raise ValueError('Invalid page layout') from None
        if not isinstance(layout,dict):raise ValueError('Invalid page layout')
        return layout
    def prune(node,wid):
        if isinstance(node,list):return [x for v in node if (x:=prune(v,wid)) is not None]
        if not isinstance(node,dict):return node
        if node.get('leaf')==wid:return None
        value={k:prune(v,wid) if k in ('children','collapsed') else v for k,v in node.items()}
        if 'children' in value and not value['children'] and 'leaf' not in value:return None
        return value
    linkcols=[r for r in meta if r['fields'].get('parentId')==table_id and r['fields'].get('colId')=='AssessmentLink']
    if len(linkcols)>1:raise ValueError('Assessment link column is ambiguous')
    if linkcols:
        f=linkcols[0]['fields']
        old_formula='SELF_HYPERLINK(label="Assess", page='+str(old_id)+') + "#a1.s'+str(old_grid_id)+'.r" + str($id)'
        if not old_id or not old_grid_id or f.get('type')!='Text' or not f.get('isFormula') or f.get('formula')!=old_formula:raise ValueError('Existing Assessment link requires review')
        if any(r['id']!=linkcols[0]['id'] and re.search(r'\$AssessmentLink\b',r['fields'].get('formula') or '') for r in meta):raise ValueError('Assessment link has a formula dependency')
        fields=client.records('_grist_Views_section_field')
        allowed={grid,table[0]['fields'].get('rawViewSectionRef'),table[0]['fields'].get('recordCardViewSectionRef')}
        if any(r['fields'].get('colRef')==linkcols[0]['id'] and r['fields'].get('parentId') not in allowed for r in fields):raise ValueError('Assessment link is visible in another widget')
    layout=prune(copy.deepcopy(load_layout(exp[0])),card_id) or {}
    layout.setdefault('collapsed',[]).append({'leaf':card_id})
    actions=[]
    if linkcols:actions.append(['RemoveColumn','Experiments','AssessmentLink'])
    if cf.get('title')!='Assess experiment':actions.append(['UpdateRecord','_grist_Views_section',card_id,{'title':'Assess experiment'}])
    if source!=page or cf.get('linkSrcSectionRef')!=grid or cf.get('linkSrcColRef') or cf.get('linkTargetColRef'):
        actions.append(['UpdateRecord','_grist_Views_section',card_id,{'parentId':page,'linkSrcSectionRef':grid,'linkSrcColRef':0,'linkTargetColRef':0}])
    if layout!=load_layout(exp[0]):actions.append(['UpdateRecord','_grist_Views',page,{'layoutSpec':json.dumps(layout,separators=(',',':'))}])
    if source==old_id:
        previous=load_layout(old[0]);updated=prune(copy.deepcopy(previous),card_id) or {}
        if updated!=previous:actions.append(['UpdateRecord','_grist_Views',old_id,{'layoutSpec':json.dumps(updated,separators=(',',':'))}])
    if actions:client.apply(actions)
    html=assessment_html(section_id=card_id)
    if created:call('grist_set_creator_html',widget_id=card_id,html=html)
    settings=call('grist_get_custom_widget_settings',widget_id=card_id)
    if settings.get('access')!='full':settings=call('grist_set_custom_widget_settings',widget_id=card_id,access='full')
    if settings.get('access')!='full' or settings.get('unmapped_required_columns'):raise ValueError('Assessment widget settings incomplete')
    if clean_old:call('grist_remove_page',page_id=old_id)
    return {'page':page,'grid':grid,'card':card_id,'removed_old_page':clean_old}


def apply_electrostatics_ui(client,doc_id):
    """Install or refresh the matrix and linked acquisition view."""
    from .coverage_install import install_coverage_matrix
    from .navigation import apply_experiment_navigation
    result=install_coverage_matrix(client,doc_id)
    result['experiment_navigation']=apply_experiment_navigation(client,dry_run=False)
    from .electrostatics_assessment import install_electrostatics_assessment
    result['assessment']=install_electrostatics_assessment(client,doc_id)
    return result


def _apply_registration_widget(client,doc_id,mode):
    """Refresh an existing form without resetting other page widgets."""
    from .registration_widget import registration_html
    if mode not in ('sample','experiment'): raise ValueError('Unknown registration form')
    page,title,table=('Samples','Register samples','Samples') if mode=='sample' else ('Experiments','Register experiment','Experiments')
    call=lambda name,**kwargs:_mcp(client,name,{'doc_id':doc_id,**kwargs})
    pages=[p for p in call('grist_get_pages') if p.get('name')==page]
    if len(pages)!=1: raise ValueError(page+' page missing or ambiguous')
    matches=[w for w in call('grist_get_page_widgets',page_id=pages[0]['id']) if w.get('title')==title]
    if len(matches)!=1 or matches[0].get('type')!='creator' or matches[0].get('table_id')!=table:
        raise ValueError(title+' form missing or ambiguous')
    widget=matches[0]['id']; html=registration_html(mode)
    if call('grist_get_creator_html',widget_id=widget).get('html')!=html:
        call('grist_set_creator_html',widget_id=widget,html=html)
    return widget

def apply_experiment_registration_widget(client,doc_id):
    return _apply_registration_widget(client,doc_id,'experiment')

def apply_sample_registration_widget(client,doc_id):
    return _apply_registration_widget(client,doc_id,'sample')


def _migrate_sample_dropdown(options):
    """Replace only the exact shipped condition; protect custom filters."""
    import re
    if not isinstance(options,dict): raise ValueError('Invalid sample widget options')
    if 'dropdownCondition' not in options: return options
    condition=options['dropdownCondition']
    if not isinstance(condition,dict): raise ValueError('Invalid sample dropdown condition')
    text=condition.get('text','')
    if not isinstance(text,str): raise ValueError('Invalid sample dropdown condition')
    legacy='(not $Series or choice.Series == $Series) and (not $ProteinSelector or choice.Protein == $ProteinSelector)'
    if re.sub(r'\s+','',text)==re.sub(r'\s+','',legacy):
        return {**options,'dropdownCondition':{**condition,'text':'not $Series or choice.Series == $Series'}}
    if re.search(r'\bProteinSelector\b',text):
        raise ValueError('Custom sample dropdown uses retired ProteinSelector; review the filter before proceeding')
    return options


def apply_experiments_ui(client,doc_id):
    """Compact existing experiment page; preserve records and other widgets."""
    call=lambda name,**kwargs:_mcp(client,name,{'doc_id':doc_id,**kwargs})
    pages=call('grist_get_pages')
    matches=[p for p in pages if p.get('name')=='Experiments']
    if len(matches)!=1: raise ValueError('Experiments page missing or ambiguous')
    page=matches[0]['id']; widgets=call('grist_get_page_widgets',page_id=page)
    def widget(title):
        found=[w for w in widgets if w.get('title')==title]
        if len(found)!=1: raise ValueError('Missing or ambiguous widget: '+title)
        return found[0]['id']
    grid=widget('Experiments'); detail=widget('Experiment details'); widget('Register experiment')
    tables={r['fields']['tableId']:r['id'] for r in client.records('_grist_Tables')}
    table=tables.get('Experiments')
    columns={r['fields']['colId']:r for r in client.records('_grist_Tables_column') if r['fields'].get('parentId')==table}
    for name in ('Sample','AcquisitionStatus','AcquiredAt','AcquisitionDateText','ProteinSelector'):
        if name not in columns: raise ValueError('Missing experiment column: '+name)
    sample=columns['Sample']; original=json.loads(sample['fields'].get('widgetOptions') or '{}')
    options=_migrate_sample_dropdown(original)
    sections=client.records('_grist_Views_section')
    compact=_compact_experiment_rows(sections,grid)
    actions=[compact] if compact else []
    if original != options:
        actions.append(['UpdateRecord','_grist_Tables_column',sample['id'],{'widgetOptions':json.dumps(options)}])
    for name,label,description in (
        ('AcquisitionStatus','Acquisition status','Collected from the acquisition log when available.'),
        ('AcquiredAt','Acquisition date',None)):
        row=columns[name]; update={}
        if row['fields'].get('label')!=label: update.update(label=label,untieColIdFromLabel=True)
        if description is not None and row['fields'].get('description')!=description: update['description']=description
        if update: actions.append(['UpdateRecord','_grist_Tables_column',row['id'],update])
    attention=[p['id'] for p in pages if p.get('name')=='Needs attention']
    if len(attention)>1: raise ValueError('Ambiguous Needs attention page')
    detail_sections={detail}
    if attention:
        detail_sections.update(r['id'] for r in sections
            if r['fields'].get('parentId')==attention[0] and r['fields'].get('title')=='Attention details')
    electrostatics=[p['id'] for p in pages if p.get('name')=='Electrostatics']
    if len(electrostatics)>1:raise ValueError('Ambiguous Electrostatics page')
    if electrostatics:
        detail_sections.update(r['id'] for r in sections if r['fields'].get('parentId')==electrostatics[0] and r['fields'].get('title')=='Acquisition details')
    hidden={columns[name]['id'] for name in ('ProteinSelector','AcquisitionDateText')}
    field_rows=client.records('_grist_Views_section_field')
    aliases={'EffectiveInstrument':'Instrument','EffectiveType':'ExperimentType','EffectiveAcquiredAt':'AcquiredAt','EffectiveAcquisitionStatus':'AcquisitionStatus','EffectiveField_MHz':'Field_MHz','EffectiveTemperature_C':'Temperature_C','EffectivePulseProgram':'PulseProgram','EffectiveD1_s':'D1_s','EffectiveDelayList_s':'DelayList_s','EffectiveSpinLock_Hz':'SpinLock_Hz','EffectiveSpinLockOffset_Hz':'SpinLockOffset_Hz'}
    for sid in detail_sections:
        current=[x for x in field_rows if x['fields'].get('parentId')==sid];refs={x['fields'].get('colRef') for x in current}
        for alias,target in aliases.items():
            if alias not in columns or target not in columns:raise ValueError('Missing canonical detail column')
            old=next((x for x in current if x['fields'].get('colRef')==columns[alias]['id']),None)
            if old:
                if columns[target]['id'] in refs:actions.append(['RemoveRecord','_grist_Views_section_field',old['id']])
                else:actions.append(['UpdateRecord','_grist_Views_section_field',old['id'],{'colRef':columns[target]['id']}]);refs.add(columns[target]['id'])
        if columns['Title']['id'] not in refs:
            position=max([r['fields'].get('parentPos') or 0 for r in current]+[0])+1
            actions.append(['AddRecord','_grist_Views_section_field',None,{'parentId':sid,'colRef':columns['Title']['id'],'parentPos':position,'width':200}])
    for field in field_rows:
        f=field['fields']
        if f.get('parentId') in detail_sections and f.get('colRef') in hidden:
            actions.append(['RemoveRecord','_grist_Views_section_field',field['id']])
        elif f.get('colRef')==sample['id'] and f.get('widgetOptions'):
            original_field=json.loads(f['widgetOptions'])
            field_options=_migrate_sample_dropdown(original_field)
            if original_field != field_options:
                actions.append(['UpdateRecord','_grist_Views_section_field',field['id'],{'widgetOptions':json.dumps(field_options)}])
    if actions: client.apply(actions)
    visible=['Assessment','AcquiredAt','Folder','Expno','Sample','Instrument','EffectiveAcquisitionStatus','ExperimentCategory','ExperimentType','PulseProgram','Title']
    call('grist_set_widget_fields',widget_id=grid,column_ids=visible)
    for name,width in (('AcquiredAt',95),('Folder',220),('Instrument',95),('EffectiveAcquisitionStatus',95),('ExperimentType',85),('PulseProgram',150),('Expno',65),('Sample',200),('Title',210),('Assessment',110)):
        call('grist_update_widget_field',widget_id=grid,column_id=name,options={'width':width})
    call('grist_set_page_widget_select_by',widget_id=detail,widget_select_by={
        'link_from_widget_id':grid,'link_from_column_id':None,'link_to_column_id':None})
    layout=call('grist_get_page_layout',page_id=page)
    all_ids=[w['id'] for w in layout['widgets']]
    if grid not in all_ids or detail not in all_ids: raise ValueError('Experiment widgets disappeared')
    collapsed=[{'leaf':wid} for wid in all_ids if wid not in (grid,detail)]
    call('grist_set_page_layout',page_id=page,layout={'children':[{'leaf':grid,'size':7},{'leaf':detail,'size':3}], 'collapsed':collapsed})
    from .experiments_ui import editable_experiment_status, apply_experiment_styles, display_partial_assessment
    editable_experiment_status(client,apply=True)
    apply_experiment_styles(client,doc_id)
    display_partial_assessment(client,apply=True)
    from .experiments_ui import configure_experiment_type_filter
    type_filter = configure_experiment_type_filter(client,apply=True)
    return {'page':page,'experiments':grid,'detail':detail,'collapsed':[x['leaf'] for x in collapsed], 'type_filter':type_filter}


def apply_completion_ui(client):
    """Add status to existing linked views without rebuilding their layouts."""
    tables={r['fields']['tableId']:r['id'] for r in client.records('_grist_Tables')}
    table=tables.get('Experiments')
    cols={r['fields']['colId']:r['id'] for r in client.records('_grist_Tables_column') if r['fields'].get('parentId')==table}
    for name in ('EffectiveAcquisitionStatus','AcquisitionStatus','ObservedAcquisitionEvidence','Assessment'):
        if name not in cols:raise ValueError('Missing completion column: '+name)
    views={r['id']:r['fields'].get('name') for r in client.records('_grist_Views')}
    sections=client.records('_grist_Views_section'); fields=client.records('_grist_Views_section_field')
    grid_targets={('Experiments','Experiments'),('Samples','Sample experiments'),('Needs attention','Needs attention'),('Electrostatics','Individual acquisitions')}
    card_targets={('Experiments','Experiment details'),('Needs attention','Attention details'),('Electrostatics','Acquisition details')}
    actions=[]
    for section in sections:
        f=section['fields'];key=(views.get(f.get('parentId')),f.get('title'));sid=section['id']
        if key not in grid_targets|card_targets:continue
        current=[r for r in fields if r['fields'].get('parentId')==sid]
        bycol={r['fields'].get('colRef'):r for r in current}
        retired=bycol.get(cols.get('ObservedAcquisitionStatus'))
        if retired:actions.append(['RemoveRecord','_grist_Views_section_field',retired['id']])
        if key in grid_targets:
            effective=bycol.get(cols['EffectiveAcquisitionStatus']); manual=bycol.get(cols['AcquisitionStatus'])
            if manual and effective:actions.append(['RemoveRecord','_grist_Views_section_field',effective['id']])
            elif effective:actions.append(['UpdateRecord','_grist_Views_section_field',effective['id'],{'colRef':cols['AcquisitionStatus']}])
            elif not manual:
                assessment=bycol.get(cols['Assessment']);position=(assessment['fields'].get('parentPos',100)-0.01) if assessment else max([r['fields'].get('parentPos',0) or 0 for r in current]+[0])+1
                actions.append(['AddRecord','_grist_Views_section_field',None,{'parentId':sid,'colRef':cols['AcquisitionStatus'],'parentPos':position,'width':100}])
        else:
            manual=bycol.get(cols['AcquisitionStatus'])
            position=(manual['fields'].get('parentPos',0) or 0) if manual else max([r['fields'].get('parentPos',0) or 0 for r in current]+[0])
            for offset,name in enumerate(('AcquisitionStatus','ObservedAcquisitionEvidence'),1):
                if cols[name] not in bycol:
                    actions.append(['AddRecord','_grist_Views_section_field',None,{'parentId':sid,'colRef':cols[name],'parentPos':position+offset/100,'width':150}])
    old='Records / done includes all records and those marked Done or Acquired. Tubes are distinct linked physical samples. Analyzed / CSV means marked analyzed / selected CSV. Assessment shows Use / Limited use / Do not use / blank. Unknown conditions stay separate; counts are not biological replicates.'
    new=old.replace('Records / done includes all records and those marked Done or Acquired.','Records / completed includes all records and those with completed acquisition status.')
    for section in sections:
        f=section['fields']
        if views.get(f.get('parentId'))=='Electrostatics' and f.get('title')=='Condition coverage' and f.get('description')==old:
            actions.append(['UpdateRecord','_grist_Views_section',section['id'],{'description':new}])
    if actions:client.apply(actions)
    return actions


def upgrade_experiment_taxonomy(client, apply=False):
    """Upgrade shipped experiment metadata without replacing custom formulas or layouts."""
    from .schema import columns, ATTENTION_FORMULA, PRE_CATEGORY_ATTENTION_FORMULA
    tables = [r for r in client.records('_grist_Tables') if r['fields'].get('tableId') == 'Experiments']
    if len(tables) != 1: raise ValueError('Experiments table missing or ambiguous')
    tid = tables[0]['id']
    rows = [r for r in client.records('_grist_Tables_column') if r['fields'].get('parentId') == tid]
    existing = {r['fields'].get('colId'): r for r in rows}
    wanted = {c['id']: c for c in columns()['Experiments']}
    actions = []
    category = existing.get('ExperimentCategory')
    target = wanted['ExperimentCategory']
    if not category:
        actions.append(['AddColumn', 'Experiments', 'ExperimentCategory', {k:v for k,v in target.items() if k != 'id'}])
    elif any(category['fields'].get(k) != target[k] for k in ('type', 'isFormula', 'formula')):
        raise ValueError('ExperimentCategory requires review')
    for name in ('ExperimentType', 'EffectiveType'):
        row = existing.get(name)
        if not row: raise ValueError(name + ' missing')
        f = row['fields']
        if name == 'EffectiveType' and (not f.get('isFormula') or f.get('formula') != '$ExperimentType'):
            raise ValueError('EffectiveType requires review')
        if f.get('label') in (None, '', name, 'Type', 'ExperimentType', 'Experiment'):
            update = {}
            if f.get('label') != 'Experiment': update['label'] = 'Experiment'
            if not f.get('untieColIdFromLabel'): update['untieColIdFromLabel'] = True
            if update: actions.append(['UpdateRecord', '_grist_Tables_column', row['id'], update])
    attention = existing.get('AttentionReason')
    if not attention or not attention['fields'].get('isFormula'):
        raise ValueError('AttentionReason requires review')
    current = attention['fields'].get('formula')
    if current not in (PRE_CATEGORY_ATTENTION_FORMULA, PRE_REGISTRATION_REMOVAL_ATTENTION_FORMULA, ATTENTION_FORMULA):
        raise ValueError('AttentionReason formula requires review')
    if current != ATTENTION_FORMULA:
        actions.append(['UpdateRecord', '_grist_Tables_column', attention['id'], {'formula': ATTENTION_FORMULA}])
    if apply and actions: client.apply(actions)
    return actions


def show_experiment_category(client, apply=False):
    """Add the category to existing main and attention views without resetting layouts."""
    tables = [r for r in client.records('_grist_Tables') if r['fields'].get('tableId') == 'Experiments']
    if len(tables) != 1: raise ValueError('Experiments table missing or ambiguous')
    tid = tables[0]['id']
    cols = [r for r in client.records('_grist_Tables_column') if r['fields'].get('parentId') == tid]
    categories = [r for r in cols if r['fields'].get('colId') == 'ExperimentCategory']
    if len(categories) != 1: raise ValueError('ExperimentCategory missing or ambiguous')
    category_id = categories[0]['id']
    experiment_ids = {r['id'] for r in cols if r['fields'].get('colId') in ('ExperimentType', 'EffectiveType')}
    views = {r['id']:r['fields'].get('name') for r in client.records('_grist_Views')}
    targets = {('Experiments','Experiments'), ('Experiments','Experiment details'), ('Needs attention','Needs attention'), ('Needs attention','Attention details')}
    sections = [r for r in client.records('_grist_Views_section') if r['fields'].get('tableRef') == tid and (views.get(r['fields'].get('parentId')), r['fields'].get('title')) in targets]
    fields = client.records('_grist_Views_section_field')
    actions = []
    for section in sections:
        current = [r for r in fields if r['fields'].get('parentId') == section['id']]
        if not any(r['fields'].get('colRef') == category_id for r in current):
            positions = [r['fields'].get('parentPos') or 0 for r in current if r['fields'].get('colRef') in experiment_ids]
            if positions:
                following = min(positions)
                previous = [r['fields'].get('parentPos') or 0 for r in current if (r['fields'].get('parentPos') or 0) < following]
                pos = (max(previous) + following) / 2 if previous else following - 1
            else:
                pos = max([r['fields'].get('parentPos') or 0 for r in current] + [0]) + 1
            actions.append(['AddRecord', '_grist_Views_section_field', None, {'parentId':section['id'], 'colRef':category_id, 'parentPos':pos, 'width':120}])
    if apply and actions:
        client.apply(actions)
        # New field IDs are allocated by Grist, so read them back before
        # placing them in explicit card layouts.
        fields = client.records('_grist_Views_section_field')
    layout_actions = []
    for section in sections:
        sf = section['fields']
        if sf.get('parentKey') != 'single' or not sf.get('layoutSpec'):
            continue
        current = [r for r in fields if r['fields'].get('parentId') == section['id']]
        category_fields = [r for r in current if r['fields'].get('colRef') == category_id]
        # A dry run cannot know the ID of a field that would be added.
        if not category_fields:
            continue
        if len(category_fields) != 1:
            raise ValueError('Ambiguous category field in card')
        category_field = category_fields[0]['id']
        precise_fields = {r['id'] for r in current if r['fields'].get('colRef') in experiment_ids}
        layout = json.loads(sf['layoutSpec'])
        def contains(node, ids):
            if isinstance(node, dict):
                return node.get('leaf') in ids or any(contains(v, ids) for v in node.values())
            if isinstance(node, list):
                return any(contains(v, ids) for v in node)
            return False
        if contains(layout, {category_field}):
            continue
        def insert_before(node):
            if isinstance(node, list):
                for index, child in enumerate(node):
                    if isinstance(child, dict) and child.get('leaf') in precise_fields:
                        node.insert(index, {'leaf': category_field})
                        return True
                    if insert_before(child): return True
            elif isinstance(node, dict):
                for child in node.values():
                    if insert_before(child): return True
            return False
        if not insert_before(layout):
            raise ValueError('Cannot place category: card layout has no experiment type field')
        layout_actions.append(['UpdateRecord', '_grist_Views_section', section['id'],
                               {'layoutSpec': json.dumps(layout, separators=(',', ':'))}])
    if apply and layout_actions: client.apply(layout_actions)
    return actions + layout_actions


def upgrade_spectral_schema(client, apply=False):
    """Plan or apply only spectral columns and their three acquisition cards.

    New column and field IDs require metadata rereads. A dry run reports the
    immediately actionable changes; call again after applying to verify them.
    """
    names = ('SpectralWidths_ppm', 'SpectralOffsets_ppm')
    tables = [r for r in client.records('_grist_Tables') if r['fields'].get('tableId') == 'Experiments']
    if len(tables) != 1: raise ValueError('Experiments table missing or ambiguous')
    tid = tables[0]['id']
    metadata = [r for r in client.records('_grist_Tables_column') if r['fields'].get('parentId') == tid]
    byname = {}
    for r in metadata:
        name = r['fields'].get('colId')
        if name in names:
            if name in byname: raise ValueError('Ambiguous spectral column: ' + name)
            byname[name] = r
            f = r['fields']
            if f.get('type') != 'Text' or f.get('isFormula') or f.get('formula', '') not in ('', None):
                raise ValueError('Unexpected spectral schema: ' + name)
    views = {r['id']: r['fields'].get('name') for r in client.records('_grist_Views')}
    targets = {('Experiments', 'Experiment details'), ('Needs attention', 'Attention details'), ('Electrostatics', 'Acquisition details')}
    sections = []
    for r in client.records('_grist_Views_section'):
        f = r['fields']; key = (views.get(f.get('parentId')), f.get('title'))
        if key in targets:
            if f.get('tableRef') != tid or f.get('parentKey') != 'single':
                raise ValueError('Unexpected spectral detail card: ' + str(key))
            sections.append(r)
    if len(sections) != 3 or len({(views.get(r['fields'].get('parentId')), r['fields'].get('title')) for r in sections}) != 3:
        raise ValueError('Spectral detail cards missing or ambiguous')
    fields = client.records('_grist_Views_section_field')
    # Validate existing explicit layouts before the first mutation.
    layouts = {}
    for section in sections:
        raw = section['fields'].get('layoutSpec')
        if raw:
            try: layout = json.loads(raw)
            except (TypeError, ValueError): raise ValueError('Invalid spectral card layout')
            if not isinstance(layout, dict) or not isinstance(layout.get('children'), list):
                raise ValueError('Unexpected spectral card layout')
            layouts[section['id']] = layout
    def contains(node, leaf):
        if isinstance(node, dict): return node.get('leaf') == leaf or any(contains(v, leaf) for v in node.values())
        if isinstance(node, list): return any(contains(v, leaf) for v in node)
        return False
    d1_cols = {r['id'] for r in metadata if r['fields'].get('colId') == 'D1_s'}
    for sid, layout in layouts.items():
        d1_fields = {r['id'] for r in fields if r['fields'].get('parentId') == sid and r['fields'].get('colRef') in d1_cols}
        if not any(contains(layout, field_id) for field_id in d1_fields):
            raise ValueError('Cannot place spectral fields: D1 is absent from card layout')
    actions = [['AddColumn', 'Experiments', c['id'], {k:v for k,v in c.items() if k != 'id'}]
               for c in spectral_columns() if c['id'] not in byname]
    if apply and actions:
        client.apply(actions)
        metadata = [r for r in client.records('_grist_Tables_column') if r['fields'].get('parentId') == tid]
        byname = {r['fields'].get('colId'): r for r in metadata if r['fields'].get('colId') in names}
        if set(byname) != set(names): raise ValueError('Spectral columns not allocated')
    field_actions = []
    for section in sections:
        current = [r for r in fields if r['fields'].get('parentId') == section['id']]
        position = max([r['fields'].get('parentPos') or 0 for r in current] + [0])
        for name in names:
            if name not in byname: continue
            found = [r for r in current if r['fields'].get('colRef') == byname[name]['id']]
            if len(found) > 1: raise ValueError('Ambiguous spectral card field: ' + name)
            if not found:
                position += 1
                field_actions.append(['AddRecord', '_grist_Views_section_field', None,
                                      {'parentId': section['id'], 'colRef': byname[name]['id'], 'parentPos': position, 'width': 150}])
    if apply and field_actions:
        client.apply(field_actions)
        fields = client.records('_grist_Views_section_field')
    actions += field_actions
    layout_actions = []
    for sid, layout in layouts.items():
        current = [r for r in fields if r['fields'].get('parentId') == sid]
        missing = []
        for name in names:
            if name not in byname: continue
            matches = [r['id'] for r in current if r['fields'].get('colRef') == byname[name]['id']]
            if len(matches) > 1: raise ValueError('Ambiguous spectral card field: ' + name)
            if matches and not contains(layout, matches[0]): missing.append({'leaf': matches[0]})
        if missing:
            row = {'children': missing} if len(missing) > 1 else missing[0]
            d1_cols = {r['id'] for r in metadata if r['fields'].get('colId') == 'D1_s'}
            d1_fields = {r['id'] for r in current if r['fields'].get('colRef') in d1_cols}
            index = next((i for i, child in enumerate(layout['children'])
                          if any(contains(child, field_id) for field_id in d1_fields)), None)
            if index is None: raise ValueError('Cannot place spectral fields: D1 is absent from card layout')
            layout['children'].insert(index + 1, row)
            layout_actions.append(['UpdateRecord', '_grist_Views_section', sid, {'layoutSpec': json.dumps(layout, separators=(',', ':'))}])
    if apply and layout_actions: client.apply(layout_actions)
    return actions + layout_actions
