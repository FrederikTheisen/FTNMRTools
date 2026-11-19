"""Conservative pilot migration. All inputs are read-only; application is opt-in."""
import csv
import hashlib
import re
from datetime import datetime, timezone
from pathlib import Path

INSTRUMENT_FIELDS = {'600_1300':600,'600Cryo':600,'600oxbru':600,'700':700,'850':850,'950':950}

REVIEW_MAP = {'usable':'Use', 'exclude':'Do not use', 'ok':'Use', 'bad':'Do not use', 'poor':'Limited use'}

def read_reviews(path):
    with open(path, newline='', encoding='utf-8-sig') as stream:
        for line in stream:
            if line.startswith('target\t'):
                reader=csv.DictReader([line, *stream], delimiter='\t')
                rows={}
                for row in reader:
                    target=(row.get('target') or '').strip()
                    if not target: continue
                    if target in rows: raise ValueError('Duplicate review target: '+target)
                    rows[target]=row
                return rows
    raise ValueError('Review TSV has no target header')

def _date(value):
    if isinstance(value, datetime):
        return value.replace(tzinfo=timezone.utc).timestamp()
    return None

def _str(value):
    return '' if value is None else str(value).strip()


def experiment_rows(sheet):
    """Expand A/B/C/D/F within a nonblank A:M block, preserving source cells.

    EXPNO and all per-acquisition conditions, type, settings, status and comments
    stay exactly on their own row. A fully blank A:M row resets every default.
    """
    inherited_columns=(0,1,2,3,5)
    context={};origins={}
    blank=lambda value:value is None or isinstance(value,str) and not value.strip()
    for number,raw in enumerate(sheet.iter_rows(min_row=2,values_only=True),2):
        values=list(raw[:13]);values.extend([None]*(13-len(values)))
        if all(blank(value) for value in values):
            context.clear();origins.clear();continue
        inherited={}
        for index in inherited_columns:
            if not blank(values[index]):context[index]=values[index];origins[index]=number
            elif index in context:values[index]=context[index];inherited[index]=origins[index]
        yield number,tuple(values),inherited


def _workbook_evidence(sheet,number,row,inherited):
    evidence=f'{sheet}!{number}'
    if row[6] is not None:evidence+='; condition: '+_str(row[6])
    if _str(row[2]):evidence+='; workbook user: '+_str(row[2])
    names={0:'Date',1:'Instrument',2:'User',3:'Folder',5:'Sample'}
    if inherited:evidence+='; inherited: '+', '.join(names[index]+' from '+chr(65+index)+str(origin) for index,origin in inherited.items())
    if 0 in inherited:evidence+='; group date context: '+_str(row[0])
    return evidence

def build_pilot(workbook, reviews, sample_ids, targets, source_id=""):
    """Plan explicitly selected workbook samples and acquisitions; never infer identities.

    targets are exact workbook folder/expno strings; rows without these are not
    selected. Linking requires the workbook sample cell to equal a selected ID.
    openpyxl is an optional dependency needed only for this migration.
    """
    import openpyxl
    book=openpyxl.load_workbook(workbook, data_only=True, read_only=True)
    review=read_reviews(reviews)
    sample_ids, targets=set(sample_ids),set(targets)
    samples=[]; experiments=[]; found=set()
    for n,row in enumerate(book['Samples'].iter_rows(values_only=True),1):
        if n == 1 or _str(row[1]) not in sample_ids: continue
        sample_id=_str(row[1]); found.add(sample_id)
        match=re.search(r'(?i)(?<![\w.])(\d+(?:\.\d+)?)\s*mM\s*NaCl\b', _str(row[5]))
        salt=float(match.group(1)) if match else None
        position=re.fullmatch(r'([A-H])([1-9]|1[0-2])', _str(row[0]))
        dna=re.search(r'(?i)_(\d+(?:\.\d+)?)uMDNA$',sample_id)
        samples.append({'SampleID':sample_id,'Protein':_str(row[2]),'NaCl_mM':salt,
            'DNA_uM':float(dna.group(1)) if dna else None,
            'StorageBox':'Unspecified box' if position else '',
            'StorageRow':position.group(1) if position else '',
            'StorageColumn':int(position.group(2)) if position else None,
            'Buffer':_str(row[5]),'Notes':' | '.join(x for x in (_str(row[3]),_str(row[6])) if x),
            'SourceEvidence':f'Samples!{n}; storage: {_str(row[0])}'})
    if found != sample_ids: raise ValueError('Missing sample IDs: '+', '.join(sorted(sample_ids-found)))
    if len({s['SampleID'] for s in samples}) != len(samples): raise ValueError('Duplicate selected sample ID')
    found_targets=set()
    for sheet in ('Database','ANP32','N-Tail','TDR1','CtIP'):
        for n,row,inherited in experiment_rows(book[sheet]):
            folder,expno=_str(row[3]),_str(row[4]); target=folder+'/'+expno
            if target not in targets: continue
            if target in found_targets:
                prior=next(x for x in experiments if x['fields']['Folder']+'/'+x['fields']['Expno']==target)
                prior['fields']['MigrationIssue'] += '; duplicate workbook target requires review'
                prior['fields']['MigrationIssue'] += f'; {sheet}!{n}: type={_str(row[8])}, status={_str(row[11])}, comment={_str(row[12])}'
                prior['fields']['SourceEvidence']=_append_evidence(prior['fields'].get('SourceEvidence',''),_workbook_evidence(sheet,n,row,inherited)+'; duplicate source row: '+repr(tuple(row[:13])))
                continue
            found_targets.add(target)
            status=_str(row[11]); notes=_str(row[12]); sample=_str(row[5]); r=review.get(target, {})
            old=_str(r.get('assessment')).lower(); decision=REVIEW_MAP.get(old,'')
            reason=_str(r.get('comment')); issue=[]
            if old and old not in REVIEW_MAP: issue.append('Legacy review assessment: '+old)
            if old == 'ok' and 'degrad' in reason.lower():
                issue.append('Conflict: review OK; comment '+reason); decision=''; reason=''
            if old in ('ok','usable') and any(token in (status+' '+notes).lower() for token in ('failed','poor sn','poor s/n')):
                issue.append('Conflict: review OK; workbook '+status+' '+notes); decision=''; reason=''
            if status.lower() == 'not for analysis':
                if decision and decision != 'Do not use':
                    issue.append('Conflict: workbook Not For Analysis; review '+old); decision=''; reason=''
                elif not decision: decision='Do not use'; reason=notes
            if decision in ('Limited use','Do not use') and not reason: issue.append('Assessment reason missing')
            if sample not in sample_ids: issue.append('Sample link unresolved: '+sample + ('; workbook condition: '+_str(row[6]) if row[6] is not None else ''))
            instrument=_str(row[1]); field=INSTRUMENT_FIELDS.get(instrument)
            fields={'SourceID':source_id,'Folder':folder,'Expno':expno,'Instrument':instrument,'Field_MHz':field,
                'ExperimentType':_str(row[8]),'AcquiredAt':_date(row[0]) if 0 not in inherited else None,
                'AcquisitionDateText':_str(row[0]) if 0 not in inherited and row[0] is not None and not isinstance(row[0],datetime) else '',
                'AcquisitionStatus':status,
                'Temperature_C':row[7] if isinstance(row[7],(int,float)) else None,'PulseProgram':_str(row[9]),
                'Assessment':decision,'AssessmentReason':reason,'Notes':notes,
                'SourceEvidence':_workbook_evidence(sheet,n,row,inherited),
                'MigrationKey':f'workbook:{sheet}:{n}','MigrationIssue':'; '.join(issue)}
            study = REVIEWED_PILOT_STUDIES.get((source_id,target))
            if study:
                fields['Study']=study
                fields['SourceEvidence']=_append_evidence(fields['SourceEvidence'],_study_evidence(source_id,target))
            experiments.append({'fields':fields,'sample_id':sample if sample in sample_ids else None})
    book.close()
    if found_targets != targets: raise ValueError('Missing targets: '+', '.join(sorted(targets-found_targets)))
    titration='20251017_15N_HeTDR1_DNA_titration'
    for sample in samples:
        if re.match(r'^20251017_15N_HeTDR1_\d+uMDNA$',sample['SampleID']): sample['_SeriesID']=titration
    series=[{'SeriesID':titration}] if any('_SeriesID' in x for x in samples) else []
    return {'scope':'pilot','series':series,'samples':samples,'experiments':experiments}

def apply_pilot(plan, client, dry_run=True, full_accepted=False):
    """Idempotent insert-only import; existing records are never overwritten.

    client.records(table) returns Grist API records; client.add(table, fields)
    inserts one row and returns its integer row ID. Full migration is prohibited.
    """
    if plan.get('scope') not in ('pilot','full'): raise ValueError('Unknown migration scope')
    if plan['scope'] == 'full' and not dry_run and not full_accepted:
        raise ValueError('Full migration requires recorded pilot and Mac acceptance')
    if not dry_run and any('Assessment' in entry.get('fields',{}) for entry in plan.get('experiments',[])):
        tables=[r['id'] for r in client.records('_grist_Tables') if r['fields'].get('tableId')=='Experiments']
        if len(tables)!=1:raise ValueError('Experiments table missing or ambiguous')
        assessment=[r for r in client.records('_grist_Tables_column') if r['fields'].get('parentId')==tables[0] and r['fields'].get('colId')=='Assessment']
        if len(assessment)!=1:raise ValueError('Assessment column missing or ambiguous')
        if assessment[0]['fields'].get('isFormula'):
            raise ValueError('Experiment migration requires reviewed QC mapping before writing formula Assessment')
    existing_samples=client.records('Samples'); existing_exp=client.records('Experiments')
    def unique(records, key):
        out={}
        for r in records:
            val=r['fields'].get(key)
            if val:
                if val in out: raise ValueError('Duplicate '+key+': '+str(val))
                out[val]=r['id']
        return out
    sample_map=unique(existing_samples,'SampleID'); exp_map=unique(existing_exp,'MigrationKey')
    def identity(fields):
        values=tuple(str(fields.get(k) or '').strip() for k in ('SourceID','Folder','Expno'))
        return values if all(values) else None
    identities={}
    migration_rows={r['fields'].get('MigrationKey'):r['fields'] for r in existing_exp if r['fields'].get('MigrationKey')}
    for record in existing_exp:
        key=identity(record['fields'])
        if key: identities.setdefault(key,[]).append(record)
    for entry in plan['experiments']:
        fields=entry['fields']; key=identity(fields)
        previous=migration_rows.get(fields.get('MigrationKey'))
        if previous is not None and identity(previous) != key:
            raise ValueError('Existing migration key points to a different acquisition: '+str(fields.get('MigrationKey')))
        if not key: continue
        matches=identities.get(key,[])
        if len(matches)>1 or (matches and matches[0]['fields'].get('MigrationKey') != fields.get('MigrationKey')):
            raise ValueError('Existing exact acquisition identity requires reconciliation: '+ '/'.join(key))
        if not matches:
            identities[key]=[{'fields':fields}]
    series_map=unique(client.records('Series'),'SeriesID') if plan.get('series') else {}
    for series in plan.get('series',[]):
        if series['SeriesID'] not in series_map:
            fields=dict(series)
            for name in ('pH','NaCl_mM','DNA_uM','TFE_pct','Protein_uM'): fields.setdefault(name,None)
            series_map[series['SeriesID']]=None if dry_run else client.add('Series',fields)
    actions={'samples':[], 'experiments':[]}
    if plan['scope'] == 'full':
        unclassified=unique(client.records('Unclassified'),'SourceKey')
        actions['unclassified']=[]
        for item in plan.get('unclassified_samples',[]):
            key='workbook:Samples:'+str(item['row'])
            if key in unclassified: continue
            actions['unclassified'].append(key)
            if not dry_run:
                client.add('Unclassified',{'SourceKey':key,'Label':item['label'],
                    'SourceEvidence':'Samples!'+str(item['row']),'Notes':item['evidence']})
            unclassified[key]=True
    for row in plan['samples']:
        if row['SampleID'] not in sample_map:
            actions['samples'].append(row['SampleID'])
            fields=dict(row); sid=fields.pop('_SeriesID',None)
            for name in ('pH','NaCl_mM','DNA_uM','TFE_pct','Protein_uM','StorageColumn'):
                fields.setdefault(name,None)
            if sid and not dry_run: fields['Series']=series_map[sid]
            sample_map[row['SampleID']] = None if dry_run else client.add('Samples',fields)
    for entry in plan['experiments']:
        fields=dict(entry['fields']); key=fields['MigrationKey']
        # Grist Numeric defaults to zero; unknown imported measurements remain null.
        for name in ('Temperature_C','Field_MHz','ObservedField_MHz','ObservedTemperature_C'):
            fields.setdefault(name,None)
        if key in exp_map: continue
        sample=entry['sample_id']
        if sample and not dry_run: fields['Sample']=sample_map[sample]
        actions['experiments'].append(key)
        if not dry_run:
            client.add('Experiments',fields)
        exp_map[key]=True
    return actions

# Curated pilot: sample links are explicit workbook cell values; DNA titration
# samples without an explicit experiment-level mapping remain unlinked.
PILOT_SAMPLES = (
    '20251104_15N_HeTDR1_500mMNaCl', '20251104_15N_HeTDR1_50mMNaCl',
    '20251017_15N_HeTDR1_0uMDNA', '20251017_15N_HeTDR1_100uMDNA',
    '20251017_15N_HeTDR1_200uMDNA', '20251017_15N_HeTDR1_500uMDNA',
    '20250325_15N_hANP32a_IDR',
)
PILOT_TARGETS = (
    '20251105_15N_heTDR1_NaCl/405', '20251105_15N_heTDR1_NaCl/406',
    '20251105_15N_heTDR1_NaCl/407', '20251218_15N_heDR1_50mMNaCl/5',
    '20251218_15N_heDR1_50mMNaCl/6', '20251218_15N_heDR1_50mMNaCl/7',
    '20250528_15N_hANP32a_50mMNaCl/5',
    '20251017_15N_HeTDR1_DNA/104', '20251017_15N_HeTDR1_DNA/204',
    '20251017_15N_HeTDR1_DNA/304', '20251017_15N_HeTDR1_DNA/404',
)


def build_full(workbook, reviews, source_id=""):
    """Build a reviewable full inventory, without asserting ambiguous labels are tubes.

    Ambiguous undated, titration-group and TFE-group labels are preserved in
    unclassified_samples and require individual human review. Incomplete planned
    entries are retained as unresolved experiments, not assigned acquisition keys.
    """
    import openpyxl
    book=openpyxl.load_workbook(workbook, data_only=True, read_only=True)
    physical=[]; ambiguous=[]; targets=set(); planned=[]
    for n,row in enumerate(book['Samples'].iter_rows(min_row=2,values_only=True),2):
        label=_str(row[1])
        if not label: continue
        lower=label.lower()
        if not re.match(r'^20[0-9]{4,6}_',label) or 'titrationsamples' in lower or lower.endswith('_tfe'):
            ambiguous.append({'row':n,'label':label,'evidence':repr(tuple(row[:7]))})
        else: physical.append(label)
    for sheet in ('Database','ANP32','N-Tail','TDR1','CtIP'):
        for n,row,inherited in experiment_rows(book[sheet]):
            folder,expno=_str(row[3]),_str(row[4])
            if folder and expno: targets.add(folder+'/'+expno)
            elif any(x is not None and x != '' for x in row[:13]) and (folder or row[8] or row[11] or row[12]):
                planned.append({'fields':{'SourceID':source_id,'Folder':folder,'Expno':expno,
                    'Instrument':_str(row[1]),'Field_MHz':INSTRUMENT_FIELDS.get(_str(row[1])),
                    'AcquiredAt':_date(row[0]) if 0 not in inherited else None,
                    'AcquisitionDateText':_str(row[0]) if 0 not in inherited and row[0] is not None and not isinstance(row[0],datetime) else '',
                    'ExperimentType':_str(row[8]),'AcquisitionStatus':_str(row[11]),
                    'Notes':_str(row[12]),'MigrationKey':f'workbook:{sheet}:{n}',
                    'MigrationIssue':'Planned or incomplete record; identity and sample unresolved',
                    'SourceEvidence':_workbook_evidence(sheet,n,row,inherited)+': '+repr(tuple(row[:13]))},
                    'sample_id':_str(row[5]) if _str(row[5]) in physical else None})
    book.close()
    plan=build_pilot(workbook,reviews,physical,targets,source_id)
    plan['scope']='full'; plan['experiments'].extend(planned); plan['unclassified_samples']=ambiguous
    present={x['fields'].get('Folder','')+'/'+x['fields'].get('Expno','') for x in plan['experiments']}
    for target,row in read_reviews(reviews).items():
        if target in present or '/' not in target: continue
        folder,expno=target.rsplit('/',1)
        if not folder or not expno or '..' in folder.split('/'): continue
        old=_str(row.get('assessment')).lower(); decision=REVIEW_MAP.get(old,'')
        plan['experiments'].append({'fields':{'SourceID':source_id,'Folder':folder,'Expno':expno,
            'Assessment':decision,'AssessmentReason':_str(row.get('comment')),
            'MigrationKey':'review:'+target,'SourceEvidence':'condition-reviews.tsv',
            'MigrationIssue':'Review-only record; sample and registration unresolved'},'sample_id':None})
    return plan

def reconcile_sources(plan, catalog=None, status_workbooks=(), overview=None, catalog_source_id='nmrdata'):
    """Reconcile exact-target evidence. Never infer a tube from condition matches.

    CSV observations and explicit investigator ratings in status workbooks can
    identify disagreements requiring attention; machine assessments are ignored.
    Inputs remain untouched. Results are added to the in-memory plan.
    """
    by_target={}; samples={s['SampleID']:s for s in plan['samples']}
    entries={}
    for entry in plan['experiments']:
        f=entry['fields']; target=f.get('Folder','')+'/'+f.get('Expno','')
        if f.get('Folder') and f.get('Expno'): by_target[target]=f; entries[target]=entry
    def issue(fields,message):
        fields['MigrationIssue']='; '.join(x for x in (fields.get('MigrationIssue',''),message) if x)
    def rating(fields, rating_text, comment, location):
        old=_str(rating_text).lower()
        if old not in REVIEW_MAP: return
        decision=REVIEW_MAP[old]
        if fields.get('Assessment') and fields['Assessment'] != decision:
            issue(fields,f'Conflict: review {fields["Assessment"]}; {location} {old}')
            fields['Assessment']=''; fields['AssessmentReason']=''
        elif not fields.get('Assessment') and not fields.get('MigrationIssue'):
            fields['Assessment']=decision; fields['AssessmentReason']=_str(comment)
            if decision in ('Limited use','Do not use') and not _str(comment): issue(fields,'Assessment reason missing')
    if catalog:
        analyzed,_ = _catalog_analysis(catalog,catalog_source_id)
        with open(catalog,newline='',encoding='utf-8-sig') as stream:
            for n,row in enumerate(csv.DictReader(stream),2):
                target=row.get('relative_path','')
                f=by_target.get(target)
                if not f:
                    if plan.get('scope') != 'full': continue
                    parts=target.rsplit('/',1)
                    if len(parts) != 2 or not all(parts) or '..' in parts[0].split('/'):
                        continue
                    source_id=next((x['fields'].get('SourceID','') for x in plan['experiments'] if x['fields'].get('SourceID')), '')
                    f={'SourceID':source_id,'Folder':parts[0],'Expno':parts[1],
                       'ExperimentType':_str(row.get('experiment_type')),
                       'MigrationKey':'catalog:'+target,'SourceEvidence':f'{Path(catalog).name}:{n}',
                       'MigrationIssue':'Catalog-only record; sample and registration unresolved'}
                    plan['experiments'].append({'fields':f,'sample_id':None})
                    by_target[target]=f
                identity = _exact_identity(f)
                evidence = analyzed.get(identity)
                if evidence and not _str(f.get('AnalysisStatus')) and _analysis_identity_safe(f):
                    f['AnalysisStatus']='Analyzed'
                    f['SourceEvidence']=_append_evidence(f.get('SourceEvidence'),evidence)
                for source,target in (('temperature_C','Temperature_C'),('field_MHz','Field_MHz')):
                    try: observed=float(row[source])
                    except (TypeError,ValueError,KeyError): continue
                    known=f.get(target)
                    matches = (round(observed) == known if target == 'Field_MHz' and isinstance(known,(int,float)) and float(known).is_integer() else observed == known)
                    if isinstance(known,(int,float)) and not matches:
                        issue(f,f'Conflict: workbook {target}={known}; catalog {target}={observed}')
    if status_workbooks:
        import openpyxl
        for filename in status_workbooks:
            book=openpyxl.load_workbook(filename,read_only=True,data_only=True)
            for sheet_name in ('Status','Experiments','Condition evidence'):
                if sheet_name not in book: continue
                sheet=book[sheet_name]; rows=sheet.iter_rows(values_only=True)
                try: header=[_str(x).lower() for x in next(rows)]
                except StopIteration: continue
                for n,row in enumerate(rows,2):
                    values=dict(zip(header,row))
                    target=_str(values.get('acquisition')) or _str(values.get('collection'))+'/'+_str(values.get('experiment'))
                    f=by_target.get(target)
                    if not f: continue
                    location=f'{Path(filename).name}:{sheet_name}!{n}'
                    rating(f,values.get('user_assessment') or values.get('assessment'),values.get('user_comment') or values.get('comment'),location)
                    sample=samples.get(entries.get(target,{}).get('sample_id'))
                    salt=values.get('salt_mM'.lower(),values.get('nacl_mm',values.get('salt (mm)')))
                    try: salt=float(salt)
                    except (TypeError,ValueError): continue
                    if sample and sample.get('NaCl_mM') is not None and sample['NaCl_mM'] != salt:
                        message=f'Conflict: sample NaCl={sample["NaCl_mM"]}; {location} NaCl={salt}'
                        if message not in f.get('MigrationIssue',''): issue(f,message)
            book.close()
    if overview:
        with open(overview,newline='',encoding='utf-8-sig') as stream:
            for n,row in enumerate(csv.DictReader(stream,delimiter='\t'),2):
                f=by_target.get(row.get('acquisition',''))
                if f: rating(f,row.get('user_assessment'),row.get('user_comment'),f'{Path(overview).name}:{n}')
    return plan


# These tuples were individually reviewed for the pilot. Folder substrings,
# zero concentrations, and absence of DNA/TFE are never scope classifiers.
REVIEWED_PILOT_IDENTITIES = {
    1: ('nmrdata','20250528_15N_hANP32a_50mMNaCl','5'),
    2: ('nmrdata','20251017_15N_HeTDR1_DNA','104'),
    3: ('nmrdata','20251017_15N_HeTDR1_DNA','204'),
    4: ('nmrdata','20251017_15N_HeTDR1_DNA','304'),
    5: ('nmrdata','20251017_15N_HeTDR1_DNA','404'),
    6: ('nmrdata','20251105_15N_heTDR1_NaCl','405'),
    7: ('nmrdata','20251105_15N_heTDR1_NaCl','406'),
    8: ('nmrdata','20251105_15N_heTDR1_NaCl','407'),
    9: ('nmrdata','20251218_15N_heDR1_50mMNaCl','5'),
    10: ('nmrdata','20251218_15N_heDR1_50mMNaCl','6'),
    11: ('nmrdata','20251218_15N_heDR1_50mMNaCl','7'),
}
REVIEWED_PILOT_STUDIES = {(source,folder+'/'+expno):('DNA' if row in (2,3,4,5) else 'Electrostatics')
    for row,(source,folder,expno) in REVIEWED_PILOT_IDENTITIES.items()}
_PILOT_PATCH_FIELDS = frozenset(('Study','AnalysisStatus','SourceEvidence'))
_GUARD_FIELDS = ('SourceID','Folder','Expno','AcquisitionKey','ObservedSourceID','ObservedFolder','ObservedExpno',
                 'ReconciliationReason','MigrationIssue','Study','AnalysisStatus','SourceEvidence')


def _exact_identity(fields):
    identity=tuple(_str(fields.get(name)) for name in ('SourceID','Folder','Expno'))
    return identity if all(identity) else None


def _append_evidence(existing,evidence):
    existing='' if existing is None else str(existing)
    if evidence in existing.split(' | '): return existing
    return ' | '.join(value for value in (existing,evidence) if value)


def _study_evidence(source,target):
    return 'Reviewed pilot study: '+source+':'+target


def _analysis_identity_safe(fields):
    identity=_exact_identity(fields)
    if not identity or fields.get('ReconciliationReason') or fields.get('MigrationIssue'): return False
    source,folder,expno=identity
    key=_str(fields.get('AcquisitionKey'))
    if key and key != source+':'+folder+'/'+expno: return False
    observed=tuple(_str(fields.get('Observed'+name)) for name in ('SourceID','Folder','Expno'))
    if any(observed) and observed != identity: return False
    return True


def _catalog_analysis(catalog,source_id):
    """Positive exact acquisition evidence, bound to an explicitly known source.

    This catalog has no source column; its reviewed origin is nmrdata. A caller
    importing a catalog from another source must provide its reviewed source ID.
    Duplicate catalog targets, derived records and raw copies never promote state.
    """
    if not _str(source_id): raise ValueError('Catalog source ID is required')
    seen={}; duplicates=set()
    with open(catalog,newline='',encoding='utf-8-sig') as stream:
        for row in csv.DictReader(stream):
            target=_str(row.get('relative_path'))
            declared=_str(row.get('SourceID') or row.get('source_id'))
            if declared and declared != source_id: continue
            folder,separator,expno=target.rpartition('/')
            if not separator or not folder or not expno.isdigit() or any(part in ('.','..') or '211@950' in part for part in folder.split('/')) or folder.startswith('/'):
                continue
            identity=(source_id,folder,expno)
            if identity in seen: duplicates.add(identity)
            seen[identity]=row
    positive={}
    for identity,row in seen.items():
        if identity in duplicates or _str(row.get('status')).casefold() not in ('analysed','analyzed'):continue
        if _str(row.get('primary_record_status')).casefold() != 'primary raw acquisition':continue
        if _str(row.get('related_primary_path')):continue
        target=identity[1]+'/'+identity[2]
        positive[identity]='Analysis catalog: '+str(Path(catalog).resolve())+'; exact acquisition '+source_id+':'+target+'; status='+_str(row['status'])
    return positive,duplicates


def reviewed_pilot_patch(records,catalog,source_id='nmrdata',expected_identities=None):
    """Return a reviewable, blank-only patch for the existing eleven pilot rows.

    No rows are created. IDs AND exact source tuples are checked; unrelated rows
    and nonblank human decisions are preserved. This function has no API writes.
    """
    if source_id != 'nmrdata':raise ValueError('The reviewed pilot belongs to source nmrdata')
    expected=REVIEWED_PILOT_IDENTITIES if expected_identities is None else expected_identities
    if any(tuple(identity) not in REVIEWED_PILOT_IDENTITIES.values() for identity in expected.values()):
        raise ValueError('Requested identity is outside the reviewed pilot')
    analyzed,duplicate_catalog=_catalog_analysis(catalog,source_id)
    by_id={};by_identity={}
    for record in records:
        by_id.setdefault(record['id'],[]).append(record)
        identity=_exact_identity(record['fields'])
        if identity:by_identity.setdefault(identity,[]).append(record)
    patches=[];guards={};issues=[]
    for row_id,identity in expected.items():
        identity=tuple(identity);candidates=by_id.get(row_id,[])
        if len(candidates)!=1 or _exact_identity(candidates[0]['fields'])!=identity or len(by_identity.get(identity,[]))!=1:
            issues.append('Pilot row '+str(row_id)+' has missing, duplicate, or inconsistent exact identity');continue
        fields=candidates[0]['fields'];source,folder,expno=identity;target=folder+'/'+expno
        key=_str(fields.get('AcquisitionKey'))
        observed=tuple(_str(fields.get('Observed'+name)) for name in ('SourceID','Folder','Expno'))
        if (key and key != source+':'+target) or (any(observed) and observed!=identity):
            issues.append('Pilot row '+str(row_id)+' has inconsistent collector identity');continue
        values={};provenance=fields.get('SourceEvidence') or ''
        if not _str(fields.get('Study')):
            values['Study']=REVIEWED_PILOT_STUDIES[(source,target)]
            provenance=_append_evidence(provenance,_study_evidence(source,target))
        if identity in duplicate_catalog:
            issues.append('Pilot row '+str(row_id)+' has ambiguous duplicate catalog evidence')
        elif not _str(fields.get('AnalysisStatus')) and _analysis_identity_safe(fields) and identity in analyzed:
            values['AnalysisStatus']='Analyzed';provenance=_append_evidence(provenance,analyzed[identity])
        if values:
            if provenance != (fields.get('SourceEvidence') or ''):values['SourceEvidence']=provenance
            patches.append({'id':row_id,'fields':values})
            guards[str(row_id)]={name:fields.get(name) for name in _GUARD_FIELDS}
    return {'scope':'reviewed_pilot_fields','source_id':source_id,'catalog':str(Path(catalog).resolve()),
            'catalog_sha256':hashlib.sha256(Path(catalog).read_bytes()).hexdigest(),
            'records':patches,'preconditions':guards,'issues':issues}


def apply_reviewed_pilot_patch(plan,client,dry_run=True):
    """Apply only the reviewed blank-field patch after rechecking every guard."""
    if plan.get('scope')!='reviewed_pilot_fields' or plan.get('source_id')!='nmrdata':
        raise ValueError('Not a reviewed pilot field patch')
    if plan.get('issues'):raise ValueError('Pilot patch has unresolved identity/evidence issues')
    if hashlib.sha256(Path(plan['catalog']).read_bytes()).hexdigest()!=plan.get('catalog_sha256'):
        raise ValueError('Catalog changed after pilot patch review')
    positive,_=_catalog_analysis(plan['catalog'],'nmrdata')
    current=client.records('Experiments');by_id={}
    for row in current:by_id.setdefault(row['id'],[]).append(row)
    for patch in plan['records']:
        if set(patch['fields'])-_PILOT_PATCH_FIELDS:raise ValueError('Pilot patch contains an unauthorized field')
        guard=plan['preconditions'].get(str(patch['id']))
        matching=by_id.get(patch['id'],[])
        if guard is None or len(matching)!=1 or any(matching[0]['fields'].get(name)!=value for name,value in guard.items()):
            raise ValueError('Pilot row changed after patch review: '+str(patch['id']))
        identity=_exact_identity(matching[0]['fields'])
        if identity not in REVIEWED_PILOT_IDENTITIES.values() or sum(_exact_identity(row['fields'])==identity for row in current)!=1:
            raise ValueError('Pilot identity is no longer uniquely reviewed')
        if 'Study' in patch['fields'] and patch['fields']['Study'] != REVIEWED_PILOT_STUDIES[(identity[0],identity[1]+'/'+identity[2])]:
            raise ValueError('Study value differs from reviewed pilot membership')
        if 'AnalysisStatus' in patch['fields'] and (patch['fields']['AnalysisStatus']!='Analyzed' or identity not in positive or not _analysis_identity_safe(matching[0]['fields'])):
            raise ValueError('Analysis status lacks exact positive current catalog evidence')
        for name,value in patch['fields'].items():
            if name!='SourceEvidence' and _str(matching[0]['fields'].get(name)):
                raise ValueError('Refusing to overwrite a nonblank human field: '+name)
    if not dry_run and plan['records']:client.update_records('Experiments',plan['records'])
    return {'records':plan['records'],'dry_run':dry_run}
