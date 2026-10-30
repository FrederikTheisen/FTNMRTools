"""Conservative pilot migration. All inputs are read-only; application is opt-in."""
import csv
import re
from datetime import datetime, timezone
from pathlib import Path

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
        for n,row in enumerate(book[sheet].iter_rows(min_row=2, values_only=True),2):
            folder,expno=_str(row[3]),_str(row[4]); target=folder+'/'+expno
            if target not in targets: continue
            if target in found_targets:
                prior=next(x for x in experiments if x['fields']['Folder']+'/'+x['fields']['Expno']==target)
                prior['fields']['MigrationIssue'] += '; duplicate workbook target requires review'
                prior['fields']['MigrationIssue'] += f'; {sheet}!{n}: type={_str(row[8])}, status={_str(row[11])}, comment={_str(row[12])}'
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
            instrument=_str(row[1]); field={'600_1300':600,'600Cryo':600,'600oxbru':600,'700':700,'850':850,'950':950}.get(instrument)
            fields={'SourceID':source_id,'Folder':folder,'Expno':expno,'Instrument':instrument,'Field_MHz':field,
                'ExperimentType':_str(row[8]),'AcquiredAt':_date(row[0]),
                'AcquisitionDateText':_str(row[0]) if row[0] is not None and not isinstance(row[0],datetime) else '',
                'AcquisitionStatus':status,
                'Temperature_C':row[7] if isinstance(row[7],(int,float)) else None,'PulseProgram':_str(row[9]),
                'Assessment':decision,'AssessmentReason':reason,'Notes':notes,
                'SourceEvidence':f'{sheet}!{n}' + ('; condition: '+_str(row[6]) if row[6] is not None else ''),
                'MigrationKey':f'workbook:{sheet}:{n}','MigrationIssue':'; '.join(issue)}
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
        if not dry_run: client.add('Experiments',fields)
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
        for n,row in enumerate(book[sheet].iter_rows(min_row=2,values_only=True),2):
            folder,expno=_str(row[3]),_str(row[4])
            if folder and expno: targets.add(folder+'/'+expno)
            elif any(x is not None and x != '' for x in row[:13]) and (folder or row[8] or row[11] or row[12]):
                planned.append({'fields':{'SourceID':source_id,'Folder':folder,'Expno':expno,
                    'AcquiredAt':_date(row[0]),
                    'AcquisitionDateText':_str(row[0]) if row[0] is not None and not isinstance(row[0],datetime) else '',
                    'ExperimentType':_str(row[8]),'AcquisitionStatus':_str(row[11]),
                    'Notes':_str(row[12]),'MigrationKey':f'workbook:{sheet}:{n}',
                    'MigrationIssue':'Planned or incomplete record; identity and sample unresolved',
                    'SourceEvidence':f'{sheet}!{n}: '+repr(tuple(row[:13]))},'sample_id':None})
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

def reconcile_sources(plan, catalog=None, status_workbooks=(), overview=None):
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
