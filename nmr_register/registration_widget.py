"""Document-embedded registration form and independently testable action builders.

Embed ``registration_html('sample')`` or ``registration_html('experiment')`` in
Grist Widget Creator, with full document access. No external service is used.
The SDK calls are documented at https://support.getgrist.com/code/interfaces/grist_plugin_api.GristDocAPI/ .
"""
from datetime import date, datetime, timezone
import json
import math
import re
from nmr_register.storage_widget import RACK_CSS, RACK_JS
from nmr_register.experiment_taxonomy import experiment_choices_by_category

NUMERIC_RULES = {'pH': [0,14], 'NaCl_mM':[0,None], 'DNA_uM':[0,None],
                 'TFE_pct':[0,100], 'Protein_uM':[0,None]}
TEXT_FIELDS = ('Protein','Construct','Isotope','Buffer','OtherConditions')
UNITS = {'pH':'pH','NaCl_mM':'mMNaCl','DNA_uM':'uMDNA','TFE_pct':'pctTFE','Protein_uM':'uMprotein'}


SAMPLE_STATUSES = ('Usable','Suspected degradation','Degraded','Dead')


def sample_status(value):
    """Legacy values remain readable; new records use the shared tube statuses."""
    if value is None or value in ('','Available','Unknown'):
        return 'Usable'
    if value in ('Unavailable','Retired'):
        return 'Dead'
    if value not in SAMPLE_STATUSES:
        raise ValueError('Invalid sample status')
    return value


def number(value,name,minimum=None,maximum=None):
    if value is None or str(value).strip() == '': return None
    if isinstance(value,bool) or not re.fullmatch(r'[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?',str(value).strip()): raise ValueError(name+' must be a number')
    try: result=float(value)
    except (TypeError,ValueError): raise ValueError(name+' must be a number') from None
    if not math.isfinite(result) or (minimum is not None and result<minimum) or (maximum is not None and result>maximum):
        raise ValueError(name+' is outside its allowed range')
    return result


def prepared(value,required=False):
    value=str(value or '').strip()
    if not value:
        if required: raise ValueError('Preparation date or month is required')
        return None,'',''
    if re.fullmatch(r'\d{4}-\d{2}',value):
        try: date.fromisoformat(value+'-01')
        except ValueError: raise ValueError('Invalid date or month') from None
        return None,value,value.replace('-','')
    try: parsed=date.fromisoformat(value)
    except ValueError: raise ValueError('Use YYYY-MM-DD or YYYY-MM') from None
    if not re.fullmatch(r'\d{4}-\d{2}-\d{2}',value): raise ValueError('Use YYYY-MM-DD or YYYY-MM')
    return datetime(parsed.year,parsed.month,parsed.day,tzinfo=timezone.utc).timestamp(),'',value.replace('-','')


def token(value):
    return re.sub(r'[^A-Za-z0-9.+-]+','_',str(value or '').strip()).strip('_')


LABELS = ('15N','13C','2H')


def labelling_state(value):
    """Recognize standard combinations while retaining custom text verbatim."""
    text=str(value or '').strip()
    parts=re.split(r'[+/,;\s]+',text) if text else []
    if all(part in LABELS for part in parts):
        return {'labels':[label for label in LABELS if label in parts],'custom':''}
    return {'labels':[],'custom':text}


def labelling(values,inherited=''):
    if any('Label'+label in values for label in LABELS) or 'IsotopeCustom' in values:
        custom=str(values.get('IsotopeCustom') or '').strip()
        if custom:return custom
        return '+'.join(label for label in LABELS if values.get('Label'+label))
    state=labelling_state(values.get('Isotope',inherited))
    return state['custom'] or '+'.join(state['labels'])


def condition_label(values,composition):
    explicit=token(values.get('DistinguishingCondition'))
    if explicit:return explicit
    parts=[token(values.get('Construct'))] if values.get('Construct') else []
    for name in ('NaCl_mM','DNA_uM','TFE_pct','Protein_uM','pH'):
        value=composition.get(name)
        if value is not None:
            text=format(value,'.15g')
            parts.append('pH'+text if name=='pH' else text+UNITS[name])
    return '_'.join(part for part in parts if part)


def _rows(tables,name): return tables.get(name,[])
def _next(tables,name): return max([int(row['id']) for row in _rows(tables,name)]+[0])+1


def sample_actions(values,tables):
    """Build one atomic AddRecord transaction, preserving blanks and literal zero.

    ``tables`` contains raw rows (including integer id) fetched immediately before
    submission. Explicit row IDs also make simultaneous stale form submissions
    collide safely rather than create duplicate visible sample IDs.
    """
    datevalue,month,stamp=prepared(values.get('Prepared'))
    selected=next((row for row in _rows(tables,'Series') if str(row['id'])==str(values.get('Series'))),{})
    common={name:str(values.get(name) or '').strip() for name in TEXT_FIELDS}
    common['Isotope']=labelling(values,selected.get('Isotope',''))
    common.update({name:number(values.get(name),name,*bounds) for name,bounds in NUMERIC_RULES.items()})
    effective={name:common[name] if common[name] is not None else selected.get(name) for name in NUMERIC_RULES}
    condition=condition_label(values,effective)
    if not common['Protein'] or not common['Isotope']:
        raise ValueError('Enter a protein and select labelling, or enter a custom labelling description')
    series=values.get('Series') or ''
    actions=[]
    if series == 'new':
        series_id=str(values.get('SeriesID') or '_'.join(part for part in [stamp,token(common['Isotope']),token(common['Protein']),condition,'series'] if part)).strip()
        if any(row.get('SeriesID')==series_id for row in _rows(tables,'Series')):
            raise ValueError('Series ID already exists; select that series')
        series_ref=_next(tables,'Series')
        actions.append(['AddRecord','Series',series_ref,{'SeriesID':series_id,**common}])
    else:
        series_ref=int(series) if series else 0
        if series_ref and not any(row['id']==series_ref for row in _rows(tables,'Series')):
            raise ValueError('Selected series no longer exists')
    varying=[]
    if values.get('mode') == 'bulk':
        field=values.get('VaryField','DNA_uM')
        if field not in NUMERIC_RULES: raise ValueError('Unknown varying condition')
        raw=re.split(r'[,;\s]+',str(values.get('VaryValues') or '').strip())
        if not raw or raw==[''] or len(raw)>96: raise ValueError('Enter between 1 and 96 varying values')
        varying=[number(value,field,*NUMERIC_RULES[field]) for value in raw]
        if len(varying)!=len(set(varying)): raise ValueError('Repeated varying values; add replicates separately with distinct counters')
    else: field=None;varying=[None]
    slot_text=str(values.get('StorageSlot') or '').strip().upper()
    if slot_text:
        match=re.fullmatch(r'([A-H])(1[0-2]|[1-9])',slot_text)
        if not match: raise ValueError('Storage position must be a letter A-H and numbered row 1-12')
        row,column=match.group(1),int(match.group(2))
    else:
        # Keep the action builder compatible with older callers while the form
        # presents the transposed rack view (numbered rows, letter columns).
        row=str(values.get('StorageRow') or '').strip().upper()
        raw_column=number(values.get('StorageColumn'),'Storage column',1,12)
        if bool(row) != (raw_column is not None) or (row and not re.fullmatch('[A-H]',row)) or (raw_column is not None and not raw_column.is_integer()):
            raise ValueError('Storage position needs a letter A-H and numbered row 1-12')
        column=int(raw_column) if raw_column is not None else None
    availability=sample_status(values.get('Available'))
    ids={row.get('SampleID') for row in _rows(tables,'Samples')}
    occupied=set()
    for sample in _rows(tables,'Samples'):
        try: slot_column=float(sample.get('StorageColumn'))
        except (TypeError,ValueError): continue
        occupied.add((str(sample.get('StorageRow') or '').strip().upper(),slot_column))
    created=[];next_id=_next(tables,'Samples');planned_slots=[]
    for index,value in enumerate(varying):
        physical={**effective,**({field:value} if field else {})}
        suffix=condition_label(values,physical)
        if field and values.get('DistinguishingCondition'):
            suffix+='_'+('pH'+format(value,'.15g') if field=='pH' else format(value,'.15g')+UNITS[field])
        base='_'.join(part for part in [stamp,token(common['Isotope']),token(common['Protein']),suffix] if part)
        counter=1
        while base+'_'+str(counter) in ids:counter+=1
        sample_id=base+'_'+str(counter)
        ids.add(sample_id);created.append(sample_id)
        slot_row,slot_col='',None
        if row:
            slot='ABCDEFGH'.index(row)*12+column-1+index
            if slot>=96: raise ValueError('Bulk storage slots exceed the rack')
            slot_row='ABCDEFGH'[slot//12];slot_col=slot%12+1
            if (slot_row,slot_col) in occupied: raise ValueError('A selected storage slot is already occupied')
            planned_slots.append(slot_row+str(slot_col))
        else: planned_slots.append('')
        fields={'SampleID':sample_id,'PreparedAt':datevalue,'PreparedDateText':month,
            'Series':series_ref,**common,'DistinguishingCondition':suffix,'SampleNumber':str(counter),
            'StorageRow':slot_row,'StorageColumn':slot_col,'Available':availability}
        if field: fields[field]=value
        actions.append(['AddRecord','Samples',next_id+index,fields])
    return {'actions':actions,'sample_ids':created,'storage_slots':planned_slots}


def experiment_actions(values,tables,now=None):
    series=int(values.get('Series') or 0);sample=int(values.get('Sample') or 0)
    samples={row['id']:row for row in _rows(tables,'Samples')}
    if series and not any(row['id']==series for row in _rows(tables,'Series')): raise ValueError('Selected series no longer exists')
    if sample and (sample not in samples or (series and samples[sample].get('Series')!=series)):
        raise ValueError('Sample is not in the selected series')
    if sample and sample_status(samples[sample].get('Available')) == 'Dead' and not values.get('IncludeUnavailable'):
        raise ValueError('Sample is dead; explicitly include dead samples to select it')
    if sample and values.get('ProteinFilter') and samples[sample].get('Protein') != values['ProteinFilter']:
        raise ValueError('Sample does not match the selected protein')
    study=str(values.get('Study') or '').strip()
    if study not in ('', 'Electrostatics', 'DNA', 'TFE', 'Other'): raise ValueError('Invalid study')
    analysis=str(values.get('AnalysisStatus') or '').strip()
    if analysis not in ('', 'Not started', 'In progress', 'Analyzed'): raise ValueError('Invalid analysis status')
    status=values.get('AcquisitionStatus') or 'Planned'
    if status not in ('Planned','Dummy'): raise ValueError('Invalid initial acquisition status')
    acquired,month,_=prepared(values.get('Acquired'))
    source=str(values.get('SourceID') or '').strip();folder=str(values.get('Folder') or '').strip();expno=str(values.get('Expno') or '').strip()
    if folder and (folder.startswith('/') or any(p in ('..','.') or '211@950' in p for p in folder.split('/'))):
        raise ValueError('Invalid or protected acquisition folder')
    if expno and not expno.isdigit(): raise ValueError('Experiment number must contain digits')
    if source and folder and expno:
        for row in _rows(tables,'Experiments'):
            if any(tuple(str(row.get(prefix+k) or '') for k in ('SourceID','Folder','Expno'))==(source,folder,expno) for prefix in ('','Observed')):
                raise ValueError('Acquisition already exists; edit its existing record')
    fields={name:str(values.get(name) or '').strip() for name in ('Instrument','ExperimentType','PulseProgram','DelayList_s','Notes')}
    fields.update(Series=series,Sample=sample,SourceID=source,Folder=folder,Expno=expno,
        Study=study,AnalysisStatus=analysis,AcquisitionStatus=status,AcquiredAt=acquired,AcquisitionDateText=month,RegisteredAt=now if now is not None else datetime.now(timezone.utc).timestamp(),
        Field_MHz=number(values.get('Field_MHz'),'Field (MHz)',0.000001),
        Temperature_C=number(values.get('Temperature_C'),'Temperature (C)',-273.15),
        D1_s=number(values.get('D1_s'),'D1 (s)',0),
        SpinLock_Hz=number(values.get('SpinLock_Hz'),'Spin-lock field (Hz)',0),
        SpinLockOffset_Hz=number(values.get('SpinLockOffset_Hz'),'Spin-lock offset (Hz)'))
    if not sample: fields.update({'Snapshot'+name:None for name in NUMERIC_RULES})
    return {'actions':[['AddRecord','Experiments',_next(tables,'Experiments'),fields]],'sample_ids':[]}


JS_LOGIC = r'''
const NMRRegistration = (() => {
 const rules = __RULES__, units = __UNITS__, textFields = __TEXT__,labels=['15N','13C','2H'];
 const sampleStatuses=['Usable','Suspected degradation','Degraded','Dead'];
 function sampleStatus(value){if(value==null||['','Available','Unknown'].includes(value))return 'Usable';if(['Unavailable','Retired'].includes(value))return 'Dead';if(!sampleStatuses.includes(value))throw Error('Invalid sample status');return value;}
 function selectableSample(row,includeDead=false){try{return sampleStatus(row.Available)!=='Dead'||includeDead;}catch(error){return false;}}

 function number(value,name,min=null,max=null) {
  if(value==null || String(value).trim()==='')return null;
  if(typeof value==='boolean'||!/^[-+]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][-+]?[0-9]+)?$/.test(String(value).trim()))throw Error(name+' must be a number');
  const n=Number(value);if(!Number.isFinite(n)||(min!=null&&n<min)||(max!=null&&n>max))throw Error(name+' is outside its allowed range');return n;
 }
 function prepared(raw,required=false){
  const value=String(raw||'').trim();if(!value){if(required)throw Error('Preparation date or month is required');return [null,'',''];}
  if(!/^\d{4}-\d{2}(?:-\d{2})?$/.test(value))throw Error('Use YYYY-MM-DD or YYYY-MM');
  const full=value.length===7?value+'-01':value,date=new Date(full+'T00:00:00Z');
  if(!Number.isFinite(date.getTime())||date.toISOString().slice(0,10)!==full)throw Error('Invalid date or month');
  return value.length===7?[null,value,value.replaceAll('-','')]:[date.getTime()/1000,'',value.replaceAll('-','')];
 }
 const token=v=>String(v||'').trim().replace(/[^A-Za-z0-9.+-]+/g,'_').replace(/^_+|_+$/g,'');
 const rows=(t,n)=>t[n]||[], next=(t,n)=>Math.max(0,...rows(t,n).map(r=>r.id))+1;
 const proteinSuggestions=t=>[...new Set([...rows(t,'Samples').map(r=>r.Protein),...rows(t,'Series').map(r=>r.Protein),...rows(t,'Experiments').map(r=>r.SnapshotProtein)].filter(v=>typeof v==='string'&&v.trim()))].sort();
 function proteinChoiceActions(plan,metadataTables,metadataColumns){
  const proteins=new Set();
  for(const [action,table,,fields]of plan.actions||[])if(action==='AddRecord'&&['Series','Samples'].includes(table)&&typeof fields?.Protein==='string'&&fields.Protein.trim()){
   proteins.add(fields.Protein);
  }
  if(!proteins.size)return [];
  const metadataRows=t=>Array.isArray(t)?t.map(r=>r.fields?{...r.fields,id:r.id}:r):(t.id||[]).map((id,i)=>Object.fromEntries(Object.keys(t).map(k=>[k,t[k][i]])));
  const tables=metadataRows(metadataTables),columns=metadataRows(metadataColumns),actions=[];
  for(const name of ['Samples','Series']){
   const table=tables.filter(t=>t.tableId===name);if(table.length!==1||!Number.isInteger(table[0].id)||table[0].id<=0)throw Error('Cannot identify '+name+' protein choices; refresh before saving.');
   const matches=columns.filter(c=>c.parentId===table[0].id&&c.colId==='Protein');
   if(matches.length!==1||!Number.isInteger(matches[0].id)||matches[0].id<=0||matches[0].type!=='Choice'||matches[0].isFormula)throw Error(name+' Protein must be the editable Choice column; refresh before saving.');
   const column=matches[0];let options;try{options=JSON.parse(column.widgetOptions||'{}');}catch(error){throw Error(name+' protein choice settings are invalid; nothing was saved.');}
   if(!options||typeof options!=='object'||Array.isArray(options)||('choices'in options&&(!Array.isArray(options.choices)||options.choices.some(v=>typeof v!=='string'))))throw Error(name+' protein choice settings are invalid; nothing was saved.');
   const choices=[...(options.choices||[])],missing=[...proteins].filter(value=>!choices.includes(value));
   if(missing.length)actions.push(['UpdateRecord','_grist_Tables_column',column.id,{widgetOptions:JSON.stringify({...options,choices:[...choices,...missing]})}]);
  }
  return actions;
 }

 function labellingState(value){
  const text=String(value||'').trim(),parts=text?text.split(/[+/,;\s]+/):[];
  return parts.every(part=>labels.includes(part))?{labels:labels.filter(label=>parts.includes(label)),custom:''}:{labels:[],custom:text};
 }
 function labelling(v,inherited=''){
  if(labels.some(label=>Object.hasOwn(v,'Label'+label))||Object.hasOwn(v,'IsotopeCustom')){
   const custom=String(v.IsotopeCustom||'').trim();return custom||labels.filter(label=>v['Label'+label]).join('+');
  }
  const state=labellingState(Object.hasOwn(v,'Isotope')?v.Isotope:inherited);return state.custom||state.labels.join('+');
 }
 function conditionLabel(v,composition){
  const explicit=token(v.DistinguishingCondition);if(explicit)return explicit;
  const parts=v.Construct?[token(v.Construct)]:[];
  for(const name of ['NaCl_mM','DNA_uM','TFE_pct','Protein_uM','pH'])if(composition[name]!=null)parts.push(name==='pH'?'pH'+composition[name]:composition[name]+units[name]);
  return parts.filter(Boolean).join('_');
 }
 function sampleActions(v,t){
  const [datevalue,month,stamp]=prepared(v.Prepared), common={};
  const selected=rows(t,'Series').find(row=>String(row.id)===String(v.Series))||{};
  textFields.forEach(k=>common[k]=String(v[k]||'').trim());common.Isotope=labelling(v,selected.Isotope||'');Object.entries(rules).forEach(([k,b])=>common[k]=number(v[k],k,...b));
  const effective=Object.fromEntries(Object.keys(rules).map(k=>[k,common[k]??selected[k]])),condition=conditionLabel(v,effective);
  if(!common.Protein||!common.Isotope)throw Error('Enter a protein and select labelling, or enter a custom labelling description');
  let seriesRef=0;const actions=[];
  if(v.Series==='new'){
   const seriesID=String(v.SeriesID||[stamp,token(common.Isotope),token(common.Protein),condition,'series'].filter(Boolean).join('_')).trim();
   if(rows(t,'Series').some(r=>r.SeriesID===seriesID))throw Error('Series ID already exists; select that series');
   seriesRef=next(t,'Series');actions.push(['AddRecord','Series',seriesRef,{SeriesID:seriesID,...common}]);
  }else{seriesRef=Number(v.Series||0);if(seriesRef&&!rows(t,'Series').some(r=>r.id===seriesRef))throw Error('Selected series no longer exists');}
  let varying=[null],field=null;
  if(v.mode==='bulk'){
   field=v.VaryField||'DNA_uM';if(!Object.hasOwn(rules,field))throw Error('Unknown varying condition');
   const raw=String(v.VaryValues||'').trim().split(/[,;\s]+/);if(!raw[0]||raw.length>96)throw Error('Enter between 1 and 96 varying values');
   varying=raw.map(x=>number(x,field,...rules[field]));if(new Set(varying).size!==varying.length)throw Error('Repeated varying values; add replicates separately with distinct counters');
  }
  let row='',column=null;
  const position=String(v.StorageSlot||'').trim().toUpperCase();
  if(position){const match=/^([A-H])(1[0-2]|[1-9])$/.exec(position);if(!match)throw Error('Storage position must be a letter A-H and numbered row 1-12');row=match[1];column=Number(match[2]);}
  else{row=String(v.StorageRow||'').trim().toUpperCase();column=number(v.StorageColumn,'Storage column',1,12);
   if(Boolean(row)!==(column!==null)||(row&&!/^[A-H]$/.test(row))||(column!==null&&!Number.isInteger(column)))throw Error('Storage position needs a letter A-H and numbered row 1-12');}
  const available=sampleStatus(v.Available);
  const ids=new Set(rows(t,'Samples').map(r=>r.SampleID)),created=[],plannedSlots=[],firstID=next(t,'Samples');
  varying.forEach((value,i)=>{
   const physical={...effective,...(field?{[field]:value}:{})};let suffix=conditionLabel(v,physical);
   if(field&&v.DistinguishingCondition)suffix+='_'+(field==='pH'?'pH'+value:value+units[field]);
   const base=[stamp,token(common.Isotope),token(common.Protein),suffix].filter(Boolean).join('_');let counter=1;while(ids.has(base+'_'+counter))counter++;
   const sampleID=base+'_'+counter;ids.add(sampleID);created.push(sampleID);
   let slotRow='',slotCol=null;
   if(row){const slot='ABCDEFGH'.indexOf(row)*12+column-1+i;if(slot>=96)throw Error('Bulk storage slots exceed the rack');slotRow='ABCDEFGH'[Math.floor(slot/12)];slotCol=slot%12+1;
    if(rows(t,'Samples').some(s=>String(s.StorageRow||'').trim().toUpperCase()===slotRow&&Number(s.StorageColumn)===slotCol))throw Error('A selected storage slot is already occupied');}
   plannedSlots.push(slotRow?slotRow+slotCol:'');
   const fields={SampleID:sampleID,PreparedAt:datevalue,PreparedDateText:month,Series:seriesRef,...common,DistinguishingCondition:suffix,SampleNumber:String(counter),StorageRow:slotRow,StorageColumn:slotCol,Available:available};
   if(field)fields[field]=value;actions.push(['AddRecord','Samples',firstID+i,fields]);
  });return {actions,sample_ids:created,storage_slots:plannedSlots};
 }
 function experimentActions(v,t,now=Date.now()/1000){
  const series=Number(v.Series||0),sample=Number(v.Sample||0),selected=rows(t,'Samples').find(r=>r.id===sample);
  if(series&&!rows(t,'Series').some(r=>r.id===series))throw Error('Selected series no longer exists');
  if(sample&&(!selected||(series&&selected.Series!==series)))throw Error('Sample is not in the selected series');
  if(sample&&sampleStatus(selected.Available)==='Dead'&&!v.IncludeUnavailable)throw Error('Sample is dead; explicitly include dead samples to select it');
  if(sample&&v.ProteinFilter&&selected.Protein!==v.ProteinFilter)throw Error('Sample does not match the selected protein');
  const study=String(v.Study||'').trim();if(!['','Electrostatics','DNA','TFE','Other'].includes(study))throw Error('Invalid study');
  const analysis=String(v.AnalysisStatus||'').trim();if(!['','Not started','In progress','Analyzed'].includes(analysis))throw Error('Invalid analysis status');
  const status=v.AcquisitionStatus||'Planned';if(!['Planned','Dummy'].includes(status))throw Error('Invalid initial acquisition status');
  const [acquired,month]=prepared(v.Acquired),source=String(v.SourceID||'').trim(),folder=String(v.Folder||'').trim(),expno=String(v.Expno||'').trim();
  if(folder&&(folder.startsWith('/')||folder.split('/').some(p=>p==='..'||p==='.'||p.includes('211@950'))))throw Error('Invalid or protected acquisition folder');
  if(expno&&!/^\d+$/.test(expno))throw Error('Experiment number must contain digits');
  if(source&&folder&&expno&&rows(t,'Experiments').some(r=>['','Observed'].some(p=>String(r[p+'SourceID']||'')===source&&String(r[p+'Folder']||'')===folder&&String(r[p+'Expno']||'')===expno)))throw Error('Acquisition already exists; edit its existing record');
  const fields={};['Instrument','ExperimentType','PulseProgram','DelayList_s','Notes'].forEach(k=>fields[k]=String(v[k]||'').trim());
  Object.assign(fields,{Series:series,Sample:sample,SourceID:source,Folder:folder,Expno:expno,Study:study,AnalysisStatus:analysis,AcquisitionStatus:status,AcquiredAt:acquired,AcquisitionDateText:month,RegisteredAt:now,
   Field_MHz:number(v.Field_MHz,'Field (MHz)',.000001),Temperature_C:number(v.Temperature_C,'Temperature (C)',-273.15),D1_s:number(v.D1_s,'D1 (s)',0),DelayList_s:String(v.DelayList_s||'').trim(),SpinLock_Hz:number(v.SpinLock_Hz,'Spin-lock field (Hz)',0),SpinLockOffset_Hz:number(v.SpinLockOffset_Hz,'Spin-lock offset (Hz)')});
  if(!sample)Object.keys(rules).forEach(k=>fields['Snapshot'+k]=null);
  return {actions:[['AddRecord','Experiments',next(t,'Experiments'),fields]],sample_ids:[]};
 }
 return {number,prepared,labellingState,labelling,sampleStatus,selectableSample,proteinSuggestions,proteinChoiceActions,sampleActions,experimentActions};
})();
'''.replace('__RULES__',json.dumps(NUMERIC_RULES)).replace('__UNITS__',json.dumps(UNITS)).replace('__TEXT__',json.dumps(TEXT_FIELDS))

HTML = r'''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<script src="grist-plugin-api.js"></script><style>
body{margin:0;padding:12px;font:13px var(--grist-theme-font-family,system-ui,sans-serif);color:var(--grist-theme-text,#262633);background:var(--grist-theme-widget-bg,#fff)}
form{display:grid;gap:12px;max-width:800px}h3,legend,summary{font-size:13px;font-weight:600}h3{margin:0 0 8px}section,details{padding:12px;border:1px solid var(--grist-theme-widget-border,#d9d9d9);border-radius:6px}summary{cursor:pointer}details[open]>summary{margin-bottom:12px}fieldset{margin:0;padding:8px;border:1px solid var(--grist-theme-input-border,#d9d9d9);border-radius:4px}small{font-size:11px;line-height:1.4;color:var(--grist-theme-text-light,#929299)}.wide{grid-column:1/-1}.labels{display:flex;flex-wrap:wrap;gap:16px}.labels label{display:flex;align-items:center;gap:6px;font-size:13px}h2{font-size:15px;font-weight:600;margin:0 0 4px}p{margin:0 0 8px;color:var(--grist-theme-text-medium,#494949)}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:8px}label{display:grid;grid-template-rows:minmax(28px,auto) auto;align-content:start;gap:4px;font-size:11px;color:var(--grist-theme-text-medium,#494949)}
input,select,textarea{box-sizing:border-box;min-width:0;width:100%;padding:6px 8px;border:1px solid var(--grist-theme-input-border,#d9d9d9);border-radius:4px;background:var(--grist-theme-widget-bg,#fff);color:var(--grist-theme-text,#262633);font:13px inherit}input:not([type=checkbox]),select{height:32px}textarea{min-height:52px}input[type=checkbox]{width:auto;justify-self:start}
button{min-height:32px;border:0;border-radius:4px;background:var(--grist-theme-control-primary-bg,#16b378);color:var(--grist-theme-control-primary-fg,#fff);padding:4px 12px;font:inherit;cursor:pointer}button:disabled{opacity:.6;cursor:wait}
select[data-tube-status="Suspected degradation"]{background:#FFF1CC;color:#5C4300}select[data-tube-status="Degraded"]{background:#FCE6DF;color:#743125}select[data-tube-status="Dead"]{background:#F8B4B4;color:#7F1D1D}
.secondary{background:var(--grist-theme-widget-bg,#fff);color:var(--grist-theme-control-fg,#16b378);border:1px solid var(--grist-theme-control-primary-bg,#16b378)}
.actions{display:flex;gap:8px;align-items:center}.hidden{display:none!important}#message{white-space:pre-wrap;overflow-wrap:anywhere;line-height:1.5}#message.error{color:var(--grist-theme-text-error,#d0021b)}
/*RACK_CSS*/
.rack-help{margin:8px 0 5px;font-size:11px;line-height:1.4;color:var(--grist-theme-text-light,#686875)}#rackPicker{max-width:560px}#storagePreview{min-height:18px;margin-top:7px;font-size:11px;line-height:1.4;overflow-wrap:anywhere}#storagePreview.error{color:#8d3028;font-weight:600}#storagePreview.ok{color:var(--grist-theme-text-light,#686875)}
</style></head><body><form id="form"><h2 id="title">Add samples</h2><p>Leave unknown values blank. Enter 0 when the value is zero.</p>
<div class="grid" id="top"></div><details id="seriesBlock"><summary>Series (optional)</summary><div class="grid" id="seriesChooser"></div><div class="grid hidden" id="seriesFields"></div></details>
<section id="identityBlock"><h3>1. Identity</h3><div class="grid" id="identity"></div><fieldset><legend>Labelling</legend><div class="labels" id="labels"></div></fieldset><details id="customLabelling"><summary>Custom labelling</summary><div class="grid" id="customLabels"></div></details></section>
<section id="bufferBlock"><h3>2. Buffer and composition</h3><div class="grid" id="buffer"></div></section>
<details id="compositionBlock"><summary>Additional composition (optional)</summary><div class="grid" id="composition"></div></details>
<section id="bulkBlock" class="hidden"><h3>Values for each physical sample</h3><div class="grid" id="bulk"></div></section>
<details id="storageBlock"><summary>Rack position (optional)</summary><p class="rack-help">Choose a numbered row and lettered column. Bulk placement previews consecutive slots and stops at any occupied position.</p><input type="hidden" id="field_StorageSlot" name="StorageSlot" value=""><div id="rackPicker" aria-label="Sample rack position"></div><div id="storagePreview" class="ok" role="status" aria-live="polite"></div><div class="grid" id="storage"></div></details>
<div class="grid hidden" id="experiment"></div>
<div class="actions"><button id="submit" type="submit" disabled>Add sample</button><button id="refresh" class="secondary" type="button">Refresh choices</button></div><div id="message" role="status" aria-live="polite">Loading register…</div></form>
<script>
grist.ready({requiredAccess:'full'});
/*LOGIC*/
/*RACK_JS*/
const byId=id=>document.getElementById(id),controls={},groups={},numericLabels={pH:'pH',NaCl_mM:'NaCl (mM)',DNA_uM:'DNA (µM)',TFE_pct:'TFE (%)',Protein_uM:'Protein (µM)'};
let tables={Series:[],Samples:[],Experiments:[]},busy=false;
function field(group,name,label,kind='text',options=null){
 const wrapper=document.createElement('label');wrapper.textContent=label;const input=document.createElement(kind==='select'?'select':kind==='textarea'?'textarea':'input');
 if(kind!=='select'&&kind!=='textarea')input.type=kind;input.name=name;input.id='field_'+name;
 if(kind==='number')input.step='any';if(options)for(const [value,title]of options){const option=document.createElement('option');option.value=value;option.textContent=title;input.append(option);}
 wrapper.append(input);byId(group).append(wrapper);controls[name]=input;groups[name]=wrapper;return input;
}
field('top','mode','Create','select',[['sample','One sample tube'],['bulk','Several tubes / series'],['experiment','Experiment']]).value=/*MODE*/;
field('seriesChooser','Series','Series','select',[['','No series'],['new','Create a series']]);field('seriesFields','SeriesID','New series label (optional)');
function hint(name,text){const note=document.createElement('small');note.textContent=text;groups[name].append(note);}
field('identity','Protein','Protein').placeholder='Choose an existing name or type a new one';field('identity','Construct','Construct (optional)');field('identity','Prepared','Preparation date (optional)');
hint('Prepared','Leave blank if unknown. For example 2026-10-30, or 2026-10 if only the month is known.');
field('identity','DistinguishingCondition','Sample label (optional)');hint('DistinguishingCondition','Leave blank to use the composition in the generated label.');
for(const label of ['15N','13C','2H'])field('labels','Label'+label,label,'checkbox');
field('customLabels','IsotopeCustom','Custom labelling description');hint('IsotopeCustom','Existing custom descriptions are kept. Clear this entry to use the label choices above.');
function updateLabelling(){const custom=controls.IsotopeCustom.value.trim();for(const label of ['15N','13C','2H'])controls['Label'+label].disabled=Boolean(custom);}
function setLabelling(value){const state=NMRRegistration.labellingState(value);controls.IsotopeCustom.value=state.custom;for(const label of ['15N','13C','2H'])controls['Label'+label].checked=state.labels.includes(label);byId('customLabelling').open=Boolean(state.custom);updateLabelling();}
controls.IsotopeCustom.addEventListener('input',updateLabelling);
field('buffer','Buffer','Buffer');field('buffer','pH','pH','number');field('buffer','NaCl_mM','NaCl (mM)','number');field('buffer','Protein_uM','Protein concentration (µM)','number');
field('composition','DNA_uM','DNA (µM)','number');field('composition','TFE_pct','TFE (%)','number');field('composition','OtherConditions','Composition notes','textarea');groups.OtherConditions.classList.add('wide');
field('bulk','VaryField','Vary one condition','select',Object.entries(numericLabels));field('bulk','VaryValues','Values separated by commas or lines, e.g. 0, 100, 200, 500','textarea');controls.VaryField.value='DNA_uM';
controls.StorageSlot=byId('field_StorageSlot');
field('identity','Available','Tube status','select',[['Usable','Usable'],['Suspected degradation','Suspected degradation'],['Degraded','Degraded'],['Dead','Dead']]);hint('Available','Dead means permanently unavailable: used up, discarded, or irreversibly unusable.');
function updateTubeStatus(){controls.Available.dataset.tubeStatus=NMRRegistration.sampleStatus(controls.Available.value);}
controls.Available.addEventListener('change',updateTubeStatus);updateTubeStatus();
function rackPreviewSlots(){
 const raw=String(controls.StorageSlot.value||'').trim().toUpperCase();if(!raw)return {slots:[],error:''};
 const match=/^([A-H])(1[0-2]|[1-9])$/.exec(raw);if(!match)return {slots:[],error:'Choose a valid rack position.'};
 let count=1;
 if(controls.mode.value==='bulk'){
  const values=String(controls.VaryValues.value||'').trim().split(/[,;\s]+/).filter(Boolean);
  if(!values.length||values.length>96)return {slots:[],error:'Enter 1–96 varying values to preview positions.'};
  count=values.length;
 }
 const first='ABCDEFGH'.indexOf(match[1])*12+Number(match[2])-1,slots=[];
 for(let i=0;i<count;i++){const index=first+i;if(index>=96)return {slots,error:'This batch extends beyond the last rack position.'};slots.push('ABCDEFGH'[Math.floor(index/12)]+String(index%12+1));}
 return {slots,error:''};
}
function refreshRack(){
 const selected=String(controls.StorageSlot.value||'').trim().toUpperCase();
 const result=Rack.mountRack(byId('rackPicker'),{records:tables.Samples,destination:selected,onSelect:(key)=>{controls.StorageSlot.value=key;refreshRack();}});
 const plan=rackPreviewSlots(),occupied=Rack.occupancy(tables.Samples),blocked=plan.slots.find(slot=>(occupied.get(slot)||[]).length>0);
 let text='';if(plan.error)text=plan.error;else if(!plan.slots.length)text='No rack positions will be assigned.';else if(blocked)text='Blocked: '+blocked+' is occupied. Choose a different starting position; the form will not skip occupied slots.';else text=(controls.mode.value==='bulk'?'Bulk positions: ':'Position: ')+plan.slots.join(', ');
 if(result.invalid.length)text+=(text?' ':'')+result.invalid.length+' saved position'+(result.invalid.length===1?' is':'s are')+' incomplete or invalid.';
 const preview=byId('storagePreview');preview.textContent=text;preview.className=(plan.error||blocked)?'error':'ok';
 // Keep validation reachable even when its rack preview is inside closed details.
 byId('submit').disabled=busy;
}
controls.StorageSlot.addEventListener('change',refreshRack);
controls.VaryValues.addEventListener('input',refreshRack);
controls.mode.addEventListener('change',refreshRack);
field('experiment','ProteinFilter','Filter samples by protein','select',[['','All proteins']]);field('experiment','IncludeUnavailable','Include dead samples','checkbox');
field('experiment','Sample','Physical sample','select',[['','Not assigned yet']]);field('experiment','Acquired','Acquisition: YYYY-MM-DD or YYYY-MM');
for(const [name,label]of [['SourceID','Source ID'],['Folder','Acquisition folder'],['Expno','Experiment number'],['Instrument','Instrument']])field('experiment',name,label);
const experimentChoices=/*EXPERIMENT_CHOICES*/;
field('experiment','ExperimentCategory','Experiment type','select',[['','All experiment types'],...Object.keys(experimentChoices).map(value=>[value,value])]);
field('experiment','ExperimentType','Experiment').placeholder='Choose an experiment or enter a custom name';
field('experiment','PulseProgram','Pulse program (optional)');
const experimentList=document.createElement('datalist');experimentList.id='experimentChoices';byId('experiment').append(experimentList);controls.ExperimentType.setAttribute('list',experimentList.id);
function refreshExperimentChoices(){
 const selected=controls.ExperimentCategory.value;
 const choices=selected?experimentChoices[selected]:Object.values(experimentChoices).flat();
 experimentList.replaceChildren();for(const value of choices){const option=document.createElement('option');option.value=value;experimentList.append(option);}
}
controls.ExperimentCategory.addEventListener('change',refreshExperimentChoices);refreshExperimentChoices();
hint('ExperimentType','Custom experiment names are retained. The experiment type is calculated from the saved experiment.');
hint('PulseProgram','Optional raw acquisition pulse name; no pulse is assigned automatically.');
field('experiment','AcquisitionStatus','Entry status','select',[['Planned','Planned'],['Dummy','Dummy entry']]);field('experiment','Study','Study','select',[['','Unclassified'],['Electrostatics','Electrostatics'],['DNA','DNA'],['TFE','TFE'],['Other','Other']]);field('experiment','AnalysisStatus','Analysis status','select',[['','Not recorded'],['Not started','Not started'],['In progress','In progress'],['Analyzed','Analyzed']]);controls.SourceID.value='nmrdata';field('experiment','Field_MHz','Field (MHz)','number');field('experiment','Temperature_C','Temperature (°C)','number');field('experiment','D1_s','D1 (s)','number');field('experiment','DelayList_s','Delay list (s)');field('experiment','SpinLock_Hz','Spin-lock field (Hz)','number');field('experiment','SpinLockOffset_Hz','Spin-lock offset (Hz)','number');field('experiment','Notes','Notes','textarea');
// Protein suggestions preserve existing names. Protocol defaults use only
// human-confirmed pairs; ambiguity or a nonblank destination leaves entries alone.
const suggestions={};
for(const name of ['Protein','Instrument','PulseProgram']){
 const list=document.createElement('datalist');list.id='suggest_'+name;
 byId(name==='Protein'?'identity':'experiment').append(list);controls[name].setAttribute('list',list.id);suggestions[name]=list;
}
function refreshProtocolSuggestions(){
 for(const name of ['Protein','Instrument','PulseProgram']){
  const values=name==='Protein'?NMRRegistration.proteinSuggestions(tables):[...new Set(tables.Experiments.map(r=>r[name]).filter(v=>typeof v==='string'&&v.trim()))].sort();
  const matching=name==='PulseProgram'&&controls.ExperimentType.value.trim()?new Set(tables.Experiments.filter(r=>r.ExperimentType===controls.ExperimentType.value.trim()).map(r=>r.PulseProgram)):null;
  if(matching?.size)values.sort((a,b)=>Number(matching.has(b))-Number(matching.has(a))||a.localeCompare(b));
  suggestions[name].replaceChildren();for(const value of values){const option=document.createElement('option');option.value=value;suggestions[name].append(option);}
 }
}
const automaticDefaults={};
const manuallyEdited=new Set();
for(const name of ['ExperimentType','Field_MHz']){
 controls[name].addEventListener('input',()=>{manuallyEdited.add(name);delete automaticDefaults[name];});
}
function fillUnique(source,target){
 if(manuallyEdited.has(target))return;
 const current=String(controls[target].value).trim();
 const wasAutomatic=Object.prototype.hasOwnProperty.call(automaticDefaults,target);
 if(wasAutomatic && current!==automaticDefaults[target]){manuallyEdited.add(target);delete automaticDefaults[target];return;}
 if(!wasAutomatic && current)return;
 const key=String(controls[source].value).trim();
 const values=key ? [...new Set(tables.Experiments.filter(r=>String(r[source]||'').trim()===key).map(r=>r[target]).filter(v=>target==='Field_MHz' ? typeof v==='number'&&Number.isFinite(v)&&v>0 : typeof v==='string'&&v.trim()&&v.trim().toLowerCase()!=='unknown'))] : [];
 const next=values.length===1 ? String(values[0]) : '';
 controls[target].value=next;
 if(next)automaticDefaults[target]=next;else delete automaticDefaults[target];
}
controls.PulseProgram.addEventListener('change',()=>fillUnique('PulseProgram','ExperimentType'));
controls.ExperimentType.addEventListener('input',refreshProtocolSuggestions);
controls.Instrument.addEventListener('change',()=>fillUnique('Instrument','Field_MHz'));
function message(text,error=false){const output=byId('message');output.setAttribute('role',error?'alert':'status');output.textContent=text;output.classList.toggle('error',error);if(error)output.scrollIntoView({block:'nearest'});}
byId('form').addEventListener('invalid',event=>{const name=event.target.name,label=groups[name]?.firstChild?.textContent||name||'Field';message(label+': '+event.target.validationMessage,true);},true);
function submitLabel(){return controls.mode.value==='experiment'?'Register experiment':controls.mode.value==='bulk'?'Add physical samples':'Add sample';}
function activeControls(){const experiment=controls.mode.value==='experiment';for(const group of ['identityBlock','bufferBlock','compositionBlock','storageBlock','bulkBlock','experiment','seriesFields']){const inactive=group==='experiment'?!experiment:group==='seriesFields'?experiment||controls.Series.value!=='new':group==='bulkBlock'?controls.mode.value!=='bulk':experiment;for(const input of byId(group).querySelectorAll('input,select,textarea'))input.disabled=inactive;}if(!experiment)updateLabelling();}

function options(control,items,base){const prior=control.value;control.replaceChildren();for(const [value,label]of [...base,...items]){const option=document.createElement('option');option.value=String(value);option.textContent=String(label);control.append(option);}if(Array.from(control.options).some(o=>o.value===prior))control.value=prior;}
function sampleChoices(){const series=Number(controls.Series.value||0),protein=controls.ProteinFilter.value;options(controls.Sample,tables.Samples.filter(s=>(!series||s.Series===series)&&(!protein||s.Protein===protein)&&NMRRegistration.selectableSample(s,controls.IncludeUnavailable.checked)).map(s=>[s.id,(s.SampleID||'(unnamed sample — fix ID)')+(NMRRegistration.sampleStatus(s.Available)==='Usable'?'':' — '+NMRRegistration.sampleStatus(s.Available))]),[['','Not assigned yet']]);}
controls.ProteinFilter.addEventListener('change',sampleChoices);controls.IncludeUnavailable.addEventListener('change',sampleChoices);
function modeChanged(){const experiment=controls.mode.value==='experiment';for(const group of ['identityBlock','bufferBlock','compositionBlock','storageBlock'])byId(group).classList.toggle('hidden',experiment);byId('experiment').classList.toggle('hidden',!experiment);byId('bulkBlock').classList.toggle('hidden',controls.mode.value!=='bulk');byId('title').textContent=experiment?'Register an experiment':'Add samples';if(experiment&&controls.Series.value==='new')controls.Series.value='';const option=Array.from(controls.Series.options).find(o=>o.value==='new');if(option)option.disabled=experiment;byId('seriesFields').classList.toggle('hidden',experiment||controls.Series.value!=='new');byId('submit').textContent=busy?(experiment?'Registering experiment…':'Adding samples…'):submitLabel();activeControls();sampleChoices();refreshRack();}
controls.mode.addEventListener('change',modeChanged);
controls.Series.addEventListener('change',()=>{modeChanged();const selected=tables.Series.find(s=>s.id===Number(controls.Series.value));if(selected&&controls.mode.value!=='experiment')for(const key of [...Object.keys(numericLabels),'Protein','Construct','Buffer','OtherConditions'])controls[key].value=selected[key]??'';if(selected&&controls.mode.value!=='experiment')setLabelling(selected.Isotope||'');});
function rowObjects(columns){return (columns.id||[]).map((id,index)=>Object.fromEntries(Object.entries(columns).map(([key,values])=>[key,values[index]])));}
async function refreshChoices(){const data=await Promise.all(['Series','Samples','Experiments'].map(name=>grist.docApi.fetchTable(name)));['Series','Samples','Experiments'].forEach((name,i)=>tables[name]=rowObjects(data[i]));refreshProtocolSuggestions();options(controls.Series,tables.Series.map(s=>[s.id,s.SeriesID||'(unnamed series — fix label)']),[['','No series'],['new','Create a series']]);options(controls.ProteinFilter,[...new Set(tables.Samples.map(s=>s.Protein).filter(Boolean))].sort().map(p=>[p,p]),[['','All proteins']]);refreshRack();modeChanged();}
byId('refresh').addEventListener('click',()=>refreshChoices().then(()=>message('Choices refreshed.')).catch(error=>message(error.message,true)));
byId('form').addEventListener('submit',async event=>{event.preventDefault();if(busy)return;busy=true;byId('submit').disabled=true;byId('refresh').disabled=true;byId('submit').textContent=controls.mode.value==='experiment'?'Registering experiment…':'Adding samples…';message(controls.mode.value==='experiment'?'Registering experiment…':'Adding samples…');
 try{const values=Object.fromEntries(Object.entries(controls).map(([name,input])=>[name,input.type==='checkbox'?input.checked:input.value]));await refreshChoices();const plan=values.mode==='experiment'?NMRRegistration.experimentActions(values,tables):NMRRegistration.sampleActions(values,tables);
 if(values.mode!=='experiment'){const [metadataTables,metadataColumns]=await Promise.all(['_grist_Tables','_grist_Tables_column'].map(name=>grist.docApi.fetchTable(name)));plan.actions=[...NMRRegistration.proteinChoiceActions(plan,metadataTables,metadataColumns),...plan.actions];}
 await grist.docApi.applyUserActions(plan.actions);const saved=plan.sample_ids.length?'Added '+plan.sample_ids.length+' physical sample'+(plan.sample_ids.length===1?'':'s')+':\n'+plan.sample_ids.join('\n'):'Experiment registered.';message(saved);if(values.mode==='bulk')controls.VaryValues.value='';try{await refreshChoices();const newSeries=plan.actions.find(a=>a[1]==='Series');if(newSeries){controls.Series.value=String(newSeries[2]);modeChanged();}}catch(error){message(saved+' Choices could not be refreshed: '+(error.message||String(error))+'. Use Refresh choices; the records were saved.');}
 }catch(error){message(error.message||String(error),true);}finally{busy=false;byId('refresh').disabled=false;byId('submit').textContent=submitLabel();refreshRack();}
});
modeChanged();refreshChoices().then(()=>{refreshRack();message('Ready.');}).catch(error=>message(error.message,true));
</script></body></html>'''


def registration_html(default_mode='sample'):
    if default_mode not in ('sample','bulk','experiment'): raise ValueError('Unknown registration mode')
    return (HTML.replace('/*LOGIC*/',JS_LOGIC).replace('/*RACK_JS*/',RACK_JS)
        .replace('/*EXPERIMENT_CHOICES*/',json.dumps(experiment_choices_by_category(),ensure_ascii=False)).replace('/*RACK_CSS*/',RACK_CSS).replace('/*MODE*/',json.dumps(default_mode)))
