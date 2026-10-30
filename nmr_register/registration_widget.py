"""Document-embedded registration form and independently testable action builders.

Embed ``registration_html('sample')`` or ``registration_html('experiment')`` in
Grist Widget Creator, with full document access. No external service is used.
The SDK calls are documented at https://support.getgrist.com/code/interfaces/grist_plugin_api.GristDocAPI/ .
"""
from datetime import date, datetime, timezone
import json
import math
import re

NUMERIC_RULES = {'pH': [0,14], 'NaCl_mM':[0,None], 'DNA_uM':[0,None],
                 'TFE_pct':[0,100], 'Protein_uM':[0,None]}
TEXT_FIELDS = ('Protein','Construct','Isotope','Buffer','OtherConditions')
UNITS = {'pH':'pH','NaCl_mM':'mMNaCl','DNA_uM':'uMDNA','TFE_pct':'pctTFE','Protein_uM':'uMprotein'}


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


def _rows(tables,name): return tables.get(name,[])
def _next(tables,name): return max([int(row['id']) for row in _rows(tables,name)]+[0])+1


def sample_actions(values,tables):
    """Build one atomic AddRecord transaction, preserving blanks and literal zero.

    ``tables`` contains raw rows (including integer id) fetched immediately before
    submission. Explicit row IDs also make simultaneous stale form submissions
    collide safely rather than create duplicate visible sample IDs.
    """
    datevalue,month,stamp=prepared(values.get('Prepared'),True)
    common={name:str(values.get(name) or '').strip() for name in TEXT_FIELDS}
    common.update({name:number(values.get(name),name,*bounds) for name,bounds in NUMERIC_RULES.items()})
    condition=token(values.get('DistinguishingCondition'))
    if not common['Protein'] or not common['Isotope'] or not condition:
        raise ValueError('Protein, labelling and distinguishing condition are required')
    start=number(1 if values.get('SampleNumber') in (None,'') else values['SampleNumber'],'Sample counter',1)
    if not start.is_integer(): raise ValueError('Sample counter must be an integer')
    series=values.get('Series') or ''
    actions=[]
    if series == 'new':
        series_id=str(values.get('SeriesID') or '_'.join([stamp,token(common['Isotope']),token(common['Protein']),condition,'series'])).strip()
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
    box=str(values.get('StorageBox') or '').strip();row=str(values.get('StorageRow') or '').strip().upper()
    column=number(values.get('StorageColumn'),'Storage column',1,12)
    if bool(row) != (column is not None) or (row and not re.fullmatch('[A-H]',row)) or (column is not None and not column.is_integer()):
        raise ValueError('Storage slot needs a row A-H and column 1-12')
    if row and not box: raise ValueError('Enter a box name for a storage slot')
    availability=values.get('Available') or 'Unknown'
    if availability not in ('Available','Unavailable','Unknown'): raise ValueError('Invalid availability')
    ids={row.get('SampleID') for row in _rows(tables,'Samples')}
    occupied={(str(s.get('StorageBox') or ''),str(s.get('StorageRow') or ''),s.get('StorageColumn')) for s in _rows(tables,'Samples')}
    if box:
        grid={ (r.get('StorageBox'),r.get('Row')) for r in _rows(tables,'StorageGrid') }
        grid_id=_next(tables,'StorageGrid')
        for letter in 'ABCDEFGH':
            if (box,letter) not in grid:
                actions.append(['AddRecord','StorageGrid',grid_id,{'StorageBox':box,'Row':letter}]);grid_id+=1
    created=[];next_id=_next(tables,'Samples')
    for index,value in enumerate(varying):
        suffix=condition+('_'+(format(value,'.15g')+UNITS[field] if field!='pH' else 'pH'+format(value,'.15g')) if field else '')
        sample_id='_'.join([stamp,token(common['Isotope']),token(common['Protein']),suffix,str(int(start)+index)])
        if sample_id in ids: raise ValueError('Sample ID already exists; change the human counter')
        ids.add(sample_id);created.append(sample_id)
        slot_row,slot_col='',None
        if row:
            slot='ABCDEFGH'.index(row)*12+int(column)-1+index
            if slot>=96: raise ValueError('Bulk storage slots exceed this box')
            slot_row='ABCDEFGH'[slot//12];slot_col=slot%12+1
            if (box,slot_row,slot_col) in occupied: raise ValueError('A selected storage slot is already occupied')
        fields={'SampleID':sample_id,'PreparedAt':datevalue,'PreparedDateText':month,
            'Series':series_ref,**common,'DistinguishingCondition':suffix,'SampleNumber':str(int(start)+index),
            'StorageBox':box,'StorageRow':slot_row,'StorageColumn':slot_col,'Available':availability}
        if field: fields[field]=value
        actions.append(['AddRecord','Samples',next_id+index,fields])
    return {'actions':actions,'sample_ids':created}


def experiment_actions(values,tables,now=None):
    series=int(values.get('Series') or 0);sample=int(values.get('Sample') or 0)
    samples={row['id']:row for row in _rows(tables,'Samples')}
    if series and not any(row['id']==series for row in _rows(tables,'Series')): raise ValueError('Selected series no longer exists')
    if sample and (sample not in samples or (series and samples[sample].get('Series')!=series)):
        raise ValueError('Sample is not in the selected series')
    if sample and not values.get('IncludeUnavailable') and samples[sample].get('Available') == 'Unavailable':
        raise ValueError('Sample is unavailable; explicitly include historical samples to select it')
    if sample and values.get('ProteinFilter') and samples[sample].get('Protein') != values['ProteinFilter']:
        raise ValueError('Sample does not match the selected protein')
    status=values.get('AcquisitionStatus')
    if status not in ('Planned','Acquired'): raise ValueError('Choose Planned or Acquired')
    acquired,month,_=prepared(values.get('Acquired'))
    source=str(values.get('SourceID') or '').strip();folder=str(values.get('Folder') or '').strip();expno=str(values.get('Expno') or '').strip()
    if folder and (folder.startswith('/') or any(p in ('..','.') or '211@950' in p for p in folder.split('/'))):
        raise ValueError('Invalid or protected acquisition folder')
    if expno and not expno.isdigit(): raise ValueError('Experiment number must contain digits')
    if source and folder and expno:
        for row in _rows(tables,'Experiments'):
            if any(tuple(str(row.get(prefix+k) or '') for k in ('SourceID','Folder','Expno'))==(source,folder,expno) for prefix in ('','Observed')):
                raise ValueError('Acquisition already exists; edit its existing record')
    fields={name:str(values.get(name) or '').strip() for name in ('Instrument','ExperimentType','PulseProgram','Notes')}
    fields.update(Series=series,Sample=sample,SourceID=source,Folder=folder,Expno=expno,
        AcquisitionStatus=status,AcquiredAt=acquired,AcquisitionDateText=month,RegisteredAt=now if now is not None else datetime.now(timezone.utc).timestamp(),
        Field_MHz=number(values.get('Field_MHz'),'Field (MHz)',0.000001),
        Temperature_C=number(values.get('Temperature_C'),'Temperature (C)',-273.15),ObservedField_MHz=None,RawBytes=None,
        ObservedAcquiredAt=None,LastChangedAt=None,ObservedTemperature_C=None)
    if not sample: fields.update({'Snapshot'+name:None for name in NUMERIC_RULES})
    return {'actions':[['AddRecord','Experiments',_next(tables,'Experiments'),fields]],'sample_ids':[]}


JS_LOGIC = r'''
const NMRRegistration = (() => {
 const rules = __RULES__, units = __UNITS__, textFields = __TEXT__;
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
 function sampleActions(v,t){
  const [datevalue,month,stamp]=prepared(v.Prepared,true), common={};
  textFields.forEach(k=>common[k]=String(v[k]||'').trim());Object.entries(rules).forEach(([k,b])=>common[k]=number(v[k],k,...b));
  const condition=token(v.DistinguishingCondition);if(!common.Protein||!common.Isotope||!condition)throw Error('Protein, labelling and distinguishing condition are required');
  const start=number(v.SampleNumber||1,'Sample counter',1);if(!Number.isInteger(start))throw Error('Sample counter must be an integer');
  let seriesRef=0;const actions=[];
  if(v.Series==='new'){
   const seriesID=String(v.SeriesID||[stamp,token(common.Isotope),token(common.Protein),condition,'series'].join('_')).trim();
   if(rows(t,'Series').some(r=>r.SeriesID===seriesID))throw Error('Series ID already exists; select that series');
   seriesRef=next(t,'Series');actions.push(['AddRecord','Series',seriesRef,{SeriesID:seriesID,...common}]);
  }else{seriesRef=Number(v.Series||0);if(seriesRef&&!rows(t,'Series').some(r=>r.id===seriesRef))throw Error('Selected series no longer exists');}
  let varying=[null],field=null;
  if(v.mode==='bulk'){
   field=v.VaryField||'DNA_uM';if(!Object.hasOwn(rules,field))throw Error('Unknown varying condition');
   const raw=String(v.VaryValues||'').trim().split(/[,;\s]+/);if(!raw[0]||raw.length>96)throw Error('Enter between 1 and 96 varying values');
   varying=raw.map(x=>number(x,field,...rules[field]));if(new Set(varying).size!==varying.length)throw Error('Repeated varying values; add replicates separately with distinct counters');
  }
  const box=String(v.StorageBox||'').trim(), row=String(v.StorageRow||'').trim().toUpperCase(), column=number(v.StorageColumn,'Storage column',1,12);
  if(Boolean(row)!==(column!==null)||(row&&!/^[A-H]$/.test(row))||(column!==null&&!Number.isInteger(column)))throw Error('Storage slot needs a row A-H and column 1-12');
  if(row&&!box)throw Error('Enter a box name for a storage slot');
  const available=v.Available||'Unknown';if(!['Available','Unavailable','Unknown'].includes(available))throw Error('Invalid availability');
  if(box){let gridID=next(t,'StorageGrid');for(const letter of 'ABCDEFGH')if(!rows(t,'StorageGrid').some(r=>r.StorageBox===box&&r.Row===letter))actions.push(['AddRecord','StorageGrid',gridID++,{StorageBox:box,Row:letter}]);}
  const ids=new Set(rows(t,'Samples').map(r=>r.SampleID)),created=[],firstID=next(t,'Samples');
  varying.forEach((value,i)=>{
   const suffix=condition+(field?'_'+(field==='pH'?'pH'+value:value+units[field]):''), sampleID=[stamp,token(common.Isotope),token(common.Protein),suffix,start+i].join('_');
   if(ids.has(sampleID))throw Error('Sample ID already exists; change the human counter');ids.add(sampleID);created.push(sampleID);
   let slotRow='',slotCol=null;
   if(row){const slot='ABCDEFGH'.indexOf(row)*12+column-1+i;if(slot>=96)throw Error('Bulk storage slots exceed this box');slotRow='ABCDEFGH'[Math.floor(slot/12)];slotCol=slot%12+1;
    if(rows(t,'Samples').some(s=>s.StorageBox===box&&s.StorageRow===slotRow&&s.StorageColumn===slotCol))throw Error('A selected storage slot is already occupied');}
   const fields={SampleID:sampleID,PreparedAt:datevalue,PreparedDateText:month,Series:seriesRef,...common,DistinguishingCondition:suffix,SampleNumber:String(start+i),StorageBox:box,StorageRow:slotRow,StorageColumn:slotCol,Available:available};
   if(field)fields[field]=value;actions.push(['AddRecord','Samples',firstID+i,fields]);
  });return {actions,sample_ids:created};
 }
 function experimentActions(v,t,now=Date.now()/1000){
  const series=Number(v.Series||0),sample=Number(v.Sample||0),selected=rows(t,'Samples').find(r=>r.id===sample);
  if(series&&!rows(t,'Series').some(r=>r.id===series))throw Error('Selected series no longer exists');
  if(sample&&(!selected||(series&&selected.Series!==series)))throw Error('Sample is not in the selected series');
  if(sample&&!v.IncludeUnavailable&&selected.Available==='Unavailable')throw Error('Sample is unavailable; explicitly include historical samples to select it');
  if(sample&&v.ProteinFilter&&selected.Protein!==v.ProteinFilter)throw Error('Sample does not match the selected protein');
  if(!['Planned','Acquired'].includes(v.AcquisitionStatus))throw Error('Choose Planned or Acquired');
  const [acquired,month]=prepared(v.Acquired),source=String(v.SourceID||'').trim(),folder=String(v.Folder||'').trim(),expno=String(v.Expno||'').trim();
  if(folder&&(folder.startsWith('/')||folder.split('/').some(p=>p==='..'||p==='.'||p.includes('211@950'))))throw Error('Invalid or protected acquisition folder');
  if(expno&&!/^\d+$/.test(expno))throw Error('Experiment number must contain digits');
  if(source&&folder&&expno&&rows(t,'Experiments').some(r=>['','Observed'].some(p=>String(r[p+'SourceID']||'')===source&&String(r[p+'Folder']||'')===folder&&String(r[p+'Expno']||'')===expno)))throw Error('Acquisition already exists; edit its existing record');
  const fields={};['Instrument','ExperimentType','PulseProgram','Notes'].forEach(k=>fields[k]=String(v[k]||'').trim());
  Object.assign(fields,{Series:series,Sample:sample,SourceID:source,Folder:folder,Expno:expno,AcquisitionStatus:v.AcquisitionStatus,AcquiredAt:acquired,AcquisitionDateText:month,RegisteredAt:now,
   Field_MHz:number(v.Field_MHz,'Field (MHz)',.000001),Temperature_C:number(v.Temperature_C,'Temperature (C)',-273.15),ObservedField_MHz:null,RawBytes:null,ObservedAcquiredAt:null,LastChangedAt:null,ObservedTemperature_C:null});
  if(!sample)Object.keys(rules).forEach(k=>fields['Snapshot'+k]=null);
  return {actions:[['AddRecord','Experiments',next(t,'Experiments'),fields]],sample_ids:[]};
 }
 return {number,prepared,sampleActions,experimentActions};
})();
'''.replace('__RULES__',json.dumps(NUMERIC_RULES)).replace('__UNITS__',json.dumps(UNITS)).replace('__TEXT__',json.dumps(TEXT_FIELDS))

HTML = r'''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<script src="grist-plugin-api.js"></script><style>
body{margin:0;padding:12px;font:13px var(--grist-theme-font-family,system-ui,sans-serif);color:var(--grist-theme-text,#262633);background:var(--grist-theme-widget-bg,#fff)}
form{display:grid;gap:8px;max-width:920px}h2{font-size:15px;font-weight:600;margin:0 0 4px}p{margin:0 0 8px;color:var(--grist-theme-text-medium,#494949)}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:8px}label{display:grid;gap:4px;font-size:11px;color:var(--grist-theme-text-medium,#494949)}
input,select,textarea{box-sizing:border-box;min-width:0;width:100%;padding:6px 8px;border:1px solid var(--grist-theme-input-border,#d9d9d9);border-radius:4px;background:var(--grist-theme-widget-bg,#fff);color:var(--grist-theme-text,#262633);font:13px inherit}textarea{min-height:52px}input[type=checkbox]{width:auto;justify-self:start}
button{min-height:32px;border:0;border-radius:4px;background:var(--grist-theme-control-primary-bg,#16b378);color:var(--grist-theme-control-primary-fg,#fff);padding:4px 12px;font:inherit;cursor:pointer}button:disabled{opacity:.6;cursor:wait}
.secondary{background:var(--grist-theme-widget-bg,#fff);color:var(--grist-theme-control-fg,#16b378);border:1px solid var(--grist-theme-control-primary-bg,#16b378)}
.actions{display:flex;gap:8px;align-items:center}.hidden{display:none!important}#message{white-space:pre-wrap;overflow-wrap:anywhere;line-height:1.5}#message.error{color:var(--grist-theme-text-error,#d0021b)}
</style></head><body><form id="form"><h2>Register samples or an experiment</h2><p>Leave unknown values blank. Enter 0 when the value is zero.</p>
<div class="grid" id="top"></div><div class="grid" id="seriesFields"></div><div class="grid" id="shared"></div><div class="grid" id="bulk"></div><div class="grid" id="storage"></div><div class="grid" id="experiment"></div>
<div class="actions"><button id="submit" type="submit" disabled>Add sample</button><button id="refresh" class="secondary" type="button">Refresh choices</button></div><div id="message" role="status" aria-live="polite">Loading register…</div></form>
<script>
grist.ready({requiredAccess:'full'});
/*LOGIC*/
const byId=id=>document.getElementById(id),controls={},groups={},numericLabels={pH:'pH',NaCl_mM:'NaCl (mM)',DNA_uM:'DNA (µM)',TFE_pct:'TFE (%)',Protein_uM:'Protein (µM)'};
let tables={Series:[],Samples:[],Experiments:[],StorageGrid:[]},busy=false;
function field(group,name,label,kind='text',options=null){
 const wrapper=document.createElement('label');wrapper.textContent=label;const input=document.createElement(kind==='select'?'select':kind==='textarea'?'textarea':'input');
 if(kind!=='select'&&kind!=='textarea')input.type=kind;input.name=name;input.id='field_'+name;
 if(kind==='number')input.step='any';if(options)for(const [value,title]of options){const option=document.createElement('option');option.value=value;option.textContent=title;input.append(option);}
 wrapper.append(input);byId(group).append(wrapper);controls[name]=input;groups[name]=wrapper;return input;
}
field('top','mode','Create','select',[['sample','One physical sample'],['bulk','Series / multiple physical samples'],['experiment','Experiment']]).value=/*MODE*/;
field('top','Series','Series','select',[['','No series'],['new','Create a series']]);field('seriesFields','SeriesID','New series label (optional)');
field('shared','Prepared','Prepared: YYYY-MM-DD or YYYY-MM');field('shared','Protein','Protein');field('shared','Construct','Construct');field('shared','Isotope','Labelling, e.g. 15N');field('shared','Buffer','Buffer');
for(const [name,label]of Object.entries(numericLabels))field('shared',name,label,'number');
field('shared','OtherConditions','Other conditions');field('shared','DistinguishingCondition','Distinguishing condition, e.g. DNA_titration');field('shared','SampleNumber','First human counter','number').value='1';
field('bulk','VaryField','Vary one condition','select',Object.entries(numericLabels));field('bulk','VaryValues','Values separated by commas or lines, e.g. 0, 100, 200, 500','textarea');controls.VaryField.value='DNA_uM';
field('storage','StorageBox','Storage box (optional)');field('storage','StorageRow','First row A–H');field('storage','StorageColumn','First column 1–12','number');field('storage','Available','Availability','select',[['Unknown','Unknown'],['Available','Available'],['Unavailable','Unavailable']]);
field('experiment','ProteinFilter','Protein','select',[['','All proteins']]);field('experiment','IncludeUnavailable','Include unavailable / historical samples','checkbox');
field('experiment','Sample','Physical sample','select',[['','Not assigned yet']]);field('experiment','Acquired','Acquisition: YYYY-MM-DD or YYYY-MM');
for(const [name,label]of [['SourceID','Source ID'],['Folder','Acquisition folder'],['Expno','Experiment number'],['Instrument','Instrument'],['ExperimentType','Experiment type'],['PulseProgram','Pulse program']])field('experiment',name,label);
field('experiment','AcquisitionStatus','Acquisition status','select',[['','Choose Planned or Acquired'],['Planned','Planned'],['Acquired','Acquired']]);controls.SourceID.value='nmrdata';field('experiment','Field_MHz','Field (MHz)','number');field('experiment','Temperature_C','Temperature (°C)','number');field('experiment','Notes','Notes','textarea');
function message(text,error=false){byId('message').textContent=text;byId('message').classList.toggle('error',error);}
function options(control,items,base){const prior=control.value;control.replaceChildren();for(const [value,label]of [...base,...items]){const option=document.createElement('option');option.value=String(value);option.textContent=String(label);control.append(option);}if(Array.from(control.options).some(o=>o.value===prior))control.value=prior;}
function sampleChoices(){const series=Number(controls.Series.value||0),protein=controls.ProteinFilter.value;options(controls.Sample,tables.Samples.filter(s=>(!series||s.Series===series)&&(!protein||s.Protein===protein)&&(controls.IncludeUnavailable.checked||s.Available!=='Unavailable')).map(s=>[s.id,s.SampleID||'(unnamed sample — fix ID)']),[['','Not assigned yet']]);}
controls.ProteinFilter.addEventListener('change',sampleChoices);controls.IncludeUnavailable.addEventListener('change',sampleChoices);
function modeChanged(){const experiment=controls.mode.value==='experiment';for(const group of ['shared','storage'])byId(group).classList.toggle('hidden',experiment);byId('experiment').classList.toggle('hidden',!experiment);byId('bulk').classList.toggle('hidden',controls.mode.value!=='bulk');if(experiment&&controls.Series.value==='new')controls.Series.value='';const option=Array.from(controls.Series.options).find(o=>o.value==='new');if(option)option.disabled=experiment;byId('seriesFields').classList.toggle('hidden',experiment||controls.Series.value!=='new');byId('submit').textContent=experiment?'Register experiment':controls.mode.value==='bulk'?'Add physical samples':'Add sample';sampleChoices();}
controls.mode.addEventListener('change',modeChanged);
controls.Series.addEventListener('change',()=>{modeChanged();const selected=tables.Series.find(s=>s.id===Number(controls.Series.value));if(selected&&controls.mode.value!=='experiment')for(const key of [...Object.keys(numericLabels),'Protein','Construct','Isotope','Buffer','OtherConditions'])controls[key].value=selected[key]??'';});
function rowObjects(columns){return (columns.id||[]).map((id,index)=>Object.fromEntries(Object.entries(columns).map(([key,values])=>[key,values[index]])));}
async function refreshChoices(){const data=await Promise.all(['Series','Samples','Experiments','StorageGrid'].map(name=>grist.docApi.fetchTable(name)));['Series','Samples','Experiments','StorageGrid'].forEach((name,i)=>tables[name]=rowObjects(data[i]));options(controls.Series,tables.Series.map(s=>[s.id,s.SeriesID||'(unnamed series — fix label)']),[['','No series'],['new','Create a series']]);options(controls.ProteinFilter,[...new Set(tables.Samples.map(s=>s.Protein).filter(Boolean))].sort().map(p=>[p,p]),[['','All proteins']]);modeChanged();}
byId('refresh').addEventListener('click',()=>refreshChoices().then(()=>message('Choices refreshed.')).catch(error=>message(error.message,true)));
byId('form').addEventListener('submit',async event=>{event.preventDefault();if(busy)return;busy=true;byId('submit').disabled=true;byId('refresh').disabled=true;
 try{const values=Object.fromEntries(Object.entries(controls).map(([name,input])=>[name,input.type==='checkbox'?input.checked:input.value]));await refreshChoices();const plan=values.mode==='experiment'?NMRRegistration.experimentActions(values,tables):NMRRegistration.sampleActions(values,tables);
 await grist.docApi.applyUserActions(plan.actions);message(plan.sample_ids.length?'Added '+plan.sample_ids.length+' physical sample'+(plan.sample_ids.length===1?'':'s')+':\n'+plan.sample_ids.join('\n'):'Experiment registered.');if(plan.sample_ids.length)controls.SampleNumber.value=String(Number(values.SampleNumber||1)+plan.sample_ids.length);if(values.mode==='bulk')controls.VaryValues.value='';await refreshChoices();const newSeries=plan.actions.find(a=>a[1]==='Series');if(newSeries){controls.Series.value=String(newSeries[2]);modeChanged();}
 }catch(error){message(error.message||String(error),true);}finally{busy=false;byId('submit').disabled=false;byId('refresh').disabled=false;}
});
modeChanged();refreshChoices().then(()=>{byId('submit').disabled=false;message('Ready.');}).catch(error=>message(error.message,true));
</script></body></html>'''


def registration_html(default_mode='sample'):
    if default_mode not in ('sample','bulk','experiment'): raise ValueError('Unknown registration mode')
    return HTML.replace('/*LOGIC*/',JS_LOGIC).replace('/*MODE*/',json.dumps(default_mode))
