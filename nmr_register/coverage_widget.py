"""Document-owned Electrostatics matrix. Runtime aggregation has no record writes.

SDK: https://support.getgrist.com/code/modules/grist_plugin_api/
Options layout: grist-core app/client/models/entities/ViewSectionRec.ts.
Only Save axes writes, and only this bound creator section's options metadata.
"""
import json


CORE_JS = r'''
const Coverage = (() => {
  const blank = x => x === null || x === undefined || (typeof x==='string'&&!x.trim());
  const num = x => !blank(x) && typeof x !== 'boolean' && !Array.isArray(x) && Number.isFinite(Number(x)) ? Number(x) : null;
  const rows = table => Array.isArray(table) ? table.map(r => r.fields ? {...r.fields,id:r.id} : {...r}) : (table.id||[]).map((id,i) => Object.fromEntries(Object.keys(table).map(k => [k,table[k][i]])));
  const unique = xs => [...new Set(xs)];
  const field = x => {const n=num(x);return n!==null&&n>0&&Math.round(n)>0?Math.round(n):null;};
  const norm = x => String(x||'').trim();
  const type = x => {
    const text=norm(x), key=text.toLowerCase().replace(/[\s_-]/g,'');
    return ({r1:'t1',t1:'t1','15nt1':'t1','15nt1hsqc':'t1',r1rho:'t1rho',t1rho:'t1rho','t1ρ':'t1rho','15nt1rho':'t1rho',
      r2:'t2',t2:'t2',noe:'noe',hetnoe:'noe',heteronuclearnoe:'noe','15nhetnoe':'noe',n15hetnoe:'noe',hsqc:'hsqc','15nhsqc':'hsqc','13chsqc':'13c_hsqc',cpmg:'cpmg',
      etaxy:'eta_xy','ηxy':'eta_xy',t1sat:'water_t1sat',watert1:'water_t1sat',watert1sat:'water_t1sat',nicot1sat:'water_t1sat','watert1(t1sat)':'water_t1sat'})[key] || (key==='unknown'?null:text||null);
  };
  const label = value => ({t1:'T1',t1rho:'T1ρ',noe:'Heteronuclear NOE',eta_xy:'etaXY',water_t1sat:'Water T1 (T1sat)',t2:'T2',hsqc:'HSQC','13c_hsqc':'13C HSQC',cpmg:'CPMG'})[value]||value;
  const status = r => {
    const text=norm(r.EffectiveAcquisitionStatus||r.AcquisitionStatus||r.ObservedAcquisitionStatus).toLowerCase();
    return ['completed','finished','done','acquired'].includes(text)?'Completed':['aborted','stopped'].includes(text)?'Aborted':'Pending / unknown';
  };
  const tuple = (r,prefix='') => {
    const source=norm(r[prefix+'SourceID']),rawFolder=r[prefix+'Folder'],folder=norm(rawFolder),exp=String(r[prefix+'Expno']??'').trim();
    const root=prefix==='Observed'&&rawFolder===''&&r.AcquisitionKey===source+':'+exp;
    return source&&/^[0-9]+$/.test(exp)&&(folder||root) ? [source,folder,exp] : null;
  };
  const identity = r => r.AcquisitionKey ? tuple(r,'Observed') : tuple(r);
  const indexes = records => {
    const keys=new Map(),tuples=new Map(),fingerprints=new Map();
    for(const r of records){if(r.ObservedFingerprint){const f=fingerprints.get(r.ObservedFingerprint)||new Set();f.add(JSON.stringify(identity(r)));fingerprints.set(r.ObservedFingerprint,f);}if(r.AcquisitionKey)keys.set(r.AcquisitionKey,(keys.get(r.AcquisitionKey)||0)+1);const t=identity(r);if(t){const k=JSON.stringify(t);tuples.set(k,(tuples.get(k)||0)+1);}}
    return {keys,tuples,fingerprints};
  };
  const issue = (r,index) => {
    if(norm(r.ReconciliationReason))return norm(r.ReconciliationReason);
    if(/(?:acquisition identity|duplicate acquisition|\bsplit\b|\bcopied\b|\breplacement\b|\b(?:source|folder|expno|acquisition key)(?:\s+(?:id|identity|path|key|number))?\s+(?:mismatch(?:es)?|differ(?:s|ent)?|conflict(?:s)?|ambiguous|invalid|missing|unresolved|requires?\s+(?:review|reconciliation))\b|\bfolder begins with an invalid date\b)/i.test(norm(r.MigrationIssue)))return norm(r.MigrationIssue);
    const t=identity(r);if(!t)return 'Exact source identity missing';
    if(r.ObservedFingerprint&&index.fingerprints.get(r.ObservedFingerprint)?.size>1)return 'Matching acquisition fingerprint at different sources; review possible copy';
    if(index.tuples.get(JSON.stringify(t))>1)return 'Duplicate exact source identity';
    if(r.AcquisitionKey){
      if(index.keys.get(r.AcquisitionKey)>1)return 'Duplicate collector key';
      if(r.AcquisitionKey!==t[0]+':'+(t[1]?t[1]+'/':'')+t[2])return 'Collector key differs from observed source';
      for(let i=0;i<3;i++){const human=norm(r[['SourceID','Folder','Expno'][i]]);if(human&&human!==t[i])return 'Confirmed and observed source identities differ';}
    }
    return '';
  };
  const list = (values,minimum=0) => unique((Array.isArray(values)?values:[]).map(num).filter(x=>x!==null&&x>=minimum)).sort((a,b)=>a-b);
  const axes = value => ({fields:unique((Array.isArray(value?.fields)?value.fields:[]).map(field).filter(x=>x!==null)).sort((a,b)=>a-b),salts:list(value?.salts),types:unique((Array.isArray(value?.types)?value.types:[]).map(type).filter(Boolean))});
  const csv = r => Array.isArray(r.ResultCSV) && (r.ResultCSV[0]==='L'?r.ResultCSV.slice(1):r.ResultCSV).some(x=>Number.isInteger(x)&&x>0);
  const temperature = r => num(r.EffectiveTemperature_C);
  const temperatureKey = value => value===null?'unknown':String(value);
  // Only reconcile a vanished UI selection. EffectiveTemperature_C remains the
  // sole authority for record grouping; deliberately unobserved filters survive.
  const reconcileTemperature = (value,available) => {
    const key=String(value),n=num(value);if(available.includes(key)||n===null)return key;
    const whole=Math.round(n),half=Math.round(n*2)/2;
    const normalized=Math.abs(n-whole)<=0.2+1e-9?whole:Math.abs(n-half)<=0.02+1e-9?half:n;
    const next=String(normalized);return next!==key&&available.includes(next)?next:key;
  };
  const canonical = x => Array.isArray(x) && x[0]==='L' ? x.slice(1) : x;
  const prepNames=['SnapshotConstruct','SnapshotIsotope','SnapshotBuffer','SnapshotpH','SnapshotDNA_uM','SnapshotTFE_pct','SnapshotProtein_uM','SnapshotOtherConditions'];
  const protocolNames=['EffectiveInstrument','EffectivePulseProgram','EffectiveD1_s','EffectiveDelayList_s','EffectiveSpinLock_Hz','EffectiveSpinLockOffset_Hz'];
  const groupKey = r => {
    const values=[...prepNames,...protocolNames].map(k=>canonical(r[k]));
    const uncertain=[...prepNames.slice(0,7),'EffectiveInstrument','EffectivePulseProgram'].some(k=>blank(r[k]));
    return JSON.stringify([...values,uncertain?(num(r.Sample)>0?'tube:'+r.Sample:'record:'+r.id):'known conditions']);
  };
  // Assessment is formula-owned: Use means all QC criteria passed, or an
  // imported decision is still active while every criterion is unassessed.
  const qcUse = r => r.Assessment==='Use';
  const summarize = records => {
    const trusted=records.filter(r=>!r._issue),finished=trusted.filter(r=>status(r)==='Completed');
    const groups=new Map(),goodGroups=new Map();for(const r of trusted){const k=groupKey(r);if(!groups.has(k))groups.set(k,[]);groups.get(k).push(r);if(qcUse(r)){if(!goodGroups.has(k))goodGroups.set(k,[]);goodGroups.get(k).push(r);}}
    const goodCounts=[...goodGroups.values()].map(group=>group.length).sort((a,b)=>b-a);
    return {ids:records.map(r=>r.id),good:goodCounts.reduce((a,b)=>a+b,0),goodCounts,maxCompatibleGood:Math.max(0,...goodCounts),
      unassessed:trusted.filter(r=>blank(r.Assessment)).length,limited:trusted.filter(r=>r.Assessment==='Limited use').length,
      finished:finished.length,stopped:trusted.filter(r=>status(r)==='Aborted').length,
      pending:trusted.filter(r=>status(r)==='Pending / unknown').length,review:records.filter(r=>r._issue).length,
      csv:trusted.filter(csv).length,analyzed:trusted.filter(r=>r.AnalysisStatus==='Analyzed').length,
      notAnalyzed:trusted.filter(r=>!(r.AnalysisStatus==='Analyzed'||(blank(r.AnalysisStatus)&&csv(r)))).length,
      assessments:Object.fromEntries(['Use','Limited use','Do not use','Unassessed'].map(v=>[v,trusted.filter(r=>(r.Assessment||'Unassessed')===v).length])),
      manual:finished.filter(r=>!r.AcquisitionKey).length,tubes:unique(trusted.map(r=>num(r.Sample)).filter(x=>x>0)).length,
      unknownTubes:trusted.filter(r=>!(num(r.Sample)>0)).length,groups:[...groups.values()].map(group=>({ids:group.map(r=>r.id),records:group}))};
  };
  function build(input,samples,options,protein,temp){
    const all=[...new Map(rows(input).filter(r=>Number.isInteger(r.id)&&r.id>0&&r.AcquisitionStatus!=='Dummy'&&r.EffectiveAcquisitionStatus!=='Dummy').map(r=>[r.id,r])).values()],index=indexes(all);
    const scope=all.filter(r=>r.Study==='Electrostatics');
    const counts=new Map();for(const r of scope){const name=norm(r.SnapshotProtein);counts.set(name,(counts.get(name)||0)+1);}
    const proteinNames=unique([...scope.map(r=>norm(r.SnapshotProtein)),...rows(samples).map(r=>norm(r.Protein)),...Object.keys(options?.proteins||{})].filter(Boolean)).sort((a,b)=>(counts.get(b)||0)-(counts.get(a)||0)||(a<b?-1:a>b?1:0));
    const selected=scope.filter(r=>norm(r.SnapshotProtein)===protein);
    const observed={fields:unique(selected.map(r=>field(r.EffectiveField_MHz)).filter(x=>x!==null)).sort((a,b)=>a-b),salts:list(selected.map(r=>r.SnapshotNaCl_mM)),types:unique(selected.map(r=>type(r.EffectiveType)).filter(Boolean))};
    const configured=Object.prototype.hasOwnProperty.call(options?.proteins||{},protein);
    const expected=configured?axes(options.proteins[protein]):axes({...observed,types:['t1','t1rho','noe','eta_xy']});
    const temperatures=unique(selected.map(r=>temperatureKey(temperature(r)))).sort((a,b)=>a==='unknown'?1:b==='unknown'?-1:Number(a)-Number(b));
    const CoverageField=field;
    const records=selected.filter(r=>temperatureKey(temperature(r))===String(temp)).map(r=>({...r,_issue:issue(r,index)}));
    const cells=[],placed=new Set();
    for(const field of expected.fields)for(const experiment of expected.types)for(const salt of expected.salts){
      const members=records.filter(r=>CoverageField(r.EffectiveField_MHz)===field&&type(r.EffectiveType)===experiment&&num(r.SnapshotNaCl_mM)===salt);
      members.forEach(r=>placed.add(r.id));cells.push({field,type:experiment,salt,records:members,...summarize(members)});
    }
    return {proteins:proteinNames,observed,axes:expected,temperatures,cells,records,unplaced:records.filter(r=>!placed.has(r.id)),unknownProtein:scope.filter(r=>!norm(r.SnapshotProtein)).length};
  }
  function mergeAxes(sectionRows,tableRows,sectionId,protein,value,base){
    if(!Number.isInteger(sectionId)||sectionId<=0)throw Error('This matrix is not bound to its widget; axes cannot be saved.');
    const section=rows(sectionRows).find(r=>r.id===sectionId),table=rows(tableRows).find(r=>r.tableId==='Experiments');
    if(!section||!table||section.tableRef!==table.id||!['custom','creator'].includes(section.parentKey))throw Error('Widget context changed; reload before saving axes.');
    const options=JSON.parse(section.options||'{}');
    if(typeof options.html!=='string'||!options.html.includes('data-nmr-coverage="1"'))throw Error('This section is not the expected coverage widget.');
    if(typeof options.customView!=='string')throw Error('Widget settings format changed; axes were not saved.');
    const custom=JSON.parse(options.customView);
    const object=x=>x!==null&&typeof x==='object'&&!Array.isArray(x);
    if(!object(custom)||!object(custom.widgetOptions||{}))throw Error('Invalid widget settings; axes were not saved.');
    const prefs=custom.widgetOptions||{},current=prefs.coverageAxes||{version:1,proteins:{}};
    if(!object(current)||current.version!==1||!object(current.proteins))throw Error('Invalid widget settings; axes were not saved.');
    const previous=Object.prototype.hasOwnProperty.call(current.proteins,protein)?axes(current.proteins[protein]):null;
    if(JSON.stringify(previous)!==JSON.stringify(base))throw Error('Axes for this protein changed elsewhere; reload before saving.');
    const desired=axes(value);if(!desired.fields.length||!desired.salts.length||!desired.types.length)throw Error('Choose at least one field, experiment and NaCl value.');
    const coverageAxes={...current,version:1,proteins:{...current.proteins,[protein]:{...current.proteins[protein],...desired}}};
    const updated={...options,customView:JSON.stringify({...custom,widgetOptions:{...prefs,coverageAxes}})};
    return {action:['UpdateRecord','_grist_Views_section',sectionId,{options:JSON.stringify(updated)}],coverageAxes};
  }
  return {blank,num,rows,type,label,status,axes,build,mergeAxes,temperatureKey,reconcileTemperature,groupKey,qcUse};
})();
'''


UI_JS = r'''
const SECTION_ID=__SECTION_ID__;
const $=id=>document.getElementById(id),el=(tag,text,cls)=>{const n=document.createElement(tag);if(text!==undefined)n.textContent=text;if(cls)n.className=cls;return n;};
let records=[],samples=[],saved={version:1,proteins:{}},protein='',temp='25',model=null,draft=null,base=null,selected=null,loading=false;
function message(text,error=false){$('message').textContent=text;$('message').className=error?'error':'muted';}
async function clearSelection(){selected=null;$('selection').replaceChildren();try{await grist.setSelectedRows([]);}catch(error){message('Linked acquisitions could not be cleared; reload the page before using their selection.',true);}}
function option(select,value,label){const item=el('option',label);item.value=String(value);select.append(item);}
function savedBase(){return Object.prototype.hasOwnProperty.call(saved.proteins||{},protein)?Coverage.axes(saved.proteins[protein]):null;}
function selectProtein(){const probe=Coverage.build(records,samples,saved,protein,temp);const values=probe.temperatures;
 temp=values.includes('25')?'25':values.length===1?values[0]:'25';base=savedBase();draft=null;draw();}
function checkboxGroup(id,choices,checked,key){const box=$(id);box.replaceChildren();for(const value of choices){const label=el('label'),input=el('input');input.type='checkbox';input.checked=checked.includes(value);input.value=String(value);
 input.addEventListener('change',async()=>{const values=[...box.querySelectorAll('input:checked')].map(n=>key==='types'?n.value:Number(n.value));await clearSelection();draft={...draft,[key]:values};renderMatrix();});label.append(input,el('span',key==='types'?Coverage.label(value):String(value)));box.append(label);}}
function draw(){model=Coverage.build(records,samples,saved,protein,temp);if(!protein&&model.proteins.length){protein=model.proteins[0];return selectProtein();}
 const selector=$('protein');selector.replaceChildren();for(const name of model.proteins)option(selector,name,name);selector.value=protein;
 const temperature=$('temperature');temperature.replaceChildren();for(const value of [...new Set([...model.temperatures,temp])])option(temperature,value,value==='unknown'?'Unknown temperature':value+' °C');temperature.value=temp;
 if(!draft)draft={...model.axes,fields:[...model.axes.fields],salts:[...model.axes.salts],types:[...model.axes.types]};
 checkboxGroup('fields',[...new Set([...model.observed.fields,...draft.fields])].sort((a,b)=>a-b),draft.fields,'fields');
 const sampleSalts=Coverage.rows(samples).filter(r=>r.Protein===protein).map(r=>Coverage.num(r.NaCl_mM)).filter(x=>x!==null);
 checkboxGroup('salts',[...new Set([...model.observed.salts,...draft.salts,...sampleSalts])].sort((a,b)=>a-b),draft.salts,'salts');
 checkboxGroup('types',[...new Set(['t1','t1rho','noe','eta_xy','water_t1sat','t2',...model.observed.types,...draft.types])],draft.types,'types');renderMatrix();}
function renderMatrix(){const prefs={version:1,proteins:{...saved.proteins,[protein]:draft}};model=Coverage.build(records,samples,prefs,protein,temp);
 const table=$('matrix');table.replaceChildren();table.style.minWidth=(144+model.axes.fields.length*170)+'px';
 table.append(el('caption','Experiment by nominal field; each line is a NaCl condition','sr-only'));
 const head=el('thead'),header=el('tr'),corner=el('th','Experiment');corner.scope='col';header.append(corner);
 for(const field of model.axes.fields){const th=el('th',field+' MHz');th.scope='col';header.append(th);}head.append(header);table.append(head);
 const body=el('tbody');for(const type of model.axes.types){const row=el('tr'),label=el('th',Coverage.label(type));label.scope='row';row.append(label);
 for(const field of model.axes.fields){const td=el('td');for(const salt of model.axes.salts){const cell=model.cells.find(c=>c.field===field&&c.type===type&&c.salt===salt),button=el('button',undefined,'cell');
 button.type='button';button.dataset.field=field;button.dataset.type=type;button.dataset.salt=salt;
 button.classList.toggle('good',cell.maxCompatibleGood>=2);button.classList.toggle('single',cell.good>0&&cell.maxCompatibleGood<2);
 button.title=cell.good+' acquisitions with overall Use';
 button.append(el('strong',salt+' mM ('+cell.good+'×)'));
 const flags=[];if(cell.limited)flags.push('+'+cell.limited+' limited');if(cell.assessments['Do not use'])flags.push('+'+cell.assessments['Do not use']+' discarded');
 if(flags.length)button.append(el('small',flags.join(' · '),'muted'));
 const needs=[];if(cell.unassessed)needs.push(cell.unassessed+' unassessed');if(cell.notAnalyzed)needs.push(cell.notAnalyzed+' analysis not confirmed');if(cell.review)needs.push(cell.review+' identity review');
 if(needs.length)button.append(el('small',needs.join(' · '),'needs'));
 button.setAttribute('aria-label',Coverage.label(type)+', '+field+' MHz, '+salt+' mM NaCl: '+[button.title,...flags,...needs].join(', '));
 button.addEventListener('click',()=>choose(cell));td.append(button);}row.append(td);}body.append(row);}table.append(body);
 $('empty').hidden=Boolean(model.axes.fields.length&&model.axes.salts.length&&model.axes.types.length);
 const outside=$('outside');outside.replaceChildren();if(model.unplaced.length){const b=el('button',model.unplaced.length+' records have missing or unselected axes');b.type='button';b.addEventListener('click',()=>choose({records:model.unplaced,ids:model.unplaced.map(r=>r.id),...{finished:0}},true));outside.append(b);}
 if(model.unknownProtein)outside.append(el('span',model.unknownProtein+' Electrostatics records have no frozen protein.','muted'));
 $('save').disabled=!SECTION_ID||loading;}
async function choose(cell,outside=false){$('selection').replaceChildren();selected=null;try{await grist.setSelectedRows(cell.ids);}catch(error){await clearSelection();message('Selection failed; linked acquisitions may be stale. Reload the page.',true);return;}
 selected=cell;message('');const box=$('selection');box.append(el('h3',outside?'Records outside the selected axes':cell.field+' MHz · '+Coverage.label(cell.type)+' · '+cell.salt+' mM NaCl'));
 if(!outside){box.append(el('p',cell.good+' overall Use · '+cell.finished+' completed total · '+cell.stopped+' aborted · '+cell.pending+' pending / unknown · '+cell.review+' identity review'));
 box.append(el('p',cell.csv+' selected CSV · '+cell.analyzed+' marked Analyzed · '+cell.tubes+' known tubes'+(cell.unknownTubes?' · '+cell.unknownTubes+' unknown sample':''),'muted'));
 box.append(el('p',Object.entries(cell.assessments).map(([key,value])=>key+': '+value).join(' · '),'muted'));
 if(cell.manual)box.append(el('p',cell.manual+' Completed status recorded manually; no matching collector observation.','muted'));
 if(cell.groups.length>1)box.append(el('p','Sample or acquisition settings differ or are incomplete. Green requires at least two good experiments with comparable known sample and acquisition settings.','muted'));}
 const detail=el('details');detail.append(el('summary','Condition and source details ('+cell.ids.length+' records)'));
 for(const r of cell.records){const line=el('div',undefined,'record');line.append(el('strong',(r.SnapshotSampleID||'Sample unassigned')+' · '+Coverage.status(r)+' · '+(r.Assessment||'Unassessed')),el('p',[r.SourceID||r.ObservedSourceID,r.Folder||r.ObservedFolder,r.Expno||r.ObservedExpno].filter(Boolean).join(' / '),'muted'));
 const context=['SnapshotConstruct','SnapshotIsotope','SnapshotBuffer','SnapshotpH','SnapshotProtein_uM','SnapshotDNA_uM','SnapshotTFE_pct','EffectiveInstrument','EffectivePulseProgram','EffectiveD1_s','EffectiveDelayList_s','EffectiveSpinLock_Hz','EffectiveSpinLockOffset_Hz'].filter(k=>!Coverage.blank(r[k])).map(k=>k.replace('Snapshot','').replace('Effective','')+': '+r[k]);line.append(el('p',context.join(' · '),'muted'));
 if(r._issue)line.append(el('p',r._issue,'bad'));if(r.MetadataDifference)line.append(el('p',r.MetadataDifference,'muted'));detail.append(line);}box.append(detail);
 if(!cell.ids.length)box.append(el('p','No completed acquisition recorded. The register may be incomplete.','muted'));}
async function load(){if(loading)return;loading=true;try{const data=await Promise.all([grist.docApi.fetchTable('Experiments'),grist.docApi.fetchTable('Samples')]);records=data[0];samples=data[1];const available=Coverage.build(records,samples,saved,protein,temp).temperatures;temp=Coverage.reconcileTemperature(temp,available);await clearSelection();draw();}catch(error){message('Could not read register data. Refresh to retry.',true);}finally{loading=false;$('save').disabled=!SECTION_ID;}}
async function save(){const requestedProtein=protein,requestedAxes=Coverage.axes(draft),requestedBase=base===null?null:Coverage.axes(base);loading=true;$('save').disabled=true;try{const [sections,tables]=await Promise.all([grist.docApi.fetchTable('_grist_Views_section'),grist.docApi.fetchTable('_grist_Tables')]);const change=Coverage.mergeAxes(sections,tables,SECTION_ID,requestedProtein,requestedAxes,requestedBase);
 await grist.docApi.applyUserActions([change.action]);saved=change.coverageAxes;if(protein===requestedProtein)base=requestedAxes;message('Selections saved for '+requestedProtein+'.');}catch(error){message(error.message||'Selections could not be saved; reload and retry.',true);}finally{loading=false;$('save').disabled=!SECTION_ID;}}
$('protein').addEventListener('change',async()=>{await clearSelection();protein=$('protein').value;selectProtein();});
$('temperature').addEventListener('change',async()=>{await clearSelection();temp=$('temperature').value;draw();});
$('refresh').addEventListener('click',load);$('save').addEventListener('click',save);
$('useAvailable').addEventListener('click',async()=>{if(!model||!draft)return;const available=model.observed,requestedProtein=protein;await clearSelection();if(protein!==requestedProtein)return;draft={...draft,fields:[...available.fields],salts:[...available.salts]};draw();message('Available fields and NaCl selected. Save selections to keep these choices.');});
for(const [id,key]of [['addField','fields'],['addSalt','salts']])$(id).addEventListener('click',async()=>{const input=$(id+'Value'),value=Coverage.num(input.value);if(value===null||value<0||(key==='fields'&&value===0)){message('Enter a valid '+(key==='fields'?'positive field in MHz.':'NaCl concentration in mM.'),true);return;}await clearSelection();draft={...draft,[key]:[...new Set([...draft[key],value])]};input.value='';draw();});
grist.ready({requiredAccess:'full',allowSelectBy:true,onEditOptions:()=>$('protein').focus()});
grist.onOptions(async options=>{const reference=base||Coverage.build(records,samples,saved,protein,temp).axes,dirty=draft&&JSON.stringify(Coverage.axes(draft))!==JSON.stringify(reference),previousBase=base;saved=options?.coverageAxes||{version:1,proteins:{}};const nextBase=savedBase();
 if(dirty&&JSON.stringify(nextBase)!==JSON.stringify(Coverage.axes(draft))){if(JSON.stringify(nextBase)!==JSON.stringify(previousBase))message('Selections changed elsewhere; reload before saving this draft.',true);}else{base=nextBase;if(!dirty)draft=null;}await clearSelection();if(records.id||records.length)draw();});
let refreshTimer;grist.onRecords(()=>{clearTimeout(refreshTimer);refreshTimer=setTimeout(load,80);});
'''


def coverage_html(section_id=None):
    """Bind durable configuration writes to this one creator section, if supplied."""
    if section_id is not None and (isinstance(section_id, bool) or not isinstance(section_id, int) or section_id <= 0):
        raise ValueError('section_id must be a positive integer')
    script = CORE_JS + '\n' + UI_JS.replace('__SECTION_ID__', json.dumps(section_id))
    return r'''<!doctype html><html lang="en" data-nmr-coverage="1"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><script src="https://docs.getgrist.com/grist-plugin-api.js"></script><style>
*{box-sizing:border-box}body{margin:0;padding:12px;font:13px system-ui,sans-serif;color:#273943;background:white}header{display:flex;gap:14px;align-items:end;flex-wrap:wrap;margin-bottom:10px}label{display:flex;align-items:center;gap:5px}header label{display:grid;gap:4px;font-size:12px}select,input,button{font:inherit;color:inherit}select,input{height:30px;border:1px solid #cad5db;border-radius:4px;padding:3px 7px}button{border:1px solid #cbd5da;border-radius:4px;background:white;padding:6px 9px;cursor:pointer}button:hover{background:#f2f6f8}button:disabled{opacity:.5;cursor:default}button:focus-visible{outline:2px solid #4b8cac;outline-offset:1px}h3{font-size:13px;margin:10px 0 5px}p{margin:5px 0;line-height:1.5}.muted{color:#71808a;font-size:11px}.error,.bad{color:#a0362e}.error{margin:8px 0}.matrix-scroll{overflow-x:auto;margin:10px 0;border:1px solid #dae1e5;border-radius:5px}.matrix-scroll:focus-visible{outline:2px solid #4b8cac;outline-offset:-2px}.sr-only{position:absolute;width:1px;height:1px;padding:0;margin:-1px;overflow:hidden;clip:rect(0,0,0,0);white-space:nowrap;border:0}table{border-collapse:collapse;width:100%;margin:0;table-layout:fixed}th,td{border:1px solid #dae1e5;text-align:center;padding:0;vertical-align:top}th{font-size:12px;font-weight:550;background:#fafbfc;padding:8px}th:first-child{width:144px;text-align:left}.cell{display:grid;align-content:start;justify-items:center;gap:3px;width:100%;min-height:51px;border:0;border-bottom:1px solid #e8edf0;border-radius:0;background:white;padding:8px 5px}.cell:last-child{border-bottom:0}.cell.good{background:#e5f3e9}.cell.single{background:#e7f0fc}.cell strong{font-size:13px;font-weight:550}.cell small{font-size:10px;line-height:1.35}.cell:hover{box-shadow:inset 0 0 0 2px #8ba9b7}.cell .needs{font-weight:700}details{margin:9px 0}summary{cursor:pointer;font-size:12px;color:#536975}.axisgrid{display:grid;grid-template-columns:repeat(3,minmax(130px,1fr));gap:15px;padding:9px 0}.choices{display:flex;gap:7px 11px;flex-wrap:wrap;min-height:24px}.choices label{font-size:12px}.choices input{height:auto;margin:0}.add{display:flex;gap:5px;margin-top:7px}.add input{width:90px}.record{border-top:1px solid #edf0f2;padding:6px 0}#selection{border-top:1px solid #e1e7ea;margin-top:9px}#outside{display:flex;align-items:center;gap:10px;flex-wrap:wrap;margin:7px 0}@media(max-width:650px){.axisgrid{grid-template-columns:1fr}}
</style></head><body><header><label>Protein<select id="protein"></select></label><label>Temperature group<select id="temperature"></select></label><button id="refresh" type="button">Refresh</button></header>
<div id="axes" aria-label="Coverage filters"><div class="axisgrid"><section><h3>Field (MHz)</h3><div id="fields" class="choices"></div><div class="add"><input id="addFieldValue" type="number" min="0" step="any" aria-label="Additional field in MHz"><button id="addField" type="button">Add field</button></div></section><section><h3>NaCl (mM)</h3><div id="salts" class="choices"></div><div class="add"><input id="addSaltValue" type="number" min="0" step="any" aria-label="Additional NaCl concentration in mM"><button id="addSalt" type="button">Add NaCl</button></div></section><section><h3>Experiment</h3><div id="types" class="choices"></div></section></div><button id="useAvailable" type="button">Use available</button> <button id="save" type="button">Save selections</button></div>
<p class="muted">Counts show experiments assessed Use. Blue: good experiments without a comparable pair. Green: at least two with comparable recorded sample and acquisition settings. Select a concentration for details.</p><div class="matrix-scroll" role="region" aria-label="Acquisition quality matrix" tabindex="0"><table id="matrix"></table></div><p id="empty" class="muted">Choose fields, experiment types and NaCl values above.</p><div id="outside"></div><div id="selection"></div>
<p id="message" class="muted" role="status"></p><script>''' + script + '</script></body></html>'
