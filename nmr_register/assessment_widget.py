"""Document-owned assessment card: explicit human ratings, computed decision."""
import json

CRITERIA = [
    ('QCAcquisitionCoverage', 'Acquisition coverage', ['Not assessed', 'Complete', 'Sufficient', 'Insufficient']),
    ('QCSettings', 'Acquisition settings', ['Not assessed', 'Appropriate', 'Usable with restrictions', 'Invalid']),
    ('QCFID', 'FIDs and spectrum', ['Not assessed', 'Clean', 'Usable with exclusions', 'Unusable']),
    ('QCSNR', 'Signal to noise', ['Not assessed', 'Adequate', 'Marginal', 'Inadequate']),
    ('QCSample', 'Sample condition', ['Not assessed', 'Acceptable', 'Concern', 'Unacceptable']),
    ('QCAnalysis', 'Analysis / fit', ['Not assessed', 'Reliable', 'Reliable — early points excluded', 'Limited', 'Two-phase — unresolved', 'Unreliable']),
]

CORE_JS = r'''
const AssessmentCard=(()=>{
 const criteria=__CRITERIA__,keys=criteria.map(c=>c[0]),note='AssessmentReason';
 const rows=table=>Array.isArray(table)?table.map(r=>r.fields?{...r.fields,id:r.id}:r):(table.id||[]).map((id,i)=>Object.fromEntries(Object.keys(table).map(k=>[k,table[k][i]])));
 const snapshot=record=>Object.fromEntries([...keys,note].map(k=>[k,k===note?String(record?.[k]||''):String(record?.[k]||'Not assessed')]).concat([['ImportedAssessment',String(record?.ImportedAssessment||'')]]));
 const equal=(a,b)=>JSON.stringify(snapshot(a))===JSON.stringify(snapshot(b));
 const identityKeys=['AcquisitionKey','SourceID','Folder','Expno','ObservedSourceID','ObservedFolder','ObservedExpno','ObservedFingerprint','ReconciliationReason','SampleTagIssue','ObservedSampleTagIssue','ObservedSampleID','RegisteredAt','Sample','SnapshotSampleID','Title','Instrument','ExperimentType','AcquiredAt','AcquisitionDateText','AcquisitionStatus','Field_MHz','Temperature_C','PulseProgram','D1_s','DelayList_s','SpinLock_Hz','SpinLockOffset_Hz'];
 const sameIdentity=(a,b)=>JSON.stringify(identityKeys.map(k=>a?.[k]??null))===JSON.stringify(identityKeys.map(k=>b?.[k]??null));
 const decisionTone=record=>record?.Assessment==='Do not use'?'fail':record?.Assessment==='Limited use'?'limited':record?.Assessment==='Use'?'pass':record?.AssessmentDisplay&&record.AssessmentDisplay!=='Not assessed'?'limited':'unchecked';
 const restrictedAnalysis=['Reliable — early points excluded','Two-phase — unresolved'];
 const relaxationType=value=>/^(?:15n)?(?:t1|r1)(?:rho|ρ)?$/i.test(String(value||'').replace(/[\s_-]+/g,''));
 const choicesFor=(key,record)=>{const choices=criteria.find(c=>c[0]===key)?.[2]||[];return key==='QCAnalysis'&&!relaxationType(record?.ExperimentType)?choices.filter(v=>!restrictedAnalysis.includes(v)):choices;};
 const choiceLabel=(key,value,record)=>{if(key!=='QCAnalysis')return value;const fit={'Reliable':'Good fit — all points','Reliable — early points excluded':'Good fit — early points excluded','Limited':'Usable fit — other concerns','Two-phase — unresolved':'Unusable fit — two-phase','Unreliable':'Unusable fit — other problem'},other={'Reliable':'Reliable','Limited':'Usable with reservations','Unreliable':'Unusable'};return (relaxationType(record?.ExperimentType)||restrictedAnalysis.includes(value)?fit:other)[value]||value;};
 const tone=(key,value)=>{if(value==='Not assessed'||!criteria.find(c=>c[0]===key)?.[2].includes(value))return 'unchecked';if(key==='QCAnalysis')return ({'Reliable':'pass','Reliable — early points excluded':'pass','Two-phase — unresolved':'fail','Limited':'limited','Unreliable':'fail'})[value];const choices=criteria.find(c=>c[0]===key)[2],index=choices.indexOf(value);return index===choices.length-1?'fail':key==='QCAcquisitionCoverage'||index===1?'pass':'limited';};
 function registrationTimestamp(base,draft,fresh,experimentRows,sampleRows,now){
  if(fresh.RegisteredAt!=null&&fresh.RegisteredAt!==0)return undefined;
  const before=snapshot(base),after=snapshot(draft);
  if(!criteria.some(([key,,choices])=>after[key]!==before[key]&&after[key]!=='Not assessed'&&choices.includes(after[key])))return undefined;
  if(['ReconciliationReason','SampleTagIssue','ObservedSampleTagIssue'].some(k=>String(fresh[k]||'').trim()))return undefined;
  if(!Number.isInteger(fresh.Sample)||fresh.Sample<=0)return undefined;
  const sample=rows(sampleRows).filter(r=>r.id===fresh.Sample);
  if(sample.length!==1||!String(sample[0].SampleID||'').trim())return undefined;
  if(fresh.ObservedSampleID&&fresh.ObservedSampleID!==sample[0].SampleID)return undefined;
  const text=v=>String(v??'').trim(),tuple=r=>['SourceID','Folder','Expno'].map(k=>text(r[k]));
  const observed=['ObservedSourceID','ObservedFolder','ObservedExpno'].map(k=>text(fresh[k]));
  if(!observed[0]||!/^\d+$/.test(observed[2])||observed[1].split('/').some(p=>p==='.'||p==='..'||p.includes('211@950'))||observed[1].startsWith('/'))return undefined;
  if(JSON.stringify(tuple(fresh))!==JSON.stringify(observed))return undefined;
  const key=observed[0]+':'+(observed[1]?observed[1]+'/':'')+observed[2];
  if(fresh.AcquisitionKey!==key)return undefined;
  const matches=rows(experimentRows).filter(r=>r.AcquisitionKey===key||JSON.stringify(tuple(r))===JSON.stringify(observed)||JSON.stringify(['ObservedSourceID','ObservedFolder','ObservedExpno'].map(k=>text(r[k])))===JSON.stringify(observed));
  if(matches.length!==1||matches[0].id!==fresh.id)return undefined;
  return Number.isFinite(now)&&now>0?now:undefined;
 }
 function action(id,base,draft,fresh,experimentRows=[],sampleRows=[],now=Date.now()/1000){
  if(!Number.isInteger(id)||id<=0||fresh?.id!==id)throw Error('The selected experiment is unavailable.');
  if(!equal(base,fresh)||!sameIdentity(base,fresh))throw Error('This experiment changed elsewhere. Discard this draft and review the current ratings before saving.');
  const fields=snapshot(draft),imported=fields.ImportedAssessment;
  delete fields.ImportedAssessment;
  for(const [key,label]of criteria)if(!choicesFor(key,fresh).includes(fields[key])&&!(key==='QCAnalysis'&&restrictedAnalysis.includes(fields[key])&&fields[key]===snapshot(base)[key]))throw Error('Choose a valid rating for '+label+'.');
  if(imported&&keys.some(key=>fields[key]!==base[key]))fields.ImportedAssessment='';
  const registered=registrationTimestamp(base,draft,fresh,experimentRows,sampleRows,now);
  if(registered!==undefined)fields.RegisteredAt=registered;
  return ['UpdateRecord','Experiments',id,fields];
 }
 function context(sectionRows,tableRows,id){
  const section=rows(sectionRows).find(r=>r.id===id),table=rows(tableRows).find(r=>r.tableId==='Experiments');
  if(!section||!table||section.tableRef!==table.id||!['custom','creator'].includes(section.parentKey))throw Error('Assessment widget context changed; reload before saving.');
  const options=JSON.parse(section.options||'{}');
  if(typeof options.html!=='string'||!options.html.includes('data-nmr-assessment="1"'))throw Error('This is not the expected assessment widget.');
 }
 return {criteria,keys,rows,snapshot,equal,sameIdentity,tone,decisionTone,registrationTimestamp,action,context,choicesFor,choiceLabel,relaxationType};
})();
'''.replace('__CRITERIA__', json.dumps(CRITERIA, ensure_ascii=False))

UI_JS = r'''
const SECTION_ID=__SECTION_ID__,$=id=>document.getElementById(id);
const el=(tag,text)=>{const node=document.createElement(tag);if(text!==undefined)node.textContent=text;return node;};
let current=null,base=null,draft=null,pending=undefined,busy=false,selectionSequence=0;
const dirty=()=>current&&draft&&!AssessmentCard.equal(base,draft);
function message(text,error=false){$('message').textContent=text;$('message').className=error?'error':'muted';}
function render(){
 $('empty').hidden=Boolean(current);$('card').hidden=!current;if(!current)return;
 $('experiment').textContent=[current.SnapshotSampleID||'Sample unassigned',current.ExperimentType||'Experiment'].join(' · ');
 $('source').textContent=[current.Folder,current.Expno].filter(Boolean).join(' / ')+([current.ExperimentCategory,current.Instrument,current.EffectiveAcquisitionStatus||'Unknown'].filter(Boolean).length?' · '+[current.ExperimentCategory,current.Instrument,current.EffectiveAcquisitionStatus||'Unknown'].filter(Boolean).join(' · '):'');
 const field=Number(current.Field_MHz),salt=current.SnapshotNaCl_mM,date=current.AcquiredAt;let day='';
 if(date){const parsed=new Date(typeof date==='number'?date*1000:date);if(Number.isFinite(parsed.getTime()))day=parsed.toISOString().slice(0,10);}
 if(!day)day=String(current.AcquisitionDateText||'').trim();
 $('conditions').textContent=[Number.isFinite(field)&&field>0?Math.round(field)+' MHz':'',salt!==null&&salt!==undefined&&String(salt).trim()!==''?salt+' mM NaCl':'',day,current.PulseProgram,current.Temperature_C!==null&&current.Temperature_C!==undefined?current.Temperature_C+' °C':''].filter(Boolean).join(' · ');
 const title=typeof current.Title==='string'?current.Title:'';$('titleText').textContent=title;$('titleBlock').hidden=!title;
 const evidence=['ReconciliationReason','MetadataDifference','ObservedAcquisitionEvidence','TitleSource'].filter(k=>current[k]).map(k=>String(current[k]));$('warnings').hidden=!evidence.length;$('evidence').textContent=evidence.join(' · ');
 $('decision').textContent=current.AssessmentDisplay||current.Assessment||'Assessment incomplete';$('decision').className='decision '+AssessmentCard.decisionTone(current);
 $('explanation').textContent=current.AssessmentExplanation||'Rate each criterion explicitly; the overall decision is calculated from the saved ratings.';
 for(const [key]of AssessmentCard.criteria){const select=$(key),value=draft[key],choices=AssessmentCard.choicesFor(key,current);select.replaceChildren();for(const choice of choices){const option=el('option',AssessmentCard.choiceLabel(key,choice,current));option.value=choice;select.append(option);}if(!choices.includes(value)){const option=el('option','Saved rating (unavailable for this experiment): '+AssessmentCard.choiceLabel(key,value,current));option.value=value;select.append(option);}select.value=value;select.className=AssessmentCard.tone(key,value);select.disabled=busy;}
 $('analysisHelp').hidden=!AssessmentCard.relaxationType(current.ExperimentType);
 $('note').value=draft.AssessmentReason;$('note').disabled=busy;
 $('save').disabled=busy||!SECTION_ID||!dirty();$('discard').disabled=busy||(!dirty()&&pending===undefined);
 $('switch').hidden=pending===undefined;$('switch').textContent=pending===undefined?'':'Another selection or update is waiting. Save this draft, or discard it to continue.';
 $('draft').textContent=dirty()?'Unsaved changes — the decision updates after Save.':'';
}
function adopt(record){current=record&&Number.isInteger(record.id)&&record.id>0?{...record}:null;base=current?{...current}:null;draft=base?AssessmentCard.snapshot(base):null;pending=undefined;render();}
function receiveFresh(record){if(busy||dirty()){if(record?.id===current?.id&&AssessmentCard.equal(record,base)&&AssessmentCard.sameIdentity(record,base)){current={...record};base={...record};render();return;}pending=record||null;render();return;}adopt(record);}
async function receive(record){const sequence=++selectionSequence,id=record?.id;if(!Number.isInteger(id)||id<=0){receiveFresh(null);return;}
 try{const fresh=AssessmentCard.rows(await grist.docApi.fetchTable('Experiments')).find(r=>r.id===id);if(sequence===selectionSequence)receiveFresh(fresh||null);}
 catch(error){if(sequence===selectionSequence)message('Could not read the selected experiment. Your current draft is retained.',true);}}
async function save(){if(!current||!dirty()||busy||!SECTION_ID)return;
 const requestedId=current.id,requestedBase={...base},requestedDraft={...draft};busy=true;message('Saving…');render();
 try{
  const [experiments,sections,tables,samples]=await Promise.all([grist.docApi.fetchTable('Experiments'),grist.docApi.fetchTable('_grist_Views_section'),grist.docApi.fetchTable('_grist_Tables'),grist.docApi.fetchTable('Samples')]);
  AssessmentCard.context(sections,tables,SECTION_ID);
  const fresh=AssessmentCard.rows(experiments).find(r=>r.id===requestedId);
  if(fresh&&(!AssessmentCard.equal(requestedBase,fresh)||!AssessmentCard.sameIdentity(requestedBase,fresh))&&pending===undefined)pending=fresh;
  const action=AssessmentCard.action(requestedId,requestedBase,requestedDraft,fresh,experiments,samples);
  await grist.docApi.applyUserActions([action]);
  const updated=AssessmentCard.rows(await grist.docApi.fetchTable('Experiments')).find(r=>r.id===requestedId);
  const next=pending!==undefined&&pending?.id!==requestedId?pending:updated;
  busy=false;adopt(next);message('Assessment saved.');
 }catch(error){message(error.message||'Assessment could not be saved. Your draft is retained.',true);}
 finally{busy=false;render();}
}
for(const [key,label]of AssessmentCard.criteria){const wrapper=el('label'),caption=el('span',label),select=el('select');select.id=key;
 select.addEventListener('change',()=>{draft[key]=select.value;message('');render();});wrapper.append(caption,select);$('ratings').append(wrapper);}
$('note').addEventListener('input',()=>{draft.AssessmentReason=$('note').value;render();});
$('save').addEventListener('click',save);$('discard').addEventListener('click',()=>{const next=pending!==undefined?pending:current;adopt(next);message('Draft discarded.');});
grist.ready({requiredAccess:'full'});grist.onRecord(receive);render();
'''


def assessment_html(section_id=None):
    """An unbound preview is read-only; installation binds the actual section ID."""
    if section_id is not None and (isinstance(section_id, bool) or not isinstance(section_id, int) or section_id <= 0):
        raise ValueError('section_id must be a positive integer')
    script = CORE_JS+'\n'+UI_JS.replace('__SECTION_ID__', json.dumps(section_id))
    return r'''<!doctype html><html lang="en" data-nmr-assessment="1"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><script src="https://docs.getgrist.com/grist-plugin-api.js"></script><style>
*{box-sizing:border-box}
body{max-width:800px;margin:0 auto;padding:12px;font:13px var(--grist-theme-font-family,system-ui,sans-serif);color:var(--grist-theme-text,#262633);background:var(--grist-theme-widget-bg,#fff)}
h2{font-size:15px;font-weight:600;margin:0 0 4px}h3{font-size:13px;font-weight:600;margin:8px 0 4px}
p{line-height:1.45;margin:5px 0}.muted{color:var(--grist-theme-text-medium,#666);font-size:11px}
.decision{display:inline-block;font-weight:600;border:1px solid var(--grist-theme-input-border,#d9d9d9);border-radius:4px;padding:5px 9px;margin:8px 0 2px}.pass{background:#e7f3eb;border-color:#8eb69d;color:#244d33}.limited{background:#fff4dc;border-color:#d8bb77;color:#5c4300}.fail{background:#fbeae8;border-color:#ce9994;color:#743125}.unchecked{background:#f6f8fa;color:#273943}.error{color:var(--grist-theme-text-error,#d0021b)}
summary{cursor:pointer;font-size:11px;color:var(--grist-theme-text-medium,#666)}
.ratings{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:8px;margin:12px 0}
label{display:grid;gap:5px;font-size:12px;color:var(--grist-theme-text-medium,#494949)}
select,textarea,button{font:inherit;color:inherit}select{height:32px;width:100%;min-width:0;border:1px solid var(--grist-theme-input-border,#d9d9d9);border-radius:4px;padding:6px 8px;background:var(--grist-theme-widget-bg,#fff);color:var(--grist-theme-text,#262633)}
textarea{box-sizing:border-box;width:100%;min-height:52px;resize:vertical;border:1px solid var(--grist-theme-input-border,#d9d9d9);border-radius:4px;padding:6px 8px;background:var(--grist-theme-widget-bg,#fff);color:var(--grist-theme-text,#262633)}
button{min-height:32px;border:0;border-radius:4px;background:var(--grist-theme-control-primary-bg,#16b378);color:var(--grist-theme-control-primary-fg,#fff);padding:4px 12px;cursor:pointer}button:disabled{opacity:.6;cursor:wait}#discard{background:var(--grist-theme-widget-bg,#fff);color:var(--grist-theme-control-fg,#16b378);border:1px solid var(--grist-theme-control-primary-bg,#16b378)}
.actions{display:flex;gap:8px;align-items:center;margin-top:12px}#switch{padding:8px;background:#fff4dc;border-left:3px solid #d8bb77}#draft{min-height:16px}
#titleText{white-space:pre-wrap;overflow-wrap:anywhere;max-height:10em;overflow:auto;border:1px solid var(--grist-theme-input-border,#d9d9d9);border-radius:4px;padding:7px;margin:0;font:11px/1.4 ui-monospace,monospace}
@media(max-width:480px){.ratings{grid-template-columns:1fr}}
</style></head><body><p id="empty" class="muted">Select an experiment in the acquisition list.</p><main id="card" hidden><h2 id="experiment"></h2><p id="conditions" class="muted"></p><p id="source" class="muted"></p><section id="titleBlock" hidden><h3>Acquisition title and comments</h3><pre id="titleText"></pre></section><details id="warnings" hidden><summary>Acquisition notes</summary><p id="evidence" class="muted"></p></details><span id="decision" class="decision" aria-label="Calculated overall assessment"></span><p id="explanation" class="muted"></p><div id="ratings" class="ratings"></div><p id="analysisHelp" class="muted" hidden>Both Good fit options are fully acceptable. Usable fit is a caution; Unusable fit excludes the experiment.</p><label>Assessment note (optional)<textarea id="note" rows="2"></textarea></label><p id="draft" class="muted"></p><p id="switch" hidden></p><div class="actions"><button id="save" type="button">Save assessment</button><button id="discard" type="button">Discard changes</button></div></main><p id="message" class="muted" role="status"></p><script>''' + script + '</script></body></html>'
