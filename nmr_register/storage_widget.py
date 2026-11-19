"""Document-owned recorded-position picker; storage encoding is unchanged."""
import json

RACK_JS = r'''
const Rack=(()=>{
 const rows=t=>Array.isArray(t)?t.map(r=>r.fields?{...r.fields,id:r.id}:r):(t?.id||[]).map((id,i)=>Object.fromEntries(Object.keys(t).map(k=>[k,t[k][i]])));
 const blank=v=>v===null||v===undefined||typeof v==='string'&&!v.trim();
 function normalizePosition(record){
  const row=blank(record?.StorageRow)?'':String(record.StorageRow).trim().toUpperCase(),raw=record?.StorageColumn;
  const column=blank(raw)?null:typeof raw==='boolean'?NaN:Number(raw);
  if(!row&&column===null)return {key:'',row:'',column:null,issue:''};
  if(!/^[A-H]$/.test(row)||!Number.isInteger(column)||column<1||column>12)return {key:'',row,column:Number.isFinite(column)?column:null,issue:'Incomplete or invalid recorded position'};
  return {key:row+column,row,column,issue:''};
 }
 function occupancy(records){const result=new Map();for(const sample of rows(records)){const p=normalizePosition(sample);if(p.key){if(!result.has(p.key))result.set(p.key,[]);result.get(p.key).push(sample);}}return result;}
 const token=r=>JSON.stringify(['id','SampleID','StorageBox','StorageRow','StorageColumn'].map(k=>r?.[k]??null));
 const dead=r=>['Dead','Retired','Unavailable'].includes(r?.Available);
 const status=r=>dead(r)?'dead':r?.Available==='Degraded'?'degraded':r?.Available==='Suspected degradation'?'suspected':'occupied';
 const label=r=>String(r?.SampleID||'(unnamed sample — fix ID)');
 function validate(id,baseline,records){
  if(!Number.isInteger(id)||id<=0||baseline?.id!==id)throw Error('Select a sample before changing its recorded position.');
  const matches=rows(records).filter(r=>r.id===id);
  if(matches.length!==1||token(matches[0])!==token(baseline))throw Error('This sample or its position changed elsewhere. Reload the selected sample before saving.');
  return matches[0];
 }
 function placementActions(id,records,destination,baseline){
  const sample=validate(id,baseline,records),match=/^([A-H])([1-9]|1[0-2])$/.exec(String(destination));
  if(!match)throw Error('Choose a valid empty position.');
  const occupants=occupancy(records).get(destination)||[];
  if(occupants.length)throw Error('That recorded position is occupied. Choose another slot.');
  return [['UpdateRecord','Samples',sample.id,{StorageRow:match[1],StorageColumn:Number(match[2])}]];
 }
 function clearActions(id,records,baseline){const sample=validate(id,baseline,records);return [['UpdateRecord','Samples',sample.id,{StorageRow:'',StorageColumn:null}]];}
 function mountRack(element,{records=[],selectedId=null,destination='',onSelect=null,onInspect=null,allowOccupiedSelect=false}={}){
  const samples=rows(records),occupied=occupancy(samples),invalid=samples.filter(r=>normalizePosition(r).issue);
  const table=document.createElement('table');table.className='nmr-rack';table.setAttribute('aria-label','Recorded sample positions: columns A to H, rows 1 to 12');
  const caption=document.createElement('caption');caption.textContent='Recorded positions';table.append(caption);
  const head=document.createElement('thead'),header=document.createElement('tr'),corner=document.createElement('th');corner.setAttribute('aria-label','Numbered row');header.append(corner);
  for(const letter of 'ABCDEFGH'){const th=document.createElement('th');th.scope='col';th.textContent=letter;header.append(th);}head.append(header);table.append(head);
  const body=document.createElement('tbody');
  for(let n=1;n<=12;n++){const tr=document.createElement('tr'),number=document.createElement('th');number.scope='row';number.textContent=String(n);tr.append(number);
   for(const letter of 'ABCDEFGH'){
    const key=letter+n,occupants=occupied.get(key)||[],selected=occupants.some(r=>r.id===selectedId),ambiguous=occupants.length>1,isDead=occupants.some(dead),state=occupants.length?['dead','degraded','suspected','occupied'].find(s=>occupants.some(r=>status(r)===s)):'';
    const td=document.createElement('td'),button=document.createElement('button');button.type='button';button.dataset.slot=key;
    button.className=['rack-slot',occupants.length?'occupied':'empty',selected?'current':'',destination===key?'destination':'',state,ambiguous?'ambiguous':''].filter(Boolean).join(' ');
    button.textContent=key+(ambiguous?' !':isDead?' ×':state==='degraded'?' D':state==='suspected'?' ?':occupants.length?' •':'');button.title=key+(occupants.length?' — '+occupants.map(r=>label(r)+' ('+(r.Available||'status unrecorded')+')').join('; '):' — no sample recorded');
    if(selected)button.title+='; selected sample';if(ambiguous)button.title+='; conflicting occupancy — review required';
    button.setAttribute('aria-label',button.title);
    const canSelect=Boolean(onSelect)&&!ambiguous&&(!occupants.length||allowOccupiedSelect),canInspect=Boolean(onInspect)&&occupants.length>0;
    if(canSelect&&!canInspect)button.setAttribute('aria-pressed',String(destination===key));
    button.disabled=!canSelect&&!canInspect;
    if(!button.disabled)button.addEventListener('click',()=>{if(canInspect)onInspect(key,occupants);else if(canSelect)onSelect(key,{StorageRow:letter,StorageColumn:n});});td.append(button);tr.append(td);
   }body.append(tr);
  }table.append(body);const legend=document.createElement('p');legend.className='rack-legend';legend.textContent='White: empty · Blue •: occupied · Yellow ?: suspected degradation · Orange D: degraded · Red ×: dead · !: conflict';element.replaceChildren(table,legend);return {occupancy:occupied,invalid};
 }
 return {rows,normalizePosition,occupancy,token,dead,label,placementActions,clearActions,mountRack};
})();
'''

RACK_CSS = r'''
.nmr-rack{border-collapse:separate;border-spacing:3px;width:100%;table-layout:fixed;font:12px system-ui,sans-serif;color:var(--grist-theme-text,#262633)}
.nmr-rack caption{text-align:left;font-size:13px;font-weight:600;padding:0 0 7px}.nmr-rack th{text-align:center;font-weight:600;height:21px}.nmr-rack th:first-child{width:25px}.nmr-rack td{padding:0}
.nmr-rack .rack-slot{display:block;width:100%;height:29px;padding:0 1px;font:inherit;border:1px solid #cbd4dc;border-radius:4px;background:#fff;color:#334155;cursor:pointer;white-space:nowrap}
.rack-legend{font:11px system-ui,sans-serif;line-height:1.4;margin:5px 0}.nmr-rack .rack-slot.occupied{background:#bfdbfe;color:#1e3a8a}.nmr-rack .rack-slot.suspected{background:#fef08a;color:#713f12}.nmr-rack .rack-slot.degraded{background:#fdba74;color:#7c2d12}.nmr-rack .rack-slot.dead{background:#F8B4B4;color:#7F1D1D;border-color:#dca99b}.nmr-rack .rack-slot.ambiguous{border:2px dashed #7c2d12}.nmr-rack .rack-slot.current{outline:2px solid #426ea2;outline-offset:-2px}.nmr-rack .rack-slot.destination{background:#dceafb;border-color:#426ea2}.nmr-rack .rack-slot:focus-visible{outline:3px solid #1d4ed8;outline-offset:1px}.nmr-rack .rack-slot:disabled{cursor:default;opacity:1}
'''

UI_JS = r'''
const SECTION_ID=__SECTION_ID__,$=id=>document.getElementById(id);
let samples=[],current=null,baseline=null,selectedId=null,destination='',busy=false,sequence=0,scheduled=false,inspectedSlot='';
function message(text,error=false){$('message').textContent=text;$('message').className=error?'error':'muted';}
function renderInspector(){
 const panel=$('inspector');panel.replaceChildren();if(!inspectedSlot){panel.textContent='Click an occupied slot to inspect its samples.';return;}
 const occupants=Rack.occupancy(samples).get(inspectedSlot)||[];
 const heading=document.createElement('strong');heading.textContent=inspectedSlot+' — '+(occupants.length>1?'conflicting sample records':occupants.length?'sample details':'no sample recorded');panel.append(heading);
 for(const sample of occupants){const line=document.createElement('p');line.textContent=[Rack.label(sample),'Protein: '+(sample.Protein||'not recorded'),'Labeling: '+(sample.Isotope||'not recorded'),'Tube status: '+(sample.Available||'not recorded')].join(' · ');panel.append(line);}
}
function render(){
 $('sample').textContent=current?Rack.label(current):'Select a sample in the sample list.';
 const position=current?Rack.normalizePosition(current):null;
 $('position').textContent=current?'Recorded position: '+(position.key||position.issue||'not assigned'):'';
 $('preview').textContent=destination?'Place '+Rack.label(current)+' at '+destination+'?':'';
 const result=Rack.mountRack($('rack'),{records:samples,selectedId:current?.id,destination,onInspect:(key)=>{inspectedSlot=key;renderInspector();},onSelect:busy||!current?null:(key)=>{destination=key;message('');render();$('rack').querySelector('[data-slot="'+key+'"]').focus();}});
 renderInspector();
 const conflicts=[...result.occupancy.values()].filter(v=>v.length>1).length;
 $('issues').textContent=[conflicts?conflicts+' positions have conflicting sample records.':'',result.invalid.length?result.invalid.length+' samples have incomplete or invalid positions.':''].filter(Boolean).join(' ');
 $('place').disabled=busy||!SECTION_ID||!current||!destination;
 $('place').textContent=position?.key?'Move':'Place';
 $('clear').disabled=busy||!SECTION_ID||!current||!position||(!position.key&&!position.issue);
}
function adopt(records,id){
 samples=Rack.rows(records);const fresh=samples.find(r=>r.id===id)||null;
 if(current?.id===id&&destination&&Rack.token(fresh)!==Rack.token(baseline)){message('The sample position changed elsewhere. Reload before saving.',true);render();return;}
 current=fresh;baseline=fresh?{...fresh}:null;render();
}
async function refresh(){const requested=++sequence,id=selectedId;
 try{const records=await grist.docApi.fetchTable('Samples');if(requested!==sequence)return;adopt(records,id);}
 catch(error){if(requested===sequence){samples=[];current=null;baseline=null;destination='';render();message('Could not read all recorded positions. Reload before choosing a slot.',true);}}
}
function receive(record){const next=Number.isInteger(record?.id)&&record.id>0?record.id:null;if(next!==selectedId){selectedId=next;current=null;baseline=null;destination='';render();}refresh();}
function changed(){if(scheduled)return;scheduled=true;queueMicrotask(()=>{scheduled=false;refresh();});}
function context(sections,tables){
 const section=Rack.rows(sections).find(r=>r.id===SECTION_ID),table=Rack.rows(tables).find(r=>r.tableId==='Samples');
 if(!section||!table||section.tableRef!==table.id||!['custom','creator'].includes(section.parentKey))throw Error('Rack widget context changed; reload before saving.');
 const options=JSON.parse(section.options||'{}');if(typeof options.html!=='string'||!options.html.includes('data-nmr-storage="1"'))throw Error('Unexpected rack widget context.');
}
async function save(clear=false){
 if(busy||!SECTION_ID||!current||(!clear&&!destination))return;
 const id=current.id,base={...baseline},target=destination;busy=true;message('Checking recorded positions…');render();
 try{
  const [records,sections,tables]=await Promise.all([grist.docApi.fetchTable('Samples'),grist.docApi.fetchTable('_grist_Views_section'),grist.docApi.fetchTable('_grist_Tables')]);
  context(sections,tables);
  if(selectedId!==id||current?.id!==id)throw Error('Sample selection changed. Choose the intended sample and position again.');
  const actions=clear?Rack.clearActions(id,records,base):Rack.placementActions(id,records,target,base);
  await grist.docApi.applyUserActions(actions);
  destination='';message(clear?'Recorded position cleared.':'Recorded position saved.');
  // The action remains bound to the captured ID even if selection changes during
  // the write. A fresh full-table fetch restores the current selection safely.
  await refresh();
 }catch(error){message(error.message||'Recorded position could not be saved.',true);}
 finally{busy=false;render();}
}
$('place').addEventListener('click',()=>save(false));$('clear').addEventListener('click',()=>save(true));
$('reload').addEventListener('click',()=>{destination='';message('');refresh();});
grist.ready({requiredAccess:'full'});grist.onRecord(receive);grist.onRecords(changed);render();
'''


def storage_html(section_id=None):
    """Bind the existing Samples-linked creator section ID to enable writes."""
    if section_id is not None and (type(section_id) is not int or section_id <= 0):
        raise ValueError('section_id must be a positive integer')
    script=RACK_JS+'\n'+UI_JS.replace('__SECTION_ID__',json.dumps(section_id))
    return r'''<!doctype html><html lang="en" data-nmr-storage="1"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><script src="https://docs.getgrist.com/grist-plugin-api.js"></script><style>
*{box-sizing:border-box}body{max-width:600px;margin:0 auto;padding:10px;font:12px system-ui,sans-serif;color:var(--grist-theme-text,#262633);background:var(--grist-theme-widget-bg,#fff)}h2{font-size:14px;margin:0 0 4px;overflow-wrap:anywhere}p{margin:5px 0;line-height:1.35}.muted{font-size:11px;color:#64748b}.error{color:#a32b20}.actions{display:flex;gap:7px;flex-wrap:wrap;margin-top:7px}.actions button{min-height:30px;border:1px solid #cbd4dc;border-radius:4px;background:#fff;padding:4px 9px;font:inherit;color:inherit;cursor:pointer}.actions #place{background:#e7f0fc;border-color:#9bb8dc}.actions button:disabled{opacity:.55;cursor:default}#preview{min-height:16px;font-weight:600}#issues{color:#8d3028}
'''+RACK_CSS+r'''</style></head><body><h2 id="sample"></h2><p id="position"></p><section id="inspector" aria-label="Inspected slot details" aria-live="polite"></section><div id="rack"></div><p id="issues" class="muted"></p><p id="preview" aria-live="polite"></p><div class="actions"><button id="place" type="button">Place</button><button id="clear" type="button">Clear position</button><button id="reload" type="button">Reload</button></div><p id="message" role="status" class="muted"></p><script>'''+script+'</script></body></html>'
