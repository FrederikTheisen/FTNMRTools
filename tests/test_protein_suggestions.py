"""Free-text protein suggestions and narrow native Choice metadata updates."""
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import unittest

SOURCE=Path(__file__).resolve().parents[1]/'nmr_register'/'registration_widget.py'
SPEC=importlib.util.spec_from_file_location('protein_registration_widget',SOURCE)
MODULE=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

CHECKS_JS=r'''
const checks=[];const assert=(ok,why='assertion failed')=>{if(!ok)throw Error(why);};
const equal=(a,b)=>assert(JSON.stringify(a)===JSON.stringify(b),JSON.stringify(a)+' != '+JSON.stringify(b));
const check=(name,fn)=>{fn();checks.push(name);};const throws=fn=>{let caught=false;try{fn();}catch(e){caught=true;}assert(caught);};
const MT=[{id:41,tableId:'Samples'},{id:42,tableId:'Series'}];
const mc=()=>[{id:61,parentId:41,colId:'Protein',type:'Choice',isFormula:false,widgetOptions:JSON.stringify({choices:['Existing','Exact Case'],choiceOptions:{Existing:{fillColor:'#123456'}},alignment:'left',otherOption:{keep:true}})}, {id:62,parentId:42,colId:'Protein',type:'Choice',isFormula:false,widgetOptions:JSON.stringify({choices:['Existing','Series only'],choiceOptions:{'Series only':{textColor:'#abcdef'}},alignment:'right'})}];
const plan=protein=>({actions:[['AddRecord','Samples',3,{Protein:protein}]],sample_ids:['tube']});
check('suggest exact union of Samples Series and frozen snapshots',()=>{equal(NMRRegistration.proteinSuggestions({Samples:[{Protein:'HeTDR1'},{Protein:' HeTDR1 '},{Protein:''},{Protein:null}],Series:[{Protein:'HeTDR1'},{Protein:'Series only'}],Experiments:[{SnapshotProtein:'Frozen only'},{SnapshotProtein:'hetDR1'},{Protein:'not an experiment snapshot'}]}),[' HeTDR1 ','Frozen only','HeTDR1','Series only','hetDR1']);});
check('new arbitrary protein remains free text in sample data',()=>{const p=NMRRegistration.sampleActions({Protein:'New protein / β + variant',Prepared:'2026-11',Isotope:'15N'}, {Samples:[],Series:[],Experiments:[]});assert(p.actions[0][3].Protein==='New protein / β + variant');});
check('sample-only save updates both native Choice lists preserving options',()=>{const columns=mc(),actions=NMRRegistration.proteinChoiceActions(plan('New protein / β'),MT,columns);assert(actions.length===2);for(const [i,action]of actions.entries()){equal(action.slice(0,3),['UpdateRecord','_grist_Tables_column',61+i]);equal(Object.keys(action[3]),['widgetOptions']);const before=JSON.parse(columns[i].widgetOptions),after=JSON.parse(action[3].widgetOptions);equal(after,{...before,choices:[...before.choices,'New protein / β']});}assert(JSON.parse(columns[0].widgetOptions).choices.length===2);});
check('blank known and experiment plans never expand native choices',()=>{equal(NMRRegistration.proteinChoiceActions(plan(''),[],[]),[]);equal(NMRRegistration.proteinChoiceActions(plan('Existing'),MT,mc()),[]);equal(NMRRegistration.proteinChoiceActions({actions:[['AddRecord','Experiments',3,{Protein:'new'}]]},[],[]),[]);});
check('metadata IDs types and widget options are checked before writing',()=>{for(const bad of [{type:'Text'},{isFormula:true},{id:0},{widgetOptions:'not json'},{widgetOptions:'[]'},{widgetOptions:JSON.stringify({choices:[1]})}])throws(()=>NMRRegistration.proteinChoiceActions(plan('New'),MT,[{...mc()[0],...bad},mc()[1]]));throws(()=>NMRRegistration.proteinChoiceActions(plan('New'),MT.slice(0,1),mc()));throws(()=>NMRRegistration.proteinChoiceActions(plan('New'),MT,[...mc(),mc()[0]]));});
check('fresh choices and colors are retained without normalization',()=>{const columns=mc();columns[0].widgetOptions=JSON.stringify({choices:['Existing','new concurrent choice',' spaced '],choiceOptions:{'new concurrent choice':{fillColor:'#fedcba'}},wrap:true});const actions=NMRRegistration.proteinChoiceActions(plan('Case-sensitive / name'),MT,columns);const result=JSON.parse(actions[0][3].widgetOptions);equal(result.choices,['Existing','new concurrent choice',' spaced ','Case-sensitive / name']);assert(result.choiceOptions['new concurrent choice'].fillColor==='#fedcba'&&result.wrap===true);});
globalThis.proteinTestResults={passed:checks.length,checks};
'''


def core_check_script():
    return MODULE.JS_LOGIC+'\n'+CHECKS_JS


BROWSER_CHECKS=r'''
setTimeout(async()=>{
 const report=document.createElement('output');report.id='protein-checks';document.body.append(report);
 const assert=(ok,why)=>{if(!ok)throw Error(why);},tick=()=>new Promise(r=>setTimeout(r,0));
 const submit=async()=>{byId('form').dispatchEvent(new Event('submit',{cancelable:true,bubbles:true}));for(let i=0;i<100&&busy;i++)await tick();assert(!busy,'save completed');};
 try{
  assert(controls.Protein.tagName==='INPUT'&&controls.Protein.getAttribute('list')==='suggest_Protein','protein is free text with datalist');
  const suggestions=[...document.querySelectorAll('#suggest_Protein option')].map(o=>o.value);
  assert(suggestions.includes('Existing')&&suggestions.includes('Series only')&&suggestions.includes('Frozen only')&&suggestions.includes('<Protein & control>'),'all exact names suggested through safe options');
  const initialWrites=WRITES.length,initialMetadataReads=METADATA_READS.length;
  controls.Protein.value='New protein / β';controls.Construct.value='draft construct';controls.Prepared.value='2026-11';controls.Label15N.checked=true;controls.pH.value='0';
  await refreshChoices();assert(controls.Protein.value==='New protein / β'&&controls.Construct.value==='draft construct'&&controls.pH.value==='0','refresh retains arbitrary text and other draft fields');
  assert(WRITES.length===initialWrites&&METADATA_READS.length===initialMetadataReads,'refresh only reads register rows');
  await submit();const first=WRITES.at(-1);assert(first.length===3&&first[0][1]==='_grist_Tables_column'&&first[1][1]==='_grist_Tables_column'&&first[2][1]==='Samples','choice options and sample saved in one batch');
  assert(DATA.Samples.at(-1).Protein==='New protein / β','new name is stored without requiring a suggestion');
  for(const column of COLUMNS){const options=JSON.parse(column.widgetOptions);assert(options.choices.includes('New protein / β')&&options.choiceOptions.Existing.fillColor==='#123456'&&options.wrap===true,'both native lists retain existing styling and append new name');}
  assert([...document.querySelectorAll('#suggest_Protein option')].some(o=>o.value==='New protein / β'),'saved name appears on refreshed suggestions');
  await submit();assert(WRITES.at(-1).every(a=>a[0]==='AddRecord'&&a[1]==='Samples'),'known name does not rewrite metadata');
  const before=METADATA_READS.length;controls.mode.value='experiment';controls.mode.dispatchEvent(new Event('change'));await submit();assert(METADATA_READS.length===before&&WRITES.at(-1).every(a=>a[0]==='AddRecord'&&a[1]==='Experiments'),'experiment registration never expands protein choices');
  assert(WRITES.every(batch=>batch.every(a=>a[0]==='AddRecord'||(a[0]==='UpdateRecord'&&a[1]==='_grist_Tables_column'&&[61,62].includes(a[2])&&Object.keys(a[3]).join()==='widgetOptions'))),'metadata mutation whitelist');
  controls.mode.value='sample';controls.mode.dispatchEvent(new Event('change'));message('Protein browser checks passed.');
  report.textContent='Protein checks passed: exact suggestions, arbitrary free text, retained draft, atomic shared native choices, preserved formatting, known/experiment no expansion.';
 }catch(error){report.textContent='Protein checks FAILED: '+error.message;report.style.color='red';}
},250);
'''


def browser_harness():
    fake=r'''<script>
const DATA={Series:[{id:1,SeriesID:'Existing series',Protein:'Series only'}],Samples:[{id:1,SampleID:'Existing tube',Protein:'Existing'},{id:2,SampleID:'Special tube',Protein:'<Protein & control>'}],Experiments:[{id:1,SnapshotProtein:'Frozen only'}]};
const TABLES=[{id:41,tableId:'Samples'},{id:42,tableId:'Series'},{id:43,tableId:'Experiments'}];
const COLUMNS=[41,42].map((parentId,i)=>({id:61+i,parentId,colId:'Protein',type:'Choice',isFormula:false,widgetOptions:JSON.stringify({choices:['Existing','Series only'],choiceOptions:{Existing:{fillColor:'#123456'}},wrap:true})}));
const WRITES=[],METADATA_READS=[];const table=rs=>{const keys=[...new Set(['id',...rs.flatMap(r=>Object.keys(r))])];return Object.fromEntries(keys.map(k=>[k,rs.map(r=>r[k])]));};
window.grist={ready:()=>{},docApi:{fetchTable:async name=>{if(name.startsWith('_grist_'))METADATA_READS.push(name);return table(name==='_grist_Tables'?TABLES:name==='_grist_Tables_column'?COLUMNS:DATA[name]);},applyUserActions:async actions=>{WRITES.push(actions);for(const [kind,name,id,fields]of actions){if(kind==='UpdateRecord'){Object.assign(COLUMNS.find(c=>c.id===id),fields);}else DATA[name].push({id,...fields});}return {};}}};
</script>'''
    return MODULE.registration_html('sample').replace('<script src="grist-plugin-api.js"></script>',fake).replace('</body>','<script>'+BROWSER_CHECKS+'</script></body>')


class ProteinSuggestionTests(unittest.TestCase):
    def test_protein_stays_free_text_and_metadata_is_save_only(self):
        html=MODULE.registration_html('sample')
        self.assertIn("field('identity','Protein','Protein')",html)
        self.assertIn("['Protein','Instrument','PulseProgram']",html)
        self.assertIn('NMRRegistration.proteinSuggestions(tables)',html)
        self.assertIn("if(values.mode!=='experiment'){const [metadataTables,metadataColumns]",html)
        self.assertNotIn('innerHTML',html)
        self.assertIn('...NMRRegistration.proteinChoiceActions(plan,metadataTables,metadataColumns),...plan.actions',html)

    @unittest.skipUnless(shutil.which('node'),'Node CLI absent; use core_check_script in JS runtime or browser fixture')
    def test_core(self):
        result=subprocess.run(['node','-e',core_check_script()+'\nconsole.log(JSON.stringify(proteinTestResults))'],check=True,capture_output=True,text=True)
        self.assertEqual(json.loads(result.stdout)['passed'],6)


if __name__=='__main__':
    unittest.main()
