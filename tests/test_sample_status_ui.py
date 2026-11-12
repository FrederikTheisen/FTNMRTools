"""Shared tube status compatibility and optional Series registration controls."""
import importlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
MODULE=importlib.import_module('nmr_register.registration_widget')
PROTEIN_SPEC=importlib.util.spec_from_file_location('status_protein_tests',Path(__file__).with_name('test_protein_suggestions.py'))
PROTEIN=importlib.util.module_from_spec(PROTEIN_SPEC);PROTEIN_SPEC.loader.exec_module(PROTEIN)

CHECKS_JS=r'''
const checks=[];const assert=(ok,why='assertion failed')=>{if(!ok)throw Error(why);};
const check=(name,fn)=>{fn();checks.push(name);};const throws=fn=>{let caught=false;try{fn();}catch(e){caught=true;}assert(caught);};
const values={Protein:'P',Isotope:'15N',Prepared:'2026-11'},tables={Samples:[],Series:[],Experiments:[]};
check('new default is Usable and all four statuses are explicit',()=>{for(const status of [undefined,'Usable','Suspected degradation','Degraded','Dead']){const p=NMRRegistration.sampleActions({...values,Available:status},tables);assert(p.actions[0][3].Available===(status||'Usable'));assert(!('Assessment'in p.actions[0][3]));}});
check('legacy Retired and Unavailable values map to Dead without accepting custom status',()=>{for(const [old,expected]of [['Available','Usable'],['Unknown','Usable'],['','Usable'],[null,'Usable'],['Retired','Dead'],['Unavailable','Dead']])assert(NMRRegistration.sampleStatus(old)===expected);for(const custom of ['custom','usable',' usable ',false]){throws(()=>NMRRegistration.sampleStatus(custom));assert(!NMRRegistration.selectableSample({Available:custom},true));}});
check('default picker keeps concerns and excludes dead samples with explicit opt in',()=>{for(const status of ['Usable','Suspected degradation','Degraded','Available','Unknown',''])assert(NMRRegistration.selectableSample({Available:status}));for(const status of ['Dead','Retired','Unavailable']){assert(!NMRRegistration.selectableSample({Available:status}));assert(NMRRegistration.selectableSample({Available:status},true));}});
check('new experiment linking dead samples requires explicit opt in, custom fails',()=>{for(const status of ['Dead','Retired','Unavailable']){const data={Samples:[{id:1,Protein:'P',Available:status}],Series:[],Experiments:[]};throws(()=>NMRRegistration.experimentActions({Sample:'1'},data));assert(NMRRegistration.experimentActions({Sample:'1',IncludeUnavailable:true},data).actions[0][3].Sample===1);}throws(()=>NMRRegistration.experimentActions({Sample:'1',IncludeUnavailable:true},{Samples:[{id:1,Available:'custom'}]}));});
check('bulk creation still shares defaults and never derives assessment',()=>{const p=NMRRegistration.sampleActions({...values,mode:'bulk',Series:'new',NaCl_mM:'50',VaryField:'DNA_uM',VaryValues:'0,100'},tables);const samples=p.actions.filter(a=>a[1]==='Samples');assert(samples.length===2&&samples.every(a=>a[3].Available==='Usable'&&a[3].NaCl_mM===50&&a[3].Series===1&&!('Assessment'in a[3])));});
globalThis.statusTestResults={passed:checks.length,checks};
'''


def core_check_script():return MODULE.JS_LOGIC+'\n'+CHECKS_JS


BROWSER_CHECKS=r'''
setTimeout(async()=>{
 const output=document.createElement('output');output.id='status-checks';document.body.append(output);
 const assert=(ok,why)=>{if(!ok)throw Error(why);},tick=()=>new Promise(r=>setTimeout(r,0));
 const submit=async()=>{byId('form').dispatchEvent(new Event('submit',{bubbles:true,cancelable:true}));for(let i=0;i<100&&busy;i++)await tick();assert(!busy,'save completed');};
 const choices=()=>[...controls.Sample.options].map(o=>o.value).filter(Boolean);
 try{
  assert(controls.Available.value==='Usable','new sample default Usable');
  assert([...controls.Available.options].map(o=>o.value).join('|')==='Usable|Suspected degradation|Degraded|Dead','four status choices with no legacy labels');
  assert(!byId('seriesBlock').open&&controls.Series.closest('details')===byId('seriesBlock')&&!byId('top').contains(controls.Series),'Series is optional and collapsed');
  assert([...document.querySelectorAll('#suggest_Protein option')].some(o=>o.value==='Existing'),'protein suggestions preserved');
  controls.mode.value='experiment';controls.mode.dispatchEvent(new Event('change'));
  assert(choices().join('|')==='1|2|3|5|7|8','usable and concern/degraded picker rows retained, dead/custom excluded');
  assert(controls.Sample.querySelector('option[value="2"]').textContent.includes('Suspected degradation')&&controls.Sample.querySelector('option[value="3"]').textContent.includes('Degraded'),'concerns visible when selecting');
  controls.IncludeUnavailable.checked=true;controls.IncludeUnavailable.dispatchEvent(new Event('change'));
  assert(choices().join('|')==='1|2|3|4|5|6|7|8'&&groups.IncludeUnavailable.textContent.includes('Include dead samples'),'explicit dead opt in');
  const before=JSON.stringify(DATA.Experiments);
  controls.mode.value='bulk';controls.mode.dispatchEvent(new Event('change'));byId('seriesBlock').open=true;controls.Series.value='new';controls.Series.dispatchEvent(new Event('change'));
  controls.SeriesID.value='Optional shared series';controls.Protein.value='Existing';controls.Prepared.value='2026-11';controls.Label15N.checked=true;controls.NaCl_mM.value='50';controls.VaryValues.value='0,100';
  await submit();const batch=WRITES.at(-1),series=batch.find(a=>a[1]==='Series'),samples=batch.filter(a=>a[1]==='Samples');
  assert(series&&samples.length===2&&samples.every(a=>a[3].Available==='Usable'&&a[3].NaCl_mM===50&&a[3].Series===series[2]),'optional bulk series and defaults preserved');
  assert(samples[0][3].DNA_uM===0&&samples[1][3].DNA_uM===100,'bulk literal zero preserved');
  assert(JSON.stringify(DATA.Experiments)===before&&batch.every(a=>!('Assessment'in a[3])),'existing experiments and assessments unchanged');
  controls.Available.value='Degraded';controls.Available.dispatchEvent(new Event('change'));assert(getComputedStyle(controls.Available).backgroundColor==='rgb(252, 230, 223)','status colors match native palette');
  controls.Available.value='Dead';controls.Available.dispatchEvent(new Event('change'));assert(getComputedStyle(controls.Available).backgroundColor==='rgb(248, 180, 180)','Dead uses the same red status palette as the Samples view');
  assert(document.querySelectorAll('.nmr-rack tbody tr').length===12&&document.querySelectorAll('.nmr-rack thead th').length===9,'visual rack has 12 numbered rows and A-H columns');
  assert(!byId('storage').querySelector('#field_StorageRow')&&!byId('storage').querySelector('#field_StorageColumn'),'old row and column inputs are not visible');
  controls.mode.value='sample';controls.mode.dispatchEvent(new Event('change'));byId('storageBlock').open=true;
  const firstEmpty=document.querySelector('.nmr-rack [data-slot="A1"]');assert(firstEmpty&&!firstEmpty.disabled,'empty position is selectable');firstEmpty.click();assert(controls.StorageSlot.value==='A1','rack click captures visual position');
  controls.mode.value='bulk';controls.mode.dispatchEvent(new Event('change'));controls.VaryValues.value='0,100';controls.VaryValues.dispatchEvent(new Event('input'));controls.StorageSlot.value='A12';refreshRack();
  assert(byId('storagePreview').textContent.includes('A12, B1'),'bulk preview follows the existing A1…A12/B1 order');
  tables.Samples.push({id:99,SampleID:'Dead tube',Available:'Dead',StorageRow:'B',StorageColumn:1});refreshRack();
  assert(byId('submit').disabled&&byId('storagePreview').textContent.includes('Blocked: B1 is occupied'),'occupied Dead slot blocks bulk without skipping');
  byId('seriesBlock').open=false;message('Status browser checks passed.');
  output.textContent='Status checks passed: Usable default, four shared statuses, dead opt in, visible concerns, optional collapsed series, bulk defaults/zero, rack positions and no-skip preview, existing records unchanged, protein suggestions.';
 }catch(error){output.textContent='Status checks FAILED: '+error.message;output.style.color='red';}
},250);
'''


def browser_harness():
    html=PROTEIN.browser_harness()
    rows=[{'id':i,'SampleID':'Tube '+str(i),'Protein':'Existing','Available':status} for i,status in enumerate(['Usable','Suspected degradation','Degraded','Dead','Unknown','Unavailable','Available','','custom'],1)]
    encoded=json.dumps(rows).replace('<','\\u003c')
    html=html.replace('const TABLES=', 'DATA.Samples='+encoded+';const TABLES=')
    html=html.replace('<script>'+PROTEIN.BROWSER_CHECKS+'</script>', '<script>'+BROWSER_CHECKS+'</script>')
    return html


class SampleStatusTests(unittest.TestCase):
    def test_python_default_compatibility_and_fail_closed(self):
        base={'Protein':'P','Isotope':'15N','Prepared':'2026-11'}
        for legacy,expected in [(None,'Usable'),('Unknown','Usable'),('Available','Usable'),('Unavailable','Dead'),('Retired','Dead'),('Degraded','Degraded')]:
            fields=MODULE.sample_actions({**base,'Available':legacy},{})['actions'][0][3]
            self.assertEqual(fields['Available'],expected)
        with self.assertRaisesRegex(ValueError,'Invalid sample status'):
            MODULE.sample_actions({**base,'Available':'custom'},{})
        data={'Samples':[{'id':1,'Available':'Dead'}]}
        with self.assertRaisesRegex(ValueError,'dead'):
            MODULE.experiment_actions({'Sample':'1'},data)
        self.assertEqual(MODULE.experiment_actions({'Sample':'1','IncludeUnavailable':True},data)['actions'][0][3]['Sample'],1)
        data['Samples'][0]['Available']='custom'
        with self.assertRaisesRegex(ValueError,'Invalid sample status'):
            MODULE.experiment_actions({'Sample':'1','IncludeUnavailable':True},data)

    def test_optional_series_and_dead_label_are_in_safe_form(self):
        html=MODULE.registration_html('sample')
        self.assertIn('<details id="seriesBlock"><summary>Series (optional)</summary>',html)
        self.assertIn("field('seriesChooser','Series'",html)
        self.assertNotIn("field('top','Series'",html)
        self.assertIn("[['Usable','Usable'],['Suspected degradation','Suspected degradation'],['Degraded','Degraded'],['Dead','Dead']]",html)
        self.assertNotIn("['Retired','Retired']",html)
        self.assertIn('Include dead samples',html)
        self.assertIn("NMRRegistration.selectableSample(s,controls.IncludeUnavailable.checked)",html)

    @unittest.skipUnless(shutil.which('node'),'Node CLI absent; use JS runtime or disconnected browser fixture')
    def test_core(self):
        run=subprocess.run(['node','-e',core_check_script()+'\nconsole.log(JSON.stringify(statusTestResults))'],check=True,capture_output=True,text=True)
        self.assertEqual(json.loads(run.stdout)['passed'],5)


if __name__=='__main__':unittest.main()
