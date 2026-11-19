"""Recorded rack placement preserves encoding and restricts writes to two fields."""
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import unittest

SOURCE=Path(__file__).resolve().parents[1]/'nmr_register'/'storage_widget.py'
SPEC=importlib.util.spec_from_file_location('storage_widget_under_test',SOURCE)
MODULE=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(MODULE)

CHECKS_JS=r'''
const checks=[],assert=(ok,why='assertion failed')=>{if(!ok)throw Error(why);};
const check=(name,fn)=>{fn();checks.push(name);};
const throws=fn=>{let caught=false;try{fn();}catch(e){caught=true;}assert(caught);};
const sample={id:1,SampleID:'Tube 1',StorageRow:'A',StorageColumn:6,Available:'Usable',Notes:'preserve'},dead={id:2,SampleID:'Dead tube',StorageRow:'C',StorageColumn:2,Available:'Dead'};
check('stored A6 means letter column A and numbered row 6',()=>{const p=Rack.normalizePosition(sample);assert(p.key==='A6'&&p.row==='A'&&p.column===6);});
check('empty and malformed positions remain distinct',()=>{assert(!Rack.normalizePosition({}).issue);for(const r of [{StorageRow:'A'},{StorageColumn:6},{StorageRow:'Z',StorageColumn:1},{StorageRow:'A',StorageColumn:0},{StorageRow:'A',StorageColumn:true}])assert(Rack.normalizePosition(r).issue);});
check('all occupants include dead and ambiguous duplicates',()=>{const occupied=Rack.occupancy([sample,dead,{...sample,id:3}]);assert(occupied.get('A6').length===2&&occupied.get('C2')[0].Available==='Dead');});
check('placement writes only position fields',()=>{const before=JSON.stringify(sample),actions=Rack.placementActions(1,[sample,dead],'H11',sample);assert(JSON.stringify(actions)==='[["UpdateRecord","Samples",1,{"StorageRow":"H","StorageColumn":11}]]');assert(JSON.stringify(sample)===before);});
check('occupied dead slot never cleared by status',()=>throws(()=>Rack.placementActions(1,[sample,dead],'C2',sample)));
check('concurrent new occupant prevents placement',()=>throws(()=>Rack.placementActions(1,[sample,{...dead,StorageRow:'H',StorageColumn:11}],'H11',sample)));
check('concurrent move or identity change prevents placement',()=>{for(const fresh of [{...sample,StorageColumn:7},{...sample,SampleID:'different tube'}])throws(()=>Rack.placementActions(1,[fresh],'H11',sample));});
check('invalid or changed sample selection fails closed',()=>{throws(()=>Rack.placementActions(2,[sample,dead],'H11',sample));throws(()=>Rack.placementActions(1,[],'H11',sample));throws(()=>Rack.placementActions(1,[sample],'A0',sample));});
check('clear explicitly blanks row and nulls column only',()=>{const actions=Rack.clearActions(1,[sample],sample);assert(JSON.stringify(actions)==='[["UpdateRecord","Samples",1,{"StorageRow":"","StorageColumn":null}]]');assert(sample.Notes==='preserve'&&sample.Available==='Usable');});
check('physical geometry has exactly 88 usable positions',()=>{assert(Rack.slots.length===88&&new Set(Rack.slots).size===88);for(const key of ['A1','B1','G1','H1','A12','B12','G12','H12'])assert(Rack.unavailable(key)&&!Rack.slots.includes(key));});
check('unavailable positions reject forged placement and retain legacy occupancy',()=>{for(const key of ['A1','B1','G1','H1','A12','B12','G12','H12']){let error;try{Rack.placementActions(1,[sample],key,sample);}catch(e){error=e;}assert(error?.message.includes('box construction'));const record={...sample,StorageRow:key[0],StorageColumn:Number(key.slice(1))};assert(Rack.normalizePosition(record).issue.includes('Unavailable')&&Rack.occupancy([record]).get(key)[0]===record);}});
globalThis.rackTestResults={passed:checks.length,checks};
'''


def core_check_script():return MODULE.RACK_JS+'\n'+CHECKS_JS


BROWSER_CHECKS=r'''
setTimeout(async()=>{
 const output=document.createElement('output');output.id='rack-checks';document.body.append(output);
 const assert=(ok,why)=>{if(!ok)throw Error(why);},tick=()=>new Promise(r=>setTimeout(r,0));
 try{
  assert(rackTestResults.passed===11,'pure core cases');
  assert(current.id===1&&samples.length===7,'partial onRecord and filtered onRecords still fetch full Samples');
  assert(document.querySelectorAll('.rack-slot').length===96,'96 cells');assert(document.querySelectorAll('.rack-slot.unavailable:disabled').length===8,'8 unavailable disabled cells');
  assert([...document.querySelectorAll('thead th')].slice(1).map(x=>x.textContent).join('')==='ABCDEFGH','letters across');
  assert([...document.querySelectorAll('tbody th')].map(x=>x.textContent).join(',')==='1,2,3,4,5,6,7,8,9,10,11,12','numbers down');
  assert(document.querySelector('tbody tr:nth-child(6) td:first-of-type button').dataset.slot==='A6','A6 encoding and orientation');
  assert(document.querySelector('[data-slot="A6"]').classList.contains('current'),'current position highlighted');
  assert(!document.querySelector('[data-slot="C2"]').disabled&&document.querySelector('[data-slot="C2"]').classList.contains('dead'),'dead remains occupied and red');
  assert(!document.querySelector('[data-slot="C12"]').disabled&&document.querySelector('[data-slot="C12"]').classList.contains('ambiguous'),'ambiguous occupancy inspectable');
  assert(!$('sample').querySelector('img'),'sample labels displayed as text');
  const editing=current,base=baseline;
  for(const [slot,state,marker] of [['A6','occupied','•'],['D3','suspected','?'],['E4','degraded','D'],['C2','dead','×']]){const b=document.querySelector('[data-slot="'+slot+'"]');assert(b.classList.contains(state)&&b.textContent.includes(marker)&&b.getAttribute('aria-label').includes('('),'visible and accessible status '+state);b.click();}
  document.querySelector('[data-slot="A6"]').click();assert($('inspector').textContent.includes('Protein: Test protein')&&$('inspector').textContent.includes('Labeling: 15N')&&$('inspector').textContent.includes('Tube status: Usable')&&!$('inspector').querySelector('img'),'inspector fields use text');
  document.querySelector('[data-slot="C12"]').click();assert($('inspector').textContent.includes('Conflict 1')&&$('inspector').textContent.includes('Conflict 2'),'all conflicting occupants inspectable');
  assert(!WRITES.length&&current===editing&&baseline===base&&selectedId===1&&!destination,'inspection never changes editing selection, baseline, destination or writes');
  const fixture=document.createElement('div');let selections=0,inspections=0;
  Rack.mountRack(fixture,{records:RECORDS,onSelect:()=>selections++});assert(fixture.querySelector('[data-slot="C2"]').disabled,'placement contract blocks occupied');
  Rack.mountRack(fixture,{records:RECORDS,onSelect:()=>selections++,allowOccupiedSelect:true});fixture.querySelector('[data-slot="C2"]').click();assert(selections===1&&fixture.querySelector('[data-slot="C12"]').disabled,'registration occupied selection preserved and conflicts blocked');
  Rack.mountRack(fixture,{records:RECORDS,onSelect:()=>selections++,onInspect:()=>inspections++});fixture.querySelector('[data-slot="C12"]').click();assert(inspections===1&&selections===1,'inspection separate from placement callback');
  document.querySelector('[data-slot="B6"]').click();assert(!WRITES.length&&$('preview').textContent.includes('B6'),'click previews without autosave');
  await save(false);assert(RECORDS[0].StorageRow==='B'&&RECORDS[0].StorageColumn===6,'explicit move changes position');
  await save(true);assert(RECORDS[0].StorageRow===''&&RECORDS[0].StorageColumn===null,'explicit clear');
  document.querySelector('[data-slot="H11"]').click();const count=WRITES.length;RECORDS[1].StorageRow='H';RECORDS[1].StorageColumn=11;await save(false);assert(WRITES.length===count&&$('message').className==='error','fresh occupancy preflight');
  destination='';await refresh();document.querySelector('[data-slot="B7"]').click();
  const fetch=grist.docApi.fetchTable;let release;grist.docApi.fetchTable=name=>name==='Samples'&&!release?new Promise(resolve=>{release=()=>resolve(table(RECORDS));}):fetch(name);
  const saving=save(false);await tick();receive({id:2});await tick();release();await saving;grist.docApi.fetchTable=fetch;
  assert(WRITES.length===count,'selection changed during fresh read prevents write');
  assert(WRITES.every(batch=>batch.length===1&&batch[0][0]==='UpdateRecord'&&batch[0][1]==='Samples'&&Object.keys(batch[0][3]).sort().join(',')==='StorageColumn,StorageRow'),'exact two-field mutation whitelist');
  assert(RECORDS[0].Notes==='keep'&&RECORDS[0].Available==='Usable','status and notes unchanged');
  output.textContent='Rack checks passed: 11 core cases; 12×8 orientation; full-table occupancy; dead/conflicting slots; preview; explicit move/clear; concurrent occupancy and selection guards; two-field writes.';
 }catch(error){output.textContent='Rack checks FAILED: '+error.message;output.style.color='red';}
},150);
'''


def browser_harness():
    fake=r'''<script>
const RECORDS=[{id:1,SampleID:'Tube <img src=x> 1',StorageRow:'A',StorageColumn:6,Available:'Usable',Protein:'Test protein',Isotope:'15N',Notes:'keep'},{id:2,SampleID:'Dead tube',StorageRow:'C',StorageColumn:2,Available:'Dead'},{id:3,SampleID:'Conflict 1',StorageRow:'C',StorageColumn:12,Available:'Usable'},{id:4,SampleID:'Conflict 2',StorageRow:'C',StorageColumn:12,Available:'Usable'},{id:5,SampleID:'Incomplete position',StorageRow:'B',StorageColumn:null,Available:'Usable'},{id:6,SampleID:'Suspected tube',StorageRow:'D',StorageColumn:3,Available:'Suspected degradation'},{id:7,SampleID:'Degraded tube',StorageRow:'E',StorageColumn:4,Available:'Degraded'}];
const WRITES=[],table=rows=>Object.fromEntries([...new Set(rows.flatMap(r=>Object.keys(r)))].map(k=>[k,rows.map(r=>r[k])]));
window.grist={ready:settings=>window.readySettings=settings,onRecord:cb=>cb({id:1}),onRecords:cb=>cb([{id:1}]),docApi:{fetchTable:async name=>name==='Samples'?table(RECORDS):name==='_grist_Views_section'?{id:[28],tableRef:[9],parentKey:['custom'],options:[JSON.stringify({html:'<html data-nmr-storage="1">'})]}:{id:[9],tableId:['Samples']},applyUserActions:async actions=>{WRITES.push(actions);Object.assign(RECORDS.find(r=>r.id===actions[0][2]),actions[0][3]);return {};}}};
</script>'''
    return MODULE.storage_html(28).replace('<script src="https://docs.getgrist.com/grist-plugin-api.js"></script>',fake).replace('</body>','<script>'+CHECKS_JS+BROWSER_CHECKS+'</script></body>')


class StoragePackagingTests(unittest.TestCase):
    def test_bound_inline_widget_and_full_table_reads(self):
        html=MODULE.storage_html(28)
        self.assertIn('const SECTION_ID=28',html)
        self.assertIn('data-nmr-storage="1"',html)
        self.assertIn("grist.onRecords(changed)",html)
        self.assertIn("grist.docApi.fetchTable('Samples')",html)
        self.assertIn("if(selectedId!==id||current?.id!==id)",html)
        self.assertEqual(html.count('applyUserActions('),1)
        self.assertNotIn('innerHTML',html)
        for value in (True,0,-1,'28'):
            with self.assertRaises(ValueError):MODULE.storage_html(value)

    def test_physical_geometry(self):
        self.assertEqual(MODULE.UNAVAILABLE_SLOTS, frozenset({'A1','B1','G1','H1','A12','B12','G12','H12'}))
        self.assertEqual(len(MODULE.RACK_SLOTS),88)
        self.assertEqual(MODULE.RACK_SLOTS[0],'A2')
        self.assertNotIn('__RACK_SLOTS__',MODULE.RACK_JS)

    def test_shared_styles_and_renderer_contract(self):
        self.assertIn('mountRack',MODULE.RACK_JS)
        self.assertIn('StorageColumn:Number(match[2])',MODULE.RACK_JS)
        self.assertIn('onSelect(key,{StorageRow:letter,StorageColumn:n})',MODULE.RACK_JS)
        self.assertNotIn('applyUserActions',MODULE.RACK_JS)
        self.assertNotIn('body{',MODULE.RACK_CSS)

    @unittest.skipUnless(shutil.which('node'),'Node CLI absent; core cases also run in browser fixture')
    def test_core(self):
        result=subprocess.run(['node','-e',core_check_script()+'\nconsole.log(JSON.stringify(rackTestResults))'],check=True,text=True,capture_output=True)
        self.assertEqual(json.loads(result.stdout)['passed'],11)


if __name__=='__main__':unittest.main()
