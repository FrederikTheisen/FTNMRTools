"""Pure matrix semantics; browser harness uses the same runtime core."""
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
# Load this checkout explicitly: an installed regular package takes precedence
# over a staging directory without __init__.py, even when PYTHONPATH points here.
_SOURCE=Path(__file__).resolve().parents[1]/'nmr_register'/'coverage_widget.py'
_SPEC=importlib.util.spec_from_file_location('coverage_widget_under_test',_SOURCE)
_MODULE=importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)
CORE_JS,coverage_html=_MODULE.CORE_JS,_MODULE.coverage_html


CHECKS_JS = r'''
const checks=[];
function check(name,fn){fn();checks.push(name);}
function assert(ok,message='assertion failed'){if(!ok)throw Error(message);}
function equal(a,b){assert(JSON.stringify(a)===JSON.stringify(b),JSON.stringify(a)+' != '+JSON.stringify(b));}
function throws(fn){let failed=false;try{fn();}catch(e){failed=true;}assert(failed,'expected rejection');}
function row(id,changes={}){return {id,Study:'Electrostatics',SnapshotProtein:'P',SnapshotNaCl_mM:50,EffectiveField_MHz:600,
 EffectiveType:'T1',EffectiveTemperature_C:25,EffectiveAcquisitionStatus:'Finished',SourceID:'source',Folder:'folder',Expno:String(id),
 ObservedSourceID:'source',ObservedFolder:'folder',ObservedExpno:String(id),AcquisitionKey:'source:folder/'+id,ObservedFingerprint:'fp'+id,
 Sample:1,SnapshotSampleID:'tube',SnapshotConstruct:'C',SnapshotIsotope:'15N',SnapshotBuffer:'NaP',SnapshotpH:6.5,SnapshotDNA_uM:0,SnapshotTFE_pct:0,
 SnapshotProtein_uM:100,SnapshotOtherConditions:'',EffectiveInstrument:'NMR600',EffectivePulseProgram:'reviewed',EffectiveD1_s:2,EffectiveDelayList_s:'[0,0.1]',
 EffectiveSpinLock_Hz:null,EffectiveSpinLockOffset_Hz:null,Assessment:'Use',QCAcquisitionCoverage:'Complete',QCSettings:'Appropriate',QCFID:'Clean',QCSNR:'Adequate',QCSample:'Acceptable',QCAnalysis:'Reliable',...changes};}
function build(records,prefs={},temp='25'){return Coverage.build(records,[],prefs,'P',temp);}
function cell(model,type='t1',field=600,salt=50){return model.cells.find(c=>c.type===type&&c.field===field&&c.salt===salt);}
check('Dummy excluded before identity indexes and observed axes',()=>{const real=row(1),dummy=row(2,{AcquisitionStatus:'Dummy',SnapshotProtein:'Dummy protein',EffectiveField_MHz:950,SnapshotNaCl_mM:999,AcquisitionKey:real.AcquisitionKey});const m=build([real,dummy]);equal(m.proteins,['P']);equal(m.observed.fields,[600]);equal(m.observed.salts,[50]);equal(m.records.map(r=>r.id),[1]);assert(!m.records[0]._issue);});
check('exact four defaults and holes',()=>{const m=build([row(1)]);equal(m.axes.types,['t1','t1rho','noe','eta_xy']);assert(cell(m).finished===1);assert(cell(m,'t1rho').finished===0);assert(m.cells.length===4);});
check('nominal proton field and explicit zero',()=>{const m=build([row(1,{EffectiveField_MHz:599.8,SnapshotNaCl_mM:0})]);equal(m.axes.fields,[600]);equal(m.axes.salts,[0]);assert(cell(m,'t1',600,0).finished===1);assert(Coverage.num('   ')===null);assert(Coverage.num(0)===0);});
check('strict study includes no zero-dose DNA or TFE',()=>{const m=build([row(1),row(2,{Study:'DNA',SnapshotNaCl_mM:0}),row(3,{Study:'TFE',SnapshotNaCl_mM:150})]);equal(m.axes.salts,[50]);assert(cell(m).finished===1);});
check('observed field and salt defaults span temperatures, saved axes take precedence',()=>{const records=[row(1),row(2,{EffectiveTemperature_C:30,EffectiveField_MHz:950,SnapshotNaCl_mM:500})];const m=build(records,{},'25');equal(m.axes.fields,[600,950]);equal(m.axes.salts,[50,500]);const saved=build(records,{proteins:{P:{fields:[700],salts:[0],types:['t1']}}},'25');equal(saved.axes.fields,[700]);equal(saved.axes.salts,[0]);});
check('one effective temperature at a time',()=>{const records=[row(1),row(2,{EffectiveTemperature_C:30}),row(3,{EffectiveTemperature_C:null})];assert(cell(build(records)).finished===1);assert(cell(build(records,{},'30')).finished===1);assert(cell(build(records,{},'unknown')).finished===1);});
check('derived whole and half groups preserve actual temperatures and holes',()=>{const records=[row(1,{Temperature_C:24.8}),row(2,{Temperature_C:25.2}),row(3,{Temperature_C:5.499,EffectiveTemperature_C:5.5}),row(4,{Temperature_C:5.5005,EffectiveTemperature_C:5.5}),row(5,{Temperature_C:5.25,EffectiveTemperature_C:5.25}),row(6,{Temperature_C:0,EffectiveTemperature_C:0}),row(7,{Temperature_C:null,EffectiveTemperature_C:null})],before=JSON.stringify(records);equal(build(records).temperatures,['0','5.25','5.5','25','unknown']);assert(cell(build(records)).good===2);assert(cell(build(records,{},'5.5')).good===2);assert(cell(build(records,{},'5.25')).good===1);assert(cell(build(records,{},'0')).good===1);assert(cell(build(records,{},'unknown')).good===1);assert(cell(build(records,{},'24.8')).good===0);assert(JSON.stringify(records)===before);});
check('stale selection reconciles only to an observed policy-equivalent value',()=>{for(const value of ['24.8','25.2','24.8498','25.2000000001'])assert(Coverage.reconcileTemperature(value,['25'])==='25');for(const value of ['5.499','5.5005','5.48','5.52'])assert(Coverage.reconcileTemperature(value,['5.5'])==='5.5');for(const [value,available]of [['25',['5.5']],['5.5',['5','6']],['5.25',['5','5.5']],['24.7999',['25']],['5.4799',['5.5']],['0',['25']],['unknown',['25']],['24.8498',['24.8498','25']],['24.8498',['30']]])assert(Coverage.reconcileTemperature(value,available)===value);});
check('manual Finished and unassigned sample count without raw',()=>{const m=build([row(1,{AcquisitionKey:'',RawPresent:false,Sample:0,RegisteredAt:null})]);assert(cell(m).finished===1);assert(cell(m).manual===1);assert(cell(m).unknownTubes===1);});
check('metadata warning does not erase done or bad assessment',()=>{const c=cell(build([row(1,{MetadataDifference:'temperature differs',Assessment:'Do not use',ResultCSV:['L',5],AnalysisStatus:'Analyzed'})]));assert(c.finished===1&&c.csv===1&&c.analyzed===1&&c.assessments['Do not use']===1);});
check('new workflow states and legacy aliases',()=>{equal(['Planned','Setup','Running','Completed','Aborted','Unknown','Finished','Stopped'].map(x=>Coverage.status({AcquisitionStatus:x})),['Pending / unknown','Pending / unknown','Pending / unknown','Completed','Aborted','Pending / unknown','Completed','Aborted']);});
check('stopped and pending are not done',()=>{const c=cell(build([row(1,{EffectiveAcquisitionStatus:'Stopped'}),row(2,{EffectiveAcquisitionStatus:'Planned'})]));assert(c.finished===0&&c.stopped===1&&c.pending===1);});
check('repeated record delivery is same acquisition count',()=>{const r=row(1);assert(cell(build([r,r])).finished===1);});
check('duplicate tuple key fingerprint and replacement require review',()=>{for(const records of [[row(1),row(2,{ObservedExpno:'1',Expno:'1',AcquisitionKey:'source:folder/1'})],[row(1),row(2,{ObservedFingerprint:'fp1'})],[row(1,{ReconciliationReason:'possible replacement'})]]){const c=cell(build(records));assert(c.finished===0);assert(c.review===records.length);}});
check('source mismatch requires review but missing sample alone does not',()=>{assert(cell(build([row(1,{Folder:'other'})])).review===1);assert(cell(build([row(1,{Sample:0,MigrationIssue:'Sample link unresolved'})])).finished===1);});
check('assessment source conflict is not acquisition identity failure',()=>{const text='Conflict: review OK; comment Looks degraded; Conflicting source evidence: TSV decision Use: Looks degraded';const c=cell(build([row(1,{Assessment:'Limited use',MigrationIssue:text,AnalysisStatus:'Analyzed'})]));assert(c.limited===1&&c.review===0&&c.analyzed===1);for(const issue of ['Source identity differs','Folder mismatch','Incomplete acquisition identity','duplicate acquisition','split acquisition','copied data','replacement']){const x=cell(build([row(1,{Assessment:'Limited use',MigrationIssue:issue})]));assert(x.review===1&&x.limited===0,issue);}});
check('unrecorded analysis includes discarded and excludes selected CSV',()=>{const c=cell(build([row(1,{Assessment:'Do not use',AnalysisStatus:'Not started'}),row(2,{Assessment:'',AnalysisStatus:''}),row(3,{Assessment:'Limited use',AnalysisStatus:'In progress'}),row(4,{Assessment:'Use',AnalysisStatus:'Analyzed'}),row(5,{Assessment:'',AnalysisStatus:'',ResultCSV:['L',5]}),row(6,{Assessment:'',ReconciliationReason:'possible copy'})]));assert(c.notAnalyzed===3&&c.unassessed===2&&c.assessments['Do not use']===1&&c.limited===1&&c.review===1);});
check('explicit incomplete analysis overrides an attached CSV',()=>{const c=cell(build([row(1,{AnalysisStatus:'Not started',ResultCSV:['L',2]}),row(2,{AnalysisStatus:'In progress',ResultCSV:['L',3]}),row(3,{AnalysisStatus:'',ResultCSV:['L',4]}),row(4,{AnalysisStatus:'Analyzed',ResultCSV:null})]));assert(c.notAnalyzed===2);});
check('same prep distinct unique acquisitions count twice',()=>{const c=cell(build([row(1),row(2,{Sample:2,SnapshotSampleID:'tube2'})]));assert(c.finished===2&&c.groups.length===1&&c.tubes===2);});
check('different and uncertain preparations remain separate groups',()=>{assert(cell(build([row(1),row(2,{SnapshotConstruct:'other'})])).groups.length===2);assert(cell(build([row(1,{SnapshotConstruct:null}),row(2,{SnapshotConstruct:null,Sample:2})])).groups.length===2);});
check('snapshots are never backfilled from current samples',()=>{const m=Coverage.build([row(1,{SnapshotProtein:''}),row(2,{SnapshotNaCl_mM:null})],[{id:1,Protein:'P',NaCl_mM:50}],{},'P','25');assert(m.unknownProtein===1);assert(m.unplaced.length===1);assert(m.axes.salts.length===0);});
check('initial protein favors populated Electrostatics over empty sample-only names',()=>{const m=Coverage.build([row(1),row(2,{SnapshotProtein:'Z'})],[{id:1,Protein:'A'},{id:2,Protein:'P'}],{},'P','25');equal(m.proteins,['P','Z','A']);});
check('sample-only protein permits empty configuration without data',()=>{const m=Coverage.build([], [{id:1,Protein:'CtIP',NaCl_mM:150}],{},'CtIP','25');equal(m.proteins,['CtIP']);equal(m.axes.fields,[]);equal(m.axes.salts,[]);});
check('expected extra axes create holes without records',()=>{const prefs={version:1,proteins:{P:{fields:[850],salts:[100,250],types:['etaXY','Water T1']}}};const m=build([row(1)],prefs);assert(m.cells.length===4&&m.cells.every(c=>c.finished===0));assert(m.unplaced.length===1);});
check('sequence aliases canonicalize with water optional',()=>{equal(['T1','T1ρ','hetNOE','ηXY','eta_xy','T1sat','water T1','water_t1sat','nico t1sat','Water T1 (T1sat)'].map(Coverage.type),['t1','t1rho','noe','eta_xy','eta_xy','water_t1sat','water_t1sat','water_t1sat','water_t1sat','water_t1sat']);equal(['HSQC','15N HSQC','15N_HSQC','13C HSQC','13c_hsqc'].map(Coverage.type),['hsqc','hsqc','hsqc','13c_hsqc','13c_hsqc']);assert(Coverage.label(Coverage.type('13C HSQC'))==='13C HSQC');assert(!build([row(1),row(2,{EffectiveType:'water_t1sat'})]).axes.types.includes('water_t1sat'));});
const originalAxes={fields:[600],salts:[50],types:['t1']};
function metadata(other={}){const prefs={coverageAxes:{version:1,proteins:{P:originalAxes,Q:{fields:[950],salts:[500],types:['noe']}}},otherPreference:'keep'};return {id:88,tableRef:42,parentKey:'custom',title:'Coverage matrix',options:JSON.stringify({html:'<html data-nmr-coverage="1">preserve all HTML</html>',customView:JSON.stringify({access:'full',widgetOptions:prefs,url:'keep'}),otherOption:23,...other})};}
check('durable save merges fresh unrelated protein settings and HTML',()=>{const s=metadata(),before=JSON.parse(s.options),custom=JSON.parse(before.customView);custom.widgetOptions.coverageAxes.proteins.Q.salts=[100,500];s.options=JSON.stringify({...before,customView:JSON.stringify(custom)});const update=Coverage.mergeAxes([s],[{id:42,tableId:'Experiments'}],88,'P',{fields:[850],salts:[0,50],types:['t1rho']},originalAxes);equal(update.action.slice(0,3),['UpdateRecord','_grist_Views_section',88]);const after=JSON.parse(update.action[3].options),savedCustom=JSON.parse(after.customView);assert(after.html===before.html&&after.otherOption===23&&savedCustom.url==='keep');equal(savedCustom.widgetOptions.coverageAxes.proteins.Q.salts,[100,500]);assert(savedCustom.widgetOptions.otherPreference==='keep');});
check('same protein concurrent changes and invalid section never write',()=>{const s=metadata();throws(()=>Coverage.mergeAxes([s],[{id:42,tableId:'Experiments'}],88,'P',originalAxes,{...originalAxes,salts:[100]}));for(const change of [{tableRef:99},{parentKey:'record'},{options:'not json'},{options:JSON.stringify({html:'other',customView:'{}'})}])throws(()=>Coverage.mergeAxes([{...s,...change}],[{id:42,tableId:'Experiments'}],88,'P',originalAxes,originalAxes));});
check('only trusted structured QC Use counts good',()=>{const c=cell(build([row(1),row(2,{Assessment:''}),row(3,{Assessment:'Limited use'}),row(4,{Assessment:'Do not use'}),row(5,{EffectiveAcquisitionStatus:'Planned',QCAcquisitionCoverage:'Not assessed',Assessment:''}),row(6,{ReconciliationReason:'possible copy'}),row(7,{Assessment:'',ResultCSV:['L',8],AnalysisStatus:'Analyzed'})]));assert(c.good===1&&c.maxCompatibleGood===1&&c.unassessed===3&&c.limited===1&&c.pending===1&&c.review===1&&c.assessments['Do not use']===1);equal(c.ids,[1,2,3,4,5,6,7]);});
check('two good plus bad is two not three',()=>{const c=cell(build([row(1),row(2),row(3,{Assessment:'Do not use'})]));assert(c.good===2&&c.maxCompatibleGood===2&&c.assessments['Do not use']===1);equal(c.goodCounts,[2]);});
check('trusted Use total is independent of compatibility subgroups',()=>{const protocol=cell(build([row(1),row(2,{EffectiveD1_s:9})]));assert(protocol.good===2&&protocol.maxCompatibleGood===1);equal(protocol.goodCounts,[1,1]);const c=cell(build([row(1),row(2,{SnapshotConstruct:'other'})]));assert(c.good===2&&c.maxCompatibleGood===1);equal(c.goodCounts,[1,1]);const matched=cell(build([row(1),row(2),row(3,{SnapshotConstruct:'other'})]));assert(matched.good===3&&matched.maxCompatibleGood===2);equal(matched.goodCounts,[2,1]);});
check('duplicate good identities and empty holes never count',()=>{const c=cell(build([row(1),row(2,{ObservedFingerprint:'fp1'})]));assert(c.good===0&&c.review===2);const empty=cell(build([row(1)]),'t1rho');assert(empty.good===0&&empty.maxCompatibleGood===0);equal(empty.ids,[]);});
check('Stopped sufficient and unknown complete passing QC count good',()=>{const c=cell(build([row(1,{EffectiveAcquisitionStatus:'Stopped',QCAcquisitionCoverage:'Sufficient'}),row(2,{EffectiveAcquisitionStatus:'Unknown'})]));assert(c.good===2&&c.maxCompatibleGood===2&&c.finished===0);});
check('formula overall Use counts while its imported decision is active',()=>{const unsettled={QCAcquisitionCoverage:'Not assessed',QCSettings:'Not assessed',QCFID:'Not assessed',QCSNR:'Not assessed',QCSample:'Not assessed',QCAnalysis:'Not assessed'};const c=cell(build([row(1,{...unsettled,Assessment:'Use',ImportedAssessment:'Use',EffectiveAcquisitionStatus:'Unknown'}),row(2,{...unsettled,Assessment:'',ImportedAssessment:''}),row(3,{...unsettled,Assessment:'Limited use',ImportedAssessment:'Limited use'}),row(4,{...unsettled,Assessment:'Do not use',ImportedAssessment:'Do not use'})]));assert(c.good===1&&c.pending===1&&c.unassessed===1&&c.assessments['Limited use']===1&&c.assessments['Do not use']===1);});
check('root acquisition needs an explicit empty observed folder and exact root key',()=>{const root=row(1,{Folder:'',ObservedFolder:'',AcquisitionKey:'source:1'});assert(cell(build([root])).good===1&&cell(build([root])).review===0);for(const r of [{...root,ObservedFolder:undefined},{...root,ObservedFolder:null},{...root,ObservedFolder:' '},{...root,AcquisitionKey:'source:/1'},{...root,AcquisitionKey:''},{...root,Expno:'one',ObservedExpno:'one',AcquisitionKey:'source:one'}])assert(cell(build([r])).good===0&&cell(build([r])).review===1);});
check('physical tube totals include trusted stopped and unknown acquisitions',()=>{const c=cell(build([row(1,{EffectiveAcquisitionStatus:'Stopped',Sample:11}),row(2,{EffectiveAcquisitionStatus:'Unknown',Sample:12}),row(3,{EffectiveAcquisitionStatus:'Unknown',Sample:0,Assessment:''}),row(4,{EffectiveAcquisitionStatus:'Stopped',Sample:13,ReconciliationReason:'copy'})]));assert(c.finished===0&&c.tubes===2&&c.unknownTubes===1&&c.good===2);});
globalThis.coverageTestResults={passed:checks.length,checks};
'''


def core_check_script():
    return CORE_JS+'\n'+CHECKS_JS


def browser_harness(baseline):
    """A local fake SDK; captures selected rows and metadata writes only."""
    snapshot=json.dumps(baseline,ensure_ascii=False).replace('<','\\u003c')
    fake=r'''<script>const SNAPSHOT=__DATA__;let optionHandler;window.selectionCalls=[];window.metadataWrites=[];
const SETTINGS={id:[88],tableRef:[42],parentKey:['custom'],title:['Coverage matrix'],options:[JSON.stringify({html:'<html data-nmr-coverage="1">fake</html>',customView:JSON.stringify({widgetOptions:{},access:'full'})})]};
const toTable=rows=>{const fields=[...new Set(rows.flatMap(r=>Object.keys(r.fields)))],result={id:rows.map(r=>r.id)};for(const k of fields)result[k]=rows.map(r=>r.fields[k]);return result;};
window.grist={ready:settings=>window.readySettings=settings,onOptions:cb=>{optionHandler=cb;cb({});},onRecords:cb=>cb([]),setSelectedRows:async ids=>{window.selectionCalls.push(ids);},docApi:{fetchTable:async name=>name==='Experiments'||name==='Samples'?toTable(SNAPSHOT[name]):name==='_grist_Views_section'?SETTINGS:{id:[42],tableId:['Experiments']},applyUserActions:async actions=>{window.metadataWrites.push(actions);SETTINGS.options[0]=actions[0][3].options;const prefs=JSON.parse(JSON.parse(SETTINGS.options[0]).customView).widgetOptions;optionHandler(prefs);return {};}}};</script>'''.replace('__DATA__',snapshot)
    return coverage_html(88).replace('<script src="https://docs.getgrist.com/grist-plugin-api.js"></script>',fake)


class PackagingTests(unittest.TestCase):
    def test_official_sdk_and_bound_configuration_only(self):
        html=coverage_html(88)
        self.assertIn('https://docs.getgrist.com/grist-plugin-api.js',html)
        self.assertIn('const SECTION_ID=88;',html)
        self.assertIn('allowSelectBy:true',html)
        self.assertIn('Temperature group<select id="temperature">',html)
        self.assertNotIn('grist.setOption(',html)
        self.assertNotIn('innerHTML',html)
        self.assertEqual(html.count('applyUserActions('),1)
        self.assertIn("['UpdateRecord','_grist_Views_section',sectionId",html)
        for invalid in [True,-1,0,'88']:
            with self.assertRaises(ValueError):coverage_html(invalid)

    def test_all_axis_and_options_changes_clear_native_selection(self):
        html=coverage_html(88)
        self.assertIn('await clearSelection();draft={...draft,[key]:values}',html)
        self.assertIn('await clearSelection();draft={...draft,[key]:[...new Set',html)
        self.assertIn('requestedProtein=protein,requestedAxes=Coverage.axes(draft)',html)
        self.assertIn('Coverage.mergeAxes(sections,tables,SECTION_ID,requestedProtein,requestedAxes,requestedBase)',html)
        self.assertNotIn("flags.push(cell.goodCounts.length",html)
        self.assertIn("analysis not confirmed",html)
        self.assertIn("[button.title,...flags,...needs].join(', ')",html)
        self.assertNotIn('condition / protocol groups',html)
        self.assertIn('Blue: 1 Use. Green: 2 or more Use.',html)
        self.assertIn('Counts show experiments assessed Use.',html)
        self.assertIn("const qcUse = r => r.Assessment==='Use'",html)
        self.assertIn(".cell.good{background:#e5f3e9}",html)
        self.assertIn("cell.good>=2",html)
        self.assertIn("cell.good===1",html)
        self.assertNotIn("cell.maxCompatibleGood>=2",html)
        self.assertNotIn("alignPanelRows",html)
        self.assertNotIn("field-panel",html)

    @unittest.skipUnless(shutil.which('node'),'Node CLI not installed; run core_check_script with a JS runtime or local browser harness')
    def test_pure_aggregation_and_metadata_contract(self):
        script=core_check_script()+"\nconsole.log(JSON.stringify(coverageTestResults));"
        run=subprocess.run(['node','-e',script],check=True,capture_output=True,text=True)
        self.assertEqual(json.loads(run.stdout)['passed'],35)

UI_CHECKS_JS = r'''
setTimeout(async()=>{
 const report=document.createElement('output');report.id='coverage-checks';document.body.append(report);
 const check=(ok,why)=>{if(!ok)throw Error(why);};const tick=()=>new Promise(resolve=>setTimeout(resolve,0));
 try{
  check(protein==='HeTDR1','initial populated protein');check(model.cells.length===8,'initial expected holes');
  check(model.cells.reduce((n,c)=>n+c.finished,0)===6,'six HeTDR1 Finished');
  const first=model.cells.find(c=>c.finished);await choose(first);check(selectionCalls.at(-1).length===1,'cell drives source row');
  const checkbox=$('fields').querySelector('input[value="600"]');checkbox.checked=false;checkbox.dispatchEvent(new Event('change'));await tick();
  check(selectionCalls.at(-1).length===0&&$('selection').textContent==='','axis edit clears stale selection');
  $('addFieldValue').value='850';$('addField').click();await tick();check(draft.fields.includes(850),'additional expected field');
  const fetch=grist.docApi.fetchTable,waiting=[];grist.docApi.fetchTable=name=>name.startsWith('_grist_')?new Promise(resolve=>waiting.push({name,resolve})):fetch(name);
  const saving=save();$('protein').value='hANP32a IDR';$('protein').dispatchEvent(new Event('change'));await tick();
  for(const item of waiting)item.resolve(await fetch(item.name));grist.docApi.fetchTable=fetch;await saving;await tick();
  const persisted=JSON.parse(JSON.parse(SETTINGS.options[0]).customView).widgetOptions.coverageAxes;
  check(persisted.proteins.HeTDR1.fields[0]===850&&!persisted.proteins['hANP32a IDR'],'Save captures requested protein before await');
  check(protein==='hANP32a IDR'&&draft.fields.includes(950),'Save preserves another protein draft');
  const select=grist.setSelectedRows;grist.setSelectedRows=async ids=>{if(ids.length)throw Error('synthetic selection failure');selectionCalls.push(ids);};
  await choose(model.cells.find(c=>c.finished));check(!selected&&selectionCalls.at(-1).length===0&&$('message').className==='error','selection failure clears snapshot and reports error');grist.setSelectedRows=select;
  const before=metadataWrites.length;SETTINGS.tableRef[0]=99;await save();SETTINGS.tableRef[0]=42;
  check(metadataWrites.length===before&&$('message').className==='error','invalid configuration save displays error without writing');
  check(metadataWrites.every(batch=>batch.every(a=>a[0]==='UpdateRecord'&&a[1]==='_grist_Views_section'&&a[2]===88&&Object.keys(a[3]).join()==='options')),'metadata mutation whitelist');
  saved={version:1,proteins:{}};protein='HeTDR1';selectProtein();await clearSelection();message('');
  report.textContent='Browser checks passed: initial protein, finished counts/holes, selection, axes reset, async Save capture, safe metadata write, errors.';
 }catch(error){report.textContent='Browser checks FAILED: '+error.message;report.style.color='red';}
},400);
'''


def checked_browser_harness(baseline):
    return browser_harness(baseline).replace('</body>', '<script>'+UI_CHECKS_JS+'</script></body>')

QUALITY_CHECKS_JS = r"""
setTimeout(async()=>{
 const report=document.createElement('output');report.id='quality-matrix-checks';document.body.append(report);
 const check=(ok,why)=>{if(!ok)throw Error(why);},tick=()=>new Promise(resolve=>setTimeout(resolve,0));
 const button=(field,type,salt)=>document.querySelector(`.cell[data-field="${field}"][data-type="${type}"][data-salt="${salt}"]`);
 try{
  check([...document.querySelectorAll('#matrix thead th')].map(x=>x.textContent).join('|')==='Experiment|600 MHz|700 MHz|850 MHz|950 MHz','field columns only');
  check([...document.querySelectorAll('#matrix tbody th')].map(x=>x.textContent).join('|')==='T1|T1ρ|Heteronuclear NOE|etaXY','experiment rows only');
  check(document.querySelectorAll('.cell').length===32,'all four fields and two salt lines including empty holes');
  check(button(600,'t1',50).classList.contains('good')&&button(600,'t1',50).textContent.includes('50 mM (2×)'),'two Use green');
  check(model.cells.find(c=>c.field===600&&c.type==='t1'&&c.salt===50).finished===0,'Stopped/Unknown run status does not override explicit passing QC');
  check(button(600,'t1',500).classList.contains('good')&&button(600,'t1',500).textContent.includes('500 mM (2×)'),'formula Use fallback counts independently of run status');
  check(button(700,'t1',50).textContent.includes('50 mM (0×)')&&button(700,'t1',50).textContent.includes('+1 discarded'),'only bad zero and red');
  check(button(700,'t1',500).textContent.includes('1 unassessed')&&!button(700,'t1',500).classList.contains('good'),'CSV Analyzed unassessed not good');
  check(button(850,'t1',50).textContent.includes('0×')&&button(850,'t1',50).textContent.includes('analysis not confirmed'),'Use pending not good');
  check(button(850,'t1',500).textContent.includes('0×')&&button(850,'t1',500).textContent.includes('2 identity review'),'duplicate fingerprint not good');
  check(button(950,'t1',50).classList.contains('good')&&button(950,'t1',50).textContent.includes('2×')&&button(950,'t1',50).textContent.includes('+1 discarded'),'mixed two good and bad kept separate');
  const sample=button(950,'t1',50), sampleAria=sample.getAttribute('aria-label'), sampleNeeds=sample.querySelector('.needs');
  check(sample.children.length===3&&sampleAria.includes('+1 discarded')&&sampleAria.includes('analysis not confirmed')&&Number(getComputedStyle(sampleNeeds).fontWeight)>=700,'three lines, bold action count and accessible labels');
  check(!button(600,'eta_xy',50).querySelector('.needs'),'empty cell has no action alert');
  check(button(950,'t1',500).classList.contains('good')&&!button(950,'t1',500).classList.contains('single'),'two Use with different sample metadata green');
  check(model.cells.find(c=>c.field===950&&c.type==='t1'&&c.salt===500).maxCompatibleGood===1,'sample fixture has incompatible metadata');
  check(model.cells.find(c=>c.field===600&&c.type==='t1'&&c.salt===500).maxCompatibleGood===1,'D1 fixture has incompatible protocol metadata');
  check(button(600,'t1rho',50).classList.contains('single')&&!button(600,'t1rho',50).classList.contains('good'),'exactly one trusted Use blue');
  check(!button(850,'t1',500).classList.contains('single')&&!button(700,'t1',50).classList.contains('single'),'identity conflicts and discarded records never blue');
  // Timers also advance in headless virtual time, where animation frames may pause.
  const frames=()=>new Promise(resolve=>setTimeout(resolve,100));
  const geometry=label=>{
   const cells=[...document.querySelectorAll('.cell')],heights=cells.map(b=>b.getBoundingClientRect().height);
   check(Math.max(...heights)-Math.min(...heights)<0.6,label+' equal button heights');
   for(const type of ['t1','t1rho','noe','eta_xy'])for(const salt of [50,500]){
    const tops=[600,700,850,950].map(field=>button(field,type,salt).getBoundingClientRect().top);
    check(Math.max(...tops)-Math.min(...tops)<0.6,label+' aligned salt rows '+type+' '+salt);
   }
   check(cells.every(b=>b.scrollHeight<=b.clientHeight+1),label+' no clipped text');
  };
  await frames();geometry('initial viewport');
  const empty=button(600,'eta_xy',50),oneLine=button(600,'t1rho',50),long=button(700,'t1rho',50);
  check(empty.children.length===1&&oneLine.children.length===1&&sample.children.length===3,'empty, one-line and three-line fixtures');
  const scroll=document.querySelector('.matrix-scroll'),originalWidth=scroll.style.width;
  for(const width of [1440,1024,720,1366]){
   scroll.style.width=width+'px';await frames();geometry('resize '+width);
   if(width===720){const needs=long.querySelector('.needs');check(needs.getBoundingClientRect().height>parseFloat(getComputedStyle(needs).lineHeight)*1.5,'long action text wraps across multiple lines');}
  }
  scroll.style.width=originalWidth;await frames();geometry('restored viewport');
  for(const [field,type,salt,ids]of [[600,'t1',50,[101,102]],[700,'t1',50,[104]],[850,'t1',500,[107,108]],[950,'t1',50,[109,110,111]],[600,'eta_xy',50,[]]]){
   button(field,type,salt).click();await tick();check(JSON.stringify(selectionCalls.at(-1))===JSON.stringify(ids),'exact source selection '+field+' '+type+' '+salt);}
  check(metadataWrites.length===0,'rendering and selection do not mutate data or config');await clearSelection();
  report.textContent='Quality-matrix checks passed: field columns, salt lines, 2 Use threshold, excluded/bad/pending/review, sample/protocol-independent colors, equal heights and aligned salt rows at resized widths, source IDs and empty selection.';
 }catch(error){report.textContent='Quality-matrix checks FAILED: '+error.message;report.style.color='red';}
},500);
"""


def quality_matrix_harness():
    """Disconnected fake records exercise quality states; nothing is imported."""
    data={'Experiments':[], 'Samples':[]}
    def add(ident,field,salt,assessment='Use',**extra):
        fields=dict(Study='Electrostatics',SnapshotProtein='Fixture protein',SnapshotNaCl_mM=salt,
                    EffectiveField_MHz=field,EffectiveType='T1',EffectiveTemperature_C=25,EffectiveAcquisitionStatus='Finished',
                    SourceID='fixture',Folder='quality',Expno=str(ident),ObservedSourceID='fixture',ObservedFolder='quality',
                    ObservedExpno=str(ident),AcquisitionKey='fixture:quality/'+str(ident),ObservedFingerprint='fixture-'+str(ident),
                    Sample=1,SnapshotSampleID='Fixture tube',SnapshotConstruct='C',SnapshotIsotope='15N',SnapshotBuffer='NaP',
                    SnapshotpH=6.5,SnapshotDNA_uM=0,SnapshotTFE_pct=0,SnapshotProtein_uM=100,SnapshotOtherConditions='',
                    EffectiveInstrument='Fixture instrument',EffectivePulseProgram='reviewed',EffectiveD1_s=2,
                    EffectiveDelayList_s='[0,0.1]',EffectiveSpinLock_Hz=None,EffectiveSpinLockOffset_Hz=None,Assessment=assessment,QCAcquisitionCoverage='Complete',QCSettings='Appropriate',QCFID='Clean',QCSNR='Adequate',QCSample='Acceptable',QCAnalysis='Reliable')
        fields.update(extra);data['Experiments'].append({'id':ident,'fields':fields})
    add(101,600,50,EffectiveAcquisitionStatus='Stopped',QCAcquisitionCoverage='Sufficient');add(102,600,50,EffectiveAcquisitionStatus='Unknown');add(103,600,500)
    add(114,600,500,'Use',EffectiveAcquisitionStatus='Unknown',EffectiveD1_s=9,ImportedAssessment='Use',QCAcquisitionCoverage='Not assessed',QCSettings='Not assessed',QCFID='Not assessed',QCSNR='Not assessed',QCSample='Not assessed',QCAnalysis='Not assessed')
    add(104,700,50,'Do not use');add(105,700,500,'',ResultCSV=['L',1],AnalysisStatus='Analyzed')
    add(106,850,50,'',EffectiveAcquisitionStatus='Planned',QCAcquisitionCoverage='Not assessed')
    add(107,850,500);add(108,850,500,ObservedFingerprint='fixture-107')
    add(109,950,50);add(110,950,50);add(111,950,50,'Do not use')
    add(112,950,500);add(113,950,500,SnapshotConstruct='Other construct')
    add(115,600,50,EffectiveType='T1rho',AnalysisStatus='Analyzed')
    add(116,700,50,'Limited use',EffectiveType='T1rho')
    add(117,700,50,'Do not use',EffectiveType='T1rho')
    add(118,700,50,'',EffectiveType='T1rho')
    add(119,700,50,EffectiveType='T1rho',ReconciliationReason='possible copy')
    prefs={'coverageAxes':{'version':1,'proteins':{'Fixture protein':{'fields':[600,700,850,950],'salts':[50,500],'types':['t1','t1rho','noe','eta_xy']}}}}
    html=browser_harness(data)
    encoded=json.dumps(prefs).replace('<','\\u003c')
    html=html.replace('widgetOptions:{},access:', 'widgetOptions:'+encoded+',access:').replace('optionHandler=cb;cb({});','optionHandler=cb;cb('+encoded+');')
    return html.replace('</body>','<script>'+QUALITY_CHECKS_JS+'</script></body>')
