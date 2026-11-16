"""Assessment actions are restricted to explicit ratings and one optional note."""
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import unittest

SOURCE = Path(__file__).resolve().parents[1]/'nmr_register'/'assessment_widget.py'
SPEC = importlib.util.spec_from_file_location('assessment_widget_under_test', SOURCE)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
CORE_JS, assessment_html = MODULE.CORE_JS, MODULE.assessment_html

CHECKS_JS = r'''
const checks=[];
const assert=(ok,message='assertion failed')=>{if(!ok)throw Error(message);};
const check=(name,fn)=>{fn();checks.push(name);};
const throws=fn=>{let caught=false;try{fn();}catch(e){caught=true;}assert(caught);};
const record={id:1,QCAcquisitionCoverage:'Not assessed',QCSettings:'Not assessed',QCFID:'Not assessed',QCSNR:'Not assessed',QCSample:'Not assessed',QCAnalysis:'Not assessed',AssessmentReason:''};
check('six exact explicit criteria default unchecked',()=>{assert(AssessmentCard.criteria.length===6);assert(AssessmentCard.keys.join(',')==='QCAcquisitionCoverage,QCSettings,QCFID,QCSNR,QCSample,QCAnalysis');assert(Object.values(AssessmentCard.snapshot({})).filter(x=>x==='Not assessed').length===6);});
check('only six criteria and one note are writable',()=>{const action=AssessmentCard.action(1,record,{...record,QCSample:'Unacceptable',Assessment:'Use',Sample:55,AssessmentReason:'optional note'},record);assert(action[0]==='UpdateRecord'&&action[1]==='Experiments'&&action[2]===1);assert(Object.keys(action[3]).length===7&&!('Assessment'in action[3])&&!('Sample'in action[3]));assert(action[3].QCSample==='Unacceptable'&&action[3].QCSettings==='Not assessed');});
check('rating changes clear imported assessment while note-only edits preserve it',()=>{const imported={...record,ImportedAssessment:'Do not use'};const noteDraft={...AssessmentCard.snapshot(imported),AssessmentReason:'context only'};const noteAction=AssessmentCard.action(1,imported,noteDraft,imported);assert(!Object.hasOwn(noteAction[3],'ImportedAssessment'));const ratingDraft={...AssessmentCard.snapshot(imported),QCSNR:'Adequate'};const ratingAction=AssessmentCard.action(1,imported,ratingDraft,imported);assert(ratingAction[3].ImportedAssessment===''&&Object.keys(ratingAction[3]).length===8);throws(()=>AssessmentCard.action(1,imported,ratingDraft,{...imported,ImportedAssessment:'Use'}));});
check('partial ratings are allowed without auto passing other checks',()=>{const action=AssessmentCard.action(1,record,{...record,QCAcquisitionCoverage:'Sufficient'},record);assert(action[3].QCAnalysis==='Not assessed'&&action[3].AssessmentReason==='');});
check('fresh criteria and note conflicts fail closed',()=>{for(const fresh of [{...record,QCSNR:'Adequate'},{...record,AssessmentReason:'another editor'},null,{...record,id:2},{...record,ObservedFingerprint:'replacement'},{...record,ObservedFolder:'different source'},{...record,Sample:22},{...record,Title:'new title comment'},{...record,Instrument:'changed instrument'},{...record,AcquisitionStatus:'Finished'}])throws(()=>AssessmentCard.action(1,record,record,fresh));});
check('title and canonical acquisition facts are freshness guards',()=>{assert(!AssessmentCard.sameIdentity(record,{...record,Title:'changed'}));assert(!AssessmentCard.sameIdentity(record,{...record,PulseProgram:'changed'}));});
check('invalid choices are never saved',()=>throws(()=>AssessmentCard.action(1,record,{...record,QCSample:'Probably okay'},record)));
check('criterion tones never infer a rating',()=>{assert(AssessmentCard.tone('QCAcquisitionCoverage','Sufficient')==='pass');assert(AssessmentCard.tone('QCSettings','Not assessed')==='unchecked');assert(AssessmentCard.tone('QCFID','Usable with exclusions')==='limited');assert(AssessmentCard.tone('QCSample','Unacceptable')==='fail');});
check('relaxation aliases expose scoped analysis choices',()=>{for(const type of ['T1','T1rho','T1ρ','15N T1','r1','r1rho','15N R1ρ']){assert(AssessmentCard.relaxationType(type),type);assert(AssessmentCard.choicesFor('QCAnalysis',{ExperimentType:type}).length===6,type);}for(const type of ['T2','HSQC','',null])assert(AssessmentCard.choicesFor('QCAnalysis',{ExperimentType:type}).join('|')==='Not assessed|Reliable|Limited|Unreliable',String(type));});
check('analysis tones and type validation preserve saved values',()=>{const excluded='Reliable — early points excluded',two='Two-phase — unresolved';assert(AssessmentCard.tone('QCAnalysis',excluded)==='pass');assert(AssessmentCard.tone('QCAnalysis',two)==='fail');assert(AssessmentCard.tone('QCAnalysis','Limited')==='limited');for(const value of [excluded,two]){const t1={...record,ExperimentType:'T1'};assert(AssessmentCard.action(1,t1,{...t1,QCAnalysis:value},t1)[3].QCAnalysis===value);const hsqc={...record,ExperimentType:'HSQC'};throws(()=>AssessmentCard.action(1,hsqc,{...hsqc,QCAnalysis:value},hsqc));const saved={...hsqc,QCAnalysis:value};assert(AssessmentCard.action(1,saved,{...saved,AssessmentReason:'note'},saved)[3].QCAnalysis===value);}});
check('decision tone preserves incomplete amber and canonical veto red',()=>{assert(AssessmentCard.decisionTone({Assessment:'',AssessmentDisplay:'Not assessed'})==='unchecked');assert(AssessmentCard.decisionTone({Assessment:'',AssessmentDisplay:'Use so far · analysis missing'})==='limited');assert(AssessmentCard.decisionTone({Assessment:'Do not use',AssessmentDisplay:'Do not use (partial)'})==='fail');assert(AssessmentCard.decisionTone({Assessment:'Use',AssessmentDisplay:'Use'})==='pass');});
check('only the bound Experiments creator is eligible to write',()=>{const section={id:88,tableRef:42,parentKey:'custom',options:JSON.stringify({html:'<html data-nmr-assessment="1">'})};AssessmentCard.context([section],[{id:42,tableId:'Experiments'}],88);throws(()=>AssessmentCard.context([{...section,tableRef:9}],[{id:42,tableId:'Experiments'}],88));throws(()=>AssessmentCard.context([{...section,options:'{}'}],[{id:42,tableId:'Experiments'}],88));});

const linked={...record,Sample:30,SnapshotSampleID:'tube30',ObservedSampleID:'tube30',AcquisitionKey:'main:folder/1',SourceID:'main',Folder:'folder',Expno:'1',ObservedSourceID:'main',ObservedFolder:'folder',ObservedExpno:'1',RegisteredAt:null};
const samples=[{id:30,SampleID:'tube30'}],rating={...AssessmentCard.snapshot(linked),QCFID:'Unusable'},stamp=1794475200;
check('explicit rating save confirms verified discovered acquisition atomically',()=>{const action=AssessmentCard.action(1,linked,rating,linked,[linked],samples,stamp);assert(action[3].RegisteredAt===stamp&&action[3].QCFID==='Unusable'&&action[3].QCAnalysis==='Not assessed');assert(!('Assessment'in action[3])&&!('Sample'in action[3]));assert(linked.RegisteredAt===null);});
check('note only reset and normalized default differences do not register',()=>{for(const draft of [{...AssessmentCard.snapshot(linked),AssessmentReason:'note'},AssessmentCard.snapshot(linked),{...AssessmentCard.snapshot(linked),QCFID:'Not assessed'}])assert(AssessmentCard.registrationTimestamp(linked,draft,linked,[linked],samples,stamp)===undefined);const absent={...linked,QCFID:null};assert(AssessmentCard.registrationTimestamp(absent,AssessmentCard.snapshot(absent),absent,[absent],samples,stamp)===undefined);const rated={...linked,QCFID:'Clean'};assert(AssessmentCard.registrationTimestamp(rated,{...AssessmentCard.snapshot(rated),QCFID:'Not assessed'},rated,[rated],samples,stamp)===undefined);});
check('existing registration timestamps and unqualified ratings are preserved',()=>{for(const changes of [{RegisteredAt:1700000000},{Sample:0},{ObservedSampleID:'another tube'},{SampleTagIssue:'mismatch'},{ObservedSampleTagIssue:'bad marker'},{ReconciliationReason:'copy requires review'},{SourceID:'other'},{ObservedExpno:'2'},{AcquisitionKey:'wrong'},{ObservedSourceID:''}]){const fresh={...linked,...changes};assert(!Object.hasOwn(AssessmentCard.action(1,fresh,rating,fresh,[fresh],samples,stamp)[3],'RegisteredAt'));}assert(!Object.hasOwn(AssessmentCard.action(1,linked,rating,linked,[linked],[],stamp)[3],'RegisteredAt'));});
check('duplicate exact identities and conflicting samples block registration only',()=>{for(const other of [{...linked,id:2},{id:2,SourceID:'main',Folder:'folder',Expno:'1'},{id:2,ObservedSourceID:'main',ObservedFolder:'folder',ObservedExpno:'1'}])assert(AssessmentCard.registrationTimestamp(linked,rating,linked,[linked,other],samples,stamp)===undefined);assert(AssessmentCard.registrationTimestamp(linked,rating,linked,[linked],[...samples,...samples],stamp)===undefined);});
check('timestamp and tag changes while editing fail freshness guard',()=>{for(const changes of [{RegisteredAt:stamp},{SampleTagIssue:'new mismatch'},{ObservedSampleTagIssue:'malformed'},{ObservedSampleID:'another tube'}])throws(()=>AssessmentCard.action(1,linked,rating,{...linked,...changes},[linked],samples,stamp));});
check('root folder identity and zero missing timestamp remain supported',()=>{const root={...linked,RegisteredAt:0,Folder:'',ObservedFolder:'',AcquisitionKey:'main:1'};assert(AssessmentCard.registrationTimestamp(root,rating,root,[root],samples,stamp)===stamp);});
globalThis.assessmentTestResults={passed:checks.length,checks};
'''


def core_check_script():
    return CORE_JS+'\n'+CHECKS_JS


BROWSER_CHECKS = r'''
setTimeout(async()=>{
 const report=document.createElement('output');report.id='assessment-checks';document.body.append(report);
 const assert=(ok,why)=>{if(!ok)throw Error(why);};const tick=()=>new Promise(r=>setTimeout(r,0));
 try{
  assert(document.querySelectorAll('select').length===6&&document.querySelectorAll('textarea').length===1,'six explicit ratings and one note');
  assert(current.SnapshotSampleID===RECORDS[0].SnapshotSampleID,'partial SDK record is loaded fully by captured ID');
  assert($('titleText').textContent===RECORDS[0].Title&&!$('titleBlock').hidden&&!$('titleText').querySelector('img'),'acquisition title and comments are shown verbatim as text');
  assert($('source').textContent.includes('quality / 1')&&$('source').textContent.includes('Bruker 600')&&$('source').textContent.includes('Stopped')&&$('conditions').textContent.includes('t1rho')&&$('conditions').textContent.includes('2023-11-14'),'assessment header uses canonical acquisition facts');
  assert(AssessmentCard.keys.every(k=>$(k).value==='Not assessed'),'Finished raw metadata does not auto pass');
  assert($('decision').textContent==='Do not use'&&!document.querySelector('input[name="Assessment"]'),'imported decision is read only and participates as the current fallback');
  assert($('explanation').textContent==='Imported assessment: Do not use; checklist not yet completed','formula-provided imported-decision origin is visible in the dialog');
  $('QCSample').value='Unacceptable';$('QCSample').dispatchEvent(new Event('change'));await save();
  assert(Number.isFinite(current.RegisteredAt)&&current.RegisteredAt>0,'first verified rating save confirms registration at actual save time');
  assert($('decision').textContent==='Do not use'&&current.QCSettings==='Not assessed'&&current.ImportedAssessment===''&&$('explanation').textContent==='Sample during acquisition: Unacceptable','saved checklist rating supersedes and clears the imported decision');
  const positives=['Sufficient','Appropriate','Clean','Adequate','Acceptable','Reliable'];
  AssessmentCard.keys.forEach((k,i)=>{$(k).value=positives[i];$(k).dispatchEvent(new Event('change'));});await save();
  assert($('decision').textContent==='Use','stored computed decision after explicit six positive ratings');
  $('QCAnalysis').value='Reliable — early points excluded';$('QCAnalysis').dispatchEvent(new Event('change'));await save();
  assert(current.QCAnalysis==='Reliable — early points excluded'&&$('decision').textContent==='Use','rescued reliable fit saves its distinct rating and remains Use');
  $('QCAnalysis').value='Two-phase — unresolved';$('QCAnalysis').dispatchEvent(new Event('change'));await save();
  assert(current.QCAnalysis==='Two-phase — unresolved'&&$('decision').textContent==='Do not use','unresolved two-phase fit is a veto');
  $('QCAnalysis').value='Reliable';$('QCAnalysis').dispatchEvent(new Event('change'));await save();
  $('note').value='draft retained';$('note').dispatchEvent(new Event('input'));await receive(RECORDS[1]);
  assert(current.id===1&&pending.id===2&&draft.AssessmentReason==='draft retained'&&!$('switch').hidden,'selection switch retains draft and makes choice explicit');
  $('discard').click();assert(current.id===2&&!dirty(),'explicit discard switches');
  assert($('conditions').textContent.includes('2026-10'),'month-precision manual acquisition date remains visible');
  await receive(RECORDS[0]);$('note').value='captured record one';$('note').dispatchEvent(new Event('input'));
  const fetch=grist.docApi.fetchTable;let release;grist.docApi.fetchTable=name=>name==='Experiments'&&!release?new Promise(resolve=>{release=()=>resolve(table(RECORDS));}):fetch(name);
  const saving=save();await tick();await receive(RECORDS[1]);assert($('note').disabled,'controls disabled while saving');release();grist.docApi.fetchTable=fetch;await saving;
  assert(WRITES.at(-1)[0][2]===1&&RECORDS[0].AssessmentReason==='captured record one'&&current.id===2,'save captures source row across asynchronous selection switch');
  $('note').value='local conflicting draft';$('note').dispatchEvent(new Event('input'));const before=WRITES.length;RECORDS[1].QCSNR='Marginal';await save();
  assert(WRITES.length===before&&draft.AssessmentReason==='local conflicting draft'&&$('message').className==='error','concurrent change blocks write and retains draft');
  assert(WRITES.every(batch=>batch.length===1&&batch[0][0]==='UpdateRecord'&&batch[0][1]==='Experiments'&&Object.keys(batch[0][3]).sort().join()===(Object.hasOwn(batch[0][3],'ImportedAssessment')?['AssessmentReason','ImportedAssessment',...AssessmentCard.keys].sort():['AssessmentReason',...AssessmentCard.keys].sort()).concat(Object.hasOwn(batch[0][3],'RegisteredAt')?['RegisteredAt']:[]).sort().join()),'exact rating/note mutation whitelist with conditional imported clear');
  assert(WRITES[0][0][3].ImportedAssessment===''&&!Object.hasOwn(WRITES.at(-1)[0][3],'ImportedAssessment'),'rating edit clears imported fallback; later note-only edit does not touch it');
  $('discard').click();
  const fetchLatest=grist.docApi.fetchTable,waiting=[];grist.docApi.fetchTable=name=>name==='Experiments'?new Promise(resolve=>waiting.push(resolve)):fetchLatest(name);
  const older=receive({id:1}),newer=receive({id:2});waiting[1](table(RECORDS));await newer;waiting[0](table(RECORDS));await older;
  assert(current.id===2,'late selection response cannot restore a stale source row');grist.docApi.fetchTable=fetchLatest;
  RECORDS[1].ExperimentType='HSQC';RECORDS[1].QCAnalysis='Reliable — early points excluded';await receive(RECORDS[1]);
  assert($('QCAnalysis').value===RECORDS[1].QCAnalysis&&$('QCAnalysis').selectedOptions[0].textContent.includes('Saved rating'),'restricted saved value stays visible');
  assert([...$('QCAnalysis').options].length===5&&!([...$('QCAnalysis').options].some(o=>o.value==='Two-phase — unresolved')),'restricted choices do not leak');
  $('note').value='retained rating';$('note').dispatchEvent(new Event('input'));await save();assert(RECORDS[1].QCAnalysis==='Reliable — early points excluded','note save preserves human rating');
  await receive(RECORDS[0]);assert([...$('QCAnalysis').options].length===6,'relaxation selection restores choices');message('');
  report.textContent='Assessment checks passed: explicit six ratings, read-only decision, failed-with-unchecked, retained draft, captured row ID, fresh-data conflict, exact write whitelist.';
 }catch(error){report.textContent='Assessment checks FAILED: '+error.message;report.style.color='red';}
},200);
'''


def browser_harness():
    fake = r'''<script>
const EMPTY={QCAcquisitionCoverage:'Not assessed',QCSettings:'Not assessed',QCFID:'Not assessed',QCSNR:'Not assessed',QCSample:'Not assessed',QCAnalysis:'Not assessed',AssessmentReason:'',Assessment:'',AssessmentExplanation:'Unchecked: acquisition coverage, settings, FID, SNR, sample, analysis.'};
const RECORDS=[{...EMPTY,id:1,Sample:1,AcquisitionKey:'fixture:quality/1',ObservedSourceID:'fixture',ObservedFolder:'quality',ObservedExpno:'1',ImportedAssessment:'Do not use',Assessment:'Do not use',AssessmentExplanation:'Imported assessment: Do not use; checklist not yet completed',SnapshotSampleID:'HeTDR1 · 50 mM NaCl',SourceID:'fixture',Folder:'quality',Expno:'1',AcquisitionStatus:'Stopped',EffectiveAcquisitionStatus:'Stopped',ExperimentType:'T1ρ',Instrument:'Bruker 600',Field_MHz:600.2,AcquiredAt:1700000000,AcquisitionDateText:'',PulseProgram:'t1rho',Temperature_C:25,Title:['SAMPLE_ID=HeTDR1','Comment: <img src=x onerror=alert(1)> text'].join(String.fromCharCode(10)),TitleSource:'pdata/1/title',SnapshotNaCl_mM:50}, {...EMPTY,id:2,SnapshotSampleID:'Second sample',ExperimentType:'T1',SourceID:'fixture',Folder:'quality',Expno:'2',AcquisitionDateText:'2026-10'}];
const WRITES=[];const table=rs=>{const keys=[...new Set(rs.flatMap(r=>Object.keys(r)))];return Object.fromEntries(keys.map(k=>[k,rs.map(r=>r[k])]));};
const recalc=r=>{const keys=['QCAcquisitionCoverage','QCSettings','QCFID','QCSNR','QCSample','QCAnalysis'];const fail=['Insufficient','Invalid','Unusable','Inadequate','Unacceptable','Two-phase — unresolved','Unreliable'],limited=['Usable with restrictions','Usable with exclusions','Marginal','Concern','Limited'];r.Assessment=keys.some(k=>fail.includes(r[k]))?'Do not use':keys.every(k=>!r[k]||r[k]==='Not assessed')?(r.ImportedAssessment||''):keys.some(k=>r[k]==='Not assessed')?'':keys.filter(k=>limited.includes(r[k])).length>=2?'Limited use':'Use';r.AssessmentExplanation=keys.every(k=>!r[k]||r[k]==='Not assessed')&&r.ImportedAssessment?'Imported assessment: '+r.ImportedAssessment+'; checklist not yet completed':r.Assessment==='Use'?'All criteria acceptable':r.Assessment==='Do not use'?'Sample during acquisition: '+r.QCSample:r.Assessment?'Limited criteria remain.':'Assessment incomplete.';};
window.grist={ready:settings=>window.readySettings=settings,onRecord:cb=>cb({id:1}),docApi:{fetchTable:async name=>name==='Experiments'?table(RECORDS):name==='Samples'?{id:[1],SampleID:['HeTDR1']}:name==='_grist_Views_section'?{id:[88],tableRef:[42],parentKey:['custom'],options:[JSON.stringify({html:'<html data-nmr-assessment="1">'})]}:{id:[42],tableId:['Experiments']},applyUserActions:async actions=>{WRITES.push(actions);const action=actions[0],record=RECORDS.find(r=>r.id===action[2]);Object.assign(record,action[3]);recalc(record);return {};}}};
</script>'''
    return assessment_html(88).replace('<script src="https://docs.getgrist.com/grist-plugin-api.js"></script>', fake).replace('</body>', '<script>'+BROWSER_CHECKS+'</script></body>')


class AssessmentPackagingTests(unittest.TestCase):
    def test_title_comments_are_displayed_as_text_with_long_content_scroll(self):
        html=assessment_html(88)
        self.assertIn("$('titleText').textContent=title",html)
        self.assertIn('id="titleBlock" hidden',html)
        self.assertIn('max-height:10em;overflow:auto',html)
        self.assertIn("'TitleSource'",html)
        self.assertNotIn('innerHTML',html)

    def test_only_explicit_ratings_and_note_are_writable(self):
        html = assessment_html(88)
        self.assertIn('const SECTION_ID=88', html)
        self.assertIn("['UpdateRecord','Experiments',id,fields]", html)
        self.assertEqual(html.count('applyUserActions('), 1)
        self.assertNotIn('innerHTML', html)
        self.assertIn('grist.onRecord(receive)', html)
        self.assertIn('requestedId=current.id,requestedBase={...base},requestedDraft={...draft}', html)
        self.assertIn('pending=record||null', html)
        self.assertIn('AssessmentCard.equal(record,base)&&AssessmentCard.sameIdentity(record,base)', html)
        self.assertIn('Calculated overall assessment', html)
        for invalid in [True, 0, -1, '88']:
            with self.assertRaises(ValueError):
                assessment_html(invalid)

    @unittest.skipUnless(shutil.which('node'), 'Node CLI absent; run core_check_script in JS runtime or browser fixture')
    def test_core(self):
        result = subprocess.run(['node', '-e', core_check_script()+'\nconsole.log(JSON.stringify(assessmentTestResults))'], check=True, capture_output=True, text=True)
        self.assertEqual(json.loads(result.stdout)['passed'], 18)


if __name__ == '__main__':
    unittest.main()
