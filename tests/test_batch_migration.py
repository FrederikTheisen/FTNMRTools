import copy,hashlib,json,tempfile,unittest
from pathlib import Path
from nmr_register.batch_migration import plan_batch,apply_batch,BatchUncertainError

def observation(folder='tube',exp='1',finger='a'*64):
 return {'AcquisitionKey':'src:'+folder+'/'+exp,'ObservedSourceID':'src','ObservedFolder':folder,'ObservedExpno':exp,
         'ObservedFingerprint':finger,'SourcePath':'/safe/'+folder+'/'+exp,'ObservedPulseProgram':'pulse',
         'ObservedType':'r1','ObservedAcquisitionStatus':'Finished','ObservedField_MHz':600.13,
         'ObservedTemperature_C':0,'ObservedD1_s':0,'ObservedSpinLockOffset_Hz':0,'RawPresent':True,'RawBytes':12,
         'ScanError':'','ObservedSampleID':'','ObservedSampleTagIssue':'','ObservedSampleTagEvidence':'No tag'}
def manifest_row(o):
 return {'AcquisitionKey':o['AcquisitionKey'],'SourceID':o['ObservedSourceID'],'Folder':o['ObservedFolder'],
         'Expno':o['ObservedExpno'],'ObservedFingerprint':o['ObservedFingerprint'],'SourcePath':o['SourcePath'],
         'ObservedPulseProgram':o['ObservedPulseProgram'],'ObservedType':o['ObservedType'],
         'ObservedAcquisitionStatus':o['ObservedAcquisitionStatus'],'MigrationKey':'workbook:A:3',
         'SourceEvidence':'A!3','legacy_source_context':{'Notes':'note','MigrationIssue':'Sample link unresolved'},'exclusions':[]}
class Client:
 def __init__(self,rows=None,mode='ok'):self.rows=copy.deepcopy(rows or []);self.writes=[];self.mode=mode
 def records(self,t):assert t=='Experiments';return copy.deepcopy(self.rows)
 def add_records(self,t,fields):
  self.writes.append(copy.deepcopy(fields));n=len(fields) if self.mode!='partial' else 1
  self.rows.extend({'id':len(self.rows)+i+1,'fields':copy.deepcopy(f)} for i,f in enumerate(fields[:n]))
  if self.mode in ('timeout','partial'):raise TimeoutError('ambiguous POST')
class Tests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.path=Path(self.tmp.name)/'manifest.json';self.source=Path(self.tmp.name)/'source.json';self.input=Path(self.tmp.name)/'workbook';self.input.write_bytes(b'original');self.obs=[observation()];self.rows=[manifest_row(self.obs[0])];self.write()
 def tearDown(self):self.tmp.cleanup()
 def write(self):
  entries=[]
  for r in self.rows:
   if r.get('origin')=='inventory':continue
   fields={k:r[k] for k in ('SourceID','Folder','Expno','MigrationKey','SourceEvidence')};fields.update(r['legacy_source_context']);fields['Instrument']='600Cryo';entries.append({'fields':fields})
  source=json.dumps({'source_files':[{'path':str(self.input),'sha256':hashlib.sha256(self.input.read_bytes()).hexdigest()}],'planned_targets':entries}).encode();self.source.write_bytes(source)
  b=json.dumps({'records':self.rows,'source_sha256':{str(self.source):hashlib.sha256(source).hexdigest()}}).encode();self.path.write_bytes(b);self.digest=hashlib.sha256(b).hexdigest()
 def plan(self,c=None,observations=None):return plan_batch(self.path,self.obs if observations is None else observations,c or Client(),expected_manifest_sha256=self.digest,source_reconciliation_path=self.source)
 def test_fields_and_zero_preservation(self):
  p=self.plan();f=p['insert'][0]['fields'];self.assertEqual((f['Temperature_C'],f['D1_s'],f['SpinLockOffset_Hz']),(0,0,0));self.assertEqual(f['Field_MHz'],600)
  self.assertEqual(f['ExperimentType'],'T1');self.assertEqual(f['Instrument'],'600Cryo');self.assertNotIn('RegisteredAt',f);self.assertEqual(f['SnapshotNaCl_mM'],None);self.assertEqual(f['Notes'],'note');self.assertEqual(f['MigrationIssue'],'Sample link unresolved')
  for name in ('Assessment','AssessmentReason','Sample','Series','Study','QCSNR','ResultCSV'):self.assertNotIn(name,f)
 def test_hash_tampering_and_protected_path_fail(self):
  with self.assertRaisesRegex(ValueError,'hash'):plan_batch(self.path,self.obs,Client(),expected_manifest_sha256='bad',source_reconciliation_path=self.source)
  o=observation('prefix211@950suffix','1');self.obs=[o];self.rows=[manifest_row(o)];self.write()
  with self.assertRaisesRegex(ValueError,'identity'):self.plan()
 def test_source_hash_change(self):
  self.input.write_bytes(b'changed')
  with self.assertRaisesRegex(ValueError,'hash changed'):self.plan()
 def test_unknown_acquisition_status_does_not_fall_back(self):
  self.obs[0]['ObservedAcquisitionStatus']='Unknown';self.rows[0]['ObservedAcquisitionStatus']='Unknown';self.write()
  f=self.plan()['insert'][0]['fields'];self.assertEqual(f['AcquisitionStatus'],'Unknown');self.assertEqual(f['Instrument'],'600Cryo')
 def test_protected_original_source_rejected_before_access(self):
  data=json.loads(self.source.read_text());data['source_files'][0]['path']='/tmp/prefix211@950suffix/never-read'
  blob=json.dumps(data).encode();self.source.write_bytes(blob)
  manifest=json.loads(self.path.read_text());manifest['source_sha256'][str(self.source)]=hashlib.sha256(blob).hexdigest()
  raw=json.dumps(manifest).encode();self.path.write_bytes(raw);self.digest=hashlib.sha256(raw).hexdigest()
  with self.assertRaisesRegex(ValueError,'protected'):self.plan()
 def test_invalid_source_and_folder_syntax(self):
  for source,folder in [('bad:id','tube'),('bad id','tube'),('src','tube//child'),('src','tube\\child')]:
   self.rows=[manifest_row(self.obs[0])];self.rows[0]['SourceID']=source;self.rows[0]['Folder']=folder;self.write()
   with self.subTest(source=source,folder=folder),self.assertRaises(ValueError):self.plan()
 def test_duplicate_fingerprint_outside_batch_fails(self):
  with self.assertRaisesRegex(ValueError,'fingerprint'):self.plan(observations=self.obs+[observation('else','2')])
 def test_known_split_stale_inventory_cannot_import(self):
  from nmr_register.split import SPLIT_EVIDENCE
  obs=copy.deepcopy(self.obs);obs[0]['ObservedAcquisitionEvidence']=SPLIT_EVIDENCE
  with self.assertRaisesRegex(ValueError,'Fresh observation missing'):self.plan(observations=obs)
 def test_changed_or_invalid_observation_fails(self):
  for key,value in [('ObservedPulseProgram','new'),('ScanError','split dataset'),('ObservedSampleID','sample'),('ObservedSampleTagIssue','unsafe'),('ResultCSV',['L',1])]:
   obs=copy.deepcopy(self.obs);obs[0][key]=value
   with self.subTest(key=key),self.assertRaises(ValueError):self.plan(observations=obs)
 def test_existing_human_data_is_never_overwritten(self):
  client=Client();apply_batch(self.path,self.obs,client,expected_manifest_sha256=self.digest,source_reconciliation_path=self.source,dry_run=False)
  client.rows[0]['fields']['Notes']='human edit';before=copy.deepcopy(client.rows)
  again=apply_batch(self.path,self.obs,client,expected_manifest_sha256=self.digest,source_reconciliation_path=self.source,dry_run=False)
  self.assertEqual(client.rows,before);self.assertEqual(len(client.writes),1);self.assertEqual(again['skipped'],[self.rows[0]['AcquisitionKey']])
 def test_partial_manual_identity_and_colliding_migration_fail(self):
  for row in ({'Folder':'tube','SourceID':'src'}, {'MigrationKey':'workbook:A:3','SourceID':'else','Folder':'x','Expno':'2'}):
   c=Client([{'id':8,'fields':row}])
   with self.assertRaises(ValueError):self.plan(c)
   self.assertEqual(c.writes,[])
 def test_timeout_after_full_commit_reconciles_without_retry(self):
  client=Client(mode='timeout');result=apply_batch(self.path,self.obs,client,expected_manifest_sha256=self.digest,source_reconciliation_path=self.source,dry_run=False)
  self.assertEqual(result['inserted'],[self.rows[0]['AcquisitionKey']]);self.assertEqual(len(client.writes),1)
 def test_nonfinite_historical_numeric_rejected(self):
  source=json.loads(self.source.read_text());source['planned_targets'][0]['fields']['Field_MHz']=float('nan')
  b=json.dumps(source).encode();self.source.write_bytes(b);m=json.loads(self.path.read_text());m['source_sha256'][str(self.source)]=hashlib.sha256(b).hexdigest();b=json.dumps(m).encode();self.path.write_bytes(b);self.digest=hashlib.sha256(b).hexdigest()
  self.obs[0]['ObservedField_MHz']=None
  with self.assertRaisesRegex(ValueError,'Non-finite'):self.plan()
 def test_inventory_only_retry_and_human_preservation(self):
  r=self.rows[0];r.update(origin='inventory',MigrationKey='discovery:'+r['AcquisitionKey'],SourceEvidence='Observed acquisition: '+r['AcquisitionKey'],legacy_source_context={});self.write()
  c=Client(mode='timeout');result=apply_batch(self.path,self.obs,c,expected_manifest_sha256=self.digest,source_reconciliation_path=self.source,dry_run=False)
  self.assertEqual(len(result['inserted']),1);f=c.rows[0]['fields'];self.assertEqual(f['Notes'],'');self.assertNotIn('Sample',f);self.assertNotIn('Study',f);self.assertNotIn('Assessment',f)
  c.rows[0]['fields']['Notes']='human edit';before=copy.deepcopy(c.rows)
  apply_batch(self.path,self.obs,c,expected_manifest_sha256=self.digest,source_reconciliation_path=self.source,dry_run=False)
  self.assertEqual(c.rows,before);self.assertEqual(len(c.writes),1)
 def test_inventory_origin_cannot_hide_legacy_or_fake_evidence(self):
  r=self.rows[0];r.update(origin='inventory',MigrationKey='discovery:'+r['AcquisitionKey'],SourceEvidence='Observed acquisition: '+r['AcquisitionKey'],legacy_source_context={});self.write()
  for key,value in [('SourceEvidence','fake'),('MigrationKey','other'),('legacy_source_context',{'Notes':'fake'})]:
   old=r[key];r[key]=value;self.write()
   with self.subTest(key=key),self.assertRaises(ValueError):self.plan()
   r[key]=old
  self.write();source=json.loads(self.source.read_text());source['planned_targets'].append({'fields':{'SourceID':'src','Folder':'tube','Expno':'1','MigrationKey':'workbook:A:3','Assessment':'Do not use'}})
  b=json.dumps(source).encode();self.source.write_bytes(b);m=json.loads(self.path.read_text());m['source_sha256'][str(self.source)]=hashlib.sha256(b).hexdigest();b=json.dumps(m).encode();self.path.write_bytes(b);self.digest=hashlib.sha256(b).hexdigest()
  with self.assertRaisesRegex(ValueError,'legacy evidence'):self.plan()
 def test_partial_commit_stops_and_resumes_only_remaining(self):
  second=observation('tube','2','b'*64);self.obs.append(second);self.rows.append({**manifest_row(second),'MigrationKey':'workbook:A:4'});self.write()
  client=Client(mode='partial')
  with self.assertRaises(BatchUncertainError) as caught:apply_batch(self.path,self.obs,client,expected_manifest_sha256=self.digest,source_reconciliation_path=self.source,dry_run=False)
  self.assertEqual((len(caught.exception.committed),len(caught.exception.remaining)),(1,1));self.assertEqual(len(client.writes),1)
  client.mode='ok';res=apply_batch(self.path,self.obs,client,expected_manifest_sha256=self.digest,source_reconciliation_path=self.source,dry_run=False)
  self.assertEqual(len(res['inserted']),1);self.assertEqual(len(client.rows),2)
if __name__=='__main__':unittest.main()
