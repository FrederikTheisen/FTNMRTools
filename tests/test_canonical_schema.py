import copy,re,unittest
from nmr_register.schema import columns,canonical_experiment_type,CANONICAL_ACQUISITION_FIELDS,MACHINE_FIELDS
from nmr_register.setup import canonical_experiment_migration,LEGACY_OBSERVED_FACTS,LEGACY_EFFECTIVE_FORMULAS
class Fake:
 def __init__(self):
  self.tables=[{'id':4,'fields':{'tableId':'Experiments'}}]
  self.cols=[{'id':100+i,'fields':{'parentId':4,'colId':c['id'],**{k:v for k,v in c.items() if k!='id'}}} for i,c in enumerate(columns()['Experiments']) if c['id'] not in ('Title','TitleSource')]
  for row in self.cols:
   if row['fields']['colId'] in LEGACY_EFFECTIVE_FORMULAS:row['fields']['formula']=LEGACY_EFFECTIVE_FORMULAS[row['fields']['colId']]
  for i,name in enumerate(LEGACY_OBSERVED_FACTS):self.cols.append({'id':1000+i,'fields':{'parentId':4,'colId':name,'type':'Text','formula':'','isFormula':False}})
  self.rows=[{'id':1,'fields':{'AcquisitionKey':'src:folder/1','ObservedSourceID':'src','ObservedFolder':'folder','ObservedExpno':'1','SourceID':'src','Folder':'folder','Expno':'1','Instrument':'600_1300','ObservedInstrument':'','ExperimentType':'T1','ObservedType':'r1','AcquisitionStatus':'Done','ObservedAcquisitionStatus':'Unknown','Temperature_C':7,'ObservedTemperature_C':0,'Field_MHz':600,'ObservedField_MHz':600.4,'SpinLockOffset_Hz':99,'ObservedSpinLockOffset_Hz':0,'MetadataDifference':'observed Temperature_C differs from confirmed value','ReconciliationReason':'','Assessment':'Do not use','Notes':'keep'}},
             {'id':2,'fields':{'Instrument':'950','AcquisitionStatus':'Planned','Assessment':'Use','Notes':'keep'}}];self.batches=[]
 def records(self,t):return {'_grist_Tables':self.tables,'_grist_Tables_column':self.cols,'Experiments':self.rows}[t]
 def apply(self,a):self.batches.append(a)
class Tests(unittest.TestCase):
 def test_audit_preserves_recorded_and_changes_only_coherent(self):
  c=Fake();before=copy.deepcopy(c.rows);plan=canonical_experiment_migration(c)
  self.assertEqual(c.rows,before);self.assertEqual(plan['audit'][0]['changes']['Temperature_C'],0)
  self.assertEqual(plan['audit'][0]['changes']['SpinLockOffset_Hz'],0)
  self.assertEqual(plan['audit'][0]['changes']['AcquisitionStatus'],'Unknown')
  self.assertNotIn('Instrument',plan['audit'][0]['changes']);self.assertEqual(plan['audit'][1]['changes'],{})
  self.assertFalse(any(a[0] in ('UpdateRecord','BulkUpdateRecord') and a[1]=='Experiments' and any(k in (a[-1] or {}) for k in ('Assessment','Sample','Notes')) for a in plan['actions']))
  with self.assertRaisesRegex(ValueError,'not approved'):canonical_experiment_migration(c,approved_digest='bad',dry_run=False)
  self.assertEqual(c.batches,[]);canonical_experiment_migration(c,approved_digest=plan['digest'],dry_run=False);self.assertEqual(len(c.batches),1)
 def test_changed_value_blocks_stale_digest_and_identity_conflict_blocks_removal(self):
  c=Fake();p=canonical_experiment_migration(c);c.rows[0]['fields']['ObservedTemperature_C']=22
  with self.assertRaisesRegex(ValueError,'changed'):canonical_experiment_migration(c,approved_digest=p['digest'],dry_run=False)
  c=Fake();c.rows[0]['fields']['ReconciliationReason']='copy';
  with self.assertRaisesRegex(ValueError,'reconciliation'):canonical_experiment_migration(c)
  self.assertEqual(c.batches,[])
 def test_partial_human_identity_and_duplicate_fail_closed(self):
  c=Fake();c.rows[0]['fields']['Folder']='other'
  with self.assertRaisesRegex(ValueError,'reconciliation'):canonical_experiment_migration(c)
  c=Fake();c.rows.append(copy.deepcopy(c.rows[0]));c.rows[-1]['id']=3
  with self.assertRaisesRegex(ValueError,'reconciliation'):canonical_experiment_migration(c)
 def test_separate_manual_identity_and_unresolved_warning(self):
  c=Fake();c.rows[0]['fields']['SourceID']='';c.rows[0]['fields']['Folder']='';c.rows[0]['fields']['Expno']=''
  c.rows[1]['fields'].update(SourceID='src',Folder='folder',Expno='1')
  with self.assertRaisesRegex(ValueError,'reconciliation'):canonical_experiment_migration(c)
  c=Fake();c.rows[0]['fields']['ObservedTemperature_C']=None
  self.assertNotIn('MetadataDifference',canonical_experiment_migration(c)['audit'][0]['changes'])
 def test_root_folder_identity(self):
  c=Fake();c.rows[0]['fields'].update(AcquisitionKey='src:1',ObservedFolder='',Folder='')
  self.assertTrue(canonical_experiment_migration(c)['audit'][0]['coherent'])
 def test_source_types_and_nominal_temperature(self):
  self.assertEqual([canonical_experiment_type(x) for x in ('r1','15N T1','15N_T1','r1rho','noe','hsqc','etaXY','ηXY','EtaXY','eta_xy','nico t1sat','nico_t1sat','Water T1 (T1sat)','15N HSQC','15N_HSQC','13C HSQC','13c_hsqc','unknown')],['T1','T1','T1','T1rho','HetNOE','HSQC','etaXY','etaXY','etaXY','etaXY','Water T1 (T1sat)','Water T1 (T1sat)','Water T1 (T1sat)','HSQC','HSQC','13C HSQC','13C HSQC',None])
  self.assertTrue(CANONICAL_ACQUISITION_FIELDS <= MACHINE_FIELDS['Experiments'])
  formula=next(c['formula'] for c in columns()['Experiments'] if c['id']=='EffectiveTemperature_C')
  evaluate=lambda v:eval(formula.replace('$Temperature_C','v'),{'v':v})
  self.assertEqual([evaluate(v) for v in (None,0,24.8498,24.8,25.2,24.7,5.499,5.5005,5.48,5.52,5.4799,24.5)],
                   [None,0,25,25,25,24.7,5.5,5.5,5.5,5.5,5.4799,24.5])
 def test_temperature_migration_is_guarded_and_preserves_raw_data(self):
  from nmr_register.setup import upgrade_temperature_grouping,LEGACY_TEMPERATURE_FORMULA
  class Client:
   def __init__(self):
    self.temp={'id':152,'fields':{'parentId':4,'colId':'EffectiveTemperature_C','type':'Numeric','isFormula':True,'formula':LEGACY_TEMPERATURE_FORMULA,'label':'Temperature (°C)','untieColIdFromLabel':False}}
    self.summary={'id':161,'fields':{'parentId':9,'summarySourceCol':152,'label':'Temperature (°C)','untieColIdFromLabel':False}}
    self.custom={'id':162,'fields':{'parentId':10,'summarySourceCol':152,'label':'Custom temperature'}}
    self.rows=[{'id':1,'fields':{'Temperature_C':24.8498}},{'id':2,'fields':{'Temperature_C':0}}];self.writes=[]
   def records(self,t):return {'_grist_Tables':[{'id':4,'fields':{'tableId':'Experiments'}}], '_grist_Tables_column':[self.temp,self.summary,self.custom]}[t]
   def apply(self,actions):
    self.writes.append(actions)
    for action in actions:next(r for r in (self.temp,self.summary,self.custom) if r['id']==action[2])['fields'].update(action[3])
  c=Client();old=copy.deepcopy(c.rows);preview=upgrade_temperature_grouping(c)
  self.assertEqual(len(preview),2);self.assertEqual(c.writes,[])
  self.assertEqual(preview[0][3]['label'],'Temperature group (°C)');self.assertTrue(preview[0][3]['untieColIdFromLabel'])
  self.assertEqual(upgrade_temperature_grouping(c,dry_run=False),preview);self.assertEqual(c.rows,old)
  self.assertEqual(c.summary['fields']['label'],'Temperature group (°C)');self.assertEqual(c.custom['fields']['label'],'Custom temperature')
  self.assertEqual(upgrade_temperature_grouping(c,dry_run=False),[]);self.assertEqual(len(c.writes),1)
  self.assertEqual(upgrade_temperature_grouping(c,dry_run=False,tables={'Experiments':4},metadata=[c.temp]),[]);self.assertEqual(len(c.writes),1)
  c.temp['fields']['formula']='$custom'
  with self.assertRaisesRegex(ValueError,'requires review'):upgrade_temperature_grouping(c,dry_run=False)
  self.assertEqual(len(c.writes),1)
if __name__=='__main__':unittest.main()
