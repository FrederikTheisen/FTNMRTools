import copy,json,re,unittest
from nmr_register.schema import columns,ATTENTION_FORMULA,PRE_DUMMY_ATTENTION_FORMULA
from nmr_register.setup import upgrade_acquisition_workflow_status,LEGACY_WORKFLOW_EFFECTIVE,LEGACY_WORKFLOW_RECORDS
class Fake:
 def __init__(self):
  old=PRE_DUMMY_ATTENTION_FORMULA.replace('if $EffectiveAcquisitionStatus == "Aborted": issues.append("Acquisition aborted")\nif $EffectiveAcquisitionStatus == "Unknown": issues.append("Acquisition status unknown")','if $EffectiveAcquisitionStatus == "Stopped": issues.append("Acquisition stopped")').replace('$EffectiveAcquisitionStatus == "Completed"','$EffectiveAcquisitionStatus == "Finished"')
  self.tables=[{'id':1,'fields':{'tableId':'Experiments'}},{'id':2,'fields':{'tableId':'Experiments_summary_X'}}]
  self.cols=[{'id':3,'fields':{'parentId':1,'colId':'AcquisitionStatus','type':'Text','isFormula':False,'formula':'','widgetOptions':json.dumps({'alignment':'right','choices':['Custom']})}}, {'id':4,'fields':{'parentId':1,'colId':'EffectiveAcquisitionStatus','isFormula':True,'formula':LEGACY_WORKFLOW_EFFECTIVE}}, {'id':5,'fields':{'parentId':1,'colId':'AttentionReason','isFormula':True,'formula':old}},{'id':6,'fields':{'parentId':2,'colId':'RecordsDone','isFormula':True,'formula':LEGACY_WORKFLOW_RECORDS}}]
  self.views=[{'id':7,'fields':{'name':'Experiments'}}];self.sections=[{'id':8,'fields':{'parentId':7,'title':'Experiments','parentKey':'record'}}];self.fields=[{'id':9,'fields':{'parentId':8,'colRef':3}}]
  self.rows=[{'id':1,'fields':{'AcquisitionStatus':'Failed','Notes':'human'}},{'id':2,'fields':{'AcquisitionStatus':'Finished'}},{'id':3,'fields':{'AcquisitionStatus':'Running'}}];self.calls=[]
 def records(self,t):return {'_grist_Tables':self.tables,'_grist_Tables_column':self.cols,'_grist_Views':self.views,'_grist_Views_section':self.sections,'_grist_Views_section_field':self.fields,'Experiments':self.rows}[t]
 def apply(self,actions):
  self.calls.append(copy.deepcopy(actions))
  for a in actions:
   if a[0]=='ModifyColumn':next(x for x in self.cols if x['fields']['colId']==a[2])['fields'].update(a[3])
   elif a[1]=='_grist_Tables_column':next(x for x in self.cols if x['id']==a[2])['fields'].update(a[3])
   elif a[1]=='_grist_Views_section_field' and a[0]=='RemoveRecord':self.fields[:]=[x for x in self.fields if x['id']!=a[2]]
   elif a[1]=='_grist_Views_section_field':next(x for x in self.fields if x['id']==a[2])['fields'].update(a[3])
class Tests(unittest.TestCase):
 def test_old_aliases_unknown_custom_and_five_states(self):
  formula=next(c['formula'] for c in columns()['Experiments'] if c['id']=='EffectiveAcquisitionStatus');code=re.sub(r'\$([A-Za-z_][A-Za-z_0-9]*)',r"r['\1']",formula)
  mapping={x:eval(code,{'r':{'AcquisitionStatus':x}}) for x in ('Planned','Setup','Running','Completed','Aborted','Finished','Done','Acquired','Stopped','Failed','Not For Analysis','Unknown','')}
  self.assertEqual([mapping[x] for x in ('Planned','Setup','Running','Completed','Aborted')],['Planned','Setup','Running','Completed','Aborted'])
  self.assertEqual([mapping[x] for x in ('Finished','Done','Acquired','Stopped','Failed','Not For Analysis','Unknown','')],['Completed','Completed','Completed','Aborted','Unknown','Unknown','Unknown','Unknown'])
 def test_narrow_upgrade_preserves_every_row_and_is_idempotent(self):
  c=Fake();before=copy.deepcopy(c.rows);actions=upgrade_acquisition_workflow_status(c,dry_run=False)
  self.assertEqual(c.rows,before);self.assertEqual(c.fields[0]['fields']['colRef'],4);self.assertEqual(c.cols[0]['fields']['type'],'Choice')
  options=json.loads(c.cols[0]['fields']['widgetOptions']);self.assertEqual(options['alignment'],'right');self.assertEqual(options['choices'],['Planned','Dummy','Setup','Running','Completed','Aborted','Unknown','Custom'])
  self.assertTrue(actions);self.assertEqual(upgrade_acquisition_workflow_status(c,dry_run=False),[])
 def test_visible_status_survives_existing_hidden_effective_field(self):
  c=Fake();c.fields[0]['fields']['width']=95;c.fields.append({'id':10,'fields':{'parentId':8,'colRef':4,'width':0}})
  upgrade_acquisition_workflow_status(c,dry_run=False)
  self.assertEqual(len(c.fields),1);self.assertEqual(c.fields[0]['fields']['colRef'],4);self.assertEqual(c.fields[0]['fields']['width'],95)
  self.assertEqual(upgrade_acquisition_workflow_status(c,dry_run=False),[])
 def test_unknown_formula_fails_closed(self):
  c=Fake();c.cols[1]['fields']['formula']='$custom'
  with self.assertRaises(ValueError):upgrade_acquisition_workflow_status(c,dry_run=False)
  self.assertFalse(c.calls)

class NormalizeTests(unittest.TestCase):
 def test_equivalents_only_preserve_all_other_fields_and_idempotence(self):
  from nmr_register.setup import normalize_acquisition_status
  class C:
   def __init__(self):self.rows=[{'id':i+1,'fields':{'AcquisitionStatus':v,'Assessment':'Use','Sample':9,'RawPresent':True,'Notes':'human'}} for i,v in enumerate(['Finished','Done','Acquired','Stopped','Failed','Unknown','Setup','Running','Planned','Completed','Aborted','Not For Analysis',''])];self.writes=[]
   def records(self,t):return copy.deepcopy(self.rows)
   def apply(self,actions):
    self.writes.append(copy.deepcopy(actions))
    for _,_,ids,fields in actions:
     for rid,status in zip(ids,fields['AcquisitionStatus']):self.rows[rid-1]['fields']['AcquisitionStatus']=status
  c=C();original=copy.deepcopy(c.rows);plan=normalize_acquisition_status(c,dry_run=False,batch_size=2)
  self.assertEqual(plan['changed'],4);self.assertEqual(len(plan['actions']),2)
  self.assertEqual([r['fields']['AcquisitionStatus'] for r in c.rows][:5],['Completed','Completed','Completed','Aborted','Failed'])
  for a,b in zip(original,c.rows):
   self.assertEqual({k:v for k,v in a['fields'].items() if k!='AcquisitionStatus'},{k:v for k,v in b['fields'].items() if k!='AcquisitionStatus'})
  self.assertEqual(normalize_acquisition_status(c,dry_run=False)['changed'],0)
 def test_changed_status_prevents_write(self):
  from nmr_register.setup import normalize_acquisition_status
  class C:
   def __init__(self):self.calls=0;self.writes=[]
   def records(self,t):self.calls+=1;return [{'id':1,'fields':{'AcquisitionStatus':'Done' if self.calls==1 else 'Running'}}]
   def apply(self,a):self.writes.append(a)
  c=C()
  with self.assertRaises(ValueError):normalize_acquisition_status(c,dry_run=False)
  self.assertFalse(c.writes)
