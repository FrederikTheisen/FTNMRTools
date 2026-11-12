import copy,json,unittest
from nmr_register.setup import upgrade_sample_status
from nmr_register.schema import columns,SAMPLE_STATUSES
class Fake:
 def __init__(self,values):
  self.col={'id':5,'fields':{'parentId':1,'colId':'Available','type':'Choice','isFormula':False,'formula':'','recalcWhen':0,'recalcDeps':None,'widgetOptions':json.dumps({'choices':['Available','Unavailable','Unknown'],'alignment':'left'})}}
  self.rows=[{'id':i+1,'fields':{'Available':v,'Notes':'keep'}} for i,v in enumerate(values)];self.actions=[]
 def records(self,table):return {'_grist_Tables':[{'id':1,'fields':{'tableId':'Samples'}}], '_grist_Tables_column':[self.col], 'Samples':self.rows}[table]
 def apply(self,actions):
  self.actions.append(actions)
  for a in actions:
   if a[0]=='ModifyColumn':self.col['fields'].update(a[3])
   else:
    for i,rid in enumerate(a[2]):next(r for r in self.rows if r['id']==rid)['fields']['Available']=a[3]['Available'][i]
class Tests(unittest.TestCase):
 def test_migrates_known_values_and_is_idempotent(self):
  c=Fake([None,'','Unknown','Available','Unavailable','Suspected degradation','Degraded','Retired']);upgrade_sample_status(c)
  self.assertEqual([r['fields']['Available'] for r in c.rows],['Usable','Usable','Usable','Usable','Dead','Suspected degradation','Degraded','Dead'])
  self.assertTrue(all(r['fields']['Notes']=='keep' for r in c.rows));self.assertEqual(json.loads(c.col['fields']['widgetOptions'])['alignment'],'left')
  self.assertEqual(upgrade_sample_status(c),[]);self.assertEqual(len(c.actions),1)
 def test_unknown_value_or_custom_formula_blocks_without_writes(self):
  c=Fake(['experimental']);before=copy.deepcopy(c.rows)
  with self.assertRaisesRegex(ValueError,'requires review'):upgrade_sample_status(c)
  self.assertEqual(c.rows,before);self.assertEqual(c.actions,[])
  c=Fake(['Available']);c.col['fields']['formula']='"custom"'
  with self.assertRaisesRegex(ValueError,'requires review'):upgrade_sample_status(c)
  self.assertEqual(c.actions,[])
 def test_reviewed_status_requires_exact_old_value_and_evidence(self):
  c=Fake(['']);c.rows[0]['fields']['Notes']='Looking a bit degraded'
  evidence={1:{'expected_status':'','expected_notes':'Looking a bit degraded','status':'Suspected degradation'}}
  upgrade_sample_status(c,reviewed=evidence);self.assertEqual(c.rows[0]['fields']['Available'],'Suspected degradation')
  self.assertEqual(upgrade_sample_status(c,reviewed=evidence),[])
  c=Fake(['']);c.rows[0]['fields']['Notes']='Changed'
  with self.assertRaisesRegex(ValueError,'evidence changed'):upgrade_sample_status(c,reviewed=evidence)
  self.assertEqual(c.actions,[])
 def test_fresh_status_and_other_conditions(self):
  sample={c['id']:c for c in columns()['Samples']};self.assertEqual(sample['Available']['formula'],'"Usable"')
  self.assertEqual(json.loads(sample['Available']['widgetOptions'])['choices'],list(SAMPLE_STATUSES))
  self.assertEqual(sample['Protein']['type'],'Choice')
if __name__=='__main__':unittest.main()

class StyleTests(unittest.TestCase):
 def test_status_styles_append_and_preserve_custom_rules(self):
  from unittest.mock import patch
  from nmr_register.samples_ui import apply_sample_status_styles,STATUS_ROW_STYLES
  old={'rowRules':[{'formula':'$Protein == "special"','style':{'fontBold':True}}], 'columnRules':{'Notes':[{'formula':'True','style':{'fontItalic':True}}]}}
  calls=[];state=copy.deepcopy(old)
  def mcp(client,name,args):
   if name=='grist_get_table_conditional_styles':return copy.deepcopy(state)
   self.assertEqual(name,'grist_set_table_conditional_styles');calls.append(args['styles']);state.update(copy.deepcopy(args['styles']));return {}
  with patch('nmr_register.samples_ui._mcp',side_effect=mcp):
   apply_sample_status_styles(None,'doc');apply_sample_status_styles(None,'doc')
  self.assertEqual(len(calls),1);self.assertEqual(calls[0]['columnRules'],old['columnRules']);self.assertEqual(calls[0]['rowRules'][0],old['rowRules'][0]);self.assertEqual(calls[0]['rowRules'][1:],list(STATUS_ROW_STYLES))
 def test_retired_style_converts_without_duplicate_and_preserves_other_rules(self):
  from unittest.mock import patch
  from nmr_register.samples_ui import apply_sample_status_styles,STATUS_ROW_STYLES
  old={'formula':'$Available == "Retired"','style':{'fillColor':'#EDF0F2','textColor':'#52616B'}}
  state={'rowRules':[STATUS_ROW_STYLES[0],STATUS_ROW_STYLES[1],old],'columnRules':{}}
  calls=[]
  def mcp(client,name,args):
   if name=='grist_get_table_conditional_styles':return copy.deepcopy(state)
   calls.append(args['styles']);state.update(copy.deepcopy(args['styles']));return {}
  with patch('nmr_register.samples_ui._mcp',side_effect=mcp):
   apply_sample_status_styles(None,'doc');apply_sample_status_styles(None,'doc')
  self.assertEqual(len(calls),1);self.assertEqual(state['rowRules'],list(STATUS_ROW_STYLES))
 def test_status_style_conflict_requires_review(self):
  from unittest.mock import patch
  from nmr_register.samples_ui import apply_sample_status_styles
  state={'rowRules':[{'formula':'$Available == "Degraded"','style':{'fillColor':'#000000'}}],'columnRules':{}}
  with patch('nmr_register.samples_ui._mcp',return_value=state),self.assertRaisesRegex(ValueError,'requires review'):
   apply_sample_status_styles(None,'doc')
