import copy,unittest
from nmr_register.navigation import apply_experiment_navigation
class Fake:
 def __init__(self):
  self.tables={'_grist_Tables':[{'id':4,'fields':{'tableId':'Experiments'}}],
   '_grist_Views':[{'id':8,'fields':{'name':'Experiments','layoutSpec':'keep','filters':'keep'}},{'id':12,'fields':{'name':'Electrostatics','layoutSpec':'keep'}}],
   '_grist_Views_section':[{'id':22,'fields':{'parentId':8,'title':'Experiments','parentKey':'record','tableRef':4,'linkSrcSectionRef':0}},{'id':40,'fields':{'parentId':12,'title':'Individual acquisitions','parentKey':'record','tableRef':4,'linkSrcSectionRef':43}}],
   '_grist_Tables_column':[], '_grist_Views_section_field':[{'id':400,'fields':{'parentId':40,'colRef':77,'parentPos':3,'width':200}}]};self.calls=[]
 def records(self,t):return copy.deepcopy(self.tables[t])
 def apply(self,actions):
  self.calls+=copy.deepcopy(actions)
  for a in actions:
   if a[0]=='AddColumn':self.tables['_grist_Tables_column'].append({'id':99,'fields':{'parentId':4,'colId':a[2],**a[3]}})
   elif a[0]=='AddRecord':self.tables['_grist_Views_section_field'].append({'id':401,'fields':copy.deepcopy(a[3])})
class Tests(unittest.TestCase):
 def test_link_routes_current_row_and_preserves_layout_filter_and_existing_fields(self):
  c=Fake();before=copy.deepcopy(c.tables);p=apply_experiment_navigation(c,dry_run=False)
  column=c.tables['_grist_Tables_column'][0]['fields'];self.assertEqual(column['formula'],'SELF_HYPERLINK(label="Open", page=8) + "#a1.s22.r" + str($id)')
  self.assertEqual(c.tables['_grist_Views'],before['_grist_Views']);self.assertEqual(c.tables['_grist_Views_section'],before['_grist_Views_section'])
  self.assertEqual(c.tables['_grist_Views_section_field'][0],before['_grist_Views_section_field'][0]);self.assertEqual(c.tables['_grist_Views_section_field'][1]['fields']['colRef'],99)
  self.assertEqual(apply_experiment_navigation(c,dry_run=False)['actions'],[]);self.assertEqual(len(c.calls),2)
 def test_restores_field_after_setup_field_reset_and_uses_actual_page_ids(self):
  c=Fake();c.tables['_grist_Views'][0]['id']=15;c.tables['_grist_Views_section'][0]['fields']['parentId']=15;c.tables['_grist_Views_section'][0]['id']=66
  apply_experiment_navigation(c,dry_run=False)
  self.assertIn('page=15',c.tables['_grist_Tables_column'][0]['fields']['formula']);self.assertIn('#a1.s66.r',c.tables['_grist_Tables_column'][0]['fields']['formula'])
  c.tables['_grist_Views_section_field']=[r for r in c.tables['_grist_Views_section_field'] if r['fields']['colRef']!=99]
  apply_experiment_navigation(c,dry_run=False)
  self.assertEqual(sum(r['fields']['colRef']==99 for r in c.tables['_grist_Views_section_field']),1)
  self.assertEqual(sum(a[0]=='AddColumn' for a in c.calls),1)
 def test_unfamiliar_column_and_linked_destination_fail_closed(self):
  c=Fake();c.tables['_grist_Tables_column'].append({'id':99,'fields':{'parentId':4,'colId':'OpenExperiment','type':'Text','isFormula':True,'formula':'$other','widgetOptions':'{}'}})
  with self.assertRaises(ValueError):apply_experiment_navigation(c,dry_run=False)
  self.assertFalse(c.calls)
  c=Fake();c.tables['_grist_Views_section'][0]['fields']['linkSrcSectionRef']=7
  with self.assertRaises(ValueError):apply_experiment_navigation(c,dry_run=False)
  self.assertFalse(c.calls)
 def test_partial_install_resumes_and_user_hidden_field_is_preserved(self):
  c=Fake();plan=apply_experiment_navigation(c);c.apply(plan['actions']);apply_experiment_navigation(c,dry_run=False)
  c.tables['_grist_Views_section_field'][-1]['fields']['width']=0;self.assertEqual(apply_experiment_navigation(c,dry_run=False)['actions'],[])
