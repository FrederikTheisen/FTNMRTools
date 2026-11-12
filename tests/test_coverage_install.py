import json,sys,types,unittest
from unittest.mock import patch
from nmr_register.coverage_install import install_coverage_matrix

class Fake:
 def __init__(self):
  self.widgets=[{'id':39,'title':'Condition coverage','table_id':'Summary','type':'table'}, {'id':40,'title':'Individual acquisitions','table_id':'Experiments','type':'table','selected_by_widget_id':39},{'id':41,'title':'Acquisition details','table_id':'Experiments','type':'card','selected_by_widget_id':40}]
  self.sections=[{'id':40,'fields':{'parentId':12,'parentKey':'record','title':'Individual acquisitions','tableRef':4,'linkSrcSectionRef':39,'linkSrcColRef':0,'linkTargetColRef':0}},{'id':41,'fields':{'linkSrcSectionRef':40}}]; self.html='';self.access='read table';self.deleted=[];self.layout=None
 def records(self,name):return ([{'id':4,'fields':{'tableId':'Experiments'}}] if name=='_grist_Tables' else self.sections)
 def apply(self,actions):
  for a in actions:next(r for r in self.sections if r['id']==a[2])['fields'].update(a[3])
 def mcp_call(self,name,args):
  a=args
  if name=='grist_get_pages':v=[{'name':'Electrostatics','id':12}]
  elif name=='grist_get_page_widgets':v=self.widgets
  elif name=='grist_add_page_widget':
   v={'id':42,'title':'Electrostatics matrix','table_id':'Experiments','type':'creator'};self.widgets.append(v);self.sections.append({'id':42,'fields':{'parentId':12,'parentKey':'creator','title':'Electrostatics matrix','tableRef':4,'linkSrcSectionRef':0}})
  elif name=='grist_get_creator_html':v={'html':self.html}
  elif name=='grist_set_creator_html':self.html=a['html'];v={}
  elif name=='grist_get_custom_widget_settings':v={'access':self.access}
  elif name=='grist_set_custom_widget_settings':self.access=a['access'];v={'access':self.access}
  elif name=='grist_get_page_widget_select_by_options':v=[{'link_from_widget_id':42,'link_from_column_id':None,'link_to_column_id':None}]
  elif name=='grist_set_page_widget_select_by':
   self.sections[0]['fields']['linkSrcSectionRef']=42;self.widgets[1]['selected_by_widget_id']=42;v={'matches_by':'row id'}
  elif name=='grist_remove_page_widget':self.deleted.append(a['widget_id']);self.widgets=[w for w in self.widgets if w['id']!=a['widget_id']];v={}
  elif name=='grist_set_page_layout':self.layout=a['layout'];v={}
  else:raise AssertionError(name)
  return {'content':[{'type':'text','text':json.dumps(v)}]}
class InstallTests(unittest.TestCase):
 def test_replaces_then_reuses(self):
  fake=Fake();mod=types.ModuleType('nmr_register.coverage_widget');mod.coverage_html=lambda section_id=None:'html:%s' % section_id
  with patch.dict(sys.modules,{'nmr_register.coverage_widget':mod}):
   r=install_coverage_matrix(fake,'doc');self.assertEqual(r['coverage'],42)
   self.assertEqual(fake.deleted,[39]);self.assertEqual(fake.layout['collapsed'],[{'leaf':41}]);self.assertEqual(fake.layout['children'],[{'leaf':42,'size':7},{'leaf':40,'size':3}])
   install_coverage_matrix(fake,'doc');self.assertEqual(fake.deleted,[39]);self.assertEqual(len(fake.widgets),3)
 def test_absent_option_uses_verified_metadata(self):
  fake=Fake();orig=fake.mcp_call
  def call(name,args):
   if name=='grist_get_page_widget_select_by_options':return {'content':[{'type':'text','text':'[]'}]}
   return orig(name,args)
  fake.mcp_call=call;mod=types.ModuleType('nmr_register.coverage_widget');mod.coverage_html=lambda section_id=None:'html'
  with patch.dict(sys.modules,{'nmr_register.coverage_widget':mod}):install_coverage_matrix(fake,'doc')
  self.assertEqual(fake.sections[0]['fields']['linkSrcSectionRef'],42);self.assertEqual(fake.deleted,[39])
 def test_wrong_table_fails_before_relink(self):
  fake=Fake();orig=fake.mcp_call
  def call(name,args):
   result=orig(name,args)
   if name=='grist_add_page_widget':fake.sections[-1]['fields']['tableRef']=77
   return result
  fake.mcp_call=call;mod=types.ModuleType('nmr_register.coverage_widget');mod.coverage_html=lambda section_id=None:'html'
  with patch.dict(sys.modules,{'nmr_register.coverage_widget':mod}):
   with self.assertRaisesRegex(ValueError,'section metadata'):install_coverage_matrix(fake,'doc')
  self.assertEqual(fake.sections[0]['fields']['linkSrcSectionRef'],39);self.assertEqual(fake.deleted,[])
 def test_dependency_fails_before_removal(self):
  fake=Fake();fake.sections.append({'id':55,'fields':{'linkSrcSectionRef':39}});mod=types.ModuleType('nmr_register.coverage_widget');mod.coverage_html=lambda section_id=None:'html:%s' % section_id
  with patch.dict(sys.modules,{'nmr_register.coverage_widget':mod}):
   with self.assertRaisesRegex(ValueError,'still has linked'):install_coverage_matrix(fake,'doc')
  self.assertEqual(fake.deleted,[])
