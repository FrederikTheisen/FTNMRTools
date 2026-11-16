import importlib.util
from pathlib import Path
import shutil
import subprocess
import unittest

path=Path(__file__).resolve().parents[1]/'nmr_register'/'experiment_visibility.py'
spec=importlib.util.spec_from_file_location('experiment_visibility',path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
CHECKS=r'''
const assert=x=>{if(!x)throw Error('assertion failed')};
const rows=[{id:1,ExperimentType:'Calibration',ReconciliationReason:''},{id:2,ExperimentType:'1D',ReconciliationReason:'identity issue'},{id:3,ExperimentType:'T1',ReconciliationReason:''},{id:4,ExperimentType:'1D',ReconciliationReason:''},{id:5,ExperimentType:'Calibration',ScanError:'acquisition mismatch'}];
let x=ExperimentVisibility.select(rows,false);assert(JSON.stringify(x.ids)==='[2,3,5]'&&x.hidden===2&&x.conflicted===2);
x=ExperimentVisibility.select(rows,true);assert(JSON.stringify(x.ids)==='[1,2,3,4,5]'&&x.hidden===0&&x.conflicted===2);
for(const invalid of [null,[{id:1},{id:1}],[{id:'1'}]]){let caught=false;try{ExperimentVisibility.select(invalid,false)}catch(e){caught=true}assert(caught)}
globalThis.visibilityChecks='passed';
'''
class VisibilityTests(unittest.TestCase):
 def test_html(self):
  html=module.experiment_visibility_html(3,22)
  for text in ('Loading experiment list','Could not load experiment list','Other filters still apply','grist.setSelectedRows([])','Show calibration / 1D'):
   self.assertIn(text,html)
  self.assertNotIn('applyUserActions',html)
  for value in (0,True,'3'):
   with self.assertRaises(ValueError):module.experiment_visibility_html(value,22)
 @unittest.skipUnless(shutil.which('node'),'Node CLI unavailable')
 def test_selection(self):subprocess.run(['node','-e',module.CORE_JS+CHECKS],check=True)

class InstallerTests(unittest.TestCase):
 def test_preflight_preserves_filter_and_layout_without_writes(self):
  import json
  from nmr_register.experiments_ui import install_experiment_visibility
  class Client:
   data={'_grist_Views':[{'id':8,'fields':{'name':'Experiments','layoutSpec':json.dumps({'children':[{'leaf':22,'size':6},{'leaf':23,'size':4}],'collapsed':[{'leaf':45}]})}}],
    '_grist_Tables':[{'id':1,'fields':{'tableId':'Experiments'}}],
    '_grist_Tables_column':[{'id':58,'fields':{'parentId':1,'colId':'ExperimentType'}}],
    '_grist_Views_section':[{'id':22,'fields':{'parentId':8,'parentKey':'record','title':'Experiments','tableRef':1}}],
    '_grist_Filters':[{'id':37,'fields':{'viewSectionRef':22,'colRef':58,'filter':json.dumps({'excluded':['Calibration','Custom','1D'],'flag':'keep'})}}]}
   def records(self,name):return self.data[name]
   def apply(self,actions):raise AssertionError('dry run wrote')
  c=Client();result=install_experiment_visibility(c,'doc');self.assertEqual(result['filter_value'],{'excluded':['Custom'],'flag':'keep'});self.assertEqual(result['original_layout']['collapsed'],[{'leaf':45}])
  c.data['_grist_Filters'][0]['fields']['filter']='{"included":["T1"]}'
  with self.assertRaises(ValueError):install_experiment_visibility(c,'doc')

 def test_install_and_rerun_preserve_layout_and_other_fields(self):
  import copy,json
  from unittest.mock import patch
  from nmr_register.experiments_ui import install_experiment_visibility
  from nmr_register.navigation import apply_experiment_navigation
  class Client:
   def __init__(self):
    self.data={'_grist_Views':[{'id':8,'fields':{'name':'Experiments','layoutSpec':json.dumps({'children':[{'leaf':22,'size':6},{'leaf':23,'size':4}],'collapsed':[{'leaf':36},{'leaf':45}]})}},{'id':12,'fields':{'name':'Electrostatics'}}],
     '_grist_Tables':[{'id':1,'fields':{'tableId':'Experiments'}}], '_grist_Tables_column':[{'id':58,'fields':{'parentId':1,'colId':'ExperimentType'}}],
     '_grist_Views_section':[{'id':22,'fields':{'parentId':8,'parentKey':'record','title':'Experiments','tableRef':1,'linkSrcSectionRef':0}},{'id':40,'fields':{'parentId':12,'parentKey':'record','title':'Individual acquisitions','tableRef':1}}],
     '_grist_Filters':[{'id':37,'fields':{'viewSectionRef':22,'colRef':58,'filter':json.dumps({'excluded':['Calibration','Custom','1D']})}}], '_grist_Views_section_field':[{'id':70,'fields':{'parentId':22,'colRef':58,'width':199}}]};self.calls=[]
   def records(self,name):return copy.deepcopy(self.data[name])
   def apply(self,actions):
    self.calls.extend(copy.deepcopy(actions))
    for a in actions:
     assert a[0]=='UpdateRecord'
     row=next(r for r in self.data[a[1]] if r['id']==a[2]);row['fields'].update(a[3])
  c=Client()
  def create(client,name,args):
   assert name=='grist_add_page_widget'
   c.data['_grist_Views_section'].append({'id':50,'fields':{'parentId':8,'parentKey':'custom','title':'Experiment visibility','tableRef':1,'options':'{}'}})
   return {'id':50}
  with patch('nmr_register.experiments_ui._mcp',create):
   result=install_experiment_visibility(c,'doc',apply=True)
   self.assertEqual(result['widget'],50);self.assertEqual(install_experiment_visibility(c,'doc',apply=True)['actions'],[])
  layout=json.loads(c.data['_grist_Views'][0]['fields']['layoutSpec'])
  self.assertEqual(layout['collapsed'],[{'leaf':36},{'leaf':45}]);self.assertEqual(layout['children'][1:],[{'leaf':22,'size':6},{'leaf':23,'size':4}])
  self.assertEqual(c.data['_grist_Views_section_field'][0]['fields']['width'],199)
  self.assertEqual(json.loads(c.data['_grist_Filters'][0]['fields']['filter'])['excluded'],['Custom'])
  self.assertEqual(apply_experiment_navigation(c)['grid'],22)
  c.data['_grist_Views_section'][-1]['fields']['linkSrcSectionRef']=22
  with self.assertRaises(ValueError):apply_experiment_navigation(c)
  with self.assertRaises(ValueError):install_experiment_visibility(c,'doc')

 def test_interrupted_install_and_collapsed_reapply(self):
  import copy,json
  from nmr_register.experiments_ui import install_experiment_visibility
  from nmr_register.experiment_visibility import experiment_visibility_html
  class Client:
   def __init__(self):
    self.data={'_grist_Views':[{'id':8,'fields':{'name':'Experiments','layoutSpec':json.dumps({'children':[{'leaf':22,'size':6},{'leaf':23,'size':4},{'leaf':50}],'collapsed':[{'leaf':36},{'leaf':45}]})}}],
      '_grist_Tables':[{'id':1,'fields':{'tableId':'Experiments'}}], '_grist_Tables_column':[{'id':58,'fields':{'parentId':1,'colId':'ExperimentType'}}],
      '_grist_Views_section':[{'id':22,'fields':{'parentId':8,'parentKey':'record','title':'Experiments','tableRef':1}},{'id':50,'fields':{'parentId':8,'parentKey':'custom','title':'Experiment visibility','tableRef':1,'options':'{}'}}],
      '_grist_Filters':[{'id':37,'fields':{'viewSectionRef':22,'colRef':58,'filter':'{}'}}]}
   def records(self,name):return copy.deepcopy(self.data[name])
   def apply(self,actions):
    for a in actions:next(r for r in self.data[a[1]] if r['id']==a[2])['fields'].update(a[3])
  c=Client();install_experiment_visibility(c,'doc',apply=True)
  self.assertEqual(c.data['_grist_Views_section'][0]['fields']['linkSrcSectionRef'],50)
  self.assertEqual(json.loads(c.data['_grist_Views_section'][1]['fields']['options'])['html'],experiment_visibility_html(50,22))
  layout=json.loads(c.data['_grist_Views'][0]['fields']['layoutSpec']);self.assertEqual(layout['children'][0]['leaf'],50);self.assertEqual(layout['children'][1:],[{'leaf':22,'size':6},{'leaf':23,'size':4}]);self.assertEqual(layout['collapsed'],[{'leaf':36},{'leaf':45}])
  c.data['_grist_Views'][0]['fields']['layoutSpec']=json.dumps({'children':[{'leaf':22,'size':6},{'leaf':23,'size':4}],'collapsed':[{'leaf':36},{'leaf':50},{'leaf':45}]})
  install_experiment_visibility(c,'doc',apply=True);layout=json.loads(c.data['_grist_Views'][0]['fields']['layoutSpec']);self.assertEqual(layout['children'][0]['leaf'],50);self.assertEqual(layout['collapsed'],[{'leaf':36},{'leaf':45}]);self.assertEqual(install_experiment_visibility(c,'doc',apply=True)['actions'],[])
