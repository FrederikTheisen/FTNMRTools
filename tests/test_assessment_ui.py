import copy,json,unittest
from unittest.mock import patch
from nmr_register.setup import apply_assessment_ui

class Fake:
 def __init__(self,extra=False):
  self.tables=[{'id':4,'fields':{'tableId':'Experiments','rawViewSectionRef':11,'recordCardViewSectionRef':12}}]
  self.views=[{'id':8,'fields':{'name':'Experiments','layoutSpec':json.dumps({'children':[{'leaf':22,'size':7},{'leaf':23,'size':3}],'collapsed':[{'leaf':36}]})}},
              {'id':13,'fields':{'name':'Assessment','layoutSpec':json.dumps({'children':[{'children':[{'leaf':44},{'leaf':45}]}]})}}]
  self.sections=[{'id':i,'fields':{'parentId':p,'title':title,'parentKey':kind,'tableRef':4,'linkSrcSectionRef':link,'linkSrcColRef':0,'linkTargetColRef':0,'options':'original HTML'}} for i,p,title,kind,link in [(22,8,'Experiments','record',0),(23,8,'Experiment details','single',22),(36,8,'Register experiment','creator',0),(44,13,'Experiments to assess','record',0),(45,13,'Assessment checklist','creator',44)]]
  if extra:
   self.sections.append({'id':46,'fields':{'parentId':13,'title':'Custom','parentKey':'single','tableRef':4,'linkSrcSectionRef':0}})
   self.views[1]['fields']['layoutSpec']=json.dumps({'children':[{'leaf':44},{'leaf':45},{'leaf':46}]})
  self.cols=[{'id':1,'fields':{'parentId':4,'colId':'AssessmentLink','type':'Text','isFormula':True,'formula':'SELF_HYPERLINK(label="Assess", page=13) + "#a1.s44.r" + str($id)','widgetOptions':'{"widget":"HyperLink"}'}}]
  self.actions=[];self.calls=[];self.section_fields=[{'id':i,'fields':{'parentId':i,'colRef':1}} for i in (11,12,22)]
 def records(self,t):return {'_grist_Tables':self.tables,'_grist_Views':self.views,'_grist_Views_section':self.sections,'_grist_Tables_column':self.cols,'_grist_Views_section_field':self.section_fields,'_grist_Pages':[{'id':1,'fields':{'viewRef':8,'pagePos':1,'indentation':0}},{'id':2,'fields':{'viewRef':13,'pagePos':2,'indentation':0}}]}[t]
 def apply(self,actions):
  self.actions.extend(copy.deepcopy(actions))
  for action in actions:
   if action[0]=='RemoveColumn':
    self.cols[:]=[r for r in self.cols if r['fields']['colId']!=action[2]];continue
   _,table,rid,fields=action
   collection={'_grist_Views':self.views,'_grist_Views_section':self.sections,'_grist_Tables_column':self.cols}[table]
   next(r for r in collection if r['id']==rid)['fields'].update(fields)
 def mcp(self,client,name,args):
  self.calls.append((name,args))
  if name=='grist_get_custom_widget_settings':return {'access':'full'}
  if name=='grist_remove_page':
   self.views[:]=[r for r in self.views if r['id']!=args['page_id']]
   self.sections[:]=[r for r in self.sections if r['fields']['parentId']!=args['page_id']]
   return {}
  raise AssertionError(name)

class InlineAssessmentTests(unittest.TestCase):
 def test_move_once_preserves_editor_and_links_selector_row(self):
  c=Fake();before=c.sections[-1]['fields']['options']
  with patch('nmr_register.setup._mcp',side_effect=c.mcp):
   first=apply_assessment_ui(c,'doc');action_count=len(c.actions)
   second=apply_assessment_ui(c,'doc')
  self.assertEqual(first['card'],45);self.assertEqual(second['card'],45)
  self.assertEqual(len(c.actions),action_count)
  card=next(r['fields'] for r in c.sections if r['id']==45)
  self.assertEqual((card['parentId'],card['linkSrcSectionRef'],card['linkSrcColRef'],card['linkTargetColRef']),(8,22,0,0))
  self.assertEqual(card['options'],before)
  self.assertEqual(json.loads(c.views[0]['fields']['layoutSpec'])['collapsed'],[{'leaf':36},{'leaf':45}])
  self.assertEqual(c.cols,[]);self.assertEqual(card['title'],'Assess experiment')
  self.assertEqual([v['fields']['name'] for v in c.views],['Experiments'])
 def test_custom_page_content_is_preserved(self):
  c=Fake(extra=True)
  with patch('nmr_register.setup._mcp',side_effect=c.mcp):apply_assessment_ui(c,'doc')
  self.assertTrue(any(v['id']==13 for v in c.views));self.assertTrue(any(r['id']==46 for r in c.sections))
  self.assertNotIn('grist_remove_page',[name for name,_ in c.calls])
 def test_custom_link_and_formula_block_before_move(self):
  for mode in ('link','formula','field'):
   c=Fake()
   if mode=='link':c.sections.append({'id':46,'fields':{'parentId':8,'title':'Linked','linkSrcSectionRef':45}})
   elif mode=='formula':c.cols[0]['fields']['formula']='$custom'
   else:c.section_fields.append({'id':50,'fields':{'parentId':99,'colRef':1}})
   with patch('nmr_register.setup._mcp',side_effect=c.mcp),self.assertRaises(ValueError):apply_assessment_ui(c,'doc')
   self.assertEqual(c.actions,[]);self.assertEqual(c.calls,[])
