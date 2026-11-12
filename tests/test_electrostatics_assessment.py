"""A separate same-table assessment tray preserves the Electrostatics matrix."""
import copy
import json
import unittest
from unittest.mock import patch
from nmr_register.electrostatics_assessment import install_electrostatics_assessment
from nmr_register.assessment_widget import assessment_html

class Client:
    def __init__(self):
        self.original_layout={'children':[{'leaf':43,'size':7},{'leaf':40,'size':3}],'collapsed':[{'leaf':41,'custom':'keep'}],'customLayoutOption':'keep'}
        self.data={'_grist_Views':[{'id':12,'fields':{'name':'Electrostatics','layoutSpec':json.dumps(self.original_layout)}},{'id':8,'fields':{'name':'Experiments','layoutSpec':'{"children":[{"leaf":22}],"collapsed":[{"leaf":45}]}'}}],
                   '_grist_Tables':[{'id':4,'fields':{'tableId':'Experiments'}}],
                   '_grist_Views_section':[self.section(43,12,'Electrostatics matrix','creator'),self.section(40,12,'Individual acquisitions','record',43),self.section(41,12,'Acquisition details','single',40),self.section(22,8,'Experiments','record'),self.section(45,8,'Assess experiment','creator',22)]}
        self.html={45:assessment_html(45)};self.settings={45:{'access':'full'}};self.calls=[];self.actions=[];self.fail_html=False
    def section(self,id,page,title,kind,link=0):
        return {'id':id,'fields':{'parentId':page,'title':title,'parentKey':kind,'tableRef':4,'linkSrcSectionRef':link,'linkSrcColRef':0,'linkTargetColRef':0,'options':'{"customSetting":"keep"}'}}
    def records(self,table):return copy.deepcopy(self.data[table])
    def apply(self,actions):
        self.actions.extend(copy.deepcopy(actions))
        for _,table,id,fields in actions:next(r for r in self.data[table] if r['id']==id)['fields'].update(fields)
    def call(self,client,name,args):
        self.calls.append((name,copy.deepcopy(args)))
        if name=='grist_add_page_widget':
            id=46;self.data['_grist_Views_section'].append(self.section(id,args['page_id'],'Assess experiment','creator'))
            view=self.data['_grist_Views'][0]['fields'];layout=json.loads(view['layoutSpec']);layout['children'].append({'leaf':id});view['layoutSpec']=json.dumps(layout)
            self.html[id]='';self.settings[id]={'access':'none','customSetting':'keep'};return {'id':id}
        id=args['widget_id']
        if name=='grist_get_creator_html':return {'html':self.html[id]}
        if name=='grist_set_creator_html':
            if self.fail_html:self.fail_html=False;raise RuntimeError('synthetic interrupted upload')
            self.html[id]=args['html'];return {'html':self.html[id]}
        if name=='grist_get_custom_widget_settings':return copy.deepcopy(self.settings[id])
        if name=='grist_set_custom_widget_settings':self.settings[id]['access']=args['access'];return copy.deepcopy(self.settings[id])
        raise AssertionError(name)
    def install(self):
        with patch('nmr_register.electrostatics_assessment._mcp',side_effect=self.call):return install_electrostatics_assessment(self,'doc')

class ElectrostaticsAssessmentTests(unittest.TestCase):
    def test_create_reuses_assessment_and_preserves_layout_options_other_page_and_data(self):
        c=Client();other=copy.deepcopy(c.data['_grist_Views_section'][-1]);other_layout=c.data['_grist_Views'][1]['fields']['layoutSpec']
        result=c.install();self.assertEqual(result,{'page':12,'grid':40,'card':46})
        expected=copy.deepcopy(c.original_layout);expected['collapsed'].append({'leaf':46})
        self.assertEqual(json.loads(c.data['_grist_Views'][0]['fields']['layoutSpec']),expected)
        card=next(r for r in c.data['_grist_Views_section'] if r['id']==46)['fields']
        self.assertEqual((card['linkSrcSectionRef'],card['linkSrcColRef'],card['linkTargetColRef']),(40,0,0))
        self.assertEqual(card['options'],'{"customSetting":"keep"}')
        self.assertEqual(c.html[46],assessment_html(46));self.assertEqual(c.settings[46],{'access':'full','customSetting':'keep'})
        self.assertEqual(c.data['_grist_Views_section'][4],other);self.assertEqual(c.data['_grist_Views'][1]['fields']['layoutSpec'],other_layout)
        writes=len(c.actions);calls=len(c.calls);self.assertEqual(c.install(),result)
        self.assertEqual(len(c.actions),writes)
        self.assertTrue(all(name.startswith('grist_get_') for name,_ in c.calls[calls:]))
        self.assertEqual(sum(name=='grist_add_page_widget' for name,_ in c.calls),1)

    def test_interrupted_creation_resumes_without_duplicate_and_collapses_existing_visible_card(self):
        c=Client();c.fail_html=True
        with self.assertRaisesRegex(RuntimeError,'interrupted'):c.install()
        self.assertEqual(c.install()['card'],46)
        self.assertEqual(sum(name=='grist_add_page_widget' for name,_ in c.calls),1)
        self.assertEqual(json.loads(c.data['_grist_Views'][0]['fields']['layoutSpec'])['children'],c.original_layout['children'])

    def test_existing_custom_editor_or_link_requires_review(self):
        for mutate in (lambda c:next(r for r in c.data['_grist_Views_section'] if r['id']==46)['fields'].update(tableRef=9),
                       lambda c:next(r for r in c.data['_grist_Views_section'] if r['id']==46)['fields'].update(linkSrcSectionRef=22),
                       lambda c:c.html.update({46:'<html>custom user tool</html>'}),
                       lambda c:c.settings[46].update(unmapped_required_columns=['required']),
                       lambda c:c.data['_grist_Views_section'].append(c.section(47,12,'Custom dependent','record',46))):
            c=Client();c.install();mutate(c);before=copy.deepcopy(c.data);actions=len(c.actions);calls=len(c.calls)
            with self.assertRaises(ValueError):c.install()
            self.assertEqual(c.data,before);self.assertEqual(len(c.actions),actions)
            self.assertFalse(any(name.startswith('grist_set') or name=='grist_add_page_widget' for name,_ in c.calls[calls:]))

    def test_invalid_layout_selector_or_cycle_fails_before_creation(self):
        for mutate in (lambda c:c.data['_grist_Views'][0]['fields'].update(layoutSpec='broken'),
                       lambda c:c.data['_grist_Views'][0]['fields'].update(layoutSpec='{"children":[{"leaf":999}]}'),
                       lambda c:c.data['_grist_Views'][0]['fields'].update(layoutSpec='{"children":[{"leaf":40},{"leaf":40}]}'),
                       lambda c:c.data['_grist_Views_section'][1]['fields'].update(tableRef=9),
                       lambda c:c.data['_grist_Views_section'][0]['fields'].update(linkSrcSectionRef=40),
                       lambda c:c.data['_grist_Views_section'].append(c.section(47,12,'Individual acquisitions','record'))):
            c=Client();mutate(c)
            with self.assertRaises(ValueError):c.install()
            self.assertEqual(c.actions,[]);self.assertEqual(c.calls,[])

    def test_setup_integrates_after_matrix_and_navigation(self):
        from nmr_register.setup import apply_electrostatics_ui
        order=[]
        with patch('nmr_register.coverage_install.install_coverage_matrix',side_effect=lambda *a:order.append('matrix') or {}),patch('nmr_register.navigation.apply_experiment_navigation',side_effect=lambda *a,**k:order.append('navigation') or {}),patch('nmr_register.electrostatics_assessment.install_electrostatics_assessment',side_effect=lambda *a:order.append('assessment') or {'card':46}):
            result=apply_electrostatics_ui(object(),'doc')
        self.assertEqual(order,['matrix','navigation','assessment']);self.assertEqual(result['assessment'],{'card':46})
