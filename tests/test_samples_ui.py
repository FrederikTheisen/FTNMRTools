import unittest
from unittest.mock import patch
from nmr_register.samples_ui import apply_samples_ui, apply_more_details, unpin_sample_series_filter, TABLE_FIELDS, DETAIL_FIELDS, EXTRA_FIELDS

class Fake:
    def __init__(self):
        self.cols=[]; self.actions=[]; self.widgets=[{'id':1,'title':'Series','type':'table'},{'id':2,'title':'Samples','type':'table','filters':[{'column_id':'DNA_uM','filter':'{}','pinned':True},{'column_id':'Series','filter':'{\"included\":[1]}','pinned':True}]},{'id':3,'title':'Sample details','type':'card'},{'id':4,'title':'Sample experiments','type':'table'},{'id':5,'title':'Storage grid','type':'table'},{'id':6,'title':'Register samples','type':'creator'}]
        self.visible={}; self.layouts=[]; self.card_layouts={}; self.adds=0; self.styles={'rowRules':[],'columnRules':{}}
    def records(self,table):
        if table=='_grist_Tables': return [{'id':1,'fields':{'tableId':'Samples'}}]
        if table=='_grist_Tables_column': return self.cols
        raise AssertionError(table)
    def apply(self,actions):
        self.actions.extend(actions)
        for a in actions:
            if a[0]=='AddColumn': self.cols.append({'id':100+len(self.cols),'fields':{'parentId':1,'colId':a[2],**a[3]}})
    def mcp(self,client,name,args):
        if name=='grist_get_table_conditional_styles':return self.styles
        if name=='grist_set_table_conditional_styles':self.styles=args['styles'];return {}
        if name=='grist_get_pages': return [{'id':9,'name':'Samples'}]
        if name=='grist_get_page_widgets': return self.widgets
        if name=='grist_add_page_widget':
            self.adds+=1; w={'id':7,'title':args['widget_options']['title'],'type':'card','table_id':'Samples'};self.widgets.append(w);return w
        if name=='grist_update_page_widget':
            widget=next(w for w in self.widgets if w['id']==args['widget_id'])
            widget.update(args['widget_options']);return {}
        if name=='grist_set_widget_fields': self.visible[args['widget_id']]=args['column_ids'];return {}
        if name=='grist_get_card_layout': return {'fields':[{'column_id':c,'field_id':i+1} for i,c in enumerate(self.visible[args['widget_id']])]}
        if name=='grist_set_card_layout': self.card_layouts[args['widget_id']]=args['layout'];return {}
        if name=='grist_get_page_layout': return {'widgets':self.widgets,'layout':None}
        if name=='grist_set_page_layout': self.layouts.append(args['layout']);return {}
        return {}

class SamplesLayoutTests(unittest.TestCase):
    def test_existing_old_title_is_renamed_without_duplicate(self):
        c=Fake();c.widgets.append({'id':7,'title':'Additional conditions','type':'card','table_id':'Samples'})
        with patch('nmr_register.samples_ui._mcp',side_effect=c.mcp):
            self.assertEqual(apply_more_details(c,'doc'),7)
            self.assertEqual(apply_more_details(c,'doc'),7)
        self.assertEqual(c.adds,0)
        self.assertEqual(c.widgets[-1]['title'],'More sample details')

    def test_unpin_series_filter_alone_preserves_condition(self):
        c=Fake()
        with patch('nmr_register.samples_ui._mcp',side_effect=c.mcp):
            result=unpin_sample_series_filter(c,'doc');again=unpin_sample_series_filter(c,'doc')
        self.assertEqual(result,again)
        series=next(f for f in result if f['column_id']=='Series')
        self.assertFalse(series['pinned']);self.assertEqual(series['filter'],'{\"included\":[1]}')
        dna=next(f for f in result if f['column_id']=='DNA_uM');self.assertTrue(dna['pinned'])
        self.assertEqual(c.actions,[]);self.assertEqual(c.layouts,[])

    def test_compact_and_repeatable_without_record_mutation(self):
        c=Fake()
        with patch('nmr_register.samples_ui._mcp',side_effect=c.mcp), patch('nmr_register.storage_install.apply_storage_widget',return_value={'grid':5}):
            one=apply_samples_ui(c,'doc');two=apply_samples_ui(c,'doc')
        self.assertEqual(one,two);self.assertEqual(c.adds,1)
        series=next(f for f in c.widgets[1]['filters'] if f['column_id']=='Series')
        self.assertFalse(series['pinned']);self.assertEqual(series['filter'],'{\"included\":[1]}')
        self.assertEqual(c.visible[2],TABLE_FIELDS);self.assertEqual(c.visible[3],DETAIL_FIELDS)
        self.assertNotIn(5,c.visible)
        self.assertNotIn('DNA_uM',DETAIL_FIELDS);self.assertNotIn('PreparedDateText',DETAIL_FIELDS)
        self.assertEqual(c.visible[7],EXTRA_FIELDS)
        self.assertEqual(set(EXTRA_FIELDS),{'SampleID','PreparedAt','Series','Protein','Construct','Isotope',
            'Buffer','pH','NaCl_mM','Protein_uM','DNA_uM','TFE_pct','OtherConditions',
            'StorageRow','StorageColumn','CurrentLocation','Available','Notes','SourceEvidence'})
        self.assertNotIn('PreparedDateText',EXTRA_FIELDS);self.assertNotIn('SampleNumber',EXTRA_FIELDS)
        self.assertNotIn('DistinguishingCondition',EXTRA_FIELDS)
        expected_rows=[['SampleID'],['PreparedAt','Series'],['Protein','Construct','Isotope'],
            ['Buffer','pH','NaCl_mM'],['Protein_uM','DNA_uM','TFE_pct'],['OtherConditions'],
            ['StorageColumn','StorageRow'],['CurrentLocation','Available'],['Notes'],['SourceEvidence']]
        ids={name:i+1 for i,name in enumerate(EXTRA_FIELDS)}
        expected_layout=[{'leaf':ids[row[0]]} if len(row)==1 else {'children':[{'leaf':ids[name]} for name in row]} for row in expected_rows]
        self.assertEqual(c.card_layouts[7]['children'],expected_layout)
        layout=c.layouts[-1];self.assertEqual(layout['children'],[{'leaf':2,'size':7},{'leaf':3,'size':3}])
        self.assertEqual({x['leaf'] for x in layout['collapsed']},{1,4,5,6,7})
        self.assertTrue(all(a[0]=='AddColumn' and a[1]=='Samples' for a in c.actions))

class StorageInstallTests(unittest.TestCase):
    def test_convert_in_place_preserves_options_and_is_repeatable(self):
        import json
        from nmr_register.storage_install import plan_storage_widget
        class Client:
            def __init__(self):
                self.data={'_grist_Tables':[{'id':1,'fields':{'tableId':'Samples'}},{'id':2,'fields':{'tableId':'StorageGrid'}}],
                    '_grist_Views':[{'id':9,'fields':{'name':'Samples'}}],
                    '_grist_Tables_column':[{'id':10,'fields':{'parentId':1,'colId':'StorageRow','label':'Storage row'}},{'id':11,'fields':{'parentId':1,'colId':'StorageColumn','label':'Storage column'}},{'id':12,'fields':{'parentId':2,'colId':'Row'}}],
                    '_grist_Views_section':[{'id':25,'fields':{'parentId':9,'title':'Samples','tableRef':1,'parentKey':'record'}},
                        {'id':28,'fields':{'parentId':9,'title':'Storage grid','tableRef':2,'parentKey':'record','linkSrcSectionRef':0,'options':json.dumps({'zebraStripes':True,'customView':''})}},
                        {'id':37,'fields':{'parentId':9,'title':'More sample details','tableRef':1,'parentKey':'single','layoutSpec':json.dumps({'children':[{'children':[{'leaf':30},{'leaf':31}]}]})}}],
                    '_grist_Views_section_field':[{'id':20,'fields':{'parentId':28,'colRef':12}},{'id':30,'fields':{'parentId':37,'colRef':10}},{'id':31,'fields':{'parentId':37,'colRef':11}}]}
            def records(self,name):return self.data[name]
            def apply(self,actions):
                for a in actions:
                    if a[0]=='RemoveRecord':self.data[a[1]][:]=[r for r in self.data[a[1]] if r['id']!=a[2]]
                    else:next(r for r in self.data[a[1]] if r['id']==a[2])['fields'].update(a[3])
        c=Client();html='<body data-nmr-storage="1"></body>';p=plan_storage_widget(c,html)
        self.assertEqual((p['grid'],p['samples']),(28,25));self.assertTrue(all(a[1].startswith('_grist_') for a in p['actions']))
        c.apply(p['actions']);rack=c.data['_grist_Views_section'][1]['fields'];opts=json.loads(rack['options'])
        self.assertEqual((rack['parentKey'],rack['tableRef'],rack['linkSrcSectionRef']),('creator',1,25));self.assertTrue(opts['zebraStripes']);self.assertEqual(json.loads(opts['customView'])['access'],'full')
        card=c.data['_grist_Views_section'][2]['fields'];self.assertEqual(json.loads(card['layoutSpec'])['children'][0]['children'],[{'leaf':31},{'leaf':30}])
        self.assertEqual(plan_storage_widget(c,html)['actions'],[])
