import unittest
from nmr_register.schema import MACHINE_FIELDS, columns, trigger_actions, snapshot_values, series_defaults

class SchemaTests(unittest.TestCase):
    def test_ownership_and_snapshots(self):
        self.assertFalse({'Study','Sample','Assessment'} & MACHINE_FIELDS['Experiments'])
        self.assertEqual(snapshot_values({'NaCl_mM':0})['SnapshotNaCl_mM'],0)
        self.assertEqual(series_defaults({'NaCl_mM':50},{'NaCl_mM':0})['NaCl_mM'],0)
        ex=columns()['Experiments']; ids={(table,c['id']):i+1 for table, cs in columns().items() for i,c in enumerate(cs)}
        actions=trigger_actions({},ids)
        self.assertTrue(all(a[3]['recalcDeps']==['L',ids[('Experiments','Sample')]] for a in actions[:10]))
        self.assertTrue(all(not c['isFormula'] for c in ex if c['id'].startswith('Snapshot')))
    def test_human_override_field_uses_null_check(self):
        ex={c['id']:c for c in columns()['Experiments']}
        self.assertEqual('$Field_MHz',ex['EffectiveField_MHz']['formula'])
        self.assertIn('$Assessment',ex['AttentionReason']['formula'])
        self.assertFalse(ex['AssessmentReason']['isFormula'])
    def test_overdue_accepts_acquired_status_and_excludes_planned(self):
        import re
        from datetime import datetime, timedelta, timezone
        formula=next(c['formula'] for c in columns()['Experiments'] if c['id']=='AttentionReason')
        source=re.sub(r'\$([A-Za-z_][A-Za-z_0-9]*)',r"fields['\1']",formula)
        namespace={}
        exec('def evaluate(fields,NOW):\n'+''.join('    '+line+'\n' for line in source.splitlines()),namespace)
        now=datetime(2026,1,3,tzinfo=timezone.utc)
        base=dict(Study='Electrostatics',ExperimentCategory='Other',ExperimentType='',MigrationIssue='',SampleTagIssue='',AssessmentExplanation='',ReconciliationReason='',MetadataDifference='',ScanError='',RegisteredAt=1,Sample=1,
                  Assessment='',AssessmentDisplay='Not assessed',AssessmentReason='',RawPresent=False,AcquisitionStatus='',EffectiveAcquisitionStatus='Completed',
                  EffectiveAcquiredAt=now-timedelta(days=2))
        check=lambda **updates: namespace['evaluate']({**base,**updates},lambda:now)
        self.assertIn('unassessed',check())
        self.assertEqual(check(EffectiveAcquisitionStatus='Planned',RawPresent=True),'')
        self.assertEqual(check(EffectiveAcquisitionStatus='Unknown',RawPresent=True),'Acquisition status unknown')
        self.assertEqual(check(EffectiveAcquisitionStatus='Aborted'),'Acquisition aborted')

class LabellingChoicesTests(unittest.TestCase):
    def test_upgrade_preserves_values_and_trigger_settings(self):
        import copy,json
        from nmr_register.schema import LABELLING_CHOICES
        from nmr_register.setup import upgrade_labelling_choices
        class Fake:
            def __init__(self):
                self.cols=[{'id':1,'fields':{'parentId':1,'colId':'Isotope','type':'Text','formula':'', 'widgetOptions':'{}'}},
                    {'id':2,'fields':{'parentId':2,'colId':'Isotope','type':'Text','formula':'value or $Series.Isotope','recalcWhen':0,'recalcDeps':['L',17],'widgetOptions':json.dumps({'alignment':'left','choices':['special']})}}]
                self.rows={'Series':[{'id':1,'fields':{'Isotope':'15N custom'}}],
                           'Samples':[{'id':1,'fields':{'Isotope':'15N+13C'}},{'id':2,'fields':{'Isotope':'unusual label'}}]}
                self.batches=[]
            def records(self,name):
                if name=='_grist_Tables': return [{'id':1,'fields':{'tableId':'Series'}},{'id':2,'fields':{'tableId':'Samples'}}]
                if name=='_grist_Tables_column': return self.cols
                return self.rows[name]
            def apply(self,actions):
                self.batches.append(copy.deepcopy(actions))
                for _,table,col,fields in actions:
                    row=self.cols[0 if table=='Series' else 1]
                    row['fields'].update(fields)
        c=Fake();old=copy.deepcopy(c.rows);actions=upgrade_labelling_choices(c)
        self.assertEqual(len(actions),2);self.assertEqual(c.rows,old)
        sample=c.cols[1]['fields'];self.assertEqual(sample['formula'],'value or $Series.Isotope')
        self.assertEqual(sample['recalcDeps'],['L',17]);opts=json.loads(sample['widgetOptions'])
        self.assertEqual(opts['alignment'],'left');self.assertEqual(opts['choices'][:7],list(LABELLING_CHOICES))
        self.assertIn('special',opts['choices']);self.assertIn('unusual label',opts['choices'])
        self.assertIn('15N custom',json.loads(c.cols[0]['fields']['widgetOptions'])['choices'])
        self.assertEqual(upgrade_labelling_choices(c),[]);self.assertEqual(len(c.batches),1)
    def test_schema_uses_string_choices_but_snapshot_stays_text(self):
        for table in ('Series','Samples'):
            self.assertEqual(next(c for c in columns()[table] if c['id']=='Isotope')['type'],'Choice')
        self.assertEqual(next(c for c in columns()['Experiments'] if c['id']=='SnapshotIsotope')['type'],'Text')


class AcquisitionUpgradeTests(unittest.TestCase):
    def test_new_numeric_cells_are_null_in_same_action_batch(self):
        from nmr_register.setup import upgrade_acquisition_schema, LEGACY_ACQUISITION_FORMULAS
        from nmr_register.schema import MACHINE_FIELDS
        class Fake:
            def __init__(self):
                self.batches=[]
                self.cols=[{'id':i+1,'fields':{'parentId':1,'colId':name,'formula':formula}} for i,(name,formula) in enumerate(LEGACY_ACQUISITION_FORMULAS.items())]
                self.cols.append({'id':4,'fields':{'parentId':1,'colId':'PulseProgram','label':'PulseProgram'}})
                self.cols.append({'id':3,'fields':{'parentId':1,'colId':'D1_s','type':'Numeric','isFormula':False,'formula':''}})
            def records(self,name):
                return {'_grist_Tables':[{'id':1,'fields':{'tableId':'Experiments'}}],
                    '_grist_Tables_column':self.cols,'Experiments':[{'id':3,'fields':{'D1_s':0}},{'id':7,'fields':{}}]}[name]
            def apply(self,actions):
                self.batches.append(actions)
                for action in actions:
                    if action[0]=='AddColumn':
                        info=dict(action[3]);
                        if info.get('recalcDeps')==[['L']]: info['recalcDeps']=['L']
                        self.cols.append({'id':len(self.cols)+10,'fields':{'parentId':1,'colId':action[2],**info}})
                    elif action[0]=='ModifyColumn':
                        info=dict(action[3]);
                        if info.get('recalcDeps')==[['L']]: info['recalcDeps']=['L']
                        next(c for c in self.cols if c['fields']['colId']==action[2])['fields'].update(info)
                    elif action[0]=='UpdateRecord':
                        next(c for c in self.cols if c['id']==action[2])['fields'].update(action[3])
        c=Fake();actions=upgrade_acquisition_schema(c)
        self.assertEqual(len(c.batches),1);self.assertEqual(actions,c.batches[0])
        default=next(a for a in actions if a[0]=='ModifyColumn' and a[2]=='D1_s')
        self.assertEqual(default[3],{'formula':'None','recalcWhen':0,'recalcDeps':[['L']]})
        self.assertTrue(any(a[0]=='AddColumn' for a in actions))
        self.assertTrue(all(a[3].get('recalcDeps')==[['L']] for a in actions if a[0]=='AddColumn' and a[3]['type']=='Numeric' and not a[3]['isFormula']))
        self.assertFalse(any(a[0]=='AddVisibleColumn' for a in actions))
        bulk=next(a for a in actions if a[0]=='BulkUpdateRecord')
        self.assertEqual(bulk[2],[3,7]);self.assertNotIn('D1_s',bulk[3])
        self.assertTrue(all(values==[None,None] for values in bulk[3].values()))
        self.assertEqual(set(bulk[3]),{'SpinLock_Hz','SpinLockOffset_Hz'})
        for name in ('D1_s','DelayList_s','SpinLock_Hz','SpinLockOffset_Hz','ObservedProtocolEvidence'):
            self.assertIn(name,MACHINE_FIELDS['Experiments'])
        self.assertEqual(upgrade_acquisition_schema(c),[])
        for row in c.cols:
            if row['fields'].get('recalcDeps') == ['L']: row['fields']['recalcDeps']=None
        self.assertEqual(upgrade_acquisition_schema(c),[])
        self.assertEqual(len(c.batches),1)
    def test_core_pilot_bootstraps_before_electrostatics_fields_exist(self):
        from nmr_register.setup import upgrade_acquisition_schema,LEGACY_ACQUISITION_FORMULAS
        class CorePilot:
            def __init__(self): self.actions=[]
            def records(self,name):
                return {'_grist_Tables':[{'id':1,'fields':{'tableId':'Experiments'}}],
                    '_grist_Tables_column':[{'id':2,'fields':{'parentId':1,'colId':'EffectiveField_MHz',
                        'formula':LEGACY_ACQUISITION_FORMULAS['EffectiveField_MHz']}},
                        {'id':3,'fields':{'parentId':1,'colId':'PulseProgram','label':'PulseProgram'}}],
                    'Experiments':[{'id':9,'fields':{}}]}[name]
            def apply(self,actions): self.actions=actions
        client=CorePilot();upgrade_acquisition_schema(client)
        self.assertTrue(any(a[0]=='AddColumn' and a[2]=='EffectiveD1_s' for a in client.actions))
        self.assertTrue(any(a[0]=='BulkUpdateRecord' and a[2]==[9] for a in client.actions))
        self.assertFalse(any(a[0]=='UpdateRecord' and a[2]==2 and a[3].get('formula')==LEGACY_ACQUISITION_FORMULAS.get('CoverageConditions') for a in client.actions))
        class UnexpectedFormula(CorePilot):
            def records(self,name):
                rows=super().records(name)
                if name=='_grist_Tables_column': rows[0]['fields']['formula']='$custom'
                return rows
        bad=UnexpectedFormula()
        with self.assertRaisesRegex(ValueError,'Existing formula differs'):
            upgrade_acquisition_schema(bad)
        self.assertEqual(bad.actions,[])

    def test_effective_numeric_zero_and_nominal_proton_field(self):
        experiment={c['id']:c for c in columns()['Experiments']}
        self.assertEqual(experiment['DelayList_s']['label'],'Delay list (s)')
        self.assertEqual(experiment['PulseProgram']['label'],'Pulse program')
        self.assertEqual(experiment['EffectiveD1_s']['formula'],'$D1_s')
        self.assertEqual(experiment['EffectiveField_MHz']['formula'],'$Field_MHz')


class CompletionStatusTests(unittest.TestCase):
    def test_canonical_legacy_and_unknown(self):
        import re
        experiment={c['id']:c for c in columns()['Experiments']}
        expr=re.sub(r'\$([A-Za-z_][A-Za-z_0-9]*)',r"r['\1']",experiment['EffectiveAcquisitionStatus']['formula'])
        def effective(value):return eval(expr,{'r':{'AcquisitionStatus':value}})
        self.assertEqual([effective(x) for x in ('','Stopped','Acquired','Done','Unknown')],['Unknown','Aborted','Completed','Completed','Unknown'])
        self.assertIn('AcquisitionStatus',MACHINE_FIELDS['Experiments'])
        self.assertNotIn('ObservedAcquisitionStatus',MACHINE_FIELDS['Experiments'])

    def test_existing_widgets_receive_effective_status_without_page_reset(self):
        from nmr_register.setup import apply_completion_ui
        class Fake:
            def __init__(self):
                names=['EffectiveAcquisitionStatus','AcquisitionStatus','ObservedAcquisitionStatus','ObservedAcquisitionEvidence','Assessment']
                self.cols=[{'id':i+1,'fields':{'parentId':1,'colId':name}} for i,name in enumerate(names)]
                self.sections=[{'id':sid,'fields':{'parentId':page,'title':title}} for sid,page,title in
                    [(11,1,'Experiments'),(12,1,'Experiment details'),(21,2,'Needs attention'),
                     (22,2,'Attention details'),(31,3,'Individual acquisitions'),(32,3,'Acquisition details')]]
                self.fields=[{'id':i+1,'fields':{'parentId':sid,'colRef':ref,'parentPos':1}} for i,(sid,ref) in
                    enumerate([(11,5),(12,2),(12,3),(21,5),(22,2),(31,2),(32,2)])]
                self.batches=[]
            def records(self,name):return {'_grist_Tables':[{'id':1,'fields':{'tableId':'Experiments'}}],
                '_grist_Tables_column':self.cols,'_grist_Views':[{'id':1,'fields':{'name':'Experiments'}},
                {'id':2,'fields':{'name':'Needs attention'}},{'id':3,'fields':{'name':'Electrostatics'}}],
                '_grist_Views_section':self.sections,'_grist_Views_section_field':self.fields}[name]
            def apply(self,actions):
                self.batches.append(actions)
                for a in actions:
                    if a[0]=='AddRecord':self.fields.append({'id':100+len(self.fields),'fields':a[3]})
                    elif a[0]=='UpdateRecord':next(r for r in self.fields if r['id']==a[2])['fields'].update(a[3])
                    else:self.fields[:]=[r for r in self.fields if r['id']!=a[2]]
        c=Fake();apply_completion_ui(c)
        for sid in (11,21,31):
            refs={f['fields']['colRef'] for f in c.fields if f['fields']['parentId']==sid}
            self.assertIn(2,refs);self.assertNotIn(1,refs)
        for sid in (12,22,32):
            refs={f['fields']['colRef'] for f in c.fields if f['fields']['parentId']==sid}
            self.assertTrue({2,4} <= refs);self.assertNotIn(3,refs)
        self.assertEqual(apply_completion_ui(c),[]);self.assertEqual(len(c.batches),1)
        self.assertTrue(all(a[1]=='_grist_Views_section_field' for a in c.batches[0]))

    def test_narrow_upgrade_accepts_canonical_schema(self):
        from nmr_register.setup import upgrade_completion_schema
        expected=next(c for c in columns()['Experiments'] if c['id']=='EffectiveAcquisitionStatus')
        class Fake:
            def records(self,name):
                return {'_grist_Tables':[{'id':1,'fields':{'tableId':'Experiments'}}],
                    '_grist_Tables_column':[{'id':1,'fields':{'parentId':1,'colId':'EffectiveAcquisitionStatus','formula':expected['formula']}}]}[name]
        self.assertEqual(upgrade_completion_schema(Fake()),[])


class LayoutTests(unittest.TestCase):
    def test_layout_builds_pages_links_summary_and_storage(self):
        from nmr_register.setup import build_layout
        class Fake:
            def __init__(self):
                self.actions=[]; self.added=[]; self.next=100
                self.tables=[]; self.cols=[]; self.sections=[]; self.rack_creator=False
                for i,(name,cs) in enumerate(columns().items(),1):
                    self.tables.append({'id':i,'fields':{'tableId':name}})
                    for c in cs:
                        self.next+=1; self.cols.append({'id':self.next,'fields':{'parentId':i,'colId':c['id']}})
            def records(self,name):
                return {'_grist_Views':[],'_grist_Tables':self.tables,'_grist_Tables_column':self.cols,'_grist_Views_section':self.sections,'_grist_Views_section_field':[],'_grist_Filters':[],'StorageGrid':[],'Samples':[]}[name]
            def apply(self,actions):
                self.actions.extend(actions); result=[]
                for a in actions:
                    self.next+=1
                    if a[0]=='CreateViewSection':
                        ref=a[1]
                        if a[4] is not None:
                            ref=self.next+9000; self.tables.append({'id':ref,'fields':{'tableId':'Summary'+str(ref),'summarySourceTable':a[1]}})
                            for source in a[4]:
                                original=next(x for x in self.cols if x['id']==source)
                                self.cols.append({'id':10000+len(self.cols),'fields':{'parentId':ref,'colId':original['fields']['colId'],'summarySourceCol':source}})
                            self.cols.append({'id':10000+len(self.cols),'fields':{'parentId':ref,'colId':'count'}})
                        kind=a[3]
                        if self.rack_creator and ref==next(r['id'] for r in self.tables if r['fields']['tableId']=='StorageGrid'):
                            ref=next(r['id'] for r in self.tables if r['fields']['tableId']=='Samples');kind='creator'
                        self.sections.append({'id':self.next,'fields':{'tableRef':ref,'parentId':a[2],'parentKey':kind}})
                    if a[0]=='UpdateRecord' and a[1]=='_grist_Views_section':
                        next(r for r in self.sections if r['id']==a[2])['fields'].update(a[3])
                    if a[0]=='AddVisibleColumn':
                        ref=next(r['id'] for r in self.tables if r['fields']['tableId']==a[1]); self.cols.append({'id':10000+len(self.cols),'fields':{'parentId':ref,'colId':a[2]}})
                    result.append({'id':self.next} if a[0]=='AddView' else {'sectionRef':self.next} if a[0]=='CreateViewSection' else {})
                return {'retValues':result}
            def add(self,table,fields): self.added.append((table,fields))
        from unittest.mock import patch
        c=Fake()
        with patch('nmr_register.navigation.apply_experiment_navigation') as nav:
            build_layout(c)
        nav.assert_called_once_with(c,dry_run=False)
        self.assertEqual([a[3] for a in c.actions if a[0]=='AddView'],['Experiments','Samples','Needs attention','Overview','Electrostatics'])
        self.assertEqual(len(c.added),8)
        summaries=[a for a in c.actions if a[0]=='CreateViewSection' and a[4] is not None]
        self.assertEqual(len(summaries),2)
        grouped={next(col['fields']['colId'] for col in c.cols if col['id']==ref) for ref in summaries[-1][4]}
        self.assertTrue({'EffectivePulseProgram','EffectiveD1_s','EffectiveDelayList_s','EffectiveSpinLock_Hz','EffectiveSpinLockOffset_Hz'} <= grouped)
        self.assertTrue(any(a[0]=='AddRecord' and a[1]=='_grist_Filters' for a in c.actions))
        creator=Fake();creator.rack_creator=True
        with patch('nmr_register.navigation.apply_experiment_navigation'):
            build_layout(creator)
        rack=next(r['id'] for r in creator.sections if r['fields'].get('title')=='Storage grid')
        self.assertFalse(any(a[0]=='AddRecord' and a[1]=='_grist_Views_section_field' and a[3].get('parentId')==rack for a in creator.actions))

class ExperimentPresentationTests(unittest.TestCase):
    def test_compact_layout_hides_helpers_and_preserves_dropdown_options(self):
        import json
        from unittest.mock import patch
        from nmr_register.setup import apply_experiments_ui
        class Fake:
            def __init__(self):
                self.columns=[{'id':i,'fields':{'parentId':1,'colId':name,'label':name,'widgetOptions':opts}}
                    for i,(name,opts) in enumerate([
                        ('Sample',json.dumps({'dropdownCondition':{'text':'(not $Series or choice.Series == $Series) and (not $ProteinSelector or choice.Protein == $ProteinSelector)','extra':'keep'},'alignment':'left'})),
                        ('AcquisitionStatus',''),('AcquiredAt',''),('AcquisitionDateText',''),('ProteinSelector','')],1)]
                used={r['fields']['colId'] for r in self.columns}
                self.columns.extend({'id':i+20,'fields':{'parentId':1,'colId':c['id'],'label':c['id'],'widgetOptions':''}} for i,c in enumerate(columns()['Experiments']) if c['id'] not in used)
                self.fields=[{'id':1,'fields':{'parentId':11,'colRef':4}},
                             {'id':2,'fields':{'parentId':11,'colRef':5}},
                             {'id':3,'fields':{'parentId':21,'colRef':4}},
                             {'id':4,'fields':{'parentId':21,'colRef':5}},
                             {'id':5,'fields':{'parentId':11,'colRef':1,'widgetOptions':json.dumps({'dropdownCondition':{'text':'(not $Series or choice.Series == $Series) and (not $ProteinSelector or choice.Protein == $ProteinSelector)','other':'keep'},'wrap':True})}}]
                self.sections=[{'id':10,'fields':{'parentId':1,'title':'Experiments','options':'{}'}},{'id':21,'fields':{'parentId':2,'title':'Attention details'}}]
                self.batches=[];self.layouts=[]
            def records(self,name):
                return {'_grist_Tables':[{'id':1,'fields':{'tableId':'Experiments'}}],
                    '_grist_Tables_column':self.columns,
                    '_grist_Views_section':self.sections,
                    '_grist_Views_section_field':self.fields}[name]
            def apply(self,actions):
                self.batches.append(actions)
                for action in actions:
                    if action[0]=='RemoveRecord': self.fields[:]=[f for f in self.fields if f['id']!=action[2]]
                    elif action[0]=='AddRecord': self.fields.append({'id':100+len(self.fields),'fields':action[3]})
                    else:
                        target={'_grist_Tables_column':self.columns,'_grist_Views_section':self.sections,'_grist_Views_section_field':self.fields}[action[1]]
                        next(r for r in target if r['id']==action[2])['fields'].update(action[3])
            def mcp(self,client,name,args):
                if name=='grist_get_pages': return [{'id':1,'name':'Experiments'},{'id':2,'name':'Needs attention'}]
                if name=='grist_get_page_widgets': return [{'id':10,'title':'Experiments'},{'id':11,'title':'Experiment details'},{'id':12,'title':'Register experiment'},{'id':13,'title':'Other'}]
                if name=='grist_get_page_layout': return {'widgets':[{'id':10},{'id':11},{'id':12},{'id':13}]}
                if name=='grist_set_page_layout': self.layouts.append(args['layout'])
                return {}
        c=Fake()
        with patch('nmr_register.setup._mcp',side_effect=c.mcp), patch('nmr_register.experiments_ui._mcp',side_effect=c.mcp), patch('nmr_register.experiments_ui.configure_experiment_type_filter',return_value={'filter':99}) as type_filter:
            first=apply_experiments_ui(c,'doc');second=apply_experiments_ui(c,'doc')
        self.assertEqual(type_filter.call_count,2);self.assertEqual(first['type_filter'],{'filter':99});type_filter.assert_called_with(c,apply=True)
        self.assertEqual(first,second);self.assertEqual(len(c.batches),1)
        self.assertTrue({f['id'] for f in c.fields}.issuperset({5}))
        opts=json.loads(c.columns[0]['fields']['widgetOptions'])
        self.assertEqual(opts['alignment'],'left');self.assertEqual(opts['dropdownCondition']['extra'],'keep')
        self.assertNotIn('ProteinSelector',opts['dropdownCondition']['text'])
        fieldopts=json.loads(c.fields[0]['fields']['widgetOptions'])
        self.assertEqual(fieldopts['dropdownCondition']['other'],'keep');self.assertTrue(fieldopts['wrap'])
        self.assertEqual(c.layouts[-1],{'children':[{'leaf':10,'size':7},{'leaf':11,'size':3}],
                                        'collapsed':[{'leaf':12},{'leaf':13}]})
        custom=Fake();column_custom=json.dumps({'dropdownCondition':{'text':'choice.Available != "Unavailable"'},'wrap':True})
        field_custom=json.dumps({'dropdownCondition':{'text':'choice.Protein == "X"'},'other':'keep'})
        custom.columns[0]['fields']['widgetOptions']=column_custom
        custom.fields[-1]['fields']['widgetOptions']=field_custom
        with patch('nmr_register.setup._mcp',side_effect=custom.mcp), patch('nmr_register.experiments_ui._mcp',side_effect=custom.mcp), patch('nmr_register.experiments_ui.configure_experiment_type_filter',return_value={'filter':99}): apply_experiments_ui(custom,'doc')
        self.assertEqual(custom.columns[0]['fields']['widgetOptions'],column_custom)
        self.assertEqual(custom.fields[0]['fields']['widgetOptions'],field_custom)
        blocked=Fake();blocked.fields[-1]['fields']['widgetOptions']=json.dumps({'dropdownCondition':{'text':'$ProteinSelector and choice.Available'}})
        with patch('nmr_register.setup._mcp',side_effect=blocked.mcp), patch('nmr_register.experiments_ui._mcp',side_effect=blocked.mcp), self.assertRaisesRegex(ValueError,'Custom sample dropdown'):
            apply_experiments_ui(blocked,'doc')
        self.assertEqual(blocked.batches,[])

    def test_dropdown_migration_only_rewrites_known_legacy(self):
        from nmr_register.setup import _migrate_sample_dropdown
        legacy='( not $Series or choice.Series == $Series ) and ( not $ProteinSelector or choice.Protein == $ProteinSelector )'
        known={'dropdownCondition':{'text':legacy,'extra':'retained'},'wrap':True}
        migrated=_migrate_sample_dropdown(known)
        self.assertEqual(migrated['dropdownCondition']['text'],'not $Series or choice.Series == $Series')
        self.assertEqual(migrated['dropdownCondition']['extra'],'retained')
        custom={'dropdownCondition':{'text':'choice.Available != "Unavailable"'},'wrap':True}
        self.assertEqual(_migrate_sample_dropdown(custom),custom)
        with self.assertRaisesRegex(ValueError,'Custom sample dropdown'):
            _migrate_sample_dropdown({'dropdownCondition':{'text':'$ProteinSelector and choice.Available'}})

class SingleRackTests(unittest.TestCase):
    def test_upgrade_reuses_rows_preserves_boxes_and_is_idempotent(self):
        from nmr_register.setup import upgrade_single_rack
        from nmr_register.samples_ui import LEGACY_STORAGE_FORMULA,STORAGE_FORMULA
        class Fake:
            def __init__(self):
                self.tables=[{'id':1,'fields':{'tableId':'Samples'}},{'id':2,'fields':{'tableId':'StorageGrid'}}]
                self.cols=[{'id':1,'fields':{'parentId':1,'colId':'StorageDisplay','isFormula':True,'formula':LEGACY_STORAGE_FORMULA}}]
                for n in range(1,13):
                    self.cols.append({'id':n+1,'fields':{'parentId':2,'colId':'Col'+str(n),'isFormula':True,
                      'formula':'", ".join(s.SampleID for s in Samples.lookupRecords(StorageBox=$StorageBox, StorageRow=$Row, StorageColumn='+str(n)+'))'}})
                self.grid=[{'id':i+1,'fields':{'Row':row,'StorageBox':'Unspecified box'}} for i,row in enumerate('ABCDEFGH')]
                self.batches=[]
            def records(self,table):return {'_grist_Tables':self.tables,'_grist_Tables_column':self.cols,'StorageGrid':self.grid}[table]
            def apply(self,actions):
                self.batches.append(actions)
                for a in actions:
                    if a[0]=='UpdateRecord':next(c for c in self.cols if c['id']==a[2])['fields'].update(a[3])
                    elif a[0]=='AddRecord':self.grid.append({'id':100+len(self.grid),'fields':a[3]})
        c=Fake(); original=[r['fields'].copy() for r in c.grid];actions=upgrade_single_rack(c)
        self.assertEqual(len(actions),13);self.assertTrue(all(a[0]=='UpdateRecord' for a in actions))
        self.assertEqual([r['fields'] for r in c.grid],original)
        self.assertEqual(c.cols[0]['fields']['formula'],STORAGE_FORMULA)
        self.assertTrue(all('StorageBox' not in r['fields']['formula'] for r in c.cols[1:]))
        self.assertEqual(upgrade_single_rack(c),[]);self.assertEqual(len(c.batches),1)
        c.grid.pop();new=upgrade_single_rack(c)
        self.assertEqual(new,[['AddRecord','StorageGrid',None,{'Row':'H'}]])
        c.grid.append({'id':500,'fields':{'Row':'A','StorageBox':'other'}})
        with self.assertRaisesRegex(ValueError,'duplicated'): upgrade_single_rack(c)
        self.assertEqual(len(c.batches),2)

class ProteinChoiceTests(unittest.TestCase):
    def test_upgrade_preserves_values_and_defaults(self):
        import copy,json
        from nmr_register.setup import upgrade_protein_choices
        class Fake:
            def __init__(self):
                self.tables=[{'id':1,'fields':{'tableId':'Series'}},{'id':2,'fields':{'tableId':'Samples'}},{'id':3,'fields':{'tableId':'Experiments'}}]
                self.cols=[{'id':11,'fields':{'parentId':1,'colId':'Protein','type':'Text','formula':'','widgetOptions':json.dumps({'choices':['manual choice'],'alignment':'left'})}},
                           {'id':12,'fields':{'parentId':2,'colId':'Protein','type':'Text','formula':'value if value else $Series.Protein','recalcDeps':['L',99],'recalcWhen':0,'widgetOptions':json.dumps({'choices':['Other'],'choiceOptions':{'Other':{'fillColor':'red'}}})}},
                           {'id':13,'fields':{'parentId':3,'colId':'SnapshotProtein','type':'Text','formula':'$Sample.Protein','recalcDeps':['L',42]}}]
                self.data={'Series':[{'id':1,'fields':{'Protein':'HeTDR1'}}], 'Samples':[{'id':2,'fields':{'Protein':'hetDR1'}},{'id':3,'fields':{'Protein':' HeTDR1 '}}], 'Experiments':[{'id':4,'fields':{'SnapshotProtein':'LegacyName'}}]};self.applied=[]
            def records(self,t):return {'_grist_Tables':self.tables,'_grist_Tables_column':self.cols,**self.data}[t]
            def apply(self,actions):
                self.applied.append(actions)
                for a in actions:
                    col=next(x for x in self.cols if x['fields'].get('parentId')==({'Series':1,'Samples':2}[a[1]]) and x['fields'].get('colId')==a[2])
                    col['fields'].update(a[3])
        c=Fake();before=copy.deepcopy(c.data);original=copy.deepcopy(c.cols)
        actions=upgrade_protein_choices(c)
        self.assertEqual(len(actions),2);self.assertEqual(c.data,before)
        names=['manual choice','HeTDR1','Other','hetDR1',' HeTDR1 ','LegacyName']
        for col in c.cols[:2]:self.assertEqual(json.loads(col['fields']['widgetOptions'])['choices'],names)
        self.assertEqual(c.cols[1]['fields']['formula'],original[1]['fields']['formula'])
        self.assertEqual(c.cols[1]['fields']['recalcDeps'],['L',99]);self.assertEqual(c.cols[2],original[2])
        self.assertEqual(json.loads(c.cols[1]['fields']['widgetOptions'])['choiceOptions'],{'Other':{'fillColor':'red'}})
        self.assertEqual(upgrade_protein_choices(c),[]);self.assertEqual(len(c.applied),1)
    def test_fresh_schema_keeps_snapshot_text_and_sample_trigger(self):
        ex={c['id']:c for c in columns()['Experiments']};samples={c['id']:c for c in columns()['Samples']};series={c['id']:c for c in columns()['Series']}
        self.assertEqual(ex['SnapshotProtein']['type'],'Text');self.assertEqual(samples['Protein']['type'],'Choice');self.assertEqual(series['Protein']['type'],'Choice')
        self.assertIn('$Series.Protein',samples['Protein']['formula'])

class AttentionSourceTests(unittest.TestCase):
    def test_label_handles_discovered_planned_partial_and_conflict(self):
        import re
        from nmr_register.schema import SOURCE_LABEL_FORMULA
        code=re.sub(r'\$([A-Za-z_][A-Za-z_0-9]*)',r"row['\1']",SOURCE_LABEL_FORMULA)
        namespace={};exec('def run(row):\n'+''.join('    '+line+'\n' for line in code.splitlines()),namespace)
        blank=dict(Folder='',Expno='',SourceID='',ObservedFolder='',ObservedExpno='',ObservedSourceID='',AcquisitionKey='')
        def run(**values):return namespace['run']({**blank,**values})
        self.assertEqual(run(Folder='long/folder',Expno='4'),'4 · long/folder')
        self.assertEqual(run(ObservedFolder='folder',ObservedExpno='9',ObservedSourceID='src'),'9 · folder')
        self.assertEqual(run(ObservedExpno='4',ObservedSourceID='src',AcquisitionKey='src:4'),'4 · (root)')
        self.assertEqual(run(Folder='manual',ObservedExpno='7',ObservedFolder='observed'),'Observed: 7 · observed | Recorded: ? · manual')
        self.assertEqual(run(Folder='f',Expno='1',SourceID='a',ObservedFolder='f',ObservedExpno='1',ObservedSourceID='b'),'Observed (different source): 1 · f | Recorded (different source): 1 · f')
        self.assertEqual(run(),'')
    def test_upgrade_is_guarded_idempotent_and_preserves_custom_labels(self):
        from nmr_register.setup import upgrade_source_label,LEGACY_SOURCE_LABEL_FORMULA
        class Client:
            def __init__(self,label):self.row={'id':193,'fields':{'parentId':4,'colId':'SourceLabel','type':'Text','isFormula':True,'formula':LEGACY_SOURCE_LABEL_FORMULA,'label':label}};self.actions=[]
            def records(self,t):return {'_grist_Tables':[{'id':4,'fields':{'tableId':'Experiments'}}],'_grist_Tables_column':[self.row]}[t]
            def apply(self,actions):
                self.actions.append(actions)
                for action in actions:
                    fields=self.row['fields'];update=action[3]
                    if 'label' in update and not fields.get('untieColIdFromLabel'):fields['colId']='renamed_by_grist'
                    fields.update(update)
        c=Client('Acquisition source');planned=upgrade_source_label(c);self.assertEqual(len(planned),2);self.assertEqual(c.actions,[])
        self.assertEqual(planned[0][3],{'untieColIdFromLabel':True})
        upgrade_source_label(c,dry_run=False);self.assertEqual(c.row['fields']['label'],'EXPNO · Folder');self.assertEqual(c.row['fields']['colId'],'SourceLabel');self.assertEqual(upgrade_source_label(c,dry_run=False),[])
        c=Client('Custom source');upgrade_source_label(c,dry_run=False);self.assertEqual(c.row['fields']['label'],'Custom source')
        c.row['fields']['formula']='$custom'
        with self.assertRaisesRegex(ValueError,'requires review'):upgrade_source_label(c,dry_run=False)

    def test_setup_migrates_existing_source_before_schema_validation(self):
        from unittest.mock import patch
        from contextlib import ExitStack
        from nmr_register.setup import create_schema
        class Client:
            def __init__(self,has_source):self.has_source=has_source
            def records(self,table):
                if table=='_grist_Tables':return [{'id':4,'fields':{'tableId':'Experiments'}}]
                if table=='_grist_Tables_column':return ([{'id':193,'fields':{'parentId':4,'colId':'SourceLabel'}}] if self.has_source else [])
                if table in ('_grist_Pages','_grist_TabBar'):return []
                raise AssertionError(table)
            def apply(self,actions):pass
        for has_source in (True,False):
            calls=[]
            with ExitStack() as stack:
                stack.enter_context(patch('nmr_register.setup.schema_actions',return_value=[]))
                stack.enter_context(patch('nmr_register.setup.trigger_actions',return_value=[]))
                stack.enter_context(patch('nmr_register.setup.apply_experiment_date_format'))
                for name in ('upgrade_acquisition_workflow_status','upgrade_completion_schema','upgrade_acquisition_schema','upgrade_temperature_grouping','upgrade_source_label','upgrade_electrostatics_schema','upgrade_sample_tag_schema','upgrade_assessment_schema'):
                    stack.enter_context(patch('nmr_register.setup.'+name,side_effect=lambda *a,_name=name,**k:calls.append(_name)))
                create_schema(Client(has_source))
            if has_source:self.assertLess(calls.index('upgrade_source_label'),calls.index('upgrade_electrostatics_schema'))
            else:self.assertNotIn('upgrade_source_label',calls)

class SplitAttentionIdentityTests(unittest.TestCase):
    def test_separate_cells_cover_conflicting_and_partial_identities(self):
        import re
        from nmr_register.schema import attention_identity_formula
        functions={}
        for part in ('folder','expno'):
            code=re.sub(r'\$([A-Za-z_][A-Za-z_0-9]*)',r"row['\1']",attention_identity_formula(part))
            scope={};exec('def run(row):\n'+''.join('    '+line+'\n' for line in code.splitlines()),scope);functions[part]=scope['run']
        blank=dict(Folder='',Expno='',SourceID='',ObservedFolder='',ObservedExpno='',ObservedSourceID='',AcquisitionKey='')
        def run(**values):
            row={**blank,**values};return functions['folder'](row),functions['expno'](row)
        self.assertEqual(run(Folder='planned',Expno='4'),('planned','4'))
        self.assertEqual(run(ObservedFolder='seen',ObservedExpno='9',ObservedSourceID='source'),('seen','9'))
        self.assertEqual(run(ObservedExpno='4',ObservedSourceID='source',AcquisitionKey='source:4'),('(root)','4'))
        self.assertEqual(run(Folder='planned',ObservedFolder='seen',ObservedExpno='9'),('Observed: seen | Recorded: planned','Observed: 9 | Recorded: ?'))
        out=run(Folder='f',Expno='3',SourceID='secret1',ObservedFolder='f',ObservedExpno='3',ObservedSourceID='secret2')
        self.assertEqual(out,('Observed (different source): f | Recorded (different source): f','Observed (different source): 3 | Recorded (different source): 3'))
        self.assertNotIn('secret',' '.join(out))
    def test_upgrade_adds_missing_and_rejects_custom_formula(self):
        from nmr_register.setup import upgrade_attention_identity
        class Client:
            def __init__(self):self.rows=[];self.writes=[]
            def records(self,t):return {'_grist_Tables':[{'id':4,'fields':{'tableId':'Experiments'}}], '_grist_Tables_column':self.rows}[t]
            def apply(self,actions):
                self.writes.append(actions)
                for i,a in enumerate(actions):self.rows.append({'id':i+1,'fields':{'parentId':4,'colId':a[2],**a[3]}})
        c=Client();plan=upgrade_attention_identity(c);self.assertEqual(len(plan),2);self.assertEqual(c.writes,[])
        upgrade_attention_identity(c,dry_run=False);self.assertEqual(upgrade_attention_identity(c,dry_run=False),[])
        c.rows[0]['fields']['formula']='$custom'
        with self.assertRaisesRegex(ValueError,'requires review'):upgrade_attention_identity(c,dry_run=False)
        self.assertEqual(len(c.writes),1)
