import unittest
from nmr_register.schema import MACHINE_FIELDS, columns, trigger_actions, snapshot_values, series_defaults

class SchemaTests(unittest.TestCase):
    def test_ownership_and_snapshots(self):
        self.assertFalse({'Sample','Assessment','Instrument','ExperimentType','AcquiredAt'} & MACHINE_FIELDS['Experiments'])
        self.assertEqual(snapshot_values({'NaCl_mM':0})['SnapshotNaCl_mM'],0)
        self.assertEqual(series_defaults({'NaCl_mM':50},{'NaCl_mM':0})['NaCl_mM'],0)
        ex=columns()['Experiments']; ids={(table,c['id']):i+1 for table, cs in columns().items() for i,c in enumerate(cs)}
        actions=trigger_actions({},ids)
        self.assertTrue(all(a[3]['recalcDeps']==['L',ids[('Experiments','Sample')]] for a in actions[:10]))
        self.assertTrue(all(not c['isFormula'] for c in ex if c['id'].startswith('Snapshot')))
    def test_human_override_field_uses_null_check(self):
        ex={c['id']:c for c in columns()['Experiments']}
        self.assertIn('is not None',ex['EffectiveField_MHz']['formula'])
        self.assertIn('$AssessmentReason',ex['AttentionReason']['formula'])
    def test_overdue_accepts_acquired_status_and_excludes_planned(self):
        import re
        from datetime import datetime, timedelta, timezone
        formula=next(c['formula'] for c in columns()['Experiments'] if c['id']=='AttentionReason')
        source=re.sub(r'\$([A-Za-z_][A-Za-z_0-9]*)',r"fields['\1']",formula)
        namespace={}
        exec('def evaluate(fields,NOW):\n'+''.join('    '+line+'\n' for line in source.splitlines()),namespace)
        now=datetime(2026,1,3,tzinfo=timezone.utc)
        base=dict(MigrationIssue='',ReconciliationReason='',MetadataDifference='',ScanError='',RegisteredAt=1,Sample=1,
                  Assessment='',AssessmentReason='',RawPresent=False,AcquisitionStatus='Acquired',
                  EffectiveAcquiredAt=now-timedelta(days=2))
        check=lambda **updates: namespace['evaluate']({**base,**updates},lambda:now)
        self.assertIn('unassessed',check())
        self.assertEqual(check(AcquisitionStatus='Planned'),'')
        self.assertIn('unassessed',check(AcquisitionStatus='Planned',RawPresent=True))

class LayoutTests(unittest.TestCase):
    def test_layout_builds_pages_links_summary_and_storage(self):
        from nmr_register.setup import build_layout
        class Fake:
            def __init__(self):
                self.actions=[]; self.added=[]; self.next=100
                self.tables=[]; self.cols=[]; self.sections=[]
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
                    if a[0]=='CreateViewSection' and a[4] is not None:
                        ref=9000; self.tables.append({'id':ref,'fields':{'tableId':'Summary','summarySourceTable':a[1]}})
                        self.sections.append({'id':self.next,'fields':{'tableRef':ref,'parentId':a[2]}})
                        for source in a[4]:
                            original=next(x for x in self.cols if x['id']==source)
                            self.cols.append({'id':10000+len(self.cols),'fields':{'parentId':ref,'colId':original['fields']['colId'],'summarySourceCol':source}})
                        self.cols.append({'id':10000+len(self.cols),'fields':{'parentId':ref,'colId':'count'}})
                    if a[0]=='AddVisibleColumn':
                        self.cols.append({'id':10000+len(self.cols),'fields':{'parentId':9000,'colId':a[2]}})
                    result.append({'id':self.next} if a[0]=='AddView' else {'sectionRef':self.next} if a[0]=='CreateViewSection' else {})
                return {'retValues':result}
            def add(self,table,fields): self.added.append((table,fields))
        c=Fake(); build_layout(c)
        self.assertEqual([a[3] for a in c.actions if a[0]=='AddView'],['Experiments','Samples','Needs attention','Overview'])
        self.assertEqual(len(c.added),8)
        self.assertEqual(len([a for a in c.actions if a[0]=='CreateViewSection' and a[4] is not None]),1)
        self.assertTrue(any(a[0]=='AddRecord' and a[1]=='_grist_Filters' for a in c.actions))
