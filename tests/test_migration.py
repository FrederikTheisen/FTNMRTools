import copy
import unittest
from nmr_register.migration import apply_pilot

class Client:
    def __init__(self): self.tables={'Samples':[],'Experiments':[]}
    def records(self,t): return self.tables[t]
    def add(self,t,fields):
        row={'id':len(self.tables[t])+1,'fields':dict(fields)}; self.tables[t].append(row); return row['id']
class MigrationTests(unittest.TestCase):
    def test_dry_run_idempotent_and_full_blocked(self):
        c=Client(); p={'scope':'pilot','samples':[{'SampleID':'one'}], 'experiments':[{'fields':{'MigrationKey':'w:1'},'sample_id':'one'}]}
        self.assertEqual(apply_pilot(p,c)['samples'],['one']); self.assertFalse(c.tables['Samples'])
        apply_pilot(p,c,False); self.assertEqual(c.tables['Experiments'][0]['fields']['Sample'],1)
        self.assertEqual(apply_pilot(p,c,False),{'samples':[],'experiments':[]})
        with self.assertRaises(ValueError): apply_pilot({'scope':'full','samples':[],'experiments':[]},c,False)

    def test_new_import_omits_registration_time_without_mutating_plan(self):
        for scope in ('pilot','full'):
            with self.subTest(scope=scope):
                plan = {'scope':scope,'samples':[], 'experiments':[{'fields':{'MigrationKey':'workbook:1'},'sample_id':None}]}
                before = copy.deepcopy(plan); client = Client()
                client.tables['Unclassified'] = []
                apply_pilot(plan,client,dry_run=False,full_accepted=True)
                self.assertNotIn('RegisteredAt', client.tables['Experiments'][0]['fields'])
                self.assertEqual(plan, before)

    def test_dry_run_and_repeat_preserve_existing_records(self):
        plan = {'scope':'pilot','samples':[], 'experiments':[{'fields':{'MigrationKey':'workbook:1'},'sample_id':None}]}
        before = copy.deepcopy(plan); client = Client()
        apply_pilot(plan,client)
        self.assertEqual(client.tables['Experiments'], [])
        client.tables['Experiments'] = [{'id':1,'fields':{'MigrationKey':'workbook:1','Notes':'preserve'}}]
        existing = copy.deepcopy(client.tables)
        self.assertEqual(apply_pilot(plan,client,dry_run=False), {'samples':[],'experiments':[]})
        self.assertEqual(client.tables, existing)
        self.assertEqual(plan, before)

    def test_formula_assessment_blocks_migration_before_any_write(self):
        c=Client();c.tables['_grist_Tables']=[{'id':4,'fields':{'tableId':'Experiments'}}]
        c.tables['_grist_Tables_column']=[{'id':7,'fields':{'parentId':4,'colId':'Assessment','isFormula':True}}]
        p={'scope':'pilot','samples':[{'SampleID':'one'}],'experiments':[{'fields':{'MigrationKey':'w:1','Assessment':'Use'},'sample_id':'one'}]}
        with self.assertRaisesRegex(ValueError,'reviewed QC mapping'):apply_pilot(p,c,False)
        self.assertEqual(c.tables['Samples'],[]);self.assertEqual(c.tables['Experiments'],[])

    def test_existing_manual_identity_blocks_before_any_write(self):
        c=Client()
        c.tables['Experiments'].append({'id':1,'fields':{'SourceID':'main','Folder':'f','Expno':'1'}})
        plan={'scope':'pilot','samples':[{'SampleID':'one'}], 'experiments':[{'fields':{'MigrationKey':'workbook:1','SourceID':'main','Folder':'f','Expno':'1'},'sample_id':'one'}]}
        with self.assertRaisesRegex(ValueError, 'requires reconciliation'): apply_pilot(plan,c,False)
        self.assertFalse(c.tables['Samples'])

    def test_existing_migration_key_for_different_path_blocks(self):
        c=Client(); c.tables['Experiments'].append({'id':1,'fields':{'MigrationKey':'workbook:1','SourceID':'main','Folder':'old','Expno':'1'}})
        plan={'scope':'pilot','samples':[{'SampleID':'one'}], 'experiments':[{'fields':{'MigrationKey':'workbook:1','SourceID':'main','Folder':'new','Expno':'1'},'sample_id':'one'}]}
        with self.assertRaisesRegex(ValueError, 'different acquisition'): apply_pilot(plan,c,False)
        self.assertFalse(c.tables['Samples'])

class WorkbookTests(unittest.TestCase):
    def test_month_date_and_conflict_are_preserved(self):
        try: import openpyxl
        except ImportError: self.skipTest('openpyxl optional')
        import tempfile
        from pathlib import Path
        from nmr_register.migration import build_pilot
        with tempfile.TemporaryDirectory() as d:
            p=Path(d); w=openpyxl.Workbook(); w.active.title='Samples'
            w['Samples'].append(['position','id','protein','comment','unused','buffer','comment'])
            w['Samples'].append(['A1','202501_15N_P_50mMNaCl','P','','','50 mM NaCl',''])
            for name in ('Database','ANP32','N-Tail','TDR1','CtIP'):
                sh=w.create_sheet(name); sh.append(['date','instrument','user','folder','expno','sample','buffer','temp','type','pulse','options','status','comment'])
            w['TDR1'].append(['2025-01','600_1300','','folder',1,'202501_15N_P_50mMNaCl','',25,'R1','','','Failed','Poor SN'])
            x=p/'book.xlsx'; w.save(x)
            review=p/'reviews.tsv'; review.write_text('target\tassessment\tcomment\nfolder/1\tOK\tlooks fine\n')
            plan=build_pilot(x,review,['202501_15N_P_50mMNaCl'],['folder/1'],'instrument')
            f=plan['experiments'][0]['fields']
            self.assertEqual(f['AcquisitionDateText'],'2025-01'); self.assertIsNone(f['AcquiredAt'])
            self.assertEqual(f['Field_MHz'],600); self.assertEqual(f['Assessment'],'')
            self.assertIn('Conflict',f['MigrationIssue'])
