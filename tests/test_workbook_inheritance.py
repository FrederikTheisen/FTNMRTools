from datetime import datetime
import importlib.util
from pathlib import Path
import tempfile
import unittest
import openpyxl

SOURCE=Path(__file__).resolve().parents[1]/'nmr_register/migration.py'
spec=importlib.util.spec_from_file_location('carry_migration_under_test',SOURCE)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

HEADER=['Date','Instrument','User','Folder','EXPNO','Sample','Buffer','Temperature','Type','Pulse','Options','Status','Comment']
SAMPLE='20250101_15N_P_50mMNaCl'

class Sheet:
    def __init__(self,rows):self.values=rows
    def iter_rows(self,**kwargs):return iter(self.values)

class WorkbookInheritanceTests(unittest.TestCase):
    def test_carries_only_authorized_columns_and_keeps_dates_exact(self):
        date=datetime(2025,1,1)
        rows=list(module.experiment_rows(Sheet([[date,'700','user','folder',1,SAMPLE,'50 mM NaCl',25,'T1','pulse','','Done','one'],[None,None,None,None,2,None,None,None,None,None,None,None,None]])))
        second=rows[1][1]
        self.assertEqual([second[i] for i in (0,1,2,3,5)],[date,'700','user','folder',SAMPLE])
        self.assertEqual(second[4],2)
        self.assertTrue(all(second[i] is None for i in (6,7,8,9,10,11,12)))
        self.assertEqual(rows[1][2],{0:2,1:2,2:2,3:2,5:2})
    def test_fully_blank_separator_resets_all_context(self):
        data=[['2025-01','700','user','first',1,SAMPLE],[' ',None,'',None,None,None,None,None,None,None,None,None,None,'outside A:M'],[None,None,None,None,2,None,None,None,'T1']]
        rows=list(module.experiment_rows(Sheet(data)))
        self.assertEqual([r[0] for r in rows],[2,4])
        self.assertTrue(all(rows[1][1][i] is None for i in (0,1,2,3,5)))
        self.assertEqual(rows[1][2],{})
    def test_explicit_changes_override_each_column_without_carrying_expno(self):
        rows=list(module.experiment_rows(Sheet([[None,'700','one','first',1,SAMPLE],[None,'850',None,'second',None,'explicit-other-sample'],[None,None,'two',None,3,None]])))
        self.assertIsNone(rows[1][1][4])
        self.assertEqual([rows[2][1][i] for i in (1,2,3,5)],['850','two','second','explicit-other-sample'])
    def workbook(self,root,rows):
        book=openpyxl.Workbook();book.active.title='Samples'
        book['Samples'].append(['position','id','protein','comment','unused','buffer','comment'])
        book['Samples'].append(['A1',SAMPLE,'P','','','50 mM NaCl',''])
        for name in ('Database','ANP32','N-Tail','TDR1','CtIP'):
            sheet=book.create_sheet(name);sheet.append(HEADER)
        for row in rows:book['TDR1'].append(row)
        path=root/'workbook.xlsx';book.save(path);book.close()
        review=root/'reviews.tsv';review.write_text('target\tassessment\tcomment\n')
        return path,review
    def test_pilot_and_full_share_carry_rule_without_inventing_sample(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            path,review=self.workbook(root,[['2025-01','700','person','folder',1,SAMPLE,'50 mM NaCl',25,'T1'],[None,None,None,None,2,None,None,None,'T1rho'],[None,None,None,None,3,'unknown exact label',None,None,'HetNOE']])
            pilot=module.build_pilot(path,review,[SAMPLE],['folder/2','folder/3'],'src')
            first,unknown=pilot['experiments']
            self.assertEqual(first['fields']['Instrument'],'700');self.assertEqual(first['fields']['AcquisitionDateText'],'')
            self.assertIn('group date context: 2025-01',first['fields']['SourceEvidence'])
            self.assertIsNone(first['fields']['AcquiredAt']);self.assertEqual(first['sample_id'],SAMPLE)
            self.assertIsNone(unknown['sample_id']);self.assertIn('unknown exact label',unknown['fields']['MigrationIssue'])
            self.assertIn('Folder from D2',first['fields']['SourceEvidence'])
            full=module.build_full(path,review,'src')
            self.assertEqual({(e['fields']['Folder'],e['fields']['Expno']) for e in full['experiments']},{('folder','1'),('folder','2'),('folder','3')})
    def test_conflicting_duplicate_target_keeps_both_rows_as_evidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            path,review=self.workbook(root,[['2025-01','700','person','folder',1,SAMPLE,None,None,'T1',None,None,'Done','first note'],[None,None,None,None,1,None,None,None,'T1rho',None,None,'Failed','second note']])
            plan=module.build_pilot(path,review,[SAMPLE],['folder/1'],'src')
            self.assertEqual(len(plan['experiments']),1)
            fields=plan['experiments'][0]['fields']
            self.assertEqual(fields['ExperimentType'],'T1')
            self.assertIn('duplicate workbook target',fields['MigrationIssue'])
            self.assertIn('TDR1!3',fields['SourceEvidence']);self.assertIn('second note',fields['SourceEvidence'])
    def test_blank_expno_stays_incomplete_after_inheritance(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            path,review=self.workbook(root,[[None,'700','person','folder',1,SAMPLE,None,None,'T1'],[None,None,None,None,None,None,None,None,'T1rho']])
            plan=module.build_full(path,review,'src')
            partial=next(e for e in plan['experiments'] if not e['fields']['Expno'])
            self.assertEqual(partial['fields']['Folder'],'folder');self.assertEqual(partial['fields']['Instrument'],'700')
            self.assertEqual(partial['fields']['Field_MHz'],700)
            self.assertNotIn('AcquisitionKey',partial['fields']);self.assertEqual(partial['sample_id'],SAMPLE)

    def test_inherited_date_is_context_only_explicit_dates_still_populate(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            path,review=self.workbook(root,[[datetime(2025,1,1),'700','person','folder',1,SAMPLE,None,None,'T1'],[None,None,None,None,2,None,None,None,'T1rho'],[datetime(2025,1,3),None,None,None,3,None,None,None,'HetNOE'],[None,None,None,None,None,None,None,None,'T1']])
            full=module.build_full(path,review,'src')
            rows={e['fields']['Expno']:e['fields'] for e in full['experiments']}
            self.assertIsNotNone(rows['1']['AcquiredAt']);self.assertIsNotNone(rows['3']['AcquiredAt'])
            self.assertNotEqual(rows['1']['AcquiredAt'],rows['3']['AcquiredAt'])
            for expno in ('2',''):
                self.assertIsNone(rows[expno]['AcquiredAt']);self.assertEqual(rows[expno]['AcquisitionDateText'],'')
                self.assertIn('group date context:',rows[expno]['SourceEvidence'])
            self.assertIn('2025-01-01',rows['2']['SourceEvidence']);self.assertIn('2025-01-03',rows['']['SourceEvidence'])
            self.assertEqual(rows['2']['Instrument'],'700');self.assertEqual(rows['2']['Folder'],'folder')

if __name__=='__main__':unittest.main()
