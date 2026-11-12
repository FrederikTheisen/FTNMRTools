"""Assigning a sample immediately removes only obsolete unresolved-link alerts."""
import copy
from datetime import datetime,timedelta
import unittest
from nmr_register.schema import ATTENTION_FORMULA,PRE_LIVE_SAMPLE_ATTENTION_FORMULA,PRE_ASSIGNED_SAMPLE_ATTENTION_FORMULA,columns
from nmr_register.setup import upgrade_assessment_schema
from nmr_register.experiments_ui import upgrade_dummy_schema
from tests.test_assessment_schema import Fake,formula

class LiveSampleAttentionTests(unittest.TestCase):
    def setUp(self):
        self.base=dict(RegisteredAt=1700000000,Sample=9,Assessment='Use',Study='',AcquisitionStatus='Completed',EffectiveAcquisitionStatus='Completed')

    def test_sample_assignment_and_unlink_recompute_without_changing_stored_issue(self):
        values=dict(self.base,Sample=0,MigrationIssue='Sample link unresolved: 15N_SeV_Ntail; workbook condition: 250 mM NaCl')
        before=copy.deepcopy(values)
        self.assertEqual(formula(ATTENTION_FORMULA,values),'Sample not assigned')
        values['Sample']=9
        self.assertEqual(formula(ATTENTION_FORMULA,values),'')
        values['Sample']=0
        self.assertEqual(values,before)
        self.assertEqual(formula(ATTENTION_FORMULA,values),'Sample not assigned')

    def test_clause_boundaries_empty_label_and_adjacent_condition(self):
        for separator in ('; ','\n','\r\n',';\n;; '):
            for label in ('tube',''):
                unresolved='Sample link unresolved: '+label+separator+'workbook condition: 250 mM NaCl'
                for text,wanted in [(unresolved,''),(unresolved+separator+'identity conflict','identity conflict'),
                                    ('identity conflict'+separator+unresolved,'identity conflict'),
                                    ('first'+separator+unresolved+separator+'last','first; last')]:
                    self.assertEqual(formula(ATTENTION_FORMULA,dict(self.base,MigrationIssue=text)),wanted)
                    self.assertEqual(formula(ATTENTION_FORMULA,dict(self.base,Sample=0,MigrationIssue=text)),(wanted+'; ' if wanted else '')+'Sample not assigned')
        text='workbook condition: keep; Sample link unresolved: old; genuine warning; workbook condition: keep too'
        self.assertEqual(formula(ATTENTION_FORMULA,dict(self.base,MigrationIssue=text)),
                         'workbook condition: keep; genuine warning; workbook condition: keep too')
        for text in ('Sample link conflict: salt differs','See Sample link unresolved: historical context','workbook condition: alone'):
            self.assertEqual(formula(ATTENTION_FORMULA,dict(self.base,MigrationIssue=text)),text)

    def test_other_warnings_and_study_exemptions_are_unchanged(self):
        issue='Sample link unresolved: old; workbook condition: salt; Sample link conflict: salt differs; Assessment reason missing'
        values=dict(self.base,MigrationIssue=issue,SampleTagIssue='Title/sample mismatch',ReconciliationReason='identity conflict',EffectiveAcquisitionStatus='Aborted')
        expected='Sample link conflict: salt differs; Assessment reason missing; Title/sample mismatch; identity conflict; Acquisition aborted'
        self.assertEqual(formula(ATTENTION_FORMULA,values),expected)
        now=datetime(2026,11,6)
        source=ATTENTION_FORMULA.replace('NOW()',"values.get('now')")
        values=dict(self.base,MigrationIssue='Sample link unresolved: old',Assessment='',Study='DNA',AssessmentDisplay='Use so far',EffectiveAcquiredAt=now-timedelta(days=2),now=now)
        self.assertEqual(formula(source,values),'Acquired over 24 hours ago; Use so far')
        self.assertEqual(formula(source,dict(values,Study='')),'')
        self.assertEqual(formula(source,dict(values,ExperimentType='Calibration')),'')
        self.assertEqual(formula(source,dict(values,AcquisitionStatus='Dummy',SampleTagIssue='mismatch')),'')

    def test_previous_formula_guard_upgrades_once_without_data_writes(self):
        for kind, previous in ((kind, previous) for kind in ('assessment','dummy') for previous in (PRE_LIVE_SAMPLE_ATTENTION_FORMULA,PRE_ASSIGNED_SAMPLE_ATTENTION_FORMULA)):
            upgrade=upgrade_assessment_schema if kind=='assessment' else lambda c:upgrade_dummy_schema(c,apply=True)
            client=Fake([])
            client.cols=[{'id':i+1,'fields':{'parentId':4,'colId':s['id'],**{k:v for k,v in s.items() if k!='id'}}} for i,s in enumerate(columns()['Experiments'])]
            attention=next(c for c in client.cols if c['fields']['colId']=='AttentionReason')
            attention['fields']['formula']=previous
            client.rows=[{'id':60,'fields':dict(self.base,MigrationIssue='Sample link unresolved: old',Notes='keep')}]
            before=copy.deepcopy(client.rows);upgrade(client)
            self.assertEqual(len(client.actions),1);self.assertEqual(len(client.actions[0]),1)
            self.assertEqual(client.rows,before)
            upgrade(client);self.assertEqual(len(client.actions),1)
            attention['fields']['formula']='return "custom"'
            with self.assertRaisesRegex(ValueError,'requires review'):upgrade(client)
            self.assertEqual(len(client.actions),1)
