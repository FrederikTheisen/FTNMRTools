"""Routine assessment reminders follow explicit Study assignment only."""
import copy
from datetime import datetime, timedelta
import unittest
from nmr_register.schema import ATTENTION_FORMULA, PRE_STUDY_ASSESSMENT_ATTENTION_FORMULA, columns
from nmr_register.setup import upgrade_assessment_schema
from nmr_register.experiments_ui import upgrade_dummy_schema
from tests.test_assessment_schema import Fake, formula

class StudyAssessmentAttentionTests(unittest.TestCase):
    def setUp(self):
        now=datetime(2026,11,5,12)
        self.values=dict(Study='',AcquisitionStatus='Completed',EffectiveAcquisitionStatus='Completed',
                         EffectiveAcquiredAt=now-timedelta(hours=25),RegisteredAt=now,Sample=1,
                         ExperimentType='T1',Assessment='',AssessmentDisplay='Not assessed',now=now)
        self.source=ATTENTION_FORMULA.replace('NOW()', "values.get('now')")

    def test_only_nonblank_assigned_studies_enable_routine_reminder(self):
        for study in (None,'','   ','\t\n'):
            self.assertEqual(formula(self.source,dict(self.values,Study=study)), '')
        for study in ('Electrostatics','DNA','TFE','Other','New custom study','  DNA  '):
            for display in ('Not assessed','Use so far','Use so far · analysis missing','Limited use so far'):
                self.assertEqual(formula(self.source,dict(self.values,Study=study,AssessmentDisplay=display)),
                                 'Acquired over 24 hours ago; '+('unassessed' if display=='Not assessed' else display))
        for decision in ('Use','Limited use','Do not use'):
            self.assertEqual(formula(self.source,dict(self.values,Study='DNA',Assessment=decision,AssessmentExplanation='reviewed')), '')
        self.assertEqual(formula(self.source,dict(self.values,Study='DNA',EffectiveAcquiredAt=self.values['now']-timedelta(hours=24))), '')

    def test_other_problems_remain_independent_of_study_and_exemptions(self):
        for study in (None,'DNA'):
            base=dict(self.values,Study=study,Assessment='Use')
            for changes,wanted in [({'Sample':0},'Sample not assigned'),
                                   ({'RegisteredAt':None},'Registration incomplete'),
                                   ({'ReconciliationReason':'identity conflict'},'identity conflict'),
                                   ({'EffectiveAcquisitionStatus':'Aborted'},'Acquisition aborted')]:
                self.assertEqual(formula(self.source,dict(base,**changes)),wanted)
        self.assertEqual(formula(self.source,dict(self.values,Study='DNA',ExperimentCategory='Calibration')), '')
        self.assertEqual(formula(self.source,dict(self.values,Study='DNA',AcquisitionStatus='Dummy')), '')

    def test_previous_attention_upgrades_without_data_changes(self):
        for kind in ('assessment','dummy'):
            upgrade=upgrade_assessment_schema if kind=='assessment' else lambda c:upgrade_dummy_schema(c,apply=True)
            client=Fake([])
            client.cols=[{'id':i+1,'fields':{'parentId':4,'colId':s['id'],**{k:v for k,v in s.items() if k!='id'}}} for i,s in enumerate(columns()['Experiments'])]
            attention=next(c for c in client.cols if c['fields']['colId']=='AttentionReason')
            attention['fields']['formula']=PRE_STUDY_ASSESSMENT_ATTENTION_FORMULA
            client.rows=[{'id':9,'fields':dict(Study='',Assessment='Use',Notes='keep',RegisteredAt=1700000000)}]
            before=copy.deepcopy(client.rows)
            upgrade(client)
            expected=['UpdateRecord','_grist_Tables_column',attention['id'],{'formula':ATTENTION_FORMULA}] if kind=='assessment' else ['ModifyColumn','Experiments','AttentionReason',{'formula':ATTENTION_FORMULA}]
            self.assertEqual(client.actions,[[expected]])
            self.assertEqual(client.rows,before)
            upgrade(client);self.assertEqual(len(client.actions),1)
            attention['fields']['formula']='return "custom"'
            with self.assertRaisesRegex(ValueError,'requires review'):upgrade(client)
            self.assertEqual(len(client.actions),1)
