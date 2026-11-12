"""Provisional display never changes the canonical decision contract."""
import copy
import itertools
import unittest
from nmr_register.schema import (
    QC_CRITERIA, QC_DEFAULT, ASSESSMENT_FORMULA, ASSESSMENT_DISPLAY_FORMULA,
    ASSESSMENT_EXPLANATION_FORMULA, PRE_PARTIAL_ASSESSMENT_EXPLANATION_FORMULA,
    assessment_from_values, assessment_display_from_values, ATTENTION_FORMULA, PRE_PARTIAL_ATTENTION_FORMULA,
    PRE_COMPACT_ASSESSMENT_DISPLAY_FORMULA,
)
from nmr_register.setup import upgrade_assessment_schema
from tests.test_assessment_schema import Fake, formula


class PartialAssessmentTests(unittest.TestCase):
    def setUp(self):
        self.good = {name: choices[0] for name, (_, choices) in QC_CRITERIA.items()}

    def test_all_4096_valid_combinations(self):
        names = list(QC_CRITERIA)
        count = 0
        for choices in itertools.product(*[(QC_DEFAULT, *QC_CRITERIA[n][1]) for n in names]):
            values = dict(zip(names, choices))
            self.assertEqual(formula(ASSESSMENT_DISPLAY_FORMULA, values), assessment_display_from_values(values), values)
            self.assertEqual(formula(ASSESSMENT_FORMULA, values), assessment_from_values(values), values)
            missing = any(value == QC_DEFAULT for value in choices)
            veto = any(values[name] == options[-1] for name, (_, options) in QC_CRITERIA.items())
            display = assessment_display_from_values(values)
            if veto:
                self.assertEqual(display, 'Do not use', values)
            elif missing and any(value != QC_DEFAULT for value in choices):
                cautions = sum(values[name] == options[1] for name, (_, options) in QC_CRITERIA.items()
                               if name != 'QCAcquisitionCoverage')
                expected = ('Limited use' if cautions >= 2 else 'Use') + ' so far'
                if [name for name in names if values[name] == QC_DEFAULT] == ['QCAnalysis']:
                    expected += ' · analysis missing'
                self.assertEqual(display, expected, values)
            if missing and not veto:
                self.assertEqual(assessment_from_values(values), '')
            if missing:
                explanation = formula(ASSESSMENT_EXPLANATION_FORMULA, values)
                self.assertIn('Missing checks:', explanation)
                self.assertNotIn('All criteria acceptable', explanation)
            count += 1
        self.assertEqual(count, 4096)

    def test_analysis_missing_then_complete(self):
        values = dict(self.good, QCAnalysis=QC_DEFAULT)
        self.assertEqual(assessment_display_from_values(values), 'Use so far · analysis missing')
        self.assertEqual(assessment_from_values(values), '')
        values.update(QCSettings='Usable with restrictions', QCSNR='Marginal')
        self.assertEqual(assessment_display_from_values(values), 'Limited use so far · analysis missing')
        explanation = formula(ASSESSMENT_EXPLANATION_FORMULA, values)
        for phrase in ('Limited use so far', 'Settings: Usable with restrictions', 'SNR: Marginal', 'Missing checks: Analysis quality'):
            self.assertIn(phrase, explanation)
        values['QCFID'] = 'Unusable'
        self.assertEqual(assessment_display_from_values(values), 'Do not use')
        self.assertEqual(assessment_from_values(values), 'Do not use')
        self.assertIn('Do not use: FIDs and spectrum: Unusable', formula(ASSESSMENT_EXPLANATION_FORMULA, values))
        self.assertIn('SNR: Marginal', formula(ASSESSMENT_EXPLANATION_FORMULA, values))
        for decision, modifications in [('Use', {}), ('Limited use', {'QCSettings':'Usable with restrictions','QCSNR':'Marginal'}), ('Do not use', {'QCFID':'Unusable'})]:
            values = dict(self.good, **modifications)
            self.assertEqual(assessment_display_from_values(values), decision)

    def test_unknown_defaults_and_imported(self):
        for value in (None, '', QC_DEFAULT, 'Unknown'):
            values = {name:value for name in QC_CRITERIA}
            self.assertEqual(assessment_display_from_values(values), 'Not assessed')
            self.assertEqual(formula(ASSESSMENT_DISPLAY_FORMULA, values), 'Not assessed')
        for decision in ('Use','Limited use','Do not use'):
            for blank in (None, '', QC_DEFAULT):
                values = {name:blank for name in QC_CRITERIA}
                values['ImportedAssessment'] = decision
                self.assertEqual(assessment_display_from_values(values), decision)
                self.assertEqual(formula(ASSESSMENT_DISPLAY_FORMULA, values), decision)
            values['QCSettings'] = 'Unknown'
            self.assertEqual(assessment_display_from_values(values), 'Not assessed')
        values = dict(self.good, QCFID='Unknown')
        self.assertEqual(assessment_display_from_values(values), 'Use so far')
        values['QCAnalysis'] = None
        self.assertEqual(assessment_display_from_values(values), 'Use so far')
        values = dict(self.good, QCAcquisitionCoverage='Sufficient', QCAnalysis=None)
        self.assertEqual(assessment_display_from_values(values), 'Use so far · analysis missing')

    def test_stale_attention_describes_partial_progress_and_keeps_exemptions(self):
        from datetime import datetime, timedelta
        now = datetime(2026, 11, 4, 12)
        values = dict(self.good, Study='Electrostatics', AcquisitionStatus='Completed', EffectiveAcquisitionStatus='Completed',
                      ExperimentType='Relaxation', RegisteredAt=now, Sample=True,
                      EffectiveAcquiredAt=now-timedelta(hours=25), now=now, Assessment='')
        source = ATTENTION_FORMULA.replace('NOW()', "values.get('now')")
        for display in ('Use so far · analysis missing', 'Limited use so far'):
            values['AssessmentDisplay'] = display
            self.assertEqual(formula(source, values), 'Acquired over 24 hours ago; ' + display)
        values['AssessmentDisplay'] = 'Not assessed'
        self.assertEqual(formula(source, values), 'Acquired over 24 hours ago; unassessed')
        values.update(AssessmentDisplay='Use so far · analysis missing', ExperimentType='Calibration')
        self.assertEqual(formula(source, values), '')
        values.update(ExperimentType='Relaxation', AcquisitionStatus='Dummy')
        self.assertEqual(formula(source, values), '')
        values.update(AcquisitionStatus='Completed', EffectiveAcquiredAt=now-timedelta(hours=24))
        self.assertEqual(formula(source, values), '')

    def test_compact_labels_for_each_missing_check_and_every_veto(self):
        for name in QC_CRITERIA:
            values = dict(self.good, **{name:QC_DEFAULT})
            expected = 'Use so far' + (' · analysis missing' if name == 'QCAnalysis' else '')
            self.assertEqual(assessment_display_from_values(values), expected)
            self.assertEqual(formula(ASSESSMENT_DISPLAY_FORMULA, values), expected)
            self.assertIn(QC_CRITERIA[name][0], formula(ASSESSMENT_EXPLANATION_FORMULA, values))
        for name, (_, choices) in QC_CRITERIA.items():
            values = {name:choices[-1]}
            self.assertEqual(assessment_display_from_values(values), 'Do not use')
            self.assertEqual(formula(ASSESSMENT_DISPLAY_FORMULA, values), 'Do not use')
        values = dict(self.good, QCAnalysis=QC_DEFAULT, QCFID=QC_DEFAULT,
                      QCSettings='Usable with restrictions', QCSNR='Marginal')
        self.assertEqual(assessment_display_from_values(values), 'Limited use so far')
        self.assertEqual(formula(ASSESSMENT_DISPLAY_FORMULA, values), 'Limited use so far')

    def test_previous_display_formula_upgrades_once_without_data_writes(self):
        from nmr_register.schema import columns
        client = Fake([])
        client.cols = [{'id':i+1,'fields':{'parentId':4,'colId':spec['id'],**{k:v for k,v in spec.items() if k!='id'}}}
                       for i,spec in enumerate(columns()['Experiments'])]
        display = next(c for c in client.cols if c['fields']['colId']=='AssessmentDisplay')
        display['fields']['formula'] = PRE_COMPACT_ASSESSMENT_DISPLAY_FORMULA
        client.rows = [{'id':5,'fields':dict(self.good,QCAnalysis=QC_DEFAULT,AssessmentReason='keep note')}]
        before = copy.deepcopy(client.rows)
        self.assertEqual(upgrade_assessment_schema(client),
                         [['UpdateRecord','_grist_Tables_column',display['id'],{'formula':ASSESSMENT_DISPLAY_FORMULA}]])
        self.assertEqual(client.rows, before)
        self.assertEqual(upgrade_assessment_schema(client), [])
        display['fields']['formula'] = 'return "custom"'
        with self.assertRaisesRegex(ValueError, 'Existing QC schema differs'):
            upgrade_assessment_schema(client)
        self.assertEqual(len(client.actions), 1)

    def test_guarded_upgrade_adds_display_preserving_rows_and_custom_columns(self):
        client = Fake([])
        from nmr_register.schema import columns
        client.cols = [{'id':i+1,'fields':{'parentId':4,'colId':spec['id'],**{k:v for k,v in spec.items() if k!='id'}}}
                       for i,spec in enumerate(columns()['Experiments']) if spec['id'] != 'AssessmentDisplay']
        explanation = next(c for c in client.cols if c['fields']['colId']=='AssessmentExplanation')
        explanation['fields']['formula'] = PRE_PARTIAL_ASSESSMENT_EXPLANATION_FORMULA
        next(c for c in client.cols if c['fields']['colId']=='AttentionReason')['fields']['formula'] = PRE_PARTIAL_ATTENTION_FORMULA
        client.cols.append({'id':2000,'fields':{'parentId':4,'colId':'PersonalNote','type':'Text','label':'My notes'}})
        client.rows = [{'id':5,'fields':dict(self.good,QCAnalysis=QC_DEFAULT,ImportedAssessment='Do not use',AssessmentReason='Keep note',PersonalNote='untouched')}]
        before = copy.deepcopy(client.rows)
        actions = upgrade_assessment_schema(client)
        self.assertEqual(len(actions), 3)
        add = next(a for a in actions if a[0]=='AddColumn')
        self.assertEqual(add[2], 'AssessmentDisplay')
        self.assertEqual(add[3]['label'], 'Assessment')
        self.assertTrue(add[3]['untieColIdFromLabel'])
        self.assertEqual(client.rows, before)
        self.assertEqual(upgrade_assessment_schema(client), [])
        explanation['fields']['formula'] = 'return "custom"'
        with self.assertRaisesRegex(ValueError, 'Existing QC schema differs'):
            upgrade_assessment_schema(client)
        self.assertEqual(len(client.actions), 1)


if __name__ == '__main__':
    unittest.main()
