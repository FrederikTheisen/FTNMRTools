from copy import deepcopy
import tempfile
from pathlib import Path
import unittest
from nmr_register.reconciliation import build_plan, apply_plan, ReconciliationError


class ReconciliationTests(unittest.TestCase):
    def make(self, old=None, target=None, types=None):
        a = {'AcquisitionKey': '', 'Notes': 'human', 'QCSNR': 0,
             'ExperimentType': 'R1', 'MigrationIssue': 'date conflict'}
        b = {'AcquisitionKey': 'disk:a/2', 'Notes': '', 'QCSNR': 'Not assessed',
             'ExperimentType': 'T1', 'MigrationIssue': ''}
        a.update(old or {})
        b.update(target or {})
        types = types or {}
        data = {'Experiments': [{'id': 10, 'fields': a}, {'id': 20, 'fields': b}],
                'Results': [{'id': 1, 'fields': {'Exp': 10, 'Many': ['L', 10, 20]}}]}
        cols = [{'table': 'Experiments', 'colId': k, 'type': types.get(k, 'Text')}
                for k in a.keys() | b.keys() | {'SourceEvidence'}]
        cols += [{'table': 'Results', 'colId': 'Exp', 'type': 'Ref:Experiments'},
                 {'table': 'Results', 'colId': 'Many', 'type': 'RefList:Experiments'}]
        return data, cols

    def plan(self, *args, **kwargs):
        return build_plan(*self.make(*args, **kwargs), 10, 20, 'disk:a/2')

    def test_preservation_alias_zero_qc_and_reference_relinking(self):
        plan = self.plan()
        row = plan.after['Experiments'][0]['fields']
        self.assertEqual(row['QCSNR'], 0)
        self.assertEqual(row['Notes'], 'human')
        self.assertEqual(row['MigrationIssue'], 'date conflict')
        self.assertEqual(row['ExperimentType'], 'T1')
        self.assertIn('Reviewed alias provenance', row['SourceEvidence'])
        self.assertEqual(plan.after['Results'][0]['fields'], {'Exp': 20, 'Many': ['L', 20]})
        self.assertEqual(plan.actions[-1], ['RemoveRecord', 'Experiments', 10])

    def test_decisions_conflict(self):
        for name in ('ImportedAssessment', 'AnalysisStatus', 'Sample', 'SnapshotProtein', 'QCSNR'):
            with self.subTest(name=name), self.assertRaises(ReconciliationError):
                self.plan({name: 'first'}, {name: 'second'})

    def test_unset_references_and_real_numeric_zero(self):
        row = self.plan({'Sample': 7, 'Series': 8, 'SnapshotDNA_uM': 0},
                        {'Sample': 0, 'Series': 0, 'SnapshotDNA_uM': None},
                        {'Sample': 'Ref:Samples', 'Series': 'Ref:Series'}).after['Experiments'][0]['fields']
        self.assertEqual((row['Sample'], row['Series'], row['SnapshotDNA_uM']), (7, 8, 0))

    def test_notes_reasons_dates_and_observed(self):
        row = self.plan({'ImportedAssessment': 'Good', 'AssessmentReason': 'one',
                         'AcquiredAt': 1, 'RegisteredAt': 2, 'ObservedExtra': 'old observation',
                         'SourceEvidence': 'workbook'},
                        {'ImportedAssessment': 'Good', 'AssessmentReason': 'two',
                         'AcquiredAt': 3, 'RegisteredAt': 4, 'ObservedExtra': None,
                         'SourceEvidence': 'disk'}).after['Experiments'][0]['fields']
        self.assertEqual(row['AssessmentReason'], 'two\none')
        self.assertEqual((row['AcquiredAt'], row['RegisteredAt']), (3, 4))
        self.assertEqual(row['ObservedExtra'], 'old observation')
        for text in ('workbook', 'disk', 'legacy AcquiredAt: 1', 'legacy RegisteredAt: 2'):
            self.assertIn(text, row['SourceEvidence'])

    def test_selected_results_conflict_and_source_only_bundle(self):
        with self.assertRaises(ReconciliationError):
            self.plan({'ResultCSV': ['L', 1]}, {'ResultCSV': ['L', 2]}, {'ResultCSV': 'Attachments'})
        row = self.plan({'ResultCSV': ['L', 1], 'ResultSHA256': 'hash', 'ResultSourcePath': 'path'},
                        {'ResultCSV': ['L'], 'ResultSHA256': '', 'ResultSourcePath': ''},
                        {'ResultCSV': 'Attachments'}).after['Experiments'][0]['fields']
        self.assertEqual((row['ResultCSV'], row['ResultSHA256'], row['ResultSourcePath']), (['L', 1], 'hash', 'path'))

    def test_collector_facts_primary_and_sort_position_ignored(self):
        row = self.plan({'AcquisitionStatus': 'Completed', 'Temperature_C': 25,
                         'D1_s': 2, 'Title': 'workbook title', 'manualSort': 1},
                        {'AcquisitionStatus': 'Unknown', 'Temperature_C': 24.948,
                         'D1_s': None, 'Title': 'raw title', 'manualSort': 99},
                        {'manualSort': 'ManualSortPos'}).after['Experiments'][0]['fields']
        self.assertEqual(row['AcquisitionStatus'], 'Unknown')
        self.assertEqual(row['Temperature_C'], 24.948)
        self.assertIsNone(row['D1_s'])
        self.assertEqual(row['Title'], 'raw title')
        self.assertEqual(row['manualSort'], 99)
        for text in ('legacy AcquisitionStatus: Completed', 'legacy Temperature_C: 25',
                     'legacy D1_s: 2', 'legacy Title: workbook title'):
            self.assertIn(text, row['SourceEvidence'])

    def test_computed_effective_decision_disagreement_blocks(self):
        for target_decision in ('Use', 'Do not use'):
            with self.subTest(target_decision=target_decision), self.assertRaisesRegex(ReconciliationError, 'effective Assessment'):
                self.plan({'Assessment': 'Limited use'}, {'Assessment': target_decision})

    def test_partial_checklist_cannot_mask_imported_decision(self):
        with self.assertRaisesRegex(ReconciliationError, 'mask or change'):
            self.plan({'ImportedAssessment': 'Limited use', 'QCSNR': 'Not assessed'},
                      {'ImportedAssessment': '', 'QCSNR': 'Adequate'})

    def test_only_obsolete_generated_warnings_removed(self):
        real = 'DNA differs; salt conflict; date conflict; type conflict; assessment conflict; fingerprint conflict; copy conflict'
        obsolete = ('Data not found in scanned roots; recorded from workbook; '
                    'Inherited exact workbook path; data not found in scanned roots; '
                    'Incomplete acquisition identity; source details retained; '
                    'Planned or incomplete record; identity and sample unresolved; '
                    'Review-only record; sample and registration unresolved; '
                    'Potential matching row 123 has unresolved acquisition identity; '
                    'Sample link unresolved: no exact sample; workbook condition: salt 50; '
                    'Assessment reason missing')
        row = self.plan({'MigrationIssue': obsolete + '; ' + real, 'Sample': 4,
                         'AssessmentReason': 'reviewed'}, {'MigrationIssue': ''}).after['Experiments'][0]['fields']
        self.assertEqual(row['MigrationIssue'], real)
        row = self.plan({'MigrationIssue': 'Sample link unresolved: no exact sample; Assessment reason missing'},
                        {'MigrationIssue': ''}).after['Experiments'][0]['fields']
        self.assertEqual(row['MigrationIssue'], 'Sample link unresolved: no exact sample; Assessment reason missing')
        row = self.plan({'MigrationIssue': 'Actual warning: Data not found in scanned roots; recorded from workbook'},
                        {'MigrationIssue': ''}).after['Experiments'][0]['fields']
        self.assertEqual(row['MigrationIssue'], 'Actual warning: Data not found in scanned roots; recorded from workbook')

    def test_type_disagreement(self):
        with self.assertRaises(ReconciliationError):
            self.plan(target={'ExperimentType': 'T2'})

    def test_preflight_and_timeout_after_atomic_commit(self):
        data, cols = self.make()
        plan = build_plan(data, cols, 10, 20, 'disk:a/2')
        state = deepcopy(data)
        class Client:
            calls = 0
            def apply(self, actions):
                self.calls += 1
                state.clear()
                state.update(deepcopy(plan.after))
                raise TimeoutError('response lost after commit')
        client = Client()
        with tempfile.TemporaryDirectory() as tmp:
            backup = Path(tmp) / 'backup.grist'
            backup.write_bytes(b'caller backup')
            def apply():
                return apply_plan(client, plan, read_snapshot=lambda: state,
                                  read_columns=lambda: cols, backup_path=backup, exclusive_write=True)
            state['Experiments'][0]['fields']['Notes'] = 'new edit'
            with self.assertRaises(ReconciliationError):
                apply()
            self.assertEqual(client.calls, 0)
            state.clear()
            state.update(deepcopy(data))
            with self.assertRaises(TimeoutError):
                apply()
            self.assertEqual(apply(), 'already_applied')
            self.assertEqual(client.calls, 1)


if __name__ == '__main__':
    unittest.main()
