"""Cross-module acceptance checks using synthetic data and local fakes only."""
import csv
import hashlib
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest
from datetime import datetime
from unittest.mock import patch

from nmr_register.config import Source
from nmr_register.discovery import acquisition_key, scan_source
from nmr_register.schema import CONDITIONS, DECISIONS, columns, series_defaults, trigger_actions
from nmr_register.state import State
from nmr_register.sync import flush
from nmr_register.grist import GristError
from nmr_register.grist import GristClient
import requests
from nmr_register.backup import create_backup, restore_copy
from nmr_register.results import select_results
from nmr_register.migration import build_pilot


class SyntheticEntry:
    """A scandir entry that makes any inspection fail, for pruning checks."""
    name = '211@950'
    path = '/synthetic/root/211@950'

    def is_symlink(self):
        raise AssertionError('protected child was inspected')

    def is_dir(self, *, follow_symlinks=True):
        raise AssertionError('protected child was inspected')


class SyntheticScan:
    def __enter__(self):
        return [SyntheticEntry()]

    def __exit__(self, *args):
        return False


class CommitThenTimeoutClient:
    def __init__(self, rows=None):
        self.rows = rows or []
        self.timed_out = False

    def records(self, table):
        return self.rows

    def add_records(self, table, fields):
        self.rows.extend({'id': len(self.rows) + 1, 'fields': dict(row)} for row in fields)
        if not self.timed_out:
            self.timed_out = True
            raise GristError('synthetic timeout after commit')

    def update_records(self, table, records):
        for change in records:
            row = next(item for item in self.rows if item['id'] == change['id'])
            row['fields'].update(change['fields'])


class HttpResponse:
    status_code = 200
    headers = {}

    def __init__(self, value):
        self.value = value
        self.content = b''

    def json(self):
        return self.value


class CommitTimeoutSession:
    def __init__(self):
        self.headers = {}
        self.rows = []
        self.committed_then_timed_out = False

    def request(self, method, url, **kwargs):
        if method == 'GET' and url.endswith('/tables/Experiments/records'):
            return HttpResponse({'records': self.rows})
        if method == 'POST' and url.endswith('/tables/Experiments/records'):
            body = kwargs['json']
            self.rows.extend({'id': len(self.rows) + i + 1, 'fields': dict(r['fields'])}
                             for i, r in enumerate(body['records']))
            if not self.committed_then_timed_out:
                self.committed_then_timed_out = True
                raise requests.Timeout('synthetic response lost after server commit')
            return HttpResponse({'records': []})
        if method == 'PATCH' and url.endswith('/tables/Experiments/records'):
            for change in kwargs['json']['records']:
                row = next(r for r in self.rows if r['id'] == change['id'])
                row['fields'].update(change['fields'])
            return HttpResponse({'records': []})
        raise AssertionError((method, url))


class RegisterAcceptanceTests(unittest.TestCase):
    def test_protected_path_name_is_pruned_before_any_child_stat(self):
        with patch('nmr_register.discovery.os.scandir', return_value=SyntheticScan()):
            with patch.object(Path, 'is_dir', return_value=True):
                observations, errors = scan_source(Source('test', '/synthetic/root'))
        self.assertEqual(observations, [])
        self.assertEqual(errors, [])

    def test_protected_root_is_rejected_through_symlinked_parent(self):
        with tempfile.TemporaryDirectory() as td:
            protected = Path(td) / '211@950'
            root = protected / 'source'
            root.mkdir(parents=True)
            alias = Path(td) / 'alias'
            try:
                alias.symlink_to(protected, target_is_directory=True)
            except (OSError, NotImplementedError):
                self.skipTest('directory symlinks are unavailable')
            with patch('nmr_register.discovery.os.scandir') as scandir:
                observations, errors = scan_source(Source('protected', str(alias / 'source')))
            self.assertEqual(observations, [])
            self.assertTrue(errors)
            scandir.assert_not_called()

    def test_exact_acquisition_keys_are_stable_and_reject_traversal(self):
        self.assertEqual(acquisition_key('linux', '20261030/sample', '12'),
                         'linux:20261030/sample/12')
        for folder, expno in [('../sample', '12'), ('sample', '../12'), ('sample', '12/x')]:
            with self.subTest(folder=folder, expno=expno), self.assertRaises(ValueError):
                acquisition_key('linux', folder, expno)

    def test_repeated_observation_is_idempotent_and_copy_is_reviewed(self):
        with tempfile.TemporaryDirectory() as td:
            state = State(td)
            self.addCleanup(state.close)
            record = {'AcquisitionKey': 'linux:20261030/sample/12',
                      'ObservedFingerprint': 'stable-fingerprint',
                      'ObservedPulseProgram': 'hsqc', 'RawPresent': True,
                      'RawBytes': 128, 'ScanError': ''}
            state.observe([record], now=100)
            state.acknowledge([record['AcquisitionKey']])
            state.observe([record], now=200)
            state.acknowledge([record['AcquisitionKey']])
            state.observe([record], now=300)
            self.assertEqual(state.pending(), [])
            self.assertEqual(state.status()['observations'], 1)

            copied = dict(record, AcquisitionKey='linux:20261031/copy/12')
            state.observe([copied], now=400)
            pending = dict(state.pending())
            self.assertEqual(len(pending), 2)
            for fields in pending.values():
                self.assertIn('possible copy or rename', fields['ReconciliationReason'])

    def test_commit_timeout_retry_reconciles_key_without_duplicate(self):
        with tempfile.TemporaryDirectory() as td:
            state = State(td)
            self.addCleanup(state.close)
            record = {'AcquisitionKey': 'linux:20261030/sample/12',
                      'ObservedSourceID': 'linux', 'ObservedFolder': '20261030/sample',
                      'ObservedExpno': '12', 'ObservedFingerprint': 'fp',
                      'ObservedPulseProgram': 'hsqc', 'ObservedType': 'hsqc',
                      'RawPresent': True, 'RawBytes': 128, 'ScanError': '',
                      'ObservedInstrument': '700'}
            state.observe([record], now=100)
            state.acknowledge([record['AcquisitionKey']])
            state.observe([record], now=200)  # transfer becomes stable
            client = CommitThenTimeoutClient()
            with self.assertRaises(GristError):
                flush(state, client)
            self.assertEqual(len(client.rows), 1)
            self.assertEqual(len(state.pending()), 1)
            flush(state, client)
            self.assertEqual(len(client.rows), 1)
            self.assertEqual(state.pending(), [])

    def test_http_timeout_after_commit_is_reconciled_before_retrying_mutation(self):
        with tempfile.TemporaryDirectory() as td:
            state = State(td)
            self.addCleanup(state.close)
            record = {'AcquisitionKey': 'linux:20261030/sample/12',
                      'ObservedSourceID': 'linux', 'ObservedFolder': '20261030/sample',
                      'ObservedExpno': '12', 'ObservedFingerprint': 'fp',
                      'ObservedPulseProgram': 'hsqc', 'ObservedType': 'hsqc',
                      'RawPresent': True, 'RawBytes': 128, 'ScanError': ''}
            state.observe([record], now=100)
            state.acknowledge([record['AcquisitionKey']])
            state.observe([record], now=200)
            client = GristClient(type('Config', (), {'api_base': 'https://grist.test/api',
                'timeout': 0.1})(), 'synthetic-token', session=CommitTimeoutSession())
            with self.assertRaises(GristError):
                flush(state, client)
            self.assertEqual(len(client.session.rows), 1)
            flush(state, client)
            self.assertEqual(len(client.session.rows), 1)
            self.assertEqual(state.pending(), [])
            self.assertEqual(client.calls, 3)  # POST timeout, then GET reconciliation

    def test_sync_updates_machine_fields_and_preserves_human_fields(self):
        with tempfile.TemporaryDirectory() as td:
            state = State(td)
            self.addCleanup(state.close)
            key = 'linux:20261030/sample/12'
            record = {'AcquisitionKey': key, 'ObservedSourceID': 'linux',
                      'ObservedFolder': '20261030/sample', 'ObservedExpno': '12',
                      'ObservedFingerprint': 'fp', 'ObservedPulseProgram': 'hsqc',
                      'ObservedType': 'hsqc', 'RawPresent': True, 'RawBytes': 128,
                      'ScanError': '', 'ObservedInstrument': '700'}
            state.observe([record], now=100)
            row = {'id': 7, 'fields': {'SourceID': 'linux', 'Folder': '20261030/sample',
                'Expno': '12', 'AcquisitionKey': '', 'Sample': 4,
                'Assessment': 'Limited use', 'AssessmentReason': 'check linewidth',
                'Notes': 'operator note'}}
            class Client:
                def records(self, table): return [row]
                def add_records(self, table, fields): raise AssertionError('unexpected insert')
                def update_records(self, table, changes): row['fields'].update(changes[0]['fields'])
            flush(state, Client())
            state.observe([dict(record, RawBytes=256)], now=200)
            flush(state, Client())
            for field, value in [('Sample', 4), ('Assessment', 'Limited use'),
                                 ('AssessmentReason', 'check linewidth'), ('Notes', 'operator note')]:
                self.assertEqual(row['fields'][field], value)
            self.assertEqual(row['fields']['AcquisitionKey'], key)
            self.assertTrue(row['fields']['RawPresent'])
            self.assertEqual(row['fields']['RawBytes'], 256)

    def test_zero_conditions_and_snapshot_recalculation_policy(self):
        tables = columns()
        sample_cols = {column['id']: column for column in tables['Samples']}
        for field in ('NaCl_mM', 'DNA_uM', 'TFE_pct'):
            self.assertEqual(sample_cols[field]['type'], 'Numeric')
        self.assertNotIn('Temperature_C', sample_cols)
        exp_cols = {column['id']: column for column in tables['Experiments']}
        self.assertEqual(json.loads(exp_cols['Assessment']['widgetOptions'])['choices'], list(DECISIONS))
        self.assertEqual(series_defaults({'NaCl_mM': 150}, {'NaCl_mM': 0})['NaCl_mM'], 0)
        self.assertIn('PreparedDateText', sample_cols['SampleID']['formula'])
        column_ids = {('Experiments', 'Sample'): 20}
        column_ids[('Experiments', 'Series')] = 21
        snapshot_names = ['SnapshotSampleID'] + ['Snapshot' + name for name in CONDITIONS]
        column_ids.update({('Experiments', name): i + 30 for i, name in enumerate(snapshot_names)})
        sample_names = list(dict.fromkeys(('PreparedAt', 'PreparedDateText', 'Protein', 'Isotope',
            'DistinguishingCondition', 'SampleNumber', 'SampleID', *CONDITIONS)))
        column_ids[('Samples', 'Series')] = 25
        column_ids[('Samples', 'SampleID')] = 26
        column_ids[('Series', 'SeriesID')] = 27
        column_ids.update({('Samples', name): i + 60 for i, name in enumerate(sample_names)})
        actions = trigger_actions({'Experiments': 10}, column_ids)
        update_actions = [action for action in actions if action[0] == 'UpdateRecord']
        trigger_updates = [action for action in update_actions if 'recalcWhen' in action[3]]
        display_actions = [action for action in actions if action[0] == 'SetDisplayFormula']
        self.assertEqual(len(trigger_updates), len(snapshot_names) + 1 + len(CONDITIONS))
        self.assertEqual(len(update_actions), len(trigger_updates) + 3)
        self.assertEqual(len(display_actions), 3)
        self.assertEqual({(action[1], action[4]) for action in display_actions},
                         {('Samples', '$Series.SeriesID'),
                          ('Experiments', '$Series.SeriesID'),
                          ('Experiments', '$Sample.SampleID')})
        sample_id_action = next(action for action in actions
                                if action[2] == column_ids[('Samples', 'SampleID')])
        self.assertEqual(sample_id_action[3]['recalcDeps'],
                         ['L'] + [column_ids[('Samples', name)] for name in
                                  ('PreparedAt', 'PreparedDateText', 'Protein', 'Isotope',
                                   'DistinguishingCondition', 'SampleNumber')])
        for action in trigger_updates:
            self.assertEqual(action[0:2], ['UpdateRecord', '_grist_Tables_column'])
            self.assertEqual(action[3]['recalcWhen'], 0)
            if action[2] in {column_ids[('Experiments', name)] for name in snapshot_names}:
                self.assertEqual(action[3]['recalcDeps'], ['L', 20])
            elif action is not sample_id_action:
                self.assertEqual(action[3]['recalcDeps'], ['L', 25])
        self.assertEqual(exp_cols['SnapshotDNA_uM']['type'], 'Numeric')

    def test_manifest_selects_only_hash_verified_current_csv(self):
        with tempfile.TemporaryDirectory() as td:
            result = Path(td) / 'numbers.csv'
            content = ('assignment,rate\nA1,1.25\n'
                       '# source_folder: dated/sample\n# expno: 7\n').encode()
            result.write_bytes(content)
            manifest = Path(td) / 'current.tsv'
            with manifest.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=('acquisition', 'numbers_file', 'numbers_sha256'),
                                        delimiter='\t', lineterminator='\n')
                writer.writeheader()
                writer.writerow({'acquisition': 'dated/sample/7', 'numbers_file': str(result),
                                 'numbers_sha256': hashlib.sha256(content).hexdigest()})
            chosen = select_results(manifest, 'linux')
            self.assertEqual(len(chosen), 1)
            self.assertEqual(chosen[0].content, content)
            result.write_bytes(content + b'altered')
            with self.assertRaisesRegex(ValueError, 'SHA256'):
                select_results(manifest, 'linux')

    def test_manifest_cannot_read_selected_csv_through_protected_parent_symlink(self):
        with tempfile.TemporaryDirectory() as td:
            protected = Path(td) / '211@950'
            protected.mkdir()
            selected = protected / 'numbers.csv'
            content = ('assignment,rate\nA1,1.25\n# source_folder: sample\n# expno: 7\n').encode()
            selected.write_bytes(content)
            alias = Path(td) / 'alias'
            try:
                alias.symlink_to(protected, target_is_directory=True)
            except (OSError, NotImplementedError):
                self.skipTest('directory symlinks are unavailable')
            manifest = Path(td) / 'current.tsv'
            with manifest.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=('acquisition', 'numbers_file', 'numbers_sha256'),
                                        delimiter='\t', lineterminator='\n')
                writer.writeheader()
                writer.writerow({'acquisition': 'sample/7', 'numbers_file': str(alias / 'numbers.csv'),
                                 'numbers_sha256': hashlib.sha256(content).hexdigest()})
            with self.assertRaisesRegex(ValueError, 'protected|symlink'):
                select_results(manifest, 'linux')

    def test_sqlite_grist_backup_restores_internal_attachment_bytes(self):
        with tempfile.TemporaryDirectory() as td:
            source = Path(td) / 'source.grist'
            db = sqlite3.connect(source)
            for table in ('_grist_Tables', '_grist_Tables_column', 'Experiments', 'Samples', 'Series'):
                db.execute(f'CREATE TABLE "{table}" (id INTEGER, value TEXT)')
            db.execute('CREATE TABLE _grist_Attachments (id INTEGER, fileIdent TEXT, storageId TEXT)')
            db.execute('CREATE TABLE _gristsys_Files (ident TEXT, data BLOB)')
            attachment = b'CSV result bytes including data\n'
            db.execute('INSERT INTO _grist_Attachments VALUES (1, ?, NULL)', ('file-1',))
            db.execute('INSERT INTO _gristsys_Files VALUES (?, ?)', ('file-1', attachment))
            db.execute('INSERT INTO Experiments VALUES (1, "sample experiment")')
            db.commit()
            db.close()
            class Download:
                def download(self): return source.read_bytes()
            backup = create_backup(Download(), Path(td) / 'backups', now=datetime(2026, 10, 30))
            restored = Path(td) / 'restored.grist'
            restore_copy(backup['path'], restored)
            with sqlite3.connect(restored) as restored_db:
                self.assertEqual(restored_db.execute('SELECT data FROM _gristsys_Files').fetchone()[0], attachment)
                self.assertEqual(restored_db.execute('SELECT value FROM Experiments').fetchone()[0],
                                 'sample experiment')

    def test_pilot_keeps_month_precision_and_flags_review_conflict(self):
        try:
            import openpyxl
        except ImportError:
            self.skipTest('openpyxl is optional for migration only')
        with tempfile.TemporaryDirectory() as td:
            workbook = Path(td) / 'legacy.xlsx'
            book = openpyxl.Workbook()
            samples = book.active
            samples.title = 'Samples'
            samples.append(['Storage', 'SampleID', 'Protein', 'Other', 'Unused', 'Buffer', 'Notes'])
            samples.append(['A1', 'sample-2025-10', 'ProteinX', '', '', 'buffer', ''])
            for title in ('Database', 'ANP32', 'N-Tail', 'TDR1', 'CtIP'):
                sheet = book.create_sheet(title)
                sheet.append(['date','instrument','unused','folder','expno','sample','unused2',
                              'temperature','type','pulse','unused3','status','notes'])
            book['Database'].append(['2025-10', '700', '', 'dated/sample', '7',
                                     'sample-2025-10', '', 25, 'HSQC', 'hsqc', '',
                                     'Not for analysis', 'review this run'])
            book.save(workbook)
            reviews = Path(td) / 'reviews.tsv'
            reviews.write_text('target\tassessment\tcomment\ndated/sample/7\tusable\tok\n')
            plan = build_pilot(workbook, reviews, ['sample-2025-10'], ['dated/sample/7'], 'linux')
            fields = plan['experiments'][0]['fields']
            self.assertIsNone(fields['AcquiredAt'])
            self.assertEqual(fields['AcquisitionDateText'], '2025-10')
            self.assertIn('Conflict', fields['MigrationIssue'])
            self.assertEqual(fields['Assessment'], '')


if __name__ == '__main__':
    unittest.main()
