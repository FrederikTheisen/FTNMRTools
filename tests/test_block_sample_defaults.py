"""Sample block defaults are explicit, scoped, and fail closed."""
import copy
import json
import tempfile
from types import SimpleNamespace
import unittest
from pathlib import Path
from unittest.mock import patch

from nmr_register.config import Source
from nmr_register.discovery import scan_source, resolve_folder_sample_defaults
from nmr_register.state import State


class BlockSampleDefaultTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name) / 'source'
        self.root.mkdir()
        self.source = Source('source', str(self.root))

    def acquisition(self, expno, title=None, folder='folder', root=None):
        path = (root or self.root) / folder / str(expno)
        path.mkdir(parents=True)
        (path / 'acqus').write_text('##$DATE_START= ' + str(1700000000 + int(expno)) + '\n##$PULPROG= <hsqc>\n')
        if title is not None:
            self.title(path, title)
        return path

    def title(self, path, content, procno=1):
        title = path / 'pdata' / str(procno) / 'title'
        title.parent.mkdir(parents=True, exist_ok=True)
        title.write_bytes(content.encode() if isinstance(content, str) else content)
        return title

    def scan(self, source=None):
        records, errors = scan_source(source or self.source)
        self.assertEqual(errors, [])
        return records

    def by_expno(self):
        return {r['ObservedExpno']: r for r in self.scan()}

    def test_boundaries_leading_zero_and_explicit_low_block(self):
        self.acquisition('099', 'SAMPLE_ID_BLOCK=low\n')
        self.acquisition('000')
        holder = self.acquisition('105', 'comments\nSAMPLE_ID_BLOCK=middle\n')
        self.acquisition('100')
        self.acquisition('199')
        self.acquisition('200')
        rows = self.by_expno()
        for exp in ('099', '000'):
            self.assertEqual(rows[exp]['ObservedSampleID'], 'low')
        for exp in ('105', '100', '199'):
            self.assertEqual(rows[exp]['ObservedSampleID'], 'middle')
        self.assertEqual(rows['200']['ObservedSampleID'], '')
        self.assertEqual(rows['105']['Title'], 'comments\n')
        self.assertEqual((holder / 'pdata/1/title').read_text(), 'comments\nSAMPLE_ID_BLOCK=middle\n')
        self.assertIn('105', rows['100']['ObservedSampleTagEvidence'])

    def test_exact_parent_source_and_allowlist_limit_scope(self):
        holder = self.acquisition(105, 'SAMPLE_ID_BLOCK=A\n')
        self.acquisition(106)
        self.acquisition(107, folder='folder/child')
        self.acquisition(108, folder='other')
        other = self.root.parent / 'other-source'; other.mkdir()
        self.acquisition(109, root=other)
        rows = self.by_expno()
        self.assertEqual(rows['106']['ObservedSampleID'], 'A')
        self.assertEqual(rows['107']['ObservedSampleID'], '')
        self.assertEqual(rows['108']['ObservedSampleID'], '')
        self.assertEqual(self.scan(Source('other', str(other)))[0]['ObservedSampleID'], '')
        self.source.targets = ['folder/106']
        from nmr_register import bruker
        original = bruker._read_metadata
        def guarded(path):
            self.assertFalse(Path(path).is_relative_to(holder))
            return original(path)
        with patch.object(bruker, '_read_metadata', side_effect=guarded):
            self.assertEqual(self.scan()[0]['ObservedSampleID'], '')

    def test_direct_conflict_blocks_inheritance_but_preserves_direct_ids(self):
        self.acquisition(105, 'SAMPLE_ID_BLOCK=A\n')
        direct = self.acquisition(106, 'SAMPLE_ID=B\n')
        self.acquisition(107)
        rows = self.by_expno()
        self.assertEqual(rows['106']['ObservedSampleID'], 'B')
        for exp in ('105', '107'):
            self.assertEqual(rows[exp]['ObservedSampleID'], '')
            self.assertTrue(rows[exp]['ObservedSampleTagIssue'])
        self.title(direct, 'SAMPLE_ID=\n')
        rows = self.by_expno()
        self.assertTrue(all(not r['ObservedSampleID'] for r in rows.values()))
        self.assertTrue(rows['107']['ObservedSampleTagIssue'])
        self.title(direct, 'SAMPLE_ID=A\n')
        self.assertTrue(all(r['ObservedSampleID'] == 'A' for r in self.scan()))

    def test_repeated_conflicting_and_malformed_declarations(self):
        self.acquisition(105, 'SAMPLE_ID_BLOCK=A\nSAMPLE_ID_BLOCK=A\n')
        second = self.acquisition(106, 'SAMPLE_ID_BLOCK=A\n')
        self.acquisition(107)
        self.assertTrue(all(r['ObservedSampleID'] == 'A' for r in self.scan()))
        for text in ('SAMPLE_ID_BLOCK=B\n', 'SAMPLE_ID_BLOCK=\n', 'SAMPLE_ID_BLOCK:bad\n'):
            with self.subTest(text=text):
                self.title(second, text)
                rows = self.scan()
                self.assertTrue(all(not r['ObservedSampleID'] for r in rows))
                self.assertTrue(all(r['ObservedSampleTagIssue'] for r in rows))

    def test_missing_title_is_eligible_unsafe_or_unreadable_blocks_peers(self):
        self.acquisition(105, 'SAMPLE_ID_BLOCK=A\n')
        bad = self.acquisition(106)
        self.acquisition(107)
        self.assertTrue(all(r['ObservedSampleID'] == 'A' for r in self.scan()))
        title = self.title(bad, b'invalid\xff')
        self.assertTrue(self.by_expno()['107']['ObservedSampleTagIssue'])
        title.unlink()
        outside = self.root.parent / 'outside'; outside.write_text('ordinary\n')
        title.symlink_to(outside)
        self.assertTrue(self.by_expno()['107']['ObservedSampleTagIssue'])
        title.unlink(); self.title(bad, 'ordinary\n')
        from nmr_register import bruker
        original = bruker._read_metadata
        def denied(path):
            if Path(path) == title:
                raise PermissionError('synthetic denied read')
            return original(path)
        with patch.object(bruker, '_read_metadata', side_effect=denied):
            self.assertTrue(self.by_expno()['107']['ObservedSampleTagIssue'])

    def test_folder_default_block_priority_and_broken_block(self):
        self.acquisition(1, 'SAMPLE_ID_FOLDER=folder\n')
        self.acquisition(105, 'SAMPLE_ID_BLOCK=block\n')
        self.acquisition(106)
        self.acquisition(200)
        rows = self.by_expno()
        self.assertEqual(rows['106']['ObservedSampleID'], 'block')
        self.assertEqual(rows['200']['ObservedSampleID'], 'folder')
        self.acquisition(107, 'SAMPLE_ID_BLOCK=other\n')
        rows = self.by_expno()
        self.assertEqual(rows['106']['ObservedSampleID'], '')
        self.assertTrue(rows['106']['ObservedSampleTagIssue'])
        self.assertEqual(rows['200']['ObservedSampleID'], 'folder')

    def test_ignored_processing_and_no_opt_in(self):
        holder = self.acquisition(105, 'ordinary SAMPLE_ID_BLOCK=A\n')
        self.title(holder, 'SAMPLE_ID_BLOCK=A\n', 2)
        self.acquisition(106, 'SAMPLE_ID=local\n')
        self.acquisition(107)
        rows = self.by_expno()
        self.assertEqual(rows['105']['ObservedSampleID'], '')
        self.assertEqual(rows['106']['ObservedSampleID'], 'local')
        self.assertEqual(rows['107']['ObservedSampleID'], '')

    def test_changes_and_removal_requeue_peers_without_identity_changes(self):
        holder = self.acquisition(105, 'SAMPLE_ID_BLOCK=A\n')
        self.acquisition(106)
        state = State(self.root.parent / 'state'); self.addCleanup(state.close)
        first = self.scan()
        fingerprints = {r['AcquisitionKey']: r['ObservedFingerprint'] for r in first}
        state.observe(first, now=1); state.observe(self.scan(), now=2)
        state.acknowledge(list(fingerprints))
        self.title(holder, 'SAMPLE_ID_BLOCK=B\n')
        state.observe(self.scan(), now=3)
        self.assertEqual(len(state.pending()), 2)
        for key, fields in state.pending():
            self.assertEqual(fields['ObservedSampleID'], 'B')
            self.assertEqual(fields['ObservedFingerprint'], fingerprints[key])
            self.assertNotIn('Sample', fields)
            self.assertEqual(fields['ReconciliationReason'], '')
            self.assertFalse(any(k.startswith('_Block') for k in fields))
        state.acknowledge(list(fingerprints))
        self.title(holder, 'ordinary\n'); state.observe(self.scan(), now=4)
        self.assertEqual(len(state.pending()), 2)
        for key, fields in state.pending():
            self.assertFalse(fields['ObservedSampleID'])
            self.assertEqual(fields['ObservedFingerprint'], fingerprints[key])
            self.assertEqual(fields['ReconciliationReason'], '')
            self.assertFalse(any(k.startswith('_Block') for k in fields))
        state.observe(self.scan(), now=5)  # allow transfer-stability observation to settle
        state.acknowledge(list(fingerprints)); state.observe(self.scan(), now=6)
        self.assertFalse(state.pending())

    def test_repeat_resolution_and_merge_of_independently_resolved_scans(self):
        self.acquisition(1, 'SAMPLE_ID_FOLDER=folder\n')
        self.acquisition(105, 'SAMPLE_ID_BLOCK=A\n')
        self.acquisition(106)
        self.acquisition(107, 'SAMPLE_ID_BLOCK=B\n')
        self.source.targets = ['folder/1', 'folder/105', 'folder/106']
        first = self.scan()
        before = json.dumps(first, sort_keys=True)
        resolve_folder_sample_defaults(first)
        self.assertEqual(json.dumps(first, sort_keys=True), before)
        self.source.targets = ['folder/107']
        second = self.scan()
        combined = copy.deepcopy(first + second)
        resolve_folder_sample_defaults(combined)
        rows = {r['ObservedExpno']: r for r in combined}
        for exp in ('105', '106', '107'):
            self.assertEqual(rows[exp]['ObservedSampleID'], '')
            self.assertTrue(rows[exp]['ObservedSampleTagIssue'])
        self.assertEqual(rows['1']['ObservedSampleID'], 'folder')
        before = json.dumps(combined, sort_keys=True)
        resolve_folder_sample_defaults(combined)
        self.assertEqual(json.dumps(combined, sort_keys=True), before)

    def test_actual_sample_formulas_use_block_observations_safely(self):
        from test_sample_tag_schema import evaluate
        from nmr_register.schema import SAMPLE_TAG_TRIGGER, SAMPLE_TAG_ISSUE_FORMULA, MACHINE_FIELDS
        self.acquisition(105, 'SAMPLE_ID_BLOCK=A\n')
        self.acquisition(106)
        fields = self.by_expno()['106']
        row = dict(fields, id=1, SourceID='', Folder='', Expno='', ReconciliationReason='', Sample=None)
        samples = [SimpleNamespace(id=8, SampleID='A')]
        experiments = [SimpleNamespace(**row)]
        self.assertEqual(evaluate(SAMPLE_TAG_TRIGGER, row, samples, experiments), 8)
        self.assertEqual(evaluate(SAMPLE_TAG_TRIGGER, row, samples, experiments, 7), 7)
        row['Sample'] = SimpleNamespace(id=7)
        self.assertIn('existing sample differs', evaluate(SAMPLE_TAG_ISSUE_FORMULA, row, samples, experiments))
        self.assertNotIn('Sample', MACHINE_FIELDS['Experiments'])
        self.assertNotIn('Assessment', MACHINE_FIELDS['Experiments'])
        self.acquisition(107, 'SAMPLE_ID_BLOCK=B\n')
        fields = self.by_expno()['106']
        row.update(fields, Sample=None)
        self.assertEqual(evaluate(SAMPLE_TAG_TRIGGER, row, samples, experiments), 0)
        self.assertTrue(evaluate(SAMPLE_TAG_ISSUE_FORMULA, row, samples, experiments))

    def test_protected_declaration_never_inspected(self):
        self.acquisition(105)
        protected = self.acquisition(106, 'SAMPLE_ID_BLOCK=secret\n', folder='211@950')
        from nmr_register import discovery, bruker
        original_scandir, original_read = discovery.os.scandir, bruker._read_metadata
        def guarded_scandir(path):
            self.assertNotIn('211@950', str(path))
            return original_scandir(path)
        def guarded_read(path):
            self.assertFalse(Path(path).is_relative_to(protected))
            return original_read(path)
        with patch.object(discovery.os, 'scandir', side_effect=guarded_scandir), patch.object(bruker, '_read_metadata', side_effect=guarded_read):
            rows = self.scan()
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['ObservedSampleID'], '')

    def test_split_declaration_cannot_supply_default(self):
        self.acquisition(105)
        split = self.acquisition(106, 'SAMPLE_ID_BLOCK=split\n')
        (split / 'audita.txt').write_text('created by split 4: dataset 1 of 4 from EXPNO 105')
        rows = self.scan()
        self.assertEqual([r['ObservedExpno'] for r in rows], ['105'])
        self.assertEqual(rows[0]['ObservedSampleID'], '')


if __name__ == '__main__':
    unittest.main()
