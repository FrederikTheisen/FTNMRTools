"""Explicit acquisition-local title declarations; no actual NMR source reads."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

SOURCE = Path(__file__).resolve().parents[1]/'nmr_register'/'bruker.py'
SPEC = importlib.util.spec_from_file_location('sample_tags_bruker', SOURCE)
bruker = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bruker)


class SampleTagTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.acq = self.root/'folder'/'1'
        self.acq.mkdir(parents=True)

    def title(self, procno, data):
        path = self.acq/'pdata'/str(procno)/'title'
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data.encode('utf-8') if isinstance(data, str) else data)
        return path

    def result(self):
        return bruker.sample_tag_observation(self.acq)

    def test_only_original_processing_title_is_read_without_enumeration(self):
        self.title(1, 'comments\nSAMPLE_ID=0\n')
        self.title(0, 'SAMPLE_ID=wrong-zero\n')
        self.title(2, 'SAMPLE_ID=wrong-two\nSAMPLE_ID_FOLDER=wrong-folder\n')
        (self.acq/'title').write_text('SAMPLE_ID=wrong-root\n')
        original = bruker._read_metadata
        reads = []
        def read(path):
            reads.append(Path(path))
            return original(path)
        with patch.object(bruker.os, 'scandir', side_effect=AssertionError('no processing enumeration')):
            with patch.object(bruker, '_read_metadata', side_effect=read):
                result = self.result()
        self.assertEqual(reads, [self.acq/'pdata'/'1'/'title'])
        self.assertEqual(result['ObservedSampleID'], '0')
        self.assertEqual(result['ObservedSampleTagIssue'], '')
        self.assertIn('pdata/1/title', result['ObservedSampleTagEvidence'])

    def test_adversarial_other_processing_titles_are_never_read(self):
        self.title(1, 'comments\nSAMPLE_ID=trusted\n')
        expected = bruker.title_observation(self.acq)
        for content in ('SAMPLE_ID=\nSAMPLE_ID_FOLDER=\n', b'invalid\xff', b'x'*(bruker.MAX_METADATA_BYTES+1)):
            self.title(2, content)
            self.assertEqual(bruker.title_observation(self.acq), expected)
        other = self.acq/'pdata'/'2'
        (other/'title').unlink(); other.rmdir()
        outside = self.root/'outside-processing'; outside.mkdir()
        (outside/'title').write_text('SAMPLE_ID=untrusted\n')
        other.symlink_to(outside, target_is_directory=True)
        self.assertEqual(bruker.title_observation(self.acq), expected)
        other.unlink(); other.mkdir()
        (other/'title').symlink_to(outside/'title')
        self.assertEqual(bruker.title_observation(self.acq), expected)

    def test_duplicate_tags_allowed_and_conflicts_in_one_title_fail_closed(self):
        self.title(1, 'SAMPLE_ID=A\nSAMPLE_ID=A\n')
        self.assertEqual(self.result()['ObservedSampleID'], 'A')
        self.title(1, 'SAMPLE_ID=A\nSAMPLE_ID=B\n')
        self.assertEqual(self.result()['ObservedSampleID'], '')
        self.assertIn('Conflicting SAMPLE_ID', self.result()['ObservedSampleTagIssue'])

    def test_missing_and_untagged_are_silent_optional_absence(self):
        for setup in [lambda: None, lambda: self.title(2, 'SAMPLE_ID=ignored'),
                      lambda: self.title(1, ''), lambda: self.title(1, 'ordinary comments')]:
            setup()
            self.assertEqual(self.result(), {'ObservedSampleID':'', 'ObservedSampleTagIssue':'', 'ObservedSampleTagEvidence':''})
            self.assertNotIn('_TitleIssue', bruker.title_observation(self.acq))

    def test_blank_malformed_and_invalid_utf8_fail_closed(self):
        for content in ['SAMPLE_ID=', 'SAMPLE_ID =sample', 'SAMPLE_ID:sample', ' SAMPLE_ID=sample', 'sample_id=sample',
                        'SAMPLE_ID= sample', 'SAMPLE_ID=sample ', 'SAMPLE_ID=bad\x00id', b'SAMPLE_ID=bad\xff']:
            with self.subTest(content=content):
                self.title(1, content)
                self.assertEqual(self.result()['ObservedSampleID'], '')
                self.assertTrue(self.result()['ObservedSampleTagIssue'])

    def test_failed_or_oversized_original_title_fails_closed(self):
        self.title(1, b'x'*(bruker.MAX_METADATA_BYTES+1))
        self.assertTrue(self.result()['ObservedSampleTagIssue'])
        with patch.object(bruker, '_read_metadata', side_effect=PermissionError('synthetic denied read')):
            self.assertEqual(self.result()['ObservedSampleID'], '')
            self.assertIn('PermissionError', self.result()['ObservedSampleTagIssue'])

    def test_symlink_original_title_and_ancestors_are_not_followed(self):
        outside = self.root/'outside'
        outside.mkdir()
        (outside/'title').write_text('SAMPLE_ID=outside\n')
        title = self.title(1, 'SAMPLE_ID=valid\n')
        title.unlink(); title.symlink_to(outside/'title')
        with patch.object(Path, 'open', side_effect=AssertionError('unsafe read')):
            self.assertIn('unsafe', self.result()['ObservedSampleTagIssue'])
        title.unlink(); title.parent.rmdir(); title.parent.symlink_to(outside, target_is_directory=True)
        self.assertIn('unsafe', self.result()['ObservedSampleTagIssue'])
        alias = self.root/'alias'; alias.symlink_to(self.acq, target_is_directory=True)
        self.assertIn('unsafe', bruker.sample_tag_observation(alias)['ObservedSampleTagIssue'])

    def test_protected_component_rejected_before_read(self):
        with patch.object(bruker, '_read_metadata', side_effect=AssertionError('protected read')):
            self.assertTrue(bruker.sample_tag_observation(self.root/'211@950'/'1')['ObservedSampleTagIssue'])


if __name__ == '__main__':
    unittest.main()
