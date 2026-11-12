"""Title-only edits queue observations without changing acquisition identity."""
import hashlib
import importlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

PACKAGE = Path(__file__).resolve().parents[1]/'nmr_register'
SPEC = importlib.util.spec_from_file_location('sample_tag_stage', PACKAGE/'__init__.py', submodule_search_locations=[str(PACKAGE)])
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)
bruker = importlib.import_module('sample_tag_stage.bruker')
discovery = importlib.import_module('sample_tag_stage.discovery')
State = importlib.import_module('sample_tag_stage.state').State

ACQUS = '##$DATE_START= 1700000000\n##$DATE= 1700000060\n##$PULPROG= <hsqc>\n##$BF1= 600.13\n##$SFO1= 600.13\n##$NUC1= <1H>\n##$TD= 1024\n'


class SampleTagCollectionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)/'source'
        self.root.mkdir()
        self.acq = self.experiment(1)
        self.source = SimpleNamespace(id='local', path=self.root, targets=None, exclude=[], instrument='')

    def experiment(self, expno):
        path = self.root/'folder'/str(expno)
        path.mkdir(parents=True)
        (path/'acqus').write_text(ACQUS.replace('1700000000', str(1700000000+expno)))
        (path/'fid').write_bytes(b'partial fixture')
        return path

    def title(self, path, value):
        title = path/'pdata'/'1'/'title'
        title.parent.mkdir(parents=True, exist_ok=True)
        title.write_text(value)
        return title

    def scan(self):
        observations, errors = discovery.scan_source(self.source)
        self.assertEqual(errors, [])
        return observations

    def test_title_only_edit_and_removal_queue_without_replacement(self):
        title = self.title(self.acq, 'SAMPLE_ID=sample-A\n')
        state = State(Path(self.temp.name)/'state')
        self.addCleanup(state.close)
        first = self.scan()
        observed = first[0]
        identity = {k:bruker.parameters((self.acq/'acqus').read_text()).get(k, '') for k in ('DATE_START','DATE','PULPROG','BF1','SFO1','NUC1','TD')}
        self.assertEqual(observed['ObservedFingerprint'], hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest())
        state.observe(first, now=1)
        state.observe(self.scan(), now=2)
        state.acknowledge([observed['AcquisitionKey']])
        self.assertEqual(state.observe(self.scan(), now=3), 0)
        title.write_text('SAMPLE_ID=sample-B\n')
        self.assertEqual(state.observe(self.scan(), now=4), 1)
        key, changed = state.pending()[0]
        self.assertEqual(changed['ObservedSampleID'], 'sample-B')
        self.assertEqual(changed['ObservedFingerprint'], observed['ObservedFingerprint'])
        self.assertEqual(changed['ReconciliationReason'], '')
        self.assertEqual(changed['RawBytes'], observed['RawBytes'])
        state.acknowledge([key])
        title.unlink()
        state.observe(self.scan(), now=5)
        cleared = state.pending()[0][1]
        self.assertEqual(cleared['ObservedSampleID'], '')
        self.assertEqual(cleared['ObservedSampleTagIssue'], '')
        self.assertEqual(cleared['ObservedSampleTagEvidence'], '')
        self.assertEqual(cleared['ReconciliationReason'], '')
        self.assertNotIn('Sample', cleared)

    def test_no_parent_or_sibling_inheritance(self):
        sibling = self.experiment(2)
        self.title(self.acq, 'SAMPLE_ID=only-first-tube\n')
        (self.acq.parent/'title').write_text('SAMPLE_ID=parent-must-not-bind\n')
        records = {r['ObservedExpno']:r for r in self.scan()}
        self.assertEqual(records['1']['ObservedSampleID'], 'only-first-tube')
        self.assertEqual(records['2']['ObservedSampleID'], '')
        self.assertEqual(bruker.sample_tag_observation(sibling)['ObservedSampleID'], '')

    def test_scanner_prunes_protected_name_before_inspection(self):
        (self.root/'211@950').mkdir()
        original = discovery.os.scandir
        class Protected:
            name = '211@950'
            def is_symlink(self):
                raise AssertionError('protected entry inspected')
            def is_dir(self, **kwargs):
                raise AssertionError('protected entry inspected')
        class Entries:
            def __enter__(self): return iter([Protected()])
            def __exit__(self, *args): pass
        with patch.object(discovery.os, 'scandir', return_value=Entries()):
            self.assertEqual(self.scan(), [])
        self.assertIsNotNone(original)


if __name__ == '__main__':
    unittest.main()
