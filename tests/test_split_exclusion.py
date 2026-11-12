"""Derived processing datasets never enter the acquisition inventory."""
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from nmr_register.discovery import scan_source
from nmr_register.split import is_derived_split, SPLIT_EVIDENCE, SPLIT_SCAN_ERROR

class SplitExclusionTests(unittest.TestCase):
    def test_all_confirmed_types_excluded_before_folder_defaults(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for number,pulse in enumerate(('hsqc','noe','etaXY','t1'),1):
                path=root/str(number);path.mkdir()
                (path/'acqus').write_text('##$PULPROG= <'+pulse+'>\n##$DATE= 1748545685\n')
                title=path/'pdata/1';title.mkdir(parents=True)
                (title/'title').write_text('SAMPLE_ID_FOLDER=split-only' if number>1 else 'Original')
                if number>1:
                    (path/'audita.txt').write_text('created by split 4: dataset 1 of 4 from EXPNO 1')
            rows,errors=scan_source(SimpleNamespace(id='src',path=tmp,targets=None,exclude=(),instrument=''))
            self.assertFalse(errors)
            self.assertEqual([r['ObservedExpno'] for r in rows],['1'])
            self.assertFalse(rows[0].get('ObservedSampleID'))
            self.assertTrue((root/'2/acqus').exists())

    def test_only_exact_machine_evidence_counts(self):
        for field in ('Title','Notes','ReconciliationReason','ObservedPulseProgram','SourcePath'):
            self.assertFalse(is_derived_split({field:SPLIT_EVIDENCE}))
        self.assertFalse(is_derived_split({'ScanError':'possible derived split experiment'}))
        self.assertTrue(is_derived_split({'ObservedAcquisitionEvidence':SPLIT_EVIDENCE}))

        self.assertTrue(is_derived_split({'ScanError':SPLIT_SCAN_ERROR}))
        self.assertTrue(is_derived_split({'ScanError':'acqus missing; '+SPLIT_SCAN_ERROR}))
