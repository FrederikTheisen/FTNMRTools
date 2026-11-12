import unittest, tempfile
from pathlib import Path
from nmr_register.legacy import review_snapshot
from nmr_register.migration import read_reviews
class LegacyTests(unittest.TestCase):
    def test_preserves_condition_and_blank_decisions_and_rejects_incomplete(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'in.tsv'; path.write_text('# context\ttext\ntarget\tassessment\tcomment\ncondition:one\texclude\tlegacy\na/1\tusable\told\na/2\t\t\n')
            with self.assertRaisesRegex(ValueError, 'without a Grist experiment'): review_snapshot(path,[])
            path.write_text('# context\ttext\ntarget\tassessment\tcomment\na/1\tusable\told\na/2\t\t\n')
            payload=review_snapshot(path,[{'Folder':'a','Expno':'1','Assessment':''},{'Folder':'a','Expno':'2','Assessment':'Do not use','AssessmentReason':'bad'}])
            out=Path(d)/'out.tsv'; out.write_bytes(payload); rows=read_reviews(out)
            self.assertEqual(rows['a/1']['assessment'],''); self.assertEqual(rows['a/1']['comment'],''); self.assertEqual(rows['a/2']['assessment'],'exclude')
            payload=review_snapshot(path,[{'Folder':'a','Expno':'1','Assessment':''},{'Folder':'a','Expno':'2','Assessment':'Do not use','AssessmentExplanation':'FID: Unusable'}])
            out.write_bytes(payload); self.assertEqual(read_reviews(out)['a/2']['comment'],'FID: Unusable')
            with self.assertRaises(ValueError): review_snapshot(path,[{'Folder':'a','Expno':'2','Assessment':'Do not use'}])
            self.assertIn('a/2\t\t',path.read_text())

class LegacyIntegrationTests(unittest.TestCase):
    def test_snapshot_mode_never_initialises_or_writes_review_file(self):
        import importlib.util, sys, types
        from types import SimpleNamespace
        stub=types.ModuleType('relax_tools')
        previous=sys.modules.get('relax_tools'); sys.modules['relax_tools']=stub
        try:
            module_path=Path(__file__).resolve().parents[1]/'relax'/'condition_overview.py'
            spec=importlib.util.spec_from_file_location('staged_condition_overview',module_path)
            mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
            mod.discover_reviews=lambda *a,**kw: []
            mod.initialise_reviews=lambda *a,**kw: self.fail('snapshot mode attempted a review write')
            with tempfile.TemporaryDirectory() as d:
                path=Path(d)/'snapshot.tsv'; path.write_text('target\tassessment\tcomment\n')
                args=SimpleNamespace(project_root=d,reviews_only=True,refresh_reports=False,catalog='',data_root='',fid_inventory=None,r1_qc=None,qc_decisions=None,grist_review_snapshot=str(path),user_reviews=None,update_review_context=False,dry_run=False)
                before=path.read_bytes(); mod.run(args); self.assertEqual(path.read_bytes(),before)
                args.update_review_context=True
                with self.assertRaises(ValueError): mod.run(args)
        finally:
            if previous is None: sys.modules.pop('relax_tools',None)
            else: sys.modules['relax_tools']=previous

class SourceGuardTests(unittest.TestCase):
    def test_requires_explicit_source_and_never_uses_other_source(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'reviews.tsv'; path.write_text('target\tassessment\tcomment\nf/1\t\t\n')
            data=[{'SourceID':'A','Folder':'f','Expno':'1','Assessment':'Use'}, {'SourceID':'B','Folder':'f','Expno':'1','Assessment':'Do not use','AssessmentReason':'bad'}]
            with self.assertRaisesRegex(ValueError,'source_id'): review_snapshot(path,data)
            out=Path(d)/'out.tsv'; out.write_bytes(review_snapshot(path,data,source_id='A'))
            self.assertEqual(read_reviews(out)['f/1']['assessment'],'usable')

class MetadataConflictTests(unittest.TestCase):
    def test_metadata_difference_blocks_legacy_export(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'reviews.tsv'; path.write_text('target\tassessment\tcomment\nf/1\t\t\n')
            with self.assertRaisesRegex(ValueError,'Unresolved experiment'):
                review_snapshot(path,[{'Folder':'f','Expno':'1','MetadataDifference':'date mismatch'}])
