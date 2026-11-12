"""Explicit folder defaults apply only to monitored acquisition siblings."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from nmr_register.bruker import title_observation
from nmr_register.config import Source
from nmr_register.discovery import scan_source,resolve_folder_sample_defaults
from nmr_register.state import State


class FolderSampleDefaultTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)/'source';self.root.mkdir()
        self.source=Source('source',str(self.root))

    def acquisition(self,expno,title=None,folder='folder',root=None):
        path=(root or self.root)/folder/str(expno);path.mkdir(parents=True)
        (path/'acqus').write_text('##$DATE_START= '+str(1700000000+int(expno))+'\n##$PULPROG= <hsqc>\n')
        if title is not None:self.title(path,title)
        return path

    def title(self,path,text,procno=1):
        title=path/'pdata'/str(procno)/'title';title.parent.mkdir(parents=True,exist_ok=True)
        title.write_bytes(text.encode() if isinstance(text,str) else text);return title

    def scan(self,source=None):
        records,errors=scan_source(source or self.source);self.assertEqual(errors,[])
        return records

    def by_expno(self,records):return {r['ObservedExpno']:r for r in records}

    def test_explicit_default_includes_holder_and_untagged_peers_direct_tag_wins(self):
        holder=self.acquisition(1,'comments\nSAMPLE_ID_FOLDER=folder-tube\n')
        self.acquisition(2)
        self.acquisition(3,'SAMPLE_ID=other-tube\n')
        rows=self.by_expno(self.scan())
        self.assertEqual([rows[str(i)]['ObservedSampleID'] for i in (1,2,3)],['folder-tube','folder-tube','other-tube'])
        self.assertEqual(rows['1']['Title'],'comments\n')
        self.assertIn('EXPNO 1: pdata/1/title',rows['2']['ObservedSampleTagEvidence'])
        self.assertEqual(rows['3']['ObservedSampleTagIssue'],'')
        self.assertNotIn('SAMPLE_ID_FOLDER',rows['1']['Title'])
        self.assertEqual((holder/'pdata'/'1'/'title').read_text(),'comments\nSAMPLE_ID_FOLDER=folder-tube\n')

    def test_same_default_repeats_are_allowed_conflicting_defaults_never_choose_winner(self):
        first=self.acquisition(1,'SAMPLE_ID_FOLDER=A\n')
        self.title(first,'SAMPLE_ID_FOLDER=A\n',2)
        second=self.acquisition(2,'SAMPLE_ID_FOLDER=A\n');self.acquisition(3)
        self.assertTrue(all(r['ObservedSampleID']=='A' for r in self.scan()))
        self.title(second,'SAMPLE_ID_FOLDER=B\n')
        self.acquisition(4,'SAMPLE_ID=direct\n')
        rows=self.by_expno(self.scan())
        for exp in ('1','2','3'):
            self.assertEqual(rows[exp]['ObservedSampleID'],'')
            self.assertIn('Conflicting SAMPLE_ID_FOLDER',rows[exp]['ObservedSampleTagIssue'])
        self.assertEqual(rows['4']['ObservedSampleID'],'direct');self.assertEqual(rows['4']['ObservedSampleTagIssue'],'')

    def test_scope_is_exact_source_and_immediate_folder_no_parent_or_child_default(self):
        self.acquisition(1,'SAMPLE_ID_FOLDER=A\n',folder='parent')
        self.acquisition(2,folder='parent/child');self.acquisition(3,folder='unrelated')
        rows=self.by_expno(self.scan())
        self.assertEqual(rows['1']['ObservedSampleID'],'A')
        self.assertEqual(rows['2']['ObservedSampleID'],'');self.assertEqual(rows['3']['ObservedSampleID'],'')
        other=self.root.parent/'second-source';other.mkdir();self.acquisition(4,folder='parent',root=other)
        other_rows=self.scan(Source('other',str(other)))
        resolve_folder_sample_defaults([*copy.deepcopy(list(rows.values())),*other_rows])
        self.assertEqual(other_rows[0]['ObservedSampleID'],'')

    def test_default_source_must_be_inside_allowlist_and_no_extra_title_reads(self):
        unmonitored=self.acquisition(1,'SAMPLE_ID_FOLDER=outside-scope\n');self.acquisition(2)
        self.source.targets=['folder/2']
        from nmr_register import bruker
        original=bruker._read_metadata
        def guarded(path):
            self.assertFalse(Path(path).is_relative_to(unmonitored),'unmonitored acquisition inspected')
            return original(path)
        with patch.object(bruker,'_read_metadata',side_effect=guarded):rows=self.scan()
        self.assertEqual(len(rows),1);self.assertEqual(rows[0]['ObservedSampleID'],'')
        self.source.targets=['folder/1','folder/2']
        self.assertTrue(all(r['ObservedSampleID']=='outside-scope' for r in self.scan()))

    def test_malformed_folder_defaults_block_inheritance_but_direct_sample_id_has_priority(self):
        self.acquisition(1,'SAMPLE_ID_FOLDER=A\n');bad=self.acquisition(2,'SAMPLE_ID_FOLDER=\n')
        self.acquisition(3)
        rows=self.by_expno(self.scan())
        self.assertTrue(all(not r['ObservedSampleID'] for r in rows.values()))
        self.assertTrue(all('malformed SAMPLE_ID_FOLDER' in r['ObservedSampleTagIssue'] for r in rows.values()))
        self.title(bad,'SAMPLE_ID=direct\nSAMPLE_ID_FOLDER:bad\n')
        rows=self.by_expno(self.scan());self.assertEqual(rows['2']['ObservedSampleID'],'direct')
        self.assertEqual(rows['2']['ObservedSampleTagIssue'],'');self.assertFalse(rows['3']['ObservedSampleID'])

    def test_malformed_direct_tag_never_falls_back_to_default_or_blocks_readable_siblings(self):
        self.acquisition(1,'SAMPLE_ID_FOLDER=A\n');self.acquisition(2,'SAMPLE_ID=\n');self.acquisition(3)
        rows=self.by_expno(self.scan())
        self.assertEqual(rows['2']['ObservedSampleID'],'');self.assertIn('malformed SAMPLE_ID',rows['2']['ObservedSampleTagIssue'])
        self.assertEqual(rows['3']['ObservedSampleID'],'A')

    def test_unreadable_or_unsafe_titles_block_opted_in_group_missing_titles_are_eligible(self):
        self.acquisition(1,'SAMPLE_ID_FOLDER=A\n');bad=self.acquisition(2,b'bad\xff');self.acquisition(3)
        rows=self.by_expno(self.scan());self.assertTrue(all(not r['ObservedSampleID'] for r in rows.values()))
        self.assertTrue(rows['3']['ObservedSampleTagIssue'])
        (bad/'pdata'/'1'/'title').unlink()
        self.assertTrue(all(r['ObservedSampleID']=='A' for r in self.scan()))
        outside=self.root.parent/'outside-title';outside.write_text('SAMPLE_ID_FOLDER=B\n')
        (bad/'pdata'/'1'/'title').symlink_to(outside)
        self.assertTrue(all(not r['ObservedSampleID'] for r in self.scan()))

    def test_edit_and_removal_requeue_unchanged_siblings_without_identity_change_or_public_context(self):
        holder=self.acquisition(1,'SAMPLE_ID_FOLDER=A\n');self.acquisition(2)
        state=State(self.root.parent/'state');self.addCleanup(state.close)
        first=self.scan();fingerprints={r['AcquisitionKey']:r['ObservedFingerprint'] for r in first}
        state.observe(first,now=1);state.observe(self.scan(),now=2);state.acknowledge(list(fingerprints))
        self.title(holder,'SAMPLE_ID_FOLDER=B\n');state.observe(self.scan(),now=3)
        self.assertEqual(len(state.pending()),2)
        for key,fields in state.pending():
            self.assertEqual(fields['ObservedSampleID'],'B');self.assertEqual(fields['ObservedFingerprint'],fingerprints[key])
            self.assertEqual(fields['ReconciliationReason'],'');self.assertFalse(any(k.startswith('_Folder') for k in fields))
            self.assertNotIn('Sample',fields)
        state.acknowledge(list(fingerprints));self.title(holder,'ordinary comments\n');state.observe(self.scan(),now=4)
        self.assertEqual(len(state.pending()),2)
        self.assertTrue(all(not fields['ObservedSampleID'] for _,fields in state.pending()))

    def test_legacy_no_opt_in_does_not_add_context_or_change_observations(self):
        self.acquisition(1,'ordinary comments\n');self.acquisition(2,'SAMPLE_ID=local\n')
        rows=self.scan();before=copy.deepcopy(rows)
        self.assertFalse(any(k.startswith('_Folder') for r in rows for k in r))
        resolve_folder_sample_defaults(rows);self.assertEqual(rows,before)
        state=State(self.root.parent/'state');self.addCleanup(state.close)
        state.observe(rows,now=1);state.observe(self.scan(),now=2);state.acknowledge([r['AcquisitionKey'] for r in rows])
        self.assertEqual(state.observe(self.scan(),now=3),0);self.assertFalse(state.pending())

    def test_root_folder_is_supported_but_incomplete_group_identity_is_not(self):
        self.acquisition(1,'SAMPLE_ID_FOLDER=root-tube\n',folder='');self.acquisition(2,folder='')
        rows=self.scan();self.assertTrue(all(r['ObservedSampleID']=='root-tube' for r in rows))
        invalid=[{'ObservedSourceID':None,'ObservedFolder':None,'_FolderSampleID':'A','_FolderSampleEvidence':'pdata/1/title'},
                 {'ObservedSourceID':None,'ObservedFolder':None,'ObservedSampleID':''}]
        resolve_folder_sample_defaults(invalid);self.assertEqual(invalid[1]['ObservedSampleID'],'')


if __name__=='__main__':unittest.main()
