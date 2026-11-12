"""Canonical acquisition facts and exact acquisition-local title observations."""
import copy
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from nmr_register.bruker import observe,title_observation
from nmr_register.state import State
from nmr_register.sync import flush,SOURCE_FACTS
from nmr_register.grist import GristError
from test_sync import Fake,fields


class CanonicalCollectionTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.path=self.root/'folder'/'1';self.path.mkdir(parents=True)
        (self.path/'acqus').write_text('##$DATE_START= 1700000000\n##$DATE= 1700000060\n##$PULPROG= <hsqc>\n##$BF1= 600.48\n##$NUC1= <1H>\n##$TE= 298.0982\n')
        self.state=State(self.root/'state');self.addCleanup(self.state.close)

    def title(self,proc,text):
        path=self.path/'pdata'/str(proc)/'title';path.parent.mkdir(parents=True,exist_ok=True)
        path.write_bytes(text.encode() if isinstance(text,str) else text);return path

    def observation(self):
        return {**observe(self.path),'AcquisitionKey':'local:folder/1','ObservedSourceID':'local','ObservedFolder':'folder','ObservedExpno':'1'}

    def test_original_title_preserves_exact_comments_and_ignores_other_processing(self):
        self.title(1,'first comment\nSAMPLE_ID=tube\n')
        self.title(2,'other comments\nSAMPLE_ID=wrong\n')
        result=title_observation(self.path)
        self.assertEqual(result['Title'],'first comment\n')
        self.assertEqual(result['TitleSource'],'pdata/1/title')
        self.assertEqual(result['ObservedSampleID'],'tube')
        self.title(3,b'invalid\xff')
        self.assertEqual(title_observation(self.path),result)

    def test_readable_empty_or_tag_only_clears_but_missing_unreadable_preserves(self):
        self.title(1,'original comments\n')
        incoming=self.observation();self.state.observe([incoming]);client=Fake();flush(self.state,client)
        saved=client.rows[0]['fields'];self.assertEqual(saved['Title'],'original comments\n')
        self.title(1,'SAMPLE_ID=tube\n');self.state.observe([self.observation()]);flush(self.state,client)
        self.assertEqual(saved['Title'],'');self.assertEqual(saved['ObservedSampleID'],'tube')
        self.title(1,'replacement comments\n');self.state.observe([self.observation()]);flush(self.state,client)
        source=saved['TitleSource'];self.title(1,b'invalid\xff')
        bad=self.observation();self.assertNotIn('Title',bad);self.assertIn('comments unavailable',bad['ScanError'])
        self.state.observe([bad]);flush(self.state,client)
        self.assertEqual(saved['Title'],'replacement comments\n');self.assertEqual(saved['TitleSource'],source)
        (self.path/'pdata'/'1'/'title').unlink();missing=self.observation();self.assertNotIn('Title',missing)
        self.state.observe([missing]);flush(self.state,client);self.assertEqual(saved['Title'],'replacement comments\n')

    def test_other_processing_paths_never_hide_original_comments_or_tags(self):
        self.title(1,'first comment\nSAMPLE_ID=tube\n')
        (self.path/'pdata'/'2').mkdir()
        self.title(3,'second comment\nSAMPLE_ID=ignored\n')
        with patch('nmr_register.bruker.os.scandir',side_effect=AssertionError('no processing enumeration')):
            incoming=self.observation()
        self.assertEqual(incoming['Title'],'first comment\n')
        self.assertEqual(incoming['TitleSource'],'pdata/1/title')
        self.assertEqual(incoming['ObservedSampleID'],'tube')
        self.assertFalse(incoming['ObservedSampleTagIssue'])
        self.assertFalse(incoming.get('ScanError'))

    def test_optional_missing_titles_preserve_remote_comments_without_attention_input(self):
        self.title(1,'previous comment\n')
        self.state.observe([self.observation()]);client=Fake();flush(self.state,client)
        saved=client.rows[0]['fields'];source=saved['TitleSource']
        (self.path/'pdata'/'1'/'title').unlink()
        for stage in ('missing title','empty pdata','missing pdata'):
            with self.subTest(stage=stage):
                incoming=self.observation()
                self.assertNotIn('Title',incoming);self.assertNotIn('TitleSource',incoming)
                self.assertFalse(incoming.get('ScanError'))
                self.assertFalse(incoming['ObservedSampleTagIssue'])
                self.state.observe([incoming]);flush(self.state,client)
                self.assertEqual(saved['Title'],'previous comment\n')
                self.assertEqual(saved['TitleSource'],source)
            if stage=='missing title':(self.path/'pdata'/'1').rmdir()
            elif stage=='empty pdata':(self.path/'pdata').rmdir()

    def test_readable_empty_title_clears_remote_comments(self):
        self.title(1,'previous comment\n')
        self.state.observe([self.observation()]);client=Fake();flush(self.state,client)
        self.title(1,'')
        (self.path/'pdata'/'2').mkdir()
        incoming=self.observation();self.assertEqual(incoming['Title'],'')
        self.assertFalse(incoming.get('ScanError'))
        self.state.observe([incoming]);flush(self.state,client)
        self.assertEqual(client.rows[0]['fields']['Title'],'')

    def test_unreadable_other_processing_title_is_ignored_and_does_not_requeue(self):
        self.title(1,'previous comment\n')
        self.state.observe([self.observation()]);self.state.observe([self.observation()]);client=Fake();flush(self.state,client)
        saved=client.rows[0]['fields'];source=saved['TitleSource']
        self.title(2,b'invalid\xff')
        incoming=self.observation()
        self.assertEqual(incoming['Title'],'previous comment\n')
        self.assertFalse(incoming.get('ScanError'))
        self.assertEqual(self.state.observe([incoming]),0)
        self.title(2,'SAMPLE_ID=ignored-new\nSAMPLE_ID_FOLDER=ignored-default\n')
        self.assertEqual(self.state.observe([self.observation()]),0)
        self.assertEqual(saved['Title'],'previous comment\n');self.assertEqual(saved['TitleSource'],source)

    def test_title_edits_queue_without_fingerprint_or_replacement_and_never_inherit(self):
        self.title(1,'first\n');first=self.observation();self.state.observe([first],now=1)
        self.state.observe([self.observation()],now=2);self.state.acknowledge([first['AcquisitionKey']])
        self.title(1,'edited\n');new=self.observation();self.state.observe([new],now=3)
        self.assertEqual(new['ObservedFingerprint'],first['ObservedFingerprint'])
        self.assertEqual(self.state.pending()[0][1]['ReconciliationReason'],'')
        sibling=self.root/'folder'/'2';sibling.mkdir();(sibling/'acqus').write_text('')
        self.assertNotIn('Title',observe(sibling))

    def test_actual_facts_nominal_field_exact_temperature_retained_terminal_and_human_qc(self):
        self.title(1,'actual title\n');incoming=self.observation();incoming['ObservedInstrument']=''
        human={'SourceID':'local','Folder':'folder','Expno':'1','Instrument':'Recorded spectrometer','Field_MHz':850,
               'Temperature_C':25,'ExperimentType':'old type','AcquisitionStatus':'Finished','Assessment':'Use',
               'QCSNR':'Adequate','AssessmentReason':'human note','Sample':7,'SnapshotProtein':'Protein','Notes':'keep'}
        client=Fake([{'id':7,'fields':human}]);self.state.observe([incoming]);flush(self.state,client)
        saved=client.rows[0]['fields']
        self.assertEqual(saved['Instrument'],'Recorded spectrometer');self.assertEqual(saved['Field_MHz'],600)
        self.assertEqual(saved['Temperature_C'],24.9482);self.assertEqual(saved['ExperimentType'],'HSQC')
        self.assertEqual(saved['AcquisitionStatus'],'Completed');self.assertEqual(saved['Title'],'actual title\n')
        for key in ('Assessment','QCSNR','AssessmentReason','Sample','SnapshotProtein','Notes'):
            self.assertEqual(saved[key],human[key])
        self.assertTrue(all(not set(SOURCE_FACTS.values())&set(item['fields']) for item in client.patches))

    def test_unchanged_cache_reasserts_actuals_and_title_without_signature_changes_or_extra_reads(self):
        self.title(1,'actual\n');incoming=self.observation();self.state.observe([incoming],now=1);client=Fake();flush(self.state,client)
        before=dict(self.state.db.execute('SELECT * FROM observations').fetchone())
        saved=client.rows[0]['fields'];saved.update(Field_MHz=950,Title='manual override')
        count=len(client.patches)
        with patch('nmr_register.sync.scan_source',side_effect=AssertionError('must not rescan')):
            flush(self.state,client)
        self.assertEqual(saved['Field_MHz'],600);self.assertEqual(saved['Title'],'actual\n')
        self.assertEqual(dict(self.state.db.execute('SELECT * FROM observations').fetchone()),before)
        self.assertEqual(len(client.patches),count+1);self.assertFalse(self.state.pending())
        count=len(client.patches);flush(self.state,client);self.assertEqual(len(client.patches),count)

    def test_unknown_nonstatus_values_retain_records_and_literal_zero_is_actual(self):
        incoming=fields();incoming.update(ObservedD1_s=None,ObservedDelayList_s='',ObservedSpinLock_Hz=None,
             ObservedInstrument='',ObservedField_MHz=None,ObservedTemperature_C=0,ObservedAcquisitionStatus='Unknown')
        remote={**incoming,'D1_s':2,'DelayList_s':'[0,1]','SpinLock_Hz':900,'Instrument':'NMR','Field_MHz':850,
                'Temperature_C':25,'AcquisitionStatus':'Finished'}
        client=Fake([{'id':1,'fields':remote}]);self.state.observe([incoming]);flush(self.state,client)
        saved=client.rows[0]['fields']
        for name in ('D1_s','DelayList_s','SpinLock_Hz','Instrument','Field_MHz'):self.assertEqual(saved[name],remote[name])
        self.assertEqual(saved['Temperature_C'],0);self.assertEqual(saved['AcquisitionStatus'],'Completed')

    def test_identity_or_replacement_conflict_withholds_facts_titles_and_stays_pending(self):
        incoming=self.observation();incoming.update(Title='new actual')
        for conflict in ('copy or replacement',''):
            with self.subTest(conflict=conflict):
                temp=State(self.root/('conflict' if conflict else 'mismatch'))
                try:
                    if conflict:
                        temp.observe([{**incoming,'ObservedFingerprint':'old-fingerprint'}])
                    temp.observe([incoming])
                    remote={**incoming,'SourceID':'local','Folder':'folder' if conflict else 'different','Expno':'1',
                            'ReconciliationReason':conflict,'Field_MHz':850,'Title':'recorded title','Assessment':'Use'}
                    client=Fake([{'id':1,'fields':remote}]);result=flush(temp,client)
                    saved=client.rows[0]['fields'];self.assertEqual(saved['Field_MHz'],850);self.assertEqual(saved['Title'],'recorded title')
                    self.assertEqual(saved['Assessment'],'Use');self.assertTrue(temp.pending());self.assertEqual(result['conflicts'],1)
                    count=client.posts;flush(temp,client);self.assertEqual(client.posts,count)
                finally:temp.close()

    def test_remote_fingerprint_guard_survives_fresh_local_state(self):
        incoming=self.observation();incoming['Title']='new acquisition'
        remote={**incoming,'ObservedFingerprint':'previous acquisition','Field_MHz':850,'Title':'old acquisition','Assessment':'Use'}
        self.state.observe([incoming]);client=Fake([{'id':1,'fields':remote}]);result=flush(self.state,client)
        saved=client.rows[0]['fields']
        self.assertEqual(result['conflicts'],1);self.assertIn('fingerprint differs',saved['ReconciliationReason'])
        self.assertEqual(saved['Field_MHz'],850);self.assertEqual(saved['Title'],'old acquisition')
        self.assertEqual(saved['Assessment'],'Use');self.assertTrue(self.state.pending())

    def test_partial_read_keeps_remote_fingerprint_for_later_replacement_detection(self):
        first=self.observation();self.state.observe([first]);client=Fake();flush(self.state,client)
        saved=client.rows[0]['fields'];original=saved['ObservedFingerprint'];saved['Assessment']='Use'
        partial={**first,'ObservedFingerprint':'','ObservedType':'unknown'}
        self.state.observe([partial]);flush(self.state,client)
        self.assertEqual(saved['ObservedFingerprint'],original)
        replaced={**first,'ObservedFingerprint':'new acquisition','ObservedField_MHz':950}
        self.state.observe([replaced]);flush(self.state,client)
        self.assertIn('fingerprint differs',saved['ReconciliationReason'])
        self.assertEqual(saved['Field_MHz'],600);self.assertEqual(saved['Assessment'],'Use')
        self.assertTrue(self.state.pending())

    def test_partial_confirmed_identity_disagreement_withholds_source_facts(self):
        incoming=self.observation();incoming['Title']='actual title'
        for field,value in [('SourceID','other source'),('Folder','other folder'),('Expno','99')]:
            with self.subTest(field=field):
                state=State(self.root/('partial-'+field))
                try:
                    state.observe([incoming])
                    remote={**incoming,'SourceID':'','Folder':'','Expno':'',field:value,
                            'Title':'recorded title','Field_MHz':850,'Assessment':'Use'}
                    client=Fake([{'id':1,'fields':remote}]);result=flush(state,client)
                    saved=client.rows[0]['fields']
                    self.assertEqual(result['conflicts'],1);self.assertTrue(state.pending())
                    self.assertEqual(saved['Title'],'recorded title');self.assertEqual(saved['Field_MHz'],850)
                    self.assertEqual(saved['Assessment'],'Use');self.assertIn('confirmed identity',saved['ReconciliationReason'])
                finally:state.close()

    def test_observed_identity_requires_numeric_explicit_folder_and_exact_key(self):
        incoming=self.observation();incoming['Title']='actual'
        for index,changes in enumerate([{'ObservedFolder':None},{'ObservedExpno':'one'},
                {'AcquisitionKey':'local:wrong/1'},{'ObservedSourceID':None}]):
            with self.subTest(changes=changes):
                state=State(self.root/('invalid-observed-'+str(index)))
                try:
                    payload={**incoming,**changes};state.observe([payload])
                    remote={**payload,'Field_MHz':850,'Title':'recorded','Assessment':'Use'}
                    client=Fake([{'id':1,'fields':remote}]);result=flush(state,client)
                    self.assertEqual(result['conflicts'],1);self.assertTrue(state.pending())
                    self.assertEqual(client.rows[0]['fields']['Field_MHz'],850)
                    self.assertEqual(client.rows[0]['fields']['Title'],'recorded')
                finally:state.close()
        root={**incoming,'ObservedFolder':'','AcquisitionKey':'local:1'}
        self.state.observe([root]);client=Fake();flush(self.state,client)
        self.assertEqual(client.rows[0]['fields']['Field_MHz'],600)
        self.assertEqual(client.rows[0]['fields']['Title'],'actual')
        self.assertFalse(self.state.pending())

    def test_cached_correction_timeout_after_commit_retries_without_duplicate_or_second_patch(self):
        self.title(1,'actual');self.state.observe([self.observation()]);client=Fake();flush(self.state,client)
        client.rows[0]['fields']['Title']='manual';original=client.update_records
        def timeout(table,items):original(table,items);raise GristError('timeout after commit')
        client.update_records=timeout
        with self.assertRaises(GristError):flush(self.state,client)
        self.assertTrue(self.state.pending());count=len(client.patches);client.update_records=original;flush(self.state,client)
        self.assertFalse(self.state.pending());self.assertEqual(len(client.patches),count);self.assertEqual(client.posts,1)

    def test_unsafe_or_oversized_original_title_never_replaces_known_comments(self):
        from nmr_register.bruker import MAX_METADATA_BYTES
        self.title(1,b'x'*(MAX_METADATA_BYTES+1))
        result=title_observation(self.path)
        self.assertNotIn('Title',result);self.assertTrue(result['ObservedSampleTagIssue'])
        outside=self.root/'outside';outside.write_text('outside')
        (self.path/'pdata'/'1'/'title').unlink();(self.path/'pdata'/'1'/'title').symlink_to(outside)
        self.assertNotIn('Title',title_observation(self.path))
        with patch('nmr_register.bruker._read_metadata',side_effect=AssertionError('must prune before read')):
            self.assertNotIn('Title',title_observation(self.root/'211@950'/'1'))


if __name__=='__main__':unittest.main()
