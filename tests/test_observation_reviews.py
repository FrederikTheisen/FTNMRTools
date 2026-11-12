import copy
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
import nmr_register

ROOT=Path(__file__).resolve().parents[1]
for name in ('state','sync'):
    spec=importlib.util.spec_from_file_location('nmr_register.'+name,ROOT/'nmr_register'/(name+'.py'))
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
State=sys.modules['nmr_register.state'].State
flush=sys.modules['nmr_register.sync'].flush

def observation(expno='1',fingerprint='shared',inode=10,title='source title'):
    return {'AcquisitionKey':'src:folder/'+expno,'ObservedSourceID':'src','ObservedFolder':'folder','ObservedExpno':expno,
        'ObservedFingerprint':fingerprint,'ObservedType':'r1','ObservedPulseProgram':'t1-program','ObservedField_MHz':600.48,
        'ObservedAcquisitionStatus':'Finished','RawPresent':True,'RawBytes':100,'ScanError':'','Title':title,'TitleSource':'pdata/1/title',
        '_RawFileStats':[['ser',100,1000,1,inode]]}

class Fake:
    def __init__(self):self.rows=[];self.patches=[]
    def records(self,table):return copy.deepcopy(self.rows)
    def add_records(self,table,rows):
        for row in rows:self.rows.append({'id':len(self.rows)+1,'fields':copy.deepcopy(row)})
    def update_records(self,table,rows):
        self.patches.extend(copy.deepcopy(rows))
        for row in rows:next(r for r in self.rows if r['id']==row['id'])['fields'].update(row['fields'])

class ReviewTransportTests(unittest.TestCase):
    def setUp(self):self.temp=tempfile.TemporaryDirectory();self.state=State(self.temp.name);self.client=Fake()
    def tearDown(self):self.state.close();self.temp.cleanup()
    def test_shared_metadata_is_visible_but_acknowledged_and_source_facts_sync(self):
        rows=[observation(),observation('2',inode=20)]
        self.state.observe(rows);first=flush(self.state,self.client)
        self.assertEqual(first['created'],2);self.assertEqual(first['conflicts'],0);self.assertFalse(self.state.pending())
        for row in self.client.rows:
            self.assertIn('matching acquisition metadata',row['fields']['ReconciliationReason'])
            self.assertEqual(row['fields']['PulseProgram'],'t1-program');self.assertEqual(row['fields']['Title'],'source title')
        self.client.rows[0]['fields'].update(Assessment='Use',Notes='human',Sample=22)
        rows[0]['Title']='edited source comment';self.state.observe(rows);flush(self.state,self.client)
        saved=self.client.rows[0]['fields'];self.assertEqual(saved['Title'],'edited source comment')
        self.assertEqual((saved['Assessment'],saved['Notes'],saved['Sample']),('Use','human',22))
        self.state.observe(rows);flush(self.state,self.client);self.state.observe(rows)
        stable=flush(self.state,self.client);self.assertEqual(stable['updated'],0);self.assertFalse(self.state.pending())
    def test_missing_fingerprint_remains_reviewable_but_split_is_excluded(self):
        rows=[observation(fingerprint=''),observation('2',fingerprint='unique',inode=20)]
        rows[1]['ScanError']='derived split experiment; original acquisition requires review'
        self.state.observe(rows);result=flush(self.state,self.client)
        self.assertEqual(result['conflicts'],0);self.assertFalse(self.state.pending())
        self.assertEqual([r['fields']['AcquisitionKey'] for r in self.client.rows], ['src:folder/1'])
        self.assertEqual(len(self.state.cached_observations()), 1)
        self.assertTrue(all(r['fields']['ReconciliationReason'] for r in self.client.rows))
        self.assertTrue(all(r['fields']['Title']=='source title' for r in self.client.rows))
    def test_old_split_cache_and_outbox_are_purged_before_flush(self):
        from nmr_register.state import serialized, MATCHING_METADATA_PREFIX
        from nmr_register.split import SPLIT_EVIDENCE
        rows=[observation(),observation('2',inode=20),observation('3',inode=30)]
        self.state.observe(rows)
        # Emulate the pre-exclusion durable database, including an orphan outbox.
        split=dict(self.state.cached_observations())['src:folder/2']
        split['ObservedAcquisitionEvidence']=SPLIT_EVIDENCE
        with self.state.db:
            self.state.db.execute('UPDATE observations SET payload=? WHERE key=?',(serialized(split),split['AcquisitionKey']))
            self.state.db.execute('UPDATE outbox SET payload=? WHERE key=?',(serialized(split),split['AcquisitionKey']))
            orphan=observation('4');orphan['ObservedAcquisitionEvidence']=SPLIT_EVIDENCE
            self.state.db.execute('INSERT INTO outbox VALUES (?,?,?)',(orphan['AcquisitionKey'],serialized(orphan),1))
        result=flush(self.state,self.client)
        self.assertEqual(result['created'],2)
        self.assertEqual({r['fields']['AcquisitionKey'] for r in self.client.rows},{'src:folder/1','src:folder/3'})
        for row in self.client.rows:
            self.assertNotIn('folder/2',row['fields']['ReconciliationReason'])
            self.assertIn('matching acquisition metadata',row['fields']['ReconciliationReason'])
        self.state.requeue_observations(['src:folder/2','src:folder/4'])
        self.assertFalse(self.state.pending())
        self.assertEqual(flush(self.state,self.client)['created'],0)

    def test_split_only_copy_warning_clears_and_queues_original(self):
        from nmr_register.state import serialized
        from nmr_register.split import SPLIT_SCAN_ERROR
        self.state.observe([observation(),observation('2',inode=20)])
        with self.state.db:
            split=dict(self.state.cached_observations())['src:folder/2']
            split['ScanError']=SPLIT_SCAN_ERROR
            self.state.db.execute('UPDATE observations SET payload=? WHERE key=?',(serialized(split),split['AcquisitionKey']))
            self.state.db.execute('DELETE FROM outbox')
        self.state.exclude_derived_splits()
        self.assertEqual(dict(self.state.cached_observations())['src:folder/1']['ReconciliationReason'],'')
        self.assertEqual([key for key,_ in self.state.pending()],['src:folder/1'])
        self.assertEqual(flush(self.state,self.client)['created'],1)
        self.assertEqual(flush(self.state,self.client)['created'],0)

    def test_exclusion_preserves_hard_and_ambiguous_reasons(self):
        from nmr_register.state import serialized, MATCHING_METADATA_PREFIX
        from nmr_register.split import SPLIT_EVIDENCE
        self.state.observe([observation()])
        hard='raw acquisition file replaced; '+MATCHING_METADATA_PREFIX+'src:folder/2'
        with self.state.db:
            row=dict(self.state.cached_observations())['src:folder/1'];row['ReconciliationReason']=hard
            self.state.db.execute('UPDATE observations SET payload=?',(serialized(row),))
        split=observation('2');split['ObservedAcquisitionEvidence']=SPLIT_EVIDENCE
        self.state.observe([split])
        self.assertEqual(dict(self.state.cached_observations())['src:folder/1']['ReconciliationReason'],hard)

    def test_changed_fingerprint_cannot_be_downgraded_by_another_matching_path(self):
        rows=[observation(fingerprint='old'),observation('2',fingerprint='new',inode=20)]
        self.state.observe(rows);flush(self.state,self.client)
        changed=observation(fingerprint='new',title='replacement title')
        self.state.observe([changed,rows[1]])
        reason=dict(self.state.cached_observations())['src:folder/1']['ReconciliationReason']
        self.assertIn('identity changed',reason)
        result=flush(self.state,self.client);self.assertEqual(result['conflicts'],1);self.assertTrue(self.state.pending())
        self.assertEqual(self.client.rows[0]['fields']['Title'],'source title')
    def test_replaced_raw_inode_same_fingerprint_remains_hard_blocked(self):
        rows=[observation(),observation('2',inode=20)]
        self.state.observe(rows);flush(self.state,self.client)
        changed=observation(inode=11,title='replacement title')
        self.state.observe([changed,rows[1]])
        self.assertIn('raw acquisition file replaced',dict(self.state.cached_observations())['src:folder/1']['ReconciliationReason'])
        result=flush(self.state,self.client);self.assertEqual(result['conflicts'],1);self.assertTrue(self.state.pending())
        self.assertEqual(self.client.rows[0]['fields']['Title'],'source title')
    def test_remote_serious_reason_is_not_erased_by_local_shared_metadata_warning(self):
        rows=[observation(),observation('2',inode=20)];self.state.observe(rows);flush(self.state,self.client)
        self.client.rows[0]['fields']['ReconciliationReason']='raw acquisition file replaced at existing path; reconcile before reusing assessment'
        rows[0]['Title']='do not apply';self.state.observe(rows);result=flush(self.state,self.client)
        self.assertEqual(result['conflicts'],1);self.assertTrue(self.state.pending())
        self.assertIn('raw acquisition file replaced',self.client.rows[0]['fields']['ReconciliationReason'])
        self.assertEqual(self.client.rows[0]['fields']['Title'],'source title')

    def test_only_missing_fingerprint_warning_clears_when_sufficient_evidence_arrives(self):
        from nmr_register.state import MISSING_FINGERPRINT_WARNING
        partial=observation(fingerprint='');self.state.observe([partial]);flush(self.state,self.client)
        self.assertEqual(self.client.rows[0]['fields']['ReconciliationReason'],MISSING_FINGERPRINT_WARNING)
        known=observation(fingerprint='a'*64);self.state.observe([known]);result=flush(self.state,self.client)
        self.assertEqual(result['conflicts'],0);self.assertFalse(self.state.pending())
        self.assertEqual(self.client.rows[0]['fields']['ReconciliationReason'],'')
        self.assertEqual(dict(self.state.cached_observations())['src:folder/1']['ReconciliationReason'],'')
        # A hard remote reason containing the same text is never treated as absence.
        self.client.rows[0]['fields']['ReconciliationReason']='raw acquisition replaced; '+MISSING_FINGERPRINT_WARNING
        result=flush(self.state,self.client);self.assertEqual(result['conflicts'],1)
        self.assertIn('raw acquisition replaced',self.client.rows[0]['fields']['ReconciliationReason'])

    def test_fingerprint_appearance_does_not_clear_copy_warning(self):
        rows=[observation(fingerprint=''),observation('2',fingerprint='a'*64,inode=20)]
        self.state.observe(rows);flush(self.state,self.client)
        rows[0]['ObservedFingerprint']='a'*64;self.state.observe(rows);flush(self.state,self.client)
        self.assertTrue(all('matching acquisition metadata' in r['fields']['ReconciliationReason'] for r in self.client.rows))

    def test_incomplete_workbook_identity_cannot_bind_future_root_acquisition(self):
        partial={'SourceID':'src','Folder':'','Expno':'1','MigrationIssue':'Incomplete acquisition identity; workbook source requires review','Notes':'keep partial','Sample':77}
        self.client.rows=[{'id':1,'fields':copy.deepcopy(partial)}]
        incoming=observation(fingerprint='a'*64)
        incoming.update(AcquisitionKey='src:1',ObservedFolder='')
        self.state.observe([incoming]);result=flush(self.state,self.client)
        self.assertEqual(result['created'],1);self.assertEqual(result['conflicts'],0)
        self.assertEqual(self.client.rows[0]['fields'],partial)
        self.assertEqual(self.client.rows[1]['fields']['AcquisitionKey'],'src:1')
        self.assertNotIn('Sample',self.client.rows[1]['fields'])

    def test_explicit_root_registration_still_binds_after_incomplete_marker_cleared(self):
        resolved={'SourceID':'src','Folder':'','Expno':'1','MigrationIssue':'Sample link unresolved','Notes':'true root registration'}
        self.client.rows=[{'id':1,'fields':copy.deepcopy(resolved)}]
        incoming=observation(fingerprint='a'*64);incoming.update(AcquisitionKey='src:1',ObservedFolder='')
        self.state.observe([incoming]);result=flush(self.state,self.client)
        self.assertEqual(result['created'],0);self.assertEqual(len(self.client.rows),1)
        self.assertEqual(self.client.rows[0]['fields']['AcquisitionKey'],'src:1')
        self.assertEqual(self.client.rows[0]['fields']['Notes'],'true root registration')

if __name__=='__main__':unittest.main()
