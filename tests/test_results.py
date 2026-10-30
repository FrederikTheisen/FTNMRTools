import copy
import csv
import hashlib
from pathlib import Path
import tempfile
import unittest
from nmr_register.grist import GristError
from nmr_register.results import select_results,attach_results,read_numbers_bytes
from nmr_register.state import State


class ResultClient:
    def __init__(self):
        self.rows=[{'id':9,'fields':{'AcquisitionKey':'local:sample/1','Assessment':'Use','Notes':'human'}}]
        self.files=[];self.uploads=0;self.fail_upload=False;self.fail_patch=False
    def records(self,table):return copy.deepcopy(self.rows)
    def attachments(self):return [{'id':f['id'],'fields':{'fileName':f['name']}} for f in self.files]
    def upload(self,name,content):
        self.uploads+=1
        self.files.append({'id':len(self.files)+1,'name':name,'content':content})
        if self.fail_upload:
            self.fail_upload=False;raise GristError('synthetic timeout after upload committed')
        return [self.files[-1]['id']]
    def download_attachment(self,attachment):return next(f['content'] for f in self.files if f['id']==attachment)
    def update_records(self,table,records):
        for record in records:self.rows[0]['fields'].update(record['fields'])
        if self.fail_patch:
            self.fail_patch=False;raise GristError('synthetic timeout after patch committed')


class ResultsTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name);self.state=State(self.root/'state')
        self.content=b'assignment,rate,warnings\nA1,2.5,"line1\n# warning line2"\n# source_folder: sample\n# expno: 1\n'
        self.path=self.root/'CURRENT_numbers.csv';self.path.write_bytes(self.content)
        self.manifest=self.root/'current.tsv'
        self.write_manifest()
    def tearDown(self):self.state.close();self.tmp.cleanup()
    def write_manifest(self,acquisition='sample/1',sha=None):
        with self.manifest.open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=['acquisition','numbers_file','numbers_sha256'],delimiter='\t')
            writer.writeheader();writer.writerow({'acquisition':acquisition,'numbers_file':str(self.path),
                'numbers_sha256':sha or hashlib.sha256(self.content).hexdigest()})
    def select(self):return select_results(self.manifest,'local')

    def test_original_bytes_and_multiline_metadata_preserved(self):
        selection=self.select()[0]
        self.assertEqual(selection.content,self.content)
        _,rows,metadata=read_numbers_bytes(selection.content)
        self.assertIn('# warning line2',rows[0]['warnings']);self.assertEqual(metadata['expno'],'1')
        client=ResultClient();attach_results(self.state,client,[selection])
        self.assertEqual(client.files[0]['content'],self.content)
        self.assertEqual(client.rows[0]['fields']['Assessment'],'Use')
        self.assertEqual(client.rows[0]['fields']['Notes'],'human')
        attach_results(self.state,client,[selection]);self.assertEqual(client.uploads,1)

    def test_hash_identity_and_size_are_verified(self):
        self.write_manifest(sha='0'*64)
        with self.assertRaisesRegex(ValueError,'SHA256'):self.select()
        self.write_manifest(acquisition='different/1')
        with self.assertRaisesRegex(ValueError,'source identity'):self.select()
        self.write_manifest()
        with self.assertRaisesRegex(ValueError,'size bound'):select_results(self.manifest,'local',max_bytes=10)

    def test_upload_timeout_reconciles_one_attachment(self):
        client=ResultClient();client.fail_upload=True
        with self.assertRaises(GristError):attach_results(self.state,client,self.select())
        self.assertEqual(client.uploads,1)
        result=attach_results(self.state,client)
        self.assertEqual(result['attached'],1);self.assertEqual(client.uploads,1)

    def test_patch_timeout_reconciles_hash_without_duplicate_upload(self):
        client=ResultClient();client.fail_patch=True
        with self.assertRaises(GristError):attach_results(self.state,client,self.select())
        result=attach_results(self.state,client)
        self.assertEqual(result['attached'],0);self.assertEqual(client.uploads,1)

    def test_explicit_current_selection_replaces_acknowledged_previous_only(self):
        client=ResultClient();attach_results(self.state,client,self.select())
        old_hash=client.rows[0]['fields']['ResultSHA256']
        self.content=self.content.replace(b'A1,2.5',b'A1,3.5')
        self.path.write_bytes(self.content);self.write_manifest()
        attach_results(self.state,client,self.select())
        self.assertEqual(client.uploads,2)
        self.assertNotEqual(client.rows[0]['fields']['ResultSHA256'],old_hash)
        self.assertEqual(client.rows[0]['fields']['ResultCSV'],['L',2])
        client.rows[0]['fields']['ResultSHA256']='c'*64
        self.content=self.content.replace(b'A1,3.5',b'A1,4.5')
        self.path.write_bytes(self.content);self.write_manifest()
        with self.assertRaisesRegex(GristError,'different result hash'):
            attach_results(self.state,client,self.select())
        self.assertEqual(client.rows[0]['fields']['ResultSHA256'],'c'*64)
        self.assertEqual(client.uploads,2)

    def test_state_bound_to_document(self):
        self.state.bind_document('https://example.test/api/docs/one')
        with self.assertRaisesRegex(ValueError,'different Grist document'):
            self.state.bind_document('https://example.test/api/docs/two')

    def test_metadata_difference_allows_exact_verified_result_attachment(self):
        client=ResultClient()
        client.rows[0]['fields']['MetadataDifference']='observed Temperature_C differs from confirmed value'
        result=attach_results(self.state,client,self.select())
        self.assertEqual(result['attached'],1)
        self.assertEqual(client.rows[0]['fields']['Assessment'],'Use')
        self.assertIn('MetadataDifference',client.rows[0]['fields'])

    def test_conflict_appearing_during_upload_blocks_linking(self):
        client=ResultClient();original=client.upload
        def conflicting_upload(name,content):
            result=original(name,content)
            client.rows[0]['fields']['ReconciliationReason']='raw acquisition replaced'
            return result
        client.upload=conflicting_upload
        with self.assertRaisesRegex(GristError,'became conflicted'):
            attach_results(self.state,client,self.select())
        self.assertNotIn('ResultCSV',client.rows[0]['fields'])
        self.assertEqual(self.state.db.execute('SELECT count(*) FROM result_queue').fetchone()[0],1)

    def test_existing_different_hash_refuses_overwrite(self):
        client=ResultClient();client.rows[0]['fields']['ResultSHA256']='a'*64
        with self.assertRaisesRegex(GristError,'different result hash'):attach_results(self.state,client,self.select())
        self.assertEqual(client.uploads,0)

    def test_target_filter_before_unrelated_file_reads(self):
        self.write_manifest(acquisition='unrelated/1')
        self.path.unlink()
        self.assertEqual(select_results(self.manifest,'local',targets=['sample/1']),[])
