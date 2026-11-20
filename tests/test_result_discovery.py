import copy
import csv
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from contextlib import redirect_stdout
from nmr_register.config import Config, Source
from nmr_register.discovery import scan_source
from nmr_register.result_discovery import prepare_results
from nmr_register.results import attach_results, queue_results
from nmr_register.state import State
from nmr_register.sync import run
from nmr_register.cli import main


class Client:
    def __init__(self):
        self.rows=[];self.status=[];self.files=[];self.uploads=0
    def records(self,table):return copy.deepcopy(self.status if table=='ScanStatus' else self.rows)
    def add_records(self,table,records):
        target=self.status if table=='ScanStatus' else self.rows
        target.extend({'id':len(target)+i+1,'fields':copy.deepcopy(r)} for i,r in enumerate(records))
    def update_records(self,table,records):
        target=self.status if table=='ScanStatus' else self.rows
        for r in records:next(v for v in target if v['id']==r['id'])['fields'].update(r['fields'])
    def attachments(self):return [{'id':r['id'],'fields':{'fileName':r['name']}} for r in self.files]
    def upload(self,name,content):
        self.uploads+=1;self.files.append(dict(id=len(self.files)+1,name=name,content=content));return [self.files[-1]['id']]
    def download_attachment(self,ident):return next(r['content'] for r in self.files if r['id']==ident)


class DiscoveryTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name);self.data=self.root/'data';self.data.mkdir()
        self.source=Source('local',str(self.data))
        self.config=Config(sources=[self.source],state_dir=str(self.root/'state'))
        self.state=State(self.config.home)
        self.path=self.acquisition('sample','1');self.numbers(self.path)
    def tearDown(self):self.state.close();self.tmp.cleanup()
    def acquisition(self,folder,expno):
        p=self.data/folder/expno;p.mkdir(parents=True)
        (p/'acqus').write_text(f'##$DATE_START= {100+int(expno)}\n##$PULPROG= <unfamiliar>\n##$BF1= 600\n')
        return p
    def numbers(self,path,value='2.5',name='CURRENT_numbers.csv',folder=None):
        p=path/name;p.write_text('assignment,rate\nA1,'+value+'\n# source_folder: '+(folder or path.parent.name)+'\n# expno: '+path.name+'\n');return p
    def observations(self):return scan_source(self.source)[0]
    def plan(self):return prepare_results(self.config,self.observations())
    def manifest(self,path):
        m=self.root/'manifest.tsv'
        with m.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=['acquisition','numbers_file','numbers_sha256'],delimiter='\t');w.writeheader();w.writerow(dict(acquisition='sample/1',numbers_file=str(path),numbers_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        self.config.results_manifest=str(m);self.config.result_source_id='local';return m

    def test_sync_uploads_new_changed_then_unchanged_and_keeps_human_fields(self):
        client=Client()
        first=run(self.config,self.state,client)
        self.assertEqual((first['results_discovered'],first['results_selected'],first['results_changed'],first['attached']),(1,1,1,1))
        client.rows[0]['fields'].update(Assessment='Use',Notes='user note')
        self.numbers(self.path,value='3.5')
        second=run(self.config,self.state,client)
        self.assertEqual((second['results_changed'],second['attached']),(1,1))
        third=run(self.config,self.state,client)
        self.assertEqual((third['results_changed'],third['attached'],client.uploads),(0,0,2))
        self.assertEqual(client.rows[0]['fields']['Assessment'],'Use');self.assertEqual(client.rows[0]['fields']['Notes'],'user note')

    def test_cli_dry_run_previews_results_without_network_or_queue_changes(self):
        before=self.state.status();out=io.StringIO()
        with patch('nmr_register.cli.load_config',return_value=self.config),patch('nmr_register.cli._client') as network,redirect_stdout(out):
            self.assertEqual(main(['sync','--dry-run']),0)
        report=json.loads(out.getvalue());self.assertEqual(report['results_changed'],1)
        self.assertEqual(report['results_selected'],1);network.assert_not_called()
        self.assertEqual(self.state.status(),before)

    def test_verified_local_update_supersedes_stale_manifest_hash_at_same_path(self):
        p=self.path/'CURRENT_numbers.csv';self.manifest(p);self.numbers(self.path,value='3.5')
        plan=self.plan();self.assertEqual(plan.issues,[])
        self.assertEqual(plan.selections[0].content,p.read_bytes())

    def test_archived_fallback_kept_and_new_local_supersedes_it(self):
        archived=self.root/'archive_numbers.csv';archived.write_bytes((self.path/'CURRENT_numbers.csv').read_bytes())
        self.manifest(archived)
        self.assertEqual(self.plan().selections[0].path,str(archived))
        (self.path/'CURRENT_numbers.csv').unlink()
        self.assertEqual(self.plan().selections[0].path,str(archived))
        self.numbers(self.path,value='4.5')
        self.assertEqual(self.plan().selections[0].path,str(self.path/'CURRENT_numbers.csv'))

    def test_ambiguous_local_versions_hold_stale_queue_but_other_results_attach(self):
        old=self.plan().selections[0];queue_results(self.state,[old])
        self.numbers(self.path,value='3.5');self.numbers(self.path,value='4.5',name='alternative_numbers.csv')
        other=self.acquisition('sample','2');self.numbers(other)
        plan=self.plan();self.assertEqual(len(plan.issues),1);self.assertIn('Multiple different',plan.issues[0]['reason'])
        client=Client();client.rows=[{'id':1,'fields':{'AcquisitionKey':'local:sample/1'}},{'id':2,'fields':{'AcquisitionKey':'local:sample/2'}}]
        result=attach_results(self.state,client,plan.selections,discovery_issues=plan.issues,eligible_keys=plan.eligible_keys)
        self.assertEqual((result['attached'],result['pending_results'],result['blocked_results']),(1,1,1))
        self.assertNotIn('ResultCSV',client.rows[0]['fields'])

    def test_identical_local_files_are_not_conflicting(self):
        self.numbers(self.path,name='identical_numbers.csv');plan=self.plan()
        self.assertFalse(plan.issues);self.assertEqual(len(plan.selections),1)

    def test_wrong_source_footer_is_reported_without_stale_manifest_fallback(self):
        archived=self.root/'old.csv';archived.write_bytes((self.path/'CURRENT_numbers.csv').read_bytes());self.manifest(archived)
        self.numbers(self.path,folder='different')
        plan=self.plan();self.assertFalse(plan.selections);self.assertIn('source identity',plan.issues[0]['reason'])

    def test_size_bound_and_symlink_result_are_not_uploaded(self):
        self.config.max_result_bytes=10
        self.assertIn('size',self.plan().issues[0]['reason'])
        self.config.max_result_bytes=2_000_000
        p=self.path/'CURRENT_numbers.csv';copy=self.root/'actual.csv';p.rename(copy);p.symlink_to(copy)
        self.assertFalse(self.plan().selections);self.assertTrue(self.plan().issues)

    def test_scanner_scope_protects_targets_splits_archives_and_exclusions(self):
        self.acquisition('211@950','2');self.acquisition('omitted','2');self.acquisition('other','2')
        archive=self.path/'archive';archive.mkdir();self.numbers(archive,value='9.5')
        self.source.exclude.append('omitted');self.source.targets=['sample/1']
        link=self.data/'linked';link.symlink_to(self.path.parent,target_is_directory=True)
        original=Path.iterdir
        def guarded(p):
            self.assertNotIn('211@950',str(p));self.assertNotIn('/omitted',str(p));self.assertNotIn('/linked',str(p));self.assertNotIn('/archive',str(p))
            return original(p)
        with patch('pathlib.Path.iterdir',guarded):plan=self.plan()
        self.assertEqual(len(plan.selections),1);self.assertEqual(plan.discovered,1)

    def test_derived_split_observation_is_not_a_result_source(self):
        from nmr_register.split import SPLIT_EVIDENCE
        records=self.observations();records[0]['ObservedAcquisitionEvidence']=SPLIT_EVIDENCE
        plan=prepare_results(self.config,records)
        self.assertEqual(plan.selections,[]);self.assertEqual(plan.eligible_keys,set())

    def test_invalid_archived_fallback_is_reported_per_acquisition(self):
        archived=self.root/'old.csv';archived.write_bytes((self.path/'CURRENT_numbers.csv').read_bytes());self.manifest(archived)
        (self.path/'CURRENT_numbers.csv').unlink();archived.write_bytes(b'bad replacement')
        other=self.acquisition('sample','2');self.numbers(other)
        plan=self.plan();self.assertEqual(len(plan.selections),1)
        self.assertEqual(plan.selections[0].acquisition_key,'local:sample/2')
        self.assertIn('SHA256',plan.issues[0]['reason'])

    def test_stale_bad_queued_file_does_not_prevent_other_uploads(self):
        queue_results(self.state,self.plan().selections);self.numbers(self.path,value='4.5')
        other=self.acquisition('sample','2');self.numbers(other)
        selection=next(s for s in self.plan().selections if s.acquisition_key=='local:sample/2')
        client=Client();client.rows=[{'id':1,'fields':{'AcquisitionKey':'local:sample/2'}}]
        result=attach_results(self.state,client,[selection])
        self.assertEqual((result['attached'],result['blocked_results']),(1,1))
        self.assertIn('SHA256',result['result_issues'][0]['reason'])

    def test_queued_result_outside_current_scope_is_held_before_file_read(self):
        queue_results(self.state,self.plan().selections)
        (self.path/'CURRENT_numbers.csv').unlink()
        client=Client()
        result=attach_results(self.state,client,[],eligible_keys=set())
        self.assertEqual(result['blocked_results'],1);self.assertEqual(client.uploads,0)
        self.assertIn('current scan scope',result['result_issues'][0]['reason'])

    def test_discovery_never_clears_acquisition_reconciliation(self):
        client=Client();client.rows=[{'id':1,'fields':{'AcquisitionKey':'local:sample/1','ReconciliationReason':'copy needs review'}}]
        result=attach_results(self.state,client,self.plan().selections)
        self.assertEqual(result['blocked_results'],1);self.assertEqual(client.uploads,0)
        self.assertEqual(client.rows[0]['fields']['ReconciliationReason'],'copy needs review')
