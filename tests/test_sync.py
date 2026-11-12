import copy
from http.server import BaseHTTPRequestHandler,ThreadingHTTPServer
import json
from pathlib import Path
import tempfile
import threading
import time
import unittest
from nmr_register.config import Config
from nmr_register.grist import GristClient,GristError
from nmr_register.state import State
from nmr_register.sync import flush,_conflicts


def fields(key='local:sample/1'):
    return {'AcquisitionKey':key,'ObservedSourceID':'local','ObservedFolder':'sample','ObservedExpno':'1',
            'ObservedFingerprint':'fingerprint','ObservedType':'unknown','RawPresent':False,'RawBytes':0}


class Fake:
    def __init__(self,rows=None): self.rows=copy.deepcopy(rows or []); self.posts=0;self.patches=[]
    def records(self,table): return copy.deepcopy(self.rows)
    def add_records(self,table,items):
        self.posts+=1
        for item in items: self.rows.append({'id':len(self.rows)+1,'fields':copy.deepcopy(item)})
    def update_records(self,table,items):
        self.patches+=copy.deepcopy(items)
        for item in items:
            next(r for r in self.rows if r['id']==item['id'])['fields'].update(item['fields'])


class SyncTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.state=State(Path(self.tmp.name)/'state')
    def tearDown(self): self.state.close();self.tmp.cleanup()

    def test_verified_identity_create_and_cached_blank_repair_are_idempotent(self):
        incoming=fields();self.state.observe([incoming]);client=Fake()
        flush(self.state,client)
        saved=client.rows[0]['fields']
        expected={'SourceID':'local','Folder':'sample','Expno':'1'}
        self.assertEqual({k:saved[k] for k in expected},expected)
        self.assertEqual(flush(self.state,client)['updated'],0)
        # Simulate an older collector row with unchanged, acknowledged cache.
        saved.update(SourceID='',Folder='',Expno='')
        self.assertFalse(self.state.pending())
        result=flush(self.state,client)
        self.assertEqual(result['updated'],1)
        self.assertEqual({k:saved[k] for k in expected},expected)
        self.assertEqual(flush(self.state,client)['updated'],0)
        self.assertFalse(self.state.pending())

    def test_partial_matching_identity_fills_blanks_without_overwriting_human(self):
        incoming=fields();self.state.observe([incoming])
        client=Fake([{'id':1,'fields':{**incoming,'SourceID':'local','Folder':'','Expno':None,'Notes':'keep'}}])
        flush(self.state,client)
        self.assertEqual(client.rows[0]['fields']['Folder'],'sample')
        self.assertEqual(client.rows[0]['fields']['Expno'],'1')
        self.assertEqual(client.rows[0]['fields']['Notes'],'keep')
        self.assertNotIn('SourceID',client.patches[0]['fields'])

    def test_hard_conflicts_do_not_fill_any_identity_component(self):
        for changes in ({'SourceID':'other'}, {'ObservedFingerprint':'different'},
                        {'ReconciliationReason':'possible copy or rename'},):
            with self.subTest(changes=changes):
                incoming=fields();self.state.observe([incoming])
                remote={**incoming,'SourceID':'','Folder':'','Expno':'',**changes}
                client=Fake([{'id':1,'fields':remote}]);result=flush(self.state,client)
                self.assertEqual(result['conflicts'],1)
                for name in ('SourceID','Folder','Expno'):
                    self.assertEqual(client.rows[0]['fields'][name],remote[name])
                self.assertTrue(self.state.pending())

    def test_root_identity_fills_explicit_empty_folder_and_second_sync_is_noop(self):
        incoming={**fields('local:1'),'ObservedFolder':''}
        self.state.observe([incoming]);client=Fake();flush(self.state,client)
        saved=client.rows[0]['fields']
        self.assertEqual((saved['SourceID'],saved['Folder'],saved['Expno']),('local','','1'))
        self.assertEqual(flush(self.state,client)['updated'],0)

    def test_actual_facts_replace_recorded_values_without_touching_human_assessment(self):
        human={'SourceID':'local','Folder':'sample','Expno':'1','Sample':77,'Assessment':'Use',
               'AssessmentReason':'manual','Notes':'keep','Field_MHz':850,'Instrument':'confirmed'}
        client=Fake([{'id':4,'fields':human}])
        incoming=fields();incoming['ObservedField_MHz']=600
        self.state.observe([incoming]); result=flush(self.state,client)
        self.assertEqual(result['created'],0);self.assertEqual(client.posts,0)
        for key,value in human.items():
            self.assertEqual(client.rows[0]['fields'][key],600 if key=='Field_MHz' else value)
        self.assertEqual(client.rows[0]['fields']['MetadataDifference'],'')
        self.assertFalse(client.rows[0]['fields']['ReconciliationReason'])
        self.assertTrue(all(not set(p['fields']) & (set(human)-{'Field_MHz'}) for p in client.patches))

    def test_inconclusive_copy_preserves_terminal_without_hourly_requeues(self):
        for prior in ('Completed','Aborted','Finished','Stopped'):
            for observed in ('Unknown','Setup'):
                with self.subTest(prior=prior,observed=observed):
                    incoming=fields();incoming.update(ObservedAcquisitionStatus=observed,
                        ObservedAcquisitionEvidence='audita.txt: missing audit; local copy only')
                    client=Fake([{'id':4,'fields':{**incoming,'AcquisitionStatus':prior,'Assessment':'Use'}}])
                    self.state.observe([incoming]);flush(self.state,client)
                    self.assertEqual(client.rows[0]['fields']['AcquisitionStatus'],
                        {'Finished':'Completed','Stopped':'Aborted'}.get(prior,prior))
                    self.assertEqual(client.rows[0]['fields']['Assessment'],'Use')
                    self.assertEqual(self.state.pending(),[])
                    self.assertEqual(flush(self.state,client)['updated'],0)

    def test_explicit_new_attempt_invalidates_terminal_without_assessment_write(self):
        for observed in ('Unknown','Setup'):
            incoming=fields();incoming.update(ObservedAcquisitionStatus=observed,
                ObservedAcquisitionEvidence='audita.txt: latest attempt unresolved: entry 9 records in progress; copied audit cannot establish live running')
            client=Fake([{'id':4,'fields':{**incoming,'AcquisitionStatus':'Completed','Assessment':'Use'}}])
            self.state.observe([incoming]);flush(self.state,client)
            self.assertEqual(client.rows[0]['fields']['AcquisitionStatus'],observed)
            self.assertEqual(client.rows[0]['fields']['Assessment'],'Use')
            self.assertTrue(all('Assessment' not in p['fields'] for p in client.patches))

    def test_latest_terminal_wins_and_legacy_outbox_names_translate(self):
        for observed,expected in [('Finished','Completed'),('Stopped','Aborted'),('Completed','Completed'),('Aborted','Aborted')]:
            incoming=fields();incoming['ObservedAcquisitionStatus']=observed
            client=Fake([{'id':4,'fields':{**incoming,'AcquisitionStatus':'Running'}}])
            self.state.observe([incoming]);flush(self.state,client)
            self.assertEqual(client.rows[0]['fields']['AcquisitionStatus'],expected)
            self.assertEqual(flush(self.state,client)['updated'],0)
            self.assertEqual(self.state.pending(),[])

    def test_duplicate_manual_visible_and_pending(self):
        manual={'SourceID':'local','Folder':'sample','Expno':'1','Assessment':'Use'}
        client=Fake([{'id':1,'fields':manual},{'id':2,'fields':manual}])
        self.state.observe([fields()]);result=flush(self.state,client)
        self.assertEqual(result['conflicts'],1);self.assertTrue(self.state.pending())
        self.assertTrue(all(r['fields'].get('ReconciliationReason') for r in client.rows))
        self.assertEqual(client.posts,0)

    def test_non_machine_fields_fail_closed(self):
        for name in ('Notes','SourceID','Folder','Expno'):
            incoming=fields();incoming[name]='bad'
            self.state.observe([incoming])
            with self.assertRaises(ValueError):flush(self.state,Fake())


    def test_key_and_different_manual_row_is_visible_conflict(self):
        keyrow={'id':1,'fields':{'AcquisitionKey':'local:sample/1','Assessment':'Use'}}
        manual={'id':2,'fields':{'SourceID':'local','Folder':'sample','Expno':'1','Assessment':'Use'}}
        client=Fake([keyrow,manual]);self.state.observe([fields()])
        result=flush(self.state,client)
        self.assertEqual(result['conflicts'],1);self.assertTrue(self.state.pending())
        self.assertTrue(all(r['fields'].get('ReconciliationReason') for r in client.rows))
        self.assertEqual(client.posts,0)

    def test_nominal_fields_type_synonyms_and_date_precision(self):
        from datetime import datetime,timezone
        midnight=datetime(2026,10,30,tzinfo=timezone.utc).timestamp()
        machine={'ObservedType':'r1rho','ObservedField_MHz':600.48,'ObservedAcquiredAt':midnight+3600*8}
        confirmed={'ExperimentType':'T1rho','Field_MHz':600,'AcquiredAt':midnight}
        self.assertEqual(_conflicts(machine,confirmed),[])
        confirmed['Field_MHz']=850
        self.assertIn('observed Field_MHz differs from confirmed value',_conflicts(machine,confirmed))

    def test_create_initializes_unknown_numeric_null_and_keeps_actual_zero(self):
        incoming=fields();incoming.update(ObservedField_MHz=600.48,ObservedTemperature_C=0,RawBytes=0)
        self.state.observe([incoming]);client=Fake();flush(self.state,client)
        saved=client.rows[0]['fields']
        for key in ('SnapshotpH','SnapshotNaCl_mM','SnapshotDNA_uM','SnapshotTFE_pct','SnapshotProtein_uM'):
            self.assertIn(key,saved);self.assertIsNone(saved[key])
        self.assertEqual(saved['Field_MHz'],600)
        self.assertEqual(saved['Temperature_C'],0);self.assertEqual(saved['RawBytes'],0)
        self.assertIsNone(saved['AcquiredAt'])
        self.assertNotIn('ObservedField_MHz',saved);self.assertNotIn('ObservedTemperature_C',saved)
        client.rows[0]['fields']['Temperature_C']=25
        flush(self.state,client)
        self.assertEqual(client.rows[0]['fields']['Temperature_C'],0)
        self.assertTrue(all(not {'Assessment','AssessmentReason','Sample','Notes'} & set(record['fields']) for record in client.patches))

    def test_completion_refresh_timeout_retry_replaces_status_preserves_assessment(self):
        incoming=fields()
        incoming.update(ObservedAcquisitionStatus='Unknown',
                        ObservedAcquisitionEvidence='no terminal audit record')
        human={'AcquisitionStatus':'Planned','Assessment':'Use',
               'AssessmentReason':'manual review','Notes':'keep this note'}
        remote={**incoming,**human}

        class CommitThenTimeout(Fake):
            def __init__(self,rows):
                super().__init__(rows)
                self.fail_next_patch=True
            def update_records(self,table,items):
                super().update_records(table,items)
                if self.fail_next_patch:
                    self.fail_next_patch=False
                    raise GristError('synthetic timeout after commit')

        client=CommitThenTimeout([{'id':4,'fields':remote}])
        self.state.observe([incoming],now=1000)
        with self.assertRaises(GristError):
            flush(self.state,client)
        self.assertTrue(self.state.pending())

        finished={**incoming,'ObservedAcquisitionStatus':'Finished',
                  'ObservedAcquisitionEvidence':'matching go4 terminal record'}
        self.state.observe([finished],now=1001)
        result=flush(self.state,client)

        self.assertEqual(result['created'],0)
        self.assertEqual(result['updated'],1)
        self.assertFalse(self.state.pending())
        self.assertEqual(len(client.rows),1)
        saved=client.rows[0]['fields']
        self.assertEqual(saved['AcquisitionStatus'],'Completed')
        self.assertEqual(saved['ObservedAcquisitionEvidence'],'matching go4 terminal record')
        for key,value in human.items():
            self.assertEqual(saved[key],'Completed' if key=='AcquisitionStatus' else value)
        self.assertTrue(all(not {'Assessment','AssessmentReason','Notes'}
                            & set(record['fields']) for record in client.patches))

        # The next identical scan changes only TransferStable from false to
        # true. Flush that expected machine-only stabilization before asserting
        # that a subsequent unchanged scan creates no more work.
        self.state.observe([finished],now=1002)
        result=flush(self.state,client)
        self.assertEqual(result['created'],0)
        self.assertEqual(result['updated'],1)
        self.assertEqual(client.rows[0]['fields']['AcquisitionStatus'],'Completed')
        for key,value in human.items():
            self.assertEqual(client.rows[0]['fields'][key],'Completed' if key=='AcquisitionStatus' else value)
        patch_count=len(client.patches)
        self.state.observe([finished],now=1003)
        result=flush(self.state,client)
        self.assertEqual(result['created'],0)
        self.assertEqual(result['updated'],0)
        self.assertEqual(len(client.rows),1)
        self.assertEqual(len(client.patches),patch_count)

    def test_http_timeout_after_commit_reconciles_before_post(self):
        server_state={'rows':[],'posts':0,'slow':True}
        class Handler(BaseHTTPRequestHandler):
            def log_message(self,*args): pass
            def do_GET(self):
                body=json.dumps({'records':server_state['rows']}).encode()
                self.send_response(200);self.end_headers();self.wfile.write(body)
            def do_POST(self):
                data=json.loads(self.rfile.read(int(self.headers['Content-Length'])))
                server_state['posts']+=1
                for record in data['records']:
                    server_state['rows'].append({'id':len(server_state['rows'])+1,'fields':record['fields']})
                if server_state['slow']:
                    server_state['slow']=False;time.sleep(.15)
                try:
                    self.send_response(200);self.end_headers();self.wfile.write(b'{"records":[]}')
                except (BrokenPipeError,ConnectionResetError):pass
        server=ThreadingHTTPServer(('127.0.0.1',0),Handler)
        thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
        try:
            config=Config(doc_url=f'http://127.0.0.1:{server.server_port}/doc',timeout=.05)
            client=GristClient(config,'test-local-key',sleep=lambda _:None)
            self.state.observe([fields()])
            with self.assertRaises(GristError):flush(self.state,client)
            self.assertTrue(self.state.pending())
            result=flush(self.state,client)
            self.assertEqual(server_state['posts'],1)
            self.assertEqual(len(server_state['rows']),1)
            self.assertFalse(self.state.pending());self.assertEqual(result['created'],0)
        finally:
            server.shutdown();server.server_close();thread.join()


class FieldConstraintResponse:
    headers={}
    def __init__(self,payload,status=200):self.payload=payload;self.status_code=status
    def json(self):return self.payload


class FieldConstraintSession:
    """Model the actual server requirement instead of accepting arbitrary batches."""
    def __init__(self,rows=(),timeout_method=None,timeout_number=2):
        self.headers={};self.rows=copy.deepcopy(list(rows));self.writes=[]
        self.timeout_method=timeout_method;self.timeout_number=timeout_number;self.failed=False
    def request(self,method,url,**kwargs):
        import requests
        if method=='GET':return FieldConstraintResponse({'records':copy.deepcopy(self.rows)})
        records=kwargs['json']['records'];self.writes.append((method,copy.deepcopy(records)))
        if len({frozenset(record['fields']) for record in records})>1:
            return FieldConstraintResponse({'error':method+' requires all records to have same fields'},400)
        returned=[]
        for record in records:
            if method=='POST':
                row={'id':max([r['id'] for r in self.rows]+[0])+1,'fields':copy.deepcopy(record['fields'])}
                self.rows.append(row)
            elif method=='PATCH':
                row=next(row for row in self.rows if row['id']==record['id'])
                row['fields'].update(record['fields'])
            else:raise AssertionError(method)
            returned.append({'id':row['id']})
        if method==self.timeout_method and sum(m==method for m,_ in self.writes)==self.timeout_number and not self.failed:
            self.failed=True
            raise requests.Timeout('synthetic timeout after subgroup commit')
        return FieldConstraintResponse({'records':returned})


class FieldGroupingTests(unittest.TestCase):
    def client(self,session):
        return GristClient(Config(doc_url='https://example.test/doc'),'synthetic-local-token',session=session,sleep=lambda _:None)

    def test_adapter_groups_patch_without_filling_missing_fields(self):
        session=FieldConstraintSession([{'id':1,'fields':{'RawBytes':1,'Notes':'human1'}},
                                       {'id':2,'fields':{'RawBytes':2,'Notes':'human2'}},
                                       {'id':3,'fields':{'RawBytes':3,'Notes':'human3'}}])
        result=self.client(session).update_records('Experiments',[
            {'id':1,'fields':{'RawBytes':10}},
            {'id':2,'fields':{'RawBytes':20,'RawPresent':True}},
            {'id':3,'fields':{'RawBytes':30}}])
        self.assertEqual([record['id'] for record in result['records']],[1,2,3])
        self.assertEqual(len(session.writes),2)
        self.assertEqual([r['id'] for r in session.writes[0][1]],[1,3])
        self.assertTrue(all(row['fields']['Notes']=='human'+str(row['id']) for row in session.rows))
        self.assertNotIn('RawPresent',session.rows[0]['fields'])
        self.assertNotIn('RawPresent',session.rows[2]['fields'])

    def test_post_groups_and_returns_created_ids_in_original_input_order(self):
        session=FieldConstraintSession()
        result=self.client(session).add_records('Experiments',[
            {'AcquisitionKey':'a','RawBytes':1},
            {'AcquisitionKey':'b','RawBytes':2,'RawPresent':True},
            {'AcquisitionKey':'c','RawBytes':3}])
        self.assertEqual([record['id'] for record in result['records']],[1,3,2])
        self.assertEqual([row['fields']['AcquisitionKey'] for row in session.rows],['a','c','b'])
        self.assertEqual(len(session.writes),2)

    def test_partial_post_group_commit_timeout_reconciles_without_duplicate(self):
        with tempfile.TemporaryDirectory() as directory:
            state=State(directory)
            try:
                a=fields();b=fields('local:second/2');b.update(ObservedFingerprint='second-fingerprint',ObservedFolder='second',ObservedExpno='2',ObservedPulseProgram='hsqc')
                state.observe([a,b]);session=FieldConstraintSession(timeout_method='POST');client=self.client(session)
                with self.assertRaises(GristError):flush(state,client)
                self.assertEqual(len(session.rows),2);self.assertEqual(len(state.pending()),2)
                flush(state,client)
                self.assertEqual(len(session.rows),2);self.assertFalse(state.pending())
                self.assertEqual(sum(method=='POST' for method,_ in session.writes),2)
            finally:state.close()

    def test_partial_patch_group_commit_timeout_reconciles_and_preserves_humans(self):
        with tempfile.TemporaryDirectory() as directory:
            state=State(directory)
            try:
                a=fields();b=fields('local:second/2');b.update(ObservedFingerprint='second-fingerprint',ObservedFolder='second',ObservedExpno='2')
                state.observe([a,b]);pending=state.pending()
                remote=[]
                for index,(_,payload) in enumerate(pending,1):
                    row={'id':index,'fields':{**payload,'RawBytes':-1,'Notes':'manual'}}
                    if index==2:row['fields']['RawPresent']=True
                    remote.append(row)
                session=FieldConstraintSession(remote,timeout_method='PATCH');client=self.client(session)
                with self.assertRaises(GristError):flush(state,client)
                self.assertEqual(len(state.pending()),2)
                self.assertTrue(all(row['fields']['Notes']=='manual' for row in session.rows))
                flush(state,client)
                self.assertFalse(state.pending())
                self.assertEqual(sum(method=='PATCH' for method,_ in session.writes),2)
                self.assertTrue(all('Notes' not in record['fields'] for _,records in session.writes for record in records))
            finally:state.close()


class MetadataRefreshTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.state=State(self.tmp.name)
    def tearDown(self):self.state.close();self.tmp.cleanup()

    def test_metadata_moves_from_legacy_reason_and_clears_after_human_correction_without_raw_changes(self):
        observed={**fields(),'PulseProgram':'original','ObservedPulseProgram':'original.ts3',
                  'Temperature_C':25,'ObservedTemperature_C':24.948,
                  'ReconciliationReason':'observed PulseProgram differs from confirmed value; observed Temperature_C differs from confirmed value'}
        client=Fake([{'id':1,'fields':observed}])
        self.assertFalse(self.state.pending())
        result=flush(self.state,client)
        saved=client.rows[0]['fields']
        self.assertEqual(result['metadata_warnings'],0)
        self.assertEqual(saved['ReconciliationReason'],'')
        self.assertEqual(saved.get('MetadataDifference',''),'')
        self.assertEqual(saved['Temperature_C'],25);self.assertEqual(saved['PulseProgram'],'original')
        self.state.observe([{**fields(),'ObservedTemperature_C':24.948,'ObservedPulseProgram':'original.ts3'}])
        flush(self.state,client)
        self.assertEqual(saved.get('MetadataDifference',''),'')
        writes=len(client.patches);flush(self.state,client)
        self.assertEqual(len(client.patches),writes)
        self.assertEqual(client.posts,0)

    def test_hsqc_alias_and_metadata_cleanup_never_clear_identity_or_mixed_reason(self):
        rows=[]
        for index,reason in enumerate(['observed ExperimentType differs from confirmed value',
             'matching acquisition metadata at another path; possible copy or rename: local:copy/1',
             'raw acquisition file replaced; observed ExperimentType differs from confirmed value'],1):
            rows.append({'id':index,'fields':{**fields(),'ExperimentType':'15N_HSQC','ObservedType':'hsqc',
                'ReconciliationReason':reason}})
        client=Fake(rows);flush(self.state,client)
        self.assertEqual(client.rows[0]['fields']['ReconciliationReason'],'')
        self.assertEqual(client.rows[1]['fields']['ReconciliationReason'],rows[1]['fields']['ReconciliationReason'])
        self.assertEqual(client.rows[2]['fields']['ReconciliationReason'],rows[2]['fields']['ReconciliationReason'])
        self.assertTrue(all(not row['fields'].get('MetadataDifference') for row in client.rows))

    def test_metadata_warning_does_not_become_identity_ambiguity_on_pending_update(self):
        incoming=fields();incoming['ObservedTemperature_C']=24.948
        self.state.observe([incoming])
        client=Fake([{'id':1,'fields':{**incoming,'Temperature_C':25,'Assessment':'Use'}}])
        flush(self.state,client)
        saved=client.rows[0]['fields']
        self.assertEqual(saved['ReconciliationReason'],'')
        self.assertEqual(saved.get('MetadataDifference',''),'')
        self.assertEqual(saved['Temperature_C'],24.948)
        self.assertEqual(saved['Assessment'],'Use')

    def test_cli_attachment_failure_publishes_error_and_restores_previous_full_success(self):
        from contextlib import redirect_stderr,redirect_stdout
        import io
        from unittest.mock import patch
        from nmr_register.cli import main
        self.state.set_meta('last_sync',123)
        config=Config(state_dir=self.tmp.name,sources=[])
        class Client:
            def __init__(self):self.status=[]
            def records(self,table):return copy.deepcopy(self.status) if table=='ScanStatus' else []
            def add_records(self,table,records):
                if table=='ScanStatus':self.status=[{'id':1,'fields':copy.deepcopy(records[0])}]
            def update_records(self,table,records):
                if table=='ScanStatus':self.status[0]['fields'].update(records[0]['fields'])
        client=Client()
        with patch('nmr_register.cli.load_config',return_value=config),patch('nmr_register.cli._client',return_value=client),\
             patch('nmr_register.results.attach_results',side_effect=GristError('result upload unavailable')),\
             redirect_stdout(io.StringIO()),redirect_stderr(io.StringIO()):
            result=main(['sync'])
        self.assertEqual(result,1)
        self.assertEqual(client.status[0]['fields']['LastSyncAt'],123)
        self.assertEqual(client.status[0]['fields']['LastError'],'result upload unavailable')
        self.assertEqual(self.state.get_meta('last_sync'),123)
        self.assertEqual(self.state.get_meta('last_error'),'result upload unavailable')
