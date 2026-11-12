import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch
from nmr_register.bruker import observe
from nmr_register.config import Source,Config,save_config,load_config
from nmr_register.discovery import scan_source
from nmr_register.state import State


def acquisition(root,folder='sample',expno='1',date=123456):
    path=Path(root)/folder/expno
    path.mkdir(parents=True)
    (path/'acqus').write_text(f'##$DATE_START= {date}\n##$PULPROG= <unfamiliar_custom_pulse>\n##$BF1= 600\n')
    return path


class CollectorTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.root=Path(self.tmp.name)
    def tearDown(self): self.tmp.cleanup()

    def test_prune_before_inspection_and_no_symlinks(self):
        acquisition(self.root)
        acquisition(self.root,'211@950')
        outside=self.root/'outside'
        acquisition(outside)
        (self.root/'link').symlink_to(outside,target_is_directory=True)
        original=os.scandir
        seen=[]
        def guarded(path):
            self.assertNotIn('211@950',str(path))
            self.assertNotIn('/link',str(path))
            seen.append(str(path))
            return original(path)
        with patch('nmr_register.discovery.os.scandir',side_effect=guarded):
            rows,errors=scan_source(Source('local',str(self.root)))
        self.assertFalse(errors)
        self.assertEqual({r['ObservedFolder'] for r in rows},{'sample','outside/sample'})

    def test_partial_unknown_is_observation_not_completion(self):
        path=acquisition(self.root)
        fields=observe(path)
        self.assertEqual(fields['ObservedType'],'unknown')
        self.assertFalse(fields['RawPresent'])
        (path/'ser').write_bytes(b'partial')
        self.assertTrue(observe(path)['RawPresent'])
        self.assertNotIn('AcquisitionStatus',observe(path))

    def test_incomplete_identity_does_not_mark_all_partial_as_copies(self):
        for name in ('a','b'):
            path=self.root/name/'1'; path.mkdir(parents=True)
            (path/'ser').write_bytes(b'partial')
        records,_=scan_source(Source('local',str(self.root)))
        state=State(self.root/'state')
        state.observe(records)
        self.assertTrue(all('fingerprint unavailable' in r['ReconciliationReason'] and 'copy' not in r['ReconciliationReason'] for _,r in state.pending()))
        state.close()

    def test_mutable_audit_not_replacement_and_same_size_write_not_stable(self):
        path=acquisition(self.root)
        raw=path/'ser';raw.write_bytes(b'aaaa')
        state=State(self.root/'state')
        def scan(): return scan_source(Source('local',str(self.root)))[0]
        state.observe(scan()); state.observe(scan())
        self.assertTrue(state.pending()[0][1]['TransferStable'])
        (path/'audita.txt').write_text('processing changed')
        raw.write_bytes(b'bbbb')
        state.observe(scan())
        payload=state.pending()[0][1]
        self.assertFalse(payload['TransferStable'])
        self.assertFalse(payload['ReconciliationReason'])
        self.assertNotIn('_RawFileStats',payload)
        state.close()

    def test_copy_and_replacement_flag_no_assessment_inheritance(self):
        path=acquisition(self.root); (path/'fid').write_bytes(b'raw')
        state=State(self.root/'state')
        state.observe(scan_source(Source('local',str(self.root)))[0])
        shutil.copytree(path,self.root/'copy'/'1')
        state.observe(scan_source(Source('local',str(self.root)))[0])
        self.assertTrue(all('copy or rename' in r['ReconciliationReason'] for _,r in state.pending()))
        (path/'acqus').write_text('##$DATE_START= 999999\n##$PULPROG= <newpulse>\n')
        state.observe(scan_source(Source('local',str(self.root)))[0])
        rows=dict(state.pending())
        self.assertIn('identity changed',rows['local:sample/1']['ReconciliationReason'])
        self.assertNotIn('Assessment',rows['local:sample/1'])
        state.close()

    def test_change_only_after_stability_and_reopen_durable(self):
        acquisition(self.root)
        state=State(self.root/'state'); records=scan_source(Source('local',str(self.root)))[0]
        state.observe(records,now=1);state.acknowledge([key for key,_ in state.pending()])
        state.observe(records,now=2);state.acknowledge([key for key,_ in state.pending()])
        self.assertEqual(state.observe(records,now=3),0)
        state.close();state=State(self.root/'state')
        self.assertEqual(state.status()['observations'],1)
        self.assertEqual(state.pending(),[])
        row=state.db.execute('SELECT last_seen FROM observations').fetchone()
        self.assertEqual(row[0],3)
        state.close()

    def test_source_targets_and_offline_root(self):
        acquisition(self.root,'first');acquisition(self.root,'second')
        rows,errors=scan_source(Source('local',str(self.root),targets=['first/1']))
        self.assertFalse(errors);self.assertEqual(len(rows),1)
        self.assertEqual(rows[0]['ObservedFolder'],'first')
        rows,errors=scan_source(Source('offline',str(self.root/'absent')))
        self.assertFalse(rows);self.assertTrue(errors)

    def test_process_lock(self):
        a=State(self.root/'state'); b=State(self.root/'state')
        with a.lock():
            with self.assertRaises(RuntimeError):
                with b.lock(): pass
        a.close();b.close()

    def test_existing_json_credential_reused_only_for_matching_document(self):
        from nmr_register.credentials import read_token
        home=self.root/'user';credential=home/'.config/nmr-register/grist.json'
        credential.parent.mkdir(parents=True)
        credential.write_text(json.dumps({'api_key':'test-secret-only','document_url':'https://docs.getgrist.com/doc/slug'}))
        credential.chmod(0o600)
        with patch('nmr_register.credentials.Path.home',return_value=home):
            self.assertEqual(read_token(home/'state','https://docs.getgrist.com/doc/other-slug'),'test-secret-only')
            with self.assertRaisesRegex(ValueError,'host and ID'):
                read_token(home/'state','https://docs.getgrist.com/other/slug')
            with self.assertRaisesRegex(ValueError,'host and ID'):
                read_token(home/'state','https://other.example/doc/slug')
            credential.chmod(0o644)
            with self.assertRaisesRegex(ValueError,'0600'):
                read_token(home/'state','https://docs.getgrist.com/doc/slug')
        self.assertFalse((home/'state'/'credential').exists())

    def test_existing_credential_symlink_rejected(self):
        from nmr_register.credentials import read_token
        actual=self.root/'actual.json'
        actual.write_text(json.dumps({'api_key':'test-secret-only','document_url':'https://docs.getgrist.com/doc/slug'}))
        actual.chmod(0o600)
        link=self.root/'link.json';link.symlink_to(actual)
        with self.assertRaisesRegex(ValueError,'regular file'):
            read_token(self.root/'state','https://docs.getgrist.com/doc/slug',link)


    def test_observed_temperature_kelvin_conversion_keeps_zero(self):
        path=acquisition(self.root)
        for raw,expected in [('273.15',0),('298.15',25),('0',None),('-20',None),('Infinity',None),('1001',None),('unknown',None)]:
            (path/'acqus').write_text('##$TE= '+raw+'\n')
            fields=observe(path)
            with self.subTest(raw=raw):
                if expected is None:self.assertIsNone(fields['ObservedTemperature_C'])
                else:self.assertEqual(fields['ObservedTemperature_C'],expected)
