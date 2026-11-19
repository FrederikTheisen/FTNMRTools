"""Setup promotion and mutable completion facts never authorize raw replacements."""
import copy
from datetime import datetime, timezone
import os
from pathlib import Path
import tempfile
import unittest
from nmr_register.bruker import observe, parameters
from nmr_register.state import State
from nmr_register.sync import flush,_identity_conflicts

class Fake:
    def __init__(self):self.rows=[];self.patches=[];self.fail_after_commit=False
    def records(self,table):return copy.deepcopy(self.rows)
    def add_records(self,table,rows):
        for row in rows:self.rows.append({'id':len(self.rows)+1,'fields':copy.deepcopy(row)})
    def update_records(self,table,rows):
        self.patches.extend(copy.deepcopy(rows))
        for row in rows:next(r for r in self.rows if r['id']==row['id'])['fields'].update(row['fields'])
        if self.fail_after_commit:self.fail_after_commit=False;raise RuntimeError('uncertain transport after commit')

def stamp(n):return datetime.fromtimestamp(n,timezone.utc).strftime('%Y-%m-%dT%H:%M:%S.000+0000')
def audit_record(number,when,body):return f'( {number},<{stamp(when)}>,<user>,<instrument>,<go4>,<TopSpin 4.4.1>,<{body}>)\n'
def audit(*records):return '##TITLE= Audit trail, TopSpin 4.4.1\n##AUDIT TRAIL=\n'+''.join(records)+'##END=\n'
def terminal(start,end,uuid='e1c20164-3bdc-11f0-a2ab-0030d625a2a2'):
    return "created by zg on data set with UUID '"+uuid+"'\nstarted at "+stamp(start)+',\ncompleted at '+stamp(end)

class AcquisitionTransitionTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.path=self.root/'folder'/'1';self.path.mkdir(parents=True)
        self.state=State(self.root/'state');self.addCleanup(self.state.close);self.client=Fake()
    def params(self,start=1700000000,end=1700000010,pulse='zg',td=1024):
        (self.path/'acqus').write_text(f'##$DATE_START= {start}\n##$DATE= {end}\n##$PULPROG= <{pulse}>\n##$TD= {td}\n##$BF1= 600.13\n##$SFO1= 600.13\n##$NUC1= <1H>\n')
    def completed_audit(self):
        end=int(parameters((self.path/'acqus').read_text())['DATE'])
        (self.path/'audita.txt').write_text(audit(audit_record(39,end,terminal(1700000000,end))))
    def observation(self):
        return dict(observe(self.path),AcquisitionKey='src:folder/1',ObservedSourceID='src',ObservedFolder='folder',ObservedExpno='1')
    def sync(self):self.state.observe([self.observation()]);return flush(self.state,self.client)
    def test_setup_to_raw_pulse_and_start_changes_sync_without_new_row(self):
        self.params();self.sync();row=self.client.rows[0]['fields']
        self.assertEqual(row['AcquisitionStatus'],'Setup')
        row.update(Sample=30,Study='Electrostatics',Notes='keep')
        self.params(start=1700001000,end=1700001010,pulse='hsqctretf3gpsitc3d.3')
        (self.path/'ser').write_bytes(b'fixture raw');self.sync()
        self.assertEqual(len(self.client.rows),1);self.assertEqual(row['ExperimentType'],'T1rho')
        self.assertEqual(row['PulseProgram'],'hsqctretf3gpsitc3d.3');self.assertEqual(row['Notes'],'keep')
        self.assertEqual(row['Sample'],30)
        self.assertFalse(row['ReconciliationReason']);self.assertEqual(self.state.pending(),[])
        self.sync();self.assertEqual(len(self.client.rows),1)
    def test_zero_length_template_file_may_be_replaced_by_actual_raw(self):
        self.params();(self.path/'fid').write_bytes(b'');self.sync()
        new=self.path/'new';new.write_bytes(b'actual fixture');os.replace(new,self.path/'fid')
        self.params(start=1700000100,end=1700000110,pulse='nico_nnh_tr_noSE');self.sync()
        self.assertFalse(self.client.rows[0]['fields']['ReconciliationReason'])
    def test_same_inode_end_date_update_preserves_review_and_attachment(self):
        self.params(pulse='t1');(self.path/'ser').write_bytes(b'fixture');self.completed_audit();self.sync()
        saved=self.client.rows[0]['fields'];saved.update(Assessment='Use',QCFID='Clean',ResultCSV=['L',42],Notes='keep')
        self.params(end=1700000020,pulse='t1');self.completed_audit()
        with (self.path/'ser').open('ab') as stream:stream.write(b'growing')
        self.sync()
        self.assertFalse(saved['ReconciliationReason']);self.assertEqual(saved['ResultCSV'],['L',42]);self.assertEqual(saved['Assessment'],'Use')
        self.assertEqual(saved['QCFID'],'Clean')
    def test_protocol_finalization_needs_matching_verified_progress_anchor(self):
        self.params();(self.path/'ser').write_bytes(b'fixture')
        progress=audit_record(38,1700000000,'acquisition in progress')
        (self.path/'audita.txt').write_text(audit(progress));self.sync()
        self.params(end=1700000010,pulse='hsqctretf3gpsitc3d.3',td=4096)
        (self.path/'audita.txt').write_text(audit(progress,audit_record(39,1700000010,terminal(1700000000,1700000010))))
        self.sync();saved=self.client.rows[0]['fields']
        self.assertFalse(saved['ReconciliationReason']);self.assertEqual(saved['ExperimentType'],'T1rho');self.assertEqual(saved['AcquisitionStatus'],'Completed')
    def test_protocol_change_without_progress_proof_remains_blocked_and_review_untouched(self):
        self.params();(self.path/'ser').write_bytes(b'fixture');self.sync()
        row=self.client.rows[0]['fields'];row.update(Assessment='Use',QCSNR='Adequate',ResultCSV=['L',42])
        self.params(pulse='t1',td=4096);self.sync()
        self.assertIn('identity changed',row['ReconciliationReason']);self.assertEqual(row['PulseProgram'],'zg')
        self.assertEqual(row['Assessment'],'Use');self.assertEqual(row['ResultCSV'],['L',42]);self.assertTrue(self.state.pending())
    def test_restart_and_inode_replacement_remain_hard_across_retries(self):
        self.params(pulse='t1');(self.path/'ser').write_bytes(b'fixture');self.completed_audit();self.sync();row=self.client.rows[0]['fields'];row['Assessment']='Do not use'
        self.params(start=1700001000,end=1700001010,pulse='t1');self.sync()
        self.assertIn('identity changed',row['ReconciliationReason']);self.sync();self.assertTrue(row['ReconciliationReason']);self.assertEqual(row['Assessment'],'Do not use')
        replacement=self.path/'new';replacement.write_bytes(b'other');os.replace(replacement,self.path/'ser');self.sync()
        self.assertIn('raw acquisition file replaced',row['ReconciliationReason'])
    def test_setup_with_existing_review_does_not_silently_carry_it_to_new_raw(self):
        self.params();self.sync();row=self.client.rows[0]['fields'];row.update(Assessment='Use',QCFID='Clean')
        self.params(start=1700001000,pulse='t1');(self.path/'ser').write_bytes(b'fixture');self.sync()
        self.assertIn('existing review',row['ReconciliationReason']);self.assertEqual(row['PulseProgram'],'zg');self.assertEqual(row['Assessment'],'Use')
    def test_timeout_after_transition_commit_retries_without_duplicate_or_context_leak(self):
        self.params();self.sync();self.params(start=1700001000,pulse='t1');(self.path/'ser').write_bytes(b'fixture')
        self.state.observe([self.observation()]);self.client.fail_after_commit=True
        with self.assertRaisesRegex(RuntimeError,'uncertain'):flush(self.state,self.client)
        self.assertTrue(self.state.pending());flush(self.state,self.client)
        self.assertEqual(len(self.client.rows),1);self.assertEqual(self.state.pending(),[])
        self.assertFalse(any(k.startswith('_') for k in self.client.rows[0]['fields']))
    def test_legacy_seed_and_untrusted_context_cannot_authorize_a_raw_replacement(self):
        self.params(pulse='t1');(self.path/'ser').write_bytes(b'fixture');self.completed_audit();obs=self.observation();obs.pop('_AcquisitionIdentity');self.state.observe([obs]);flush(self.state,self.client)
        self.state.observe([self.observation()]);flush(self.state,self.client) # unchanged legacy fingerprint seeds context
        self.params(end=1700000020,pulse='t1');self.completed_audit();self.sync();self.assertFalse(self.client.rows[0]['fields']['ReconciliationReason'])
        self.params(start=1700000300,pulse='t1');forged=self.observation();forged['_TransitionFingerprints']=[self.client.rows[0]['fields']['ObservedFingerprint'],forged['ObservedFingerprint']]
        self.state.observe([forged]);flush(self.state,self.client)
        self.assertIn('identity changed',self.client.rows[0]['fields']['ReconciliationReason'])

    def test_progress_anchor_requires_fractional_start_equality_not_only_same_second(self):
        self.params();(self.path/'ser').write_bytes(b'fixture')
        progress=audit_record(38,1700000000,'acquisition in progress')
        (self.path/'audita.txt').write_text(audit(progress));self.sync()
        self.params(pulse='t1')
        body=terminal(1700000000,1700000010).replace('20.000+0000,','20.500+0000,')
        (self.path/'audita.txt').write_text(audit(progress,audit_record(39,1700000010,body)))
        self.sync();self.assertIn('identity changed',self.client.rows[0]['fields']['ReconciliationReason'])

    def test_setup_birth_discards_only_obsolete_template_copy_warning(self):
        self.params();first=self.observation();other=dict(first,AcquisitionKey='src:folder/2',ObservedExpno='2')
        self.state.observe([first,other]);flush(self.state,self.client)
        self.params(start=1700000200,pulse='t1');(self.path/'ser').write_bytes(b'fixture');self.sync()
        saved=next(r['fields'] for r in self.client.rows if r['fields']['Expno']=='1')
        self.assertFalse(saved['ReconciliationReason']);self.assertEqual(saved['PulseProgram'],'t1')
        self.assertTrue(all(not r['fields']['ReconciliationReason'] for r in self.client.rows))

    def test_many_offline_end_updates_preserve_uncertain_transport_identity_with_bounded_lineage(self):
        self.params(pulse='t1');(self.path/'ser').write_bytes(b'fixture');self.completed_audit();self.sync()
        self.params(end=1700000011,pulse='t1');self.completed_audit();self.state.observe([self.observation()])
        self.client.fail_after_commit=True
        with self.assertRaises(RuntimeError):flush(self.state,self.client)
        committed=self.client.rows[0]['fields']['ObservedFingerprint']
        for end in range(1700000012,1700000040):
            self.params(end=end,pulse='t1');self.completed_audit();self.state.observe([self.observation()])
        current=self.state.cached_observations()[0][1]['ObservedFingerprint']
        lineage,_=self.state.transition_lineage('src:folder/1',current)
        self.assertIn(committed,lineage);self.assertLessEqual(len(lineage),11)
        flush(self.state,self.client);self.assertFalse(self.client.rows[0]['fields']['ReconciliationReason']);self.assertEqual(self.state.pending(),[])

    def test_raw_disappearance_cannot_reset_generation_and_disguise_replacement_as_setup(self):
        self.params(pulse='t1');(self.path/'ser').write_bytes(b'fixture');self.completed_audit();self.sync()
        (self.path/'ser').rename(self.path/'removed');self.sync()
        (self.path/'ser').write_bytes(b'new generation');self.params(start=1700001000,pulse='zg');self.sync()
        self.assertIn('raw acquisition file replaced',self.client.rows[0]['fields']['ReconciliationReason'])

    def test_local_proof_seeding_keeps_legacy_transfer_signature_stable_without_outbox(self):
        self.params(pulse='t1');(self.path/'ser').write_bytes(b'fixture')
        legacy=self.observation();legacy.pop('_AcquisitionIdentity')
        self.state.observe([legacy]);flush(self.state,self.client)
        self.state.observe([legacy]);flush(self.state,self.client)
        before=dict(self.state.db.execute('SELECT signature,payload FROM observations').fetchone())
        self.assertTrue(self.state.cached_observations()[0][1]['TransferStable'])
        patches=len(self.client.patches)
        self.assertEqual(self.state.observe([self.observation()]),0)
        after=dict(self.state.db.execute('SELECT signature,payload FROM observations').fetchone())
        self.assertEqual(before,after);self.assertEqual(self.state.pending(),[])
        flush(self.state,self.client);self.assertEqual(len(self.client.patches),patches)

    def test_raw_scientific_acqu_conflict_is_visible_without_reclassifying_calibration(self):
        self.params()
        (self.path/'acqu').write_text('##$PULPROG= <nico_nnh_tr_noSE>\n')
        self.assertNotIn('Pulse-program metadata conflict', self.observation()['ScanError'])
        (self.path/'ser').write_bytes(b'fixture')
        conflict = self.observation()
        self.assertEqual(conflict['ObservedType'], 'zg')
        self.assertIn('Pulse-program metadata conflict: acqus zg; acqu nico_nnh_tr_noSE; verify acquisition type', conflict['ScanError'])
        self.assertNotEqual(conflict['ObservedAcquisitionStatus'], 'Completed')
        (self.path/'acqu').write_text('##$PULPROG= <zg>\n')
        self.assertNotIn('Pulse-program metadata conflict', self.observation()['ScanError'])
        self.params(pulse='nico_nnh_tr_noSE')
        self.assertNotIn('Pulse-program metadata conflict', self.observation()['ScanError'])

    def test_same_protocol_fractional_restart_cannot_carry_review_or_result(self):
        self.params(pulse='t1');(self.path/'ser').write_bytes(b'fixture')
        progress=audit_record(38,1700000000,'acquisition in progress').replace('.000+', '.100+')
        (self.path/'audita.txt').write_text(audit(progress));self.sync()
        row=self.client.rows[0]['fields'];row.update(Assessment='Use', ResultCSV=['L',42])
        self.params(end=1700000020,pulse='t1')
        body=terminal(1700000000,1700000010).replace('20.000+0000,','20.700+0000,')
        (self.path/'audita.txt').write_text(audit(progress,audit_record(39,1700000010,body)))
        self.sync();self.assertIn('identity changed',row['ReconciliationReason'])
        self.assertEqual(row['Assessment'],'Use');self.assertEqual(row['ResultCSV'],['L',42])

    def test_same_protocol_exact_anchor_finalization_accepts_end_date(self):
        self.params(pulse='t1');(self.path/'ser').write_bytes(b'fixture')
        progress=audit_record(38,1700000000,'acquisition in progress')
        (self.path/'audita.txt').write_text(audit(progress));self.sync()
        self.params(end=1700000020,pulse='t1')
        (self.path/'audita.txt').write_text(audit(progress,audit_record(39,1700000020,terminal(1700000000,1700000020))))
        self.sync();self.assertFalse(self.client.rows[0]['fields']['ReconciliationReason'])

    def test_empty_or_malformed_audit_cannot_authorize_end_date_change(self):
        for content in ('', 'not a valid audit'):
            with self.subTest(content=content):
                self.state.close();self.state=State(self.root/('state'+str(len(content))))
                self.client=Fake();self.params(pulse='t1');(self.path/'ser').write_bytes(b'fixture')
                (self.path/'audita.txt').write_text(content);self.sync()
                self.params(end=1700000020,pulse='t1');self.sync()
                self.assertIn('identity changed',self.client.rows[0]['fields']['ReconciliationReason'])

    def test_retained_historical_anchor_does_not_hide_new_same_second_progress(self):
        self.params(pulse='t1');(self.path/'ser').write_bytes(b'fixture')
        first=audit_record(38,1700000000,'acquisition in progress').replace('.000+', '.100+')
        (self.path/'audita.txt').write_text(audit(first));self.sync()
        self.params(end=1700000020,pulse='t1')
        later=audit_record(39,1700000000,'acquisition in progress').replace('.000+', '.700+')
        (self.path/'audita.txt').write_text(audit(first,later));self.sync()
        self.assertIn('identity changed',self.client.rows[0]['fields']['ReconciliationReason'])

    def test_fingerprintless_clean_setup_birth_with_review_or_result_is_blocked(self):
        self.params(start=0,end=0);self.sync()
        row=self.client.rows[0]['fields'];row.update(QCFID='Clean',Assessment='Use',ResultCSV=['L',42],Notes='retain')
        self.assertFalse(row['ObservedFingerprint'])
        self.params(start=1700001000,pulse='t1');(self.path/'ser').write_bytes(b'fixture');self.sync()
        self.assertIn('setup-to-raw transition has existing review or result',row['ReconciliationReason'])
        self.assertEqual(row['Assessment'],'Use');self.assertEqual(row['ResultCSV'],['L',42]);self.assertEqual(row['Notes'],'retain')
        self.assertEqual(row['PulseProgram'],'zg')
        self.sync();self.assertIn('existing review or result',row['ReconciliationReason'])

    def test_fingerprintless_clean_setup_birth_without_review_is_allowed(self):
        self.params(start=0,end=0);self.sync();row=self.client.rows[0]['fields']
        self.params(start=1700001000,pulse='t1');(self.path/'ser').write_bytes(b'fixture');self.sync()
        self.assertFalse(row['ReconciliationReason']);self.assertEqual(row['PulseProgram'],'t1')
        self.assertTrue(row['ObservedFingerprint']);self.assertEqual(len(self.client.rows),1)

    def test_clean_placeholder_reconfiguration_updates_same_row(self):
        self.params();self.sync();row=self.client.rows[0]['fields']
        row.update(Sample=30,Notes='keep')
        self.params(start=1700000100,pulse='t1',td=4096);self.sync()
        self.assertEqual(len(self.client.rows),1)
        self.assertEqual(row['PulseProgram'],'t1');self.assertEqual(row['Notes'],'keep')
        self.assertFalse(row['ReconciliationReason']);self.assertEqual(self.state.pending(),[])

    def test_zero_byte_placeholder_recreation_is_safe(self):
        self.params();(self.path/'fid').touch();self.sync()
        replacement=self.path/'new';replacement.touch();os.replace(replacement,self.path/'fid')
        self.params(start=1700000100,pulse='t1');self.sync()
        self.assertFalse(self.client.rows[0]['fields']['ReconciliationReason'])

    def test_reviewed_placeholder_reconfiguration_is_blocked(self):
        self.params();self.sync();row=self.client.rows[0]['fields']
        row.update(Assessment='Use',ResultCSV=['L',42])
        self.params(start=1700000100,pulse='t1');self.sync()
        self.assertIn('existing review or result',row['ReconciliationReason'])
        self.assertEqual(row['PulseProgram'],'zg');self.assertEqual(row['ResultCSV'],['L',42])

    def test_clean_placeholder_peers_do_not_create_copy_warning(self):
        self.params();first=self.observation();other=dict(first,AcquisitionKey='src:folder/2',ObservedExpno='2')
        self.state.observe([first,other]);flush(self.state,self.client)
        self.assertTrue(all(not r['fields']['ReconciliationReason'] for r in self.client.rows))
        (self.path/'fid').write_bytes(b'raw')
        self.state.observe([self.observation(),other]);flush(self.state,self.client)
        self.assertTrue(all(not r['fields']['ReconciliationReason'] for r in self.client.rows))

    def test_raw_history_prevents_setup_reconfiguration_exemption(self):
        self.params();(self.path/'fid').write_bytes(b'raw');self.sync()
        (self.path/'fid').unlink();self.sync()
        self.params(start=1700000100,pulse='t1');self.sync()
        self.assertIn('identity changed',self.client.rows[0]['fields']['ReconciliationReason'])

    def test_setup_transition_requires_positive_clean_absence_on_both_sides(self):
        from nmr_register.acquisition_transition import permitted_transition
        self.params();first=self.observation()
        fields={k:v for k,v in first.items() if not k.startswith('_')}
        context={k:v for k,v in first.items() if k.startswith('_')}
        self.assertEqual(permitted_transition(fields,context,fields,context),'setup-to-setup')
        for damaged in ({'_RawFileStats':None}, {'_RawFileStats':[['fid']]},
                        {'_RawFileStats':[['fid',-1,0,1,1]]}, {'_RawGenerationFiles':[['fid',10,0,1,1]]},
                        {'_AcquisitionIdentity':dict(context['_AcquisitionIdentity'],raw_clean=False)}):
            with self.subTest(damaged=damaged):
                bad=dict(context,**damaged)
                self.assertIsNone(permitted_transition(fields,bad,fields,context))
                self.assertIsNone(permitted_transition(fields,context,fields,bad))
        (self.path/'fid').mkdir()
        malformed=self.observation()
        bad={k:v for k,v in malformed.items() if k.startswith('_')}
        self.assertIsNone(permitted_transition(fields,context,fields,bad))

    def test_remote_prior_raw_blocks_placeholder_reconfiguration(self):
        self.params();self.sync();row=self.client.rows[0]['fields']
        row.update(RawPresent=True,RawBytes=100)
        self.params(start=1700000100,pulse='t1');self.sync()
        self.assertTrue(row['ReconciliationReason'])
        self.assertEqual(row['PulseProgram'],'zg')

    def test_placeholder_lineage_preserves_review_guard_through_offline_changes(self):
        self.params();self.sync();row=self.client.rows[0]['fields'];row['Assessment']='Use'
        for start in (1700000100,1700000200):
            self.params(start=start,pulse='t1');self.state.observe([self.observation()])
        (self.path/'fid').write_bytes(b'new raw');self.state.observe([self.observation()])
        flush(self.state,self.client)
        self.assertIn('existing review or result',row['ReconciliationReason'])
        self.assertEqual(row['PulseProgram'],'zg');self.assertEqual(row['Assessment'],'Use')

    def test_same_fingerprint_zero_byte_placeholder_recreation_is_safe(self):
        self.params();(self.path/'fid').touch();self.sync()
        replacement=self.path/'new';replacement.touch();os.replace(replacement,self.path/'fid');self.sync()
        self.assertFalse(self.client.rows[0]['fields']['ReconciliationReason'])

    def test_setup_reconfiguration_retry_after_uncertain_commit(self):
        self.params();self.sync();self.params(start=1700000100,pulse='t1')
        self.state.observe([self.observation()]);self.client.fail_after_commit=True
        with self.assertRaises(RuntimeError):flush(self.state,self.client)
        flush(self.state,self.client)
        self.assertEqual(len(self.client.rows),1);self.assertEqual(self.state.pending(),[])
        self.assertFalse(self.client.rows[0]['fields']['ReconciliationReason'])

    def test_copied_terminal_audit_does_not_preclude_positive_setup_absence(self):
        self.params();self.completed_audit();self.sync()
        self.params(start=1700000100,pulse='t1');self.sync()
        self.assertFalse(self.client.rows[0]['fields']['ReconciliationReason'])
        (self.path/'fid').write_bytes(b'new raw');self.sync()
        self.assertFalse(self.client.rows[0]['fields']['ReconciliationReason'])
