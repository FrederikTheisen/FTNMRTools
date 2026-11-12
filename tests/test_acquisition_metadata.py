"""Acquisition metadata contracts, using small saved-program excerpts only."""
import hashlib
import json
import re
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from nmr_register.bruker import MAX_METADATA_BYTES, observe, protocol_parameters
from nmr_register.discovery import scan_source
from nmr_register.config import Source
from nmr_register.state import State
from nmr_register.sync import _conflicts, flush


T1 = '''define list<delay> vd_list=<$VDLIST>
"d31=vd_list[l11]"
"COUNTER1=( d31 / 40m )"
lo to 8 times COUNTER1
'''
RHO = '''define list<delay> vd_list=<$VDLIST>
"d31=vd_list[l11]"
"p25=d31"
"plw25=plw3*pow((p21*cnst27*4/1000000),2)"
d16 fq=cnst30(sfo hz):f3
2u pl25:f3
(p25 ph1):f3
'''
OLD_RHO = '''define list<delay> vd_list=<$VDLIST>
"p15=vd_list[t1loop]"
"plw8=plw3*(pow((p21*4*cnst8/1000000),2))"
(p13:sp11 ph0 pl8 p15 ph0 p13:sp12 ph0):f3
'''
VP = '''"p17=p18+1m-vp"
"d17=d1-p17"
"d8=p8+vp*0.25-p1*2-2u"
"d9=vp*0.5-p1*4-24u"
(p8:sp8 ph0 3u pl8 vp ph0 p8:sp9 ph0):f3
2m ivp
'''


def array(name, values):
    return '##$'+name+'= (0..'+str(len(values)-1)+')\n'+' '.join(map(str,values))+'\n'


class AcquisitionMetadataTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.path = self.root/'collection'/'1'
        self.path.mkdir(parents=True)

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, pulse='unknown_pulse', body='', extra='', reference='timings.vd'):
        cnst = [1]*64
        cnst[8],cnst[27],cnst[30] = 900,1500,0
        p = [1]*32
        p[21] = 25.81
        plw = [1]*32
        plw[3],plw[8] = 169.66,6.1115
        text = ('##$DATE_START= 1734600000\n##$PULPROG= <'+pulse+'>\n'
                '##$BF1= 600.48\n##$NUC1= <1H>\n##$NUC3= <15N>\n'
                +array('D',[.1,1.8,.2])+array('CNST',cnst)+array('P',p)+array('PLW',plw)
                +'##$VDLIST= <'+reference+'>\n'+extra)
        (self.path/'acqus').write_text(text)
        (self.path/'pulseprogram').write_text(body)
        return text

    def test_hsqc_carbon_classification_uses_indirect_nucleus(self):
        self.write('hsqcetgpsp', extra='##$NUC2= <13C>\n')
        # Having a carbon channel does not establish a carbon HSQC.
        self.assertEqual(observe(self.path)['ObservedType'],'hsqc')
        (self.path/'acqu2s').write_text('##$NUC1= <15N>\n')
        self.assertEqual(observe(self.path)['ObservedType'],'hsqc')
        before=observe(self.path)['ObservedFingerprint']
        (self.path/'acqu2s').write_text('##$NUC1= <13C>\n')
        observed=observe(self.path)
        self.assertEqual(observed['ObservedType'],'13c_hsqc')
        self.assertEqual(observed['ObservedFingerprint'],before)
        from nmr_register.sync import _canonical_payload
        self.assertEqual(_canonical_payload({'ObservedType':observed['ObservedType']})[1]['ExperimentType'],'13C HSQC')
        # An indirect carbon axis cannot turn a relaxation acquisition into HSQC.
        self.write('hsqct1etf3gpsitc3d.3.ts3',T1)
        self.assertEqual(observe(self.path)['ObservedType'],'r1')

    def test_multiline_d1_and_ordered_duplicate_programmed_times(self):
        self.write('hsqct1etf3gpsitc3d.3.ts3',T1)
        (self.path/'vdlist').write_text('0m\n1200m\n40ms\n.04s\n2e3us\n')
        fields=observe(self.path)
        self.assertEqual(fields['ObservedD1_s'],1.8)
        self.assertEqual(json.loads(fields['ObservedDelayList_s']),[0,1.2,.04,.04,.002])
        self.assertIn('quantization',fields['ObservedProtocolEvidence'])
        self.assertIn('programmed',fields['ObservedProtocolEvidence'])
        self.assertIn('vdlist#',fields['ObservedProtocolEvidence'])
        self.assertIsNone(fields['ObservedSpinLock_Hz'])

    def test_body_proven_spinlock_indices_and_zero_offset(self):
        (self.path/'vdlist').write_text('1m\n200m\n1m\n')
        self.write('hsqctretf3gpsitc3d.3.ts3',RHO)
        fields=observe(self.path)
        self.assertEqual(fields['ObservedSpinLock_Hz'],1500)
        self.assertEqual(fields['ObservedSpinLockOffset_Hz'],0)
        self.write('15N_T1rho',OLD_RHO)
        fields=observe(self.path)
        self.assertEqual(fields['ObservedSpinLock_Hz'],900)
        self.assertIsNone(fields['ObservedSpinLockOffset_Hz'])
        # A familiar name and comments cannot prove the active calculation.
        self.write('hsqctretf3gpsitc3d.3.ts3','\n'.join('; '+line for line in RHO.splitlines()))
        fields=observe(self.path)
        self.assertIsNone(fields['ObservedSpinLock_Hz'])
        self.assertEqual(fields['ObservedDelayList_s'],'')
        self.write('hsqctretf3gpsitc3d.3.ts3',RHO+'"plw25=999"\n')
        self.assertIsNone(observe(self.path)['ObservedSpinLock_Hz'])
        self.write('hsqctretf3gpsitc3d.3.ts3',RHO+'"cnst27=2000"\n')
        self.assertIsNone(observe(self.path)['ObservedSpinLock_Hz'])

    def test_vp_seconds_but_rf_power_is_not_hz(self):
        self.write('T1rhoN15_rk.nl',VP,extra='##$VPLIST= <durations.vp>\n')
        (self.path/'vplist').write_text('1m\n190m\n60m\n60m\n')
        fields=observe(self.path)
        self.assertEqual(json.loads(fields['ObservedDelayList_s']),[.001,.19,.06,.06])
        self.assertIsNone(fields['ObservedSpinLock_Hz'])
        self.assertIn('PLW[8]=6.1115 W',fields['ObservedProtocolEvidence'])

    def test_cycles_and_unreviewed_cpmg_vd_are_not_seconds(self):
        self.write('T1N15_rk.nl','lo to 70 times c',extra='##$VCLIST= <cycles.vc>\n')
        (self.path/'vclist').write_text('0\n2\n40\n15\n15\n')
        fields=observe(self.path)
        self.assertEqual(fields['ObservedDelayList_s'],'')
        self.assertIn('cycle counts=[0,2,40,15,15]',fields['ObservedProtocolEvidence'])
        self.write('hsqcrexetf3gpsitc3d.2')
        (self.path/'vdlist').write_text('0\n150\n50\n1000\n')
        fields=observe(self.path)
        self.assertEqual(fields['ObservedDelayList_s'],'')
        self.assertIn('uninterpreted list=[0,150,50,1000]',fields['ObservedProtocolEvidence'])

    def test_archive_resolution_exact_historical_basename_and_conflicting_copies(self):
        self.write('hsqct1etf3gpsitc3d.3.ts3',T1,reference='/old/instrument/lists/timings.vd')
        archive=self.path/'lists'/'vd'
        archive.mkdir(parents=True)
        (archive/'unrelated.vd').write_text('999\n')
        (archive/'timings.vd').write_text('0\n40m\n40m\n')
        fields=observe(self.path)
        self.assertEqual(json.loads(fields['ObservedDelayList_s']),[0,.04,.04])
        self.assertIn('lists/vd/timings.vd',fields['ObservedProtocolEvidence'])
        (self.path/'vdlist').write_text('0\n80m\n40m\n')
        fields=observe(self.path)
        self.assertEqual(fields['ObservedDelayList_s'],'')
        self.assertIn('conflicting vd',fields['ObservedProtocolEvidence'])
        (self.path/'vdlist').unlink();(archive/'timings.vd').unlink()
        self.assertEqual(observe(self.path)['ObservedDelayList_s'],'')

    def test_unsafe_symlink_and_oversized_lists_not_read(self):
        outside=self.root/'elsewhere.vd';outside.write_text('1m\n')
        self.write('hsqct1etf3gpsitc3d.3.ts3',T1,reference='../elsewhere.vd')
        self.assertIn('unsafe',observe(self.path)['ObservedProtocolEvidence'])
        self.write('hsqct1etf3gpsitc3d.3.ts3',T1)
        (self.path/'vdlist').symlink_to(outside)
        with patch.object(Path,'open',side_effect=AssertionError('symlink must not be opened')):
            from nmr_register.bruker import _selected_list
            with self.assertRaisesRegex(ValueError,'symlink'):
                _selected_list(self.path,'vd','timings.vd')
        (self.path/'vdlist').unlink()
        (self.path/'vdlist').write_bytes(b'0\n'*(MAX_METADATA_BYTES//2+1))
        self.assertIn('size limit',observe(self.path)['ObservedProtocolEvidence'])
        alias=self.root/'alias';alias.symlink_to(self.root/'collection',target_is_directory=True)
        with self.assertRaisesRegex(ValueError,'symlink'):
            observe(alias/'1')
        with self.assertRaisesRegex(ValueError,'protected'):
            observe(self.root/'211@950'/'1')

    def test_invalid_arrays_and_values_clear_unknowns(self):
        for text in ['##$D= (0..2)\n0 1\n','##$D= (0..1)\n0 1\n##$D= (0..1)\n0 2\n',
                     '##$D= (0..1)\n0 NaN\n','##$D= (0..1)\n0 -1\n']:
            with self.subTest(text=text):
                (self.path/'acqus').write_text(text)
                self.assertIsNone(observe(self.path)['ObservedD1_s'])
        self.assertEqual(protocol_parameters('##$D= (0..1)\n0 $$ inline\n0\n')['D'],['0','0'])
        self.assertEqual(protocol_parameters('##$NAME= (0..0)\n<value with spaces>\n')['NAME'],['value with spaces'])
        self.write('hsqct1etf3gpsitc3d.3.ts3',T1)
        for value in ['NaN','Infinity','-1m','1e999s']:
            (self.path/'vdlist').write_text(value+'\n')
            self.assertEqual(observe(self.path)['ObservedDelayList_s'],'')

    def test_proton_channel_identification_never_infers_instrument(self):
        (self.path/'acqus').write_text('##$BF1= 150.9\n##$NUC1= <13C>\n##$BF2= 600.48\n##$NUC2= <1H>\n')
        fields=observe(self.path)
        self.assertEqual(fields['ObservedField_MHz'],600.48)
        self.assertNotIn('ObservedInstrument',fields)
        (self.path/'acqus').write_text('##$BF1= 600.48\n')
        self.assertIsNone(observe(self.path)['ObservedField_MHz'])
        (self.path/'acqus').write_text('##$BF1= 600\n##$NUC1= <1H>\n##$BF2= 950\n##$NUC2= <1H>\n')
        self.assertIsNone(observe(self.path)['ObservedField_MHz'])

    def test_original_fingerprint_unchanged_and_protocol_changes_queue_without_replacement(self):
        original=self.write('hsqctretf3gpsitc3d.3.ts3',RHO)
        legacy={m[1]:m[2].strip() for m in re.finditer(r'^##\$([A-Za-z0-9_]+)=\s*([^\n]*)',original,re.M)}
        identity={k:legacy.get(k,'') for k in ('DATE_START','DATE','PULPROG','BF1','SFO1','NUC1','TD')}
        expected=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()
        self.assertEqual(observe(self.path)['ObservedFingerprint'],expected)
        (self.path/'vdlist').write_text('1m\n100m\n')
        state=State(self.root/'state')
        try:
            def records():return scan_source(Source('local',str(self.root),instrument='explicitly configured'))[0]
            state.observe(records());state.observe(records());state.acknowledge([key for key,_ in state.pending()])
            (self.path/'vdlist').write_text('1m\n200m\n1m\n')
            state.observe(records())
            payload=state.pending()[0][1]
            self.assertEqual(payload['ObservedFingerprint'],expected)
            self.assertFalse(payload['ReconciliationReason'])
            self.assertEqual(payload['ObservedInstrument'],'explicitly configured')
            (self.path/'acqus').write_text('##$DATE_START= 1734600000\n##$PULPROG= <hsqctretf3gpsitc3d.3.ts3>\n##$BF1= 600.48\n##$NUC1= <1H>\n')
            (self.path/'pulseprogram').unlink();(self.path/'vdlist').unlink()
            state.observe(records())
            payload=state.pending()[0][1]
            self.assertIsNone(payload['ObservedD1_s'])
            self.assertIsNone(payload['ObservedSpinLock_Hz'])
            self.assertEqual(payload['ObservedDelayList_s'],'')
            self.assertFalse(payload['ReconciliationReason'])
        finally:state.close()

    def test_protocol_differences_precision_order_and_zero(self):
        machine={'ObservedD1_s':1.8000000001,'ObservedSpinLock_Hz':1500.01,
                 'ObservedSpinLockOffset_Hz':0,'ObservedDelayList_s':'[0,0.04,0.04]'}
        human={'D1_s':1.8,'SpinLock_Hz':1500,'SpinLockOffset_Hz':0,'DelayList_s':'[0.0, 4e-2, 0.04]'}
        self.assertEqual(_conflicts(machine,human),[])
        human['D1_s']=0;human['SpinLockOffset_Hz']=10;human['DelayList_s']='[0.04,0,0.04]'
        differences=_conflicts(machine,human)
        self.assertIn('observed D1_s differs from confirmed value',differences)
        self.assertIn('observed SpinLockOffset_Hz differs from confirmed value',differences)
        self.assertIn('observed DelayList_s differs from confirmed value',differences)

    def test_actual_protocol_facts_win_and_unknowns_retain_recorded_values(self):
        from test_sync import Fake
        self.write('hsqctretf3gpsitc3d.3.ts3',RHO)
        (self.path/'vdlist').write_text('0\n1m\n')
        state=State(self.root/'state')
        try:
            records=scan_source(Source('local',str(self.root)))[0]
            state.observe(records);client=Fake();flush(state,client)
            fields=client.rows[0]['fields']
            self.assertEqual(fields['D1_s'],1.8)
            self.assertEqual(fields['SpinLock_Hz'],1500)
            self.assertEqual(fields['SpinLockOffset_Hz'],0)
            self.assertFalse(any(k.startswith('ObservedSpinLock') for k in fields))
            fields.update(D1_s=0,SpinLock_Hz=0,SpinLockOffset_Hz=0,DelayList_s='[0]',Assessment='Use',Notes='human')
            (self.path/'pulseprogram').unlink()
            state.observe(scan_source(Source('local',str(self.root)))[0]);flush(state,client)
            self.assertEqual(fields['D1_s'],1.8);self.assertEqual(fields['SpinLock_Hz'],0)
            self.assertEqual(fields['SpinLockOffset_Hz'],0)
            self.assertEqual(fields['DelayList_s'],'[0]');self.assertEqual(fields['Assessment'],'Use')
            self.assertTrue(all(not {'SpinLock_Hz','SpinLockOffset_Hz','DelayList_s','Assessment','Notes'} & record['fields'].keys()
                                for record in client.patches))
        finally:state.close()
