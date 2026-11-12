import re
import unittest
from datetime import datetime, timedelta, timezone
from nmr_register.bruker import pulse_type
from nmr_register.schema import ATTENTION_FORMULA, PRE_DUMMY_ATTENTION_FORMULA, canonical_experiment_type
from nmr_register.sync import _canonical_payload


class CalibrationTests(unittest.TestCase):
    def test_exact_normalized_names_only(self):
        for name in ('zg','zg30',' Zg30 ','15N_calib','13C_calib',' 15n_CALIB '):
            self.assertEqual(pulse_type(name),'calibration')
        for name in ('zg300','zg30_extra','zgpr','zg_calib','15N_calib_extra','13C_calibration','x15N_calib','1H_calib'):
            self.assertNotEqual(pulse_type(name),'calibration')
        self.assertEqual(canonical_experiment_type('calibration'),'Calibration')
        machine,facts=_canonical_payload({'ObservedType':'calibration','ObservedAcquisitionStatus':'Finished'})
        self.assertEqual(facts['ExperimentType'],'Calibration')
        self.assertEqual(facts['AcquisitionStatus'],'Completed')

    def test_only_age_reminder_is_exempt(self):
        now=datetime.now(timezone.utc)
        source=re.sub(r'\$(\w+)',lambda m:'fields.get('+repr(m[1])+')',ATTENTION_FORMULA)
        namespace={'NOW':lambda:now}
        exec('def run(fields):\n'+''.join('    '+line+'\n' for line in source.splitlines()),namespace)
        fields={'Study':'Electrostatics','AcquisitionStatus':'Completed','EffectiveAcquisitionStatus':'Completed','EffectiveAcquiredAt':now-timedelta(days=2),'RegisteredAt':now,'Sample':1}
        self.assertIn('unassessed',namespace['run'](fields))
        self.assertEqual(namespace['run']({**fields,'ExperimentType':'Calibration'}),'')
        reason=namespace['run']({**fields,'ExperimentType':'Calibration','Sample':0,'ReconciliationReason':'identity conflict','ScanError':'scan failed'})
        for text in ('Sample not assigned','identity conflict','scan failed'):self.assertIn(text,reason)
        self.assertNotIn('unassessed',reason)
        self.assertEqual(namespace['run']({'AcquisitionStatus':'Dummy'}),'')
        self.assertNotIn('Dummy',PRE_DUMMY_ATTENTION_FORMULA)
        self.assertNotIn('Calibration',PRE_DUMMY_ATTENTION_FORMULA)
