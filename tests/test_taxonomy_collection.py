import unittest
from nmr_register.sync import _canonical_payload


class TaxonomyCollectionTests(unittest.TestCase):
    def test_cached_coarse_calibration_retries_keep_precise_experiment(self):
        for pulse, expected in [('zg', 'zg'), ('zg30', 'zg30'), ('15N_calib', '15N calibration'), ('13C_calib', '13C calibration')]:
            with self.subTest(pulse=pulse):
                machine, facts = _canonical_payload({'ObservedType':'calibration', 'ObservedPulseProgram':pulse})
                self.assertEqual(facts['ExperimentType'], expected)
                self.assertEqual(facts['PulseProgram'], pulse)
                self.assertNotIn('ExperimentCategory', machine)
                self.assertNotIn('ExperimentCategory', facts)

    def test_cached_unknown_assignment_is_recognized_without_new_identity(self):
        fingerprint = 'existing-acquisition-fingerprint'
        machine, facts = _canonical_payload({'ObservedType':'unknown', 'ObservedPulseProgram':'BT_HNCACB+', 'ObservedFingerprint':fingerprint})
        self.assertEqual(facts['ExperimentType'], 'HNCACB')
        self.assertEqual(facts['PulseProgram'], 'BT_HNCACB+')
        self.assertEqual(machine['ObservedFingerprint'], fingerprint)

    def test_unknown_pulse_does_not_overwrite_user_experiment(self):
        _, facts = _canonical_payload({'ObservedType':'unknown', 'ObservedPulseProgram':'my_new_pulse'})
        self.assertNotIn('ExperimentType', facts)

    def test_cached_unknown_carbon_sofast_preserves_observation_identity(self):
        fingerprint = 'existing-carbon-sofast-fingerprint'
        for pulse in ('C_SOFAST', 'c_sofast'):
            with self.subTest(pulse=pulse):
                machine, facts = _canonical_payload({
                    'ObservedType': 'unknown', 'ObservedPulseProgram': pulse,
                    'ObservedFingerprint': fingerprint})
                self.assertEqual(facts, {'ExperimentType': '13C HSQC', 'PulseProgram': pulse})
                self.assertEqual(machine, {'ObservedFingerprint': fingerprint})
                for fields in (machine, facts):
                    for field in ('ExperimentCategory', 'Sample', 'Protein', 'Notes', 'Comments'):
                        self.assertNotIn(field, fields)
