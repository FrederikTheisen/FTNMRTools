"""Exact sculpting pulse classification through canonical acquisition facts."""
import tempfile
import unittest
from pathlib import Path

from nmr_register.bruker import observe, pulse_type
from nmr_register.schema import canonical_experiment_type
from nmr_register.sync import _canonical_payload


class SculptingTests(unittest.TestCase):
    def test_observation_maps_exact_normalized_sculpting_to_1d(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            for name in ('IBS_1D_SCULPTING', '  ibs_1d_ScUlPtInG  '):
                with self.subTest(name=name):
                    (path/'acqus').write_text('##$PULPROG= <'+name+'>\n')
                    observed = observe(path)
                    self.assertEqual(observed['ObservedType'], '1d')
                    self.assertEqual(_canonical_payload({'ObservedType': observed['ObservedType']})[1]['ExperimentType'], '1D')
        self.assertEqual(canonical_experiment_type('1d'), '1D')
        self.assertEqual(canonical_experiment_type('1D'), '1D')

    def test_suffixes_and_calibration_keep_existing_classification(self):
        for name in ('IBS_1D_SCULPTING_extra', 'xIBS_1D_SCULPTING', 'IBS_1D_SCULPTING.bak'):
            with self.subTest(name=name):
                self.assertEqual(pulse_type(name), 'unknown')
                self.assertNotIn('ExperimentType', _canonical_payload({'ObservedType': pulse_type(name)})[1])
        for name in ('zg', '15N_calib', '13C_calib'):
            with self.subTest(name=name):
                self.assertEqual(pulse_type(name), 'calibration')
                self.assertEqual(_canonical_payload({'ObservedType': pulse_type(name)})[1]['ExperimentType'], 'Calibration')
