import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock, patch

from nmr_register.config import Config, InstrumentProfile, Source, load_config, save_config
from nmr_register.discovery import match_instrument, scan_source
from nmr_register.sync import run


class InstrumentProfileTests(unittest.TestCase):
    def setUp(self):
        self.profiles = [InstrumentProfile('600_1300', 599.8, 600), InstrumentProfile('600Cryo', 600.48, 600)]

    def test_raw_frequency_distinguishes_two_nominal_600_instruments(self):
        for profile in self.profiles:
            self.assertEqual(match_instrument(profile.proton_mhz, self.profiles), profile)
        self.assertIsNone(match_instrument(600, self.profiles))
        self.assertIsNone(match_instrument(950.23, self.profiles))
        self.assertIsNone(match_instrument(599.8))

    def test_ambiguous_and_invalid_frequencies_stay_unresolved(self):
        profiles = self.profiles + [InstrumentProfile('overlap', 599.81, 600)]
        self.assertIsNone(match_instrument(599.8, profiles))
        for value in (None, '', '599.8', True, 0, -1, float('nan'), float('inf')):
            self.assertIsNone(match_instrument(value, self.profiles))

    def test_config_roundtrip_and_legacy_default(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'config.json'
            config = Config(instrument_profiles=self.profiles)
            save_config(config, path)
            self.assertEqual(load_config(path), config)
            path.write_text('{}')
            self.assertEqual(load_config(path).instrument_profiles, [])
            for value in (None, {}, ['bad'], [{'id':'bad','proton_mhz':'599.8','nominal_mhz':600}]):
                path.write_text(json.dumps({'instrument_profiles': value}))
                with self.assertRaises((ValueError, TypeError)):
                    load_config(path)

    def test_invalid_profile_values_and_duplicate_ids(self):
        for name in ('proton_mhz', 'nominal_mhz', 'tolerance_mhz'):
            for value in (0, -1, float('nan'), float('inf'), True, '600', None):
                values = dict(id='test', proton_mhz=600, nominal_mhz=600, tolerance_mhz=.05)
                values[name] = value
                with self.assertRaises(ValueError):
                    Config(instrument_profiles=[InstrumentProfile(**values)]).validate()
        for identifier in ('', ' ', ' padded ', None, 600):
            with self.assertRaises(ValueError):
                Config(instrument_profiles=[InstrumentProfile(identifier,600,600)]).validate()
        with self.assertRaises(ValueError):
            Config(instrument_profiles=[self.profiles[0], self.profiles[0]]).validate()

    def test_scan_inference_explicit_precedence_and_protected_pruning(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in ('sample', '211@950'):
                acquisition = root/name/'1'
                acquisition.mkdir(parents=True)
                (acquisition/'acqus').write_text('##$SFO1= 600.48\n')
            (root/'alias').symlink_to(root/'sample', target_is_directory=True)
            with patch('nmr_register.discovery.observe', side_effect=lambda path: {'ObservedField_MHz':600.48,'ObservedProtocolEvidence':'pulse evidence'}) as observe:
                rows, errors = scan_source(Source('local',directory), self.profiles)
                self.assertEqual(errors, [])
                self.assertEqual(len(rows),1)
                self.assertEqual(observe.call_args.args[0],root/'sample'/'1')
                self.assertEqual(rows[0]['ObservedInstrument'],'600Cryo')
                self.assertEqual(rows[0]['ObservedField_MHz'],600.48)
                self.assertEqual(rows[0]['ObservedProtocolEvidence'],'pulse evidence; instrument=600Cryo from proton frequency 600.48 MHz; configured match ±0.05 MHz')
                rows,_ = scan_source(Source('local',directory,instrument='explicit'),self.profiles)
                self.assertEqual(rows[0]['ObservedInstrument'],'explicit')
                self.assertEqual(rows[0]['ObservedProtocolEvidence'],'pulse evidence')
                rows,_ = scan_source(Source('local',directory))
                self.assertNotIn('ObservedInstrument',rows[0])

    def test_application_scan_passes_profiles(self):
        config = Config(sources=[Source('local','/synthetic')],instrument_profiles=self.profiles)
        scan = Mock(return_value=([],[]))
        with patch.dict(run.__globals__, {'scan_source': scan}):
            self.assertTrue(run(config,None,dry_run=True)['dry_run'])
            scan.assert_called_once_with(config.sources[0],self.profiles)
