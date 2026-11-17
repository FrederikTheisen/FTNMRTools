import unittest
from nmr_register.experiment_taxonomy import canonical_experiment, category, classify_pulse, is_calibration, category_formula


class ExperimentTaxonomyTests(unittest.TestCase):
    def test_aliases_and_unknowns(self):
        for alias, expected in [('r1', 'T1'), ('15N_HSQC', 'HSQC'), ('13c_hsqc', '13C HSQC'), ('r1rho', 'T1rho'), ('HN(CA)CO', 'HNcaCO'), ('HN(CO)CA', 'HNcoCA'), ('HN(CO)CACB', 'HNcoCACB')]:
            self.assertEqual(canonical_experiment(alias), expected)
        self.assertEqual(canonical_experiment(' Novel sequence '), 'Novel sequence')
        self.assertEqual(category('Novel sequence'), 'Other')
        self.assertIsNone(canonical_experiment('unknown'))

    def test_pulse_families(self):
        cases = {'B_HNCA':'HNCA', 'BT_HNCACO':'HNcaCO', 'B_HNCOCA':'HNcoCA',
                 'BT_HNCOCACB':'HNcoCACB', 'BT_iHNCACB':'iHNCACB', 'BT_HNCACB+':'HNCACB',
                 'N_SOFAST_TROSY':'HSQC',
                 '15N_calib':'15N calibration', 'calib13C_fast':'13C calibration', 'zg':'zg'}
        for pulse, expected in cases.items():
            self.assertEqual(classify_pulse(pulse), expected)
        self.assertIsNone(classify_pulse('unrelated_B_HNCA_extra'))
        self.assertEqual(classify_pulse('hsqcetf3gp', '13C'), '13C HSQC')
        self.assertEqual(category('HNcoCA'), 'Assignment')
        self.assertEqual(category('HNN'), 'Assignment')
        self.assertTrue(is_calibration('zg'))
        self.assertEqual(category('T1rho'), 'Relaxation')

    def test_formula_matches_categories(self):
        formula = category_formula('value', 'pulse')
        namespace = {}
        exec('def formula(value, pulse):\n' + '\n'.join('    '+line for line in formula.splitlines()), namespace)
        for value in ('r1', 'HNcoCA', '15N calibration', 'Novel sequence', '13c_hsqc'):
            self.assertEqual(namespace['formula'](value, ''), category(value))

    def test_generated_formula_all_cached_local_pulse_families(self):
        # Unique pulse names from the bounded desktop observation inventory.
        pulses = """13C_calib 15N_HETNOE 15N_T1_HSQC 15N_T1_ts4 15N_T1rho
15N_T1rho_ts4 15N_T2_HSQC 15N_calib BTROSY BT_HNCACB BT_HNCACB+
BT_HNCACO BT_HNCANH BT_HNCO BT_HNCOCACB BT_iHNCACB B_HNCA B_HNCACB
B_HNCO B_HNCOCA B_HNCOCACB C_SOFAST IBS_1D_SCULPTING M.HN.hsqc.wgfb
M.calib15N.wg M.zg M.zggpw5 N15HetNoe_rk.nl NHSQC_SE N_SOFAST_TROSY
NfHSQC PRESAT T1N15_rk.nl T1rhoN15_rk.nl calib13C_fast calib15N_fast
hsqcetf3gp hsqcetfpf3gp hsqcetfpf3gpsi hsqcetfpf3gpsi2 hsqcnoef3gpsi
hsqcrexetf3gpsitc3d hsqcrexetf3gpsitc3d.2 hsqct1etf3gpsitc3d.3
hsqct1etf3gpsitc3d.3.ts3 hsqctretf3gpsitc3d.3 hsqctretf3gpsitc3d.3.ts3
n15cpmg_nh_cw_600_v5_Topspin3.5 nico_nnh_tr_noSE nico_t1sat_v2s zg zg30
zggpw5 zgpr unknown hsqc t1 unrelated_B_HNCA_extra""".split()
        namespace = {}
        body = category_formula('value', 'pulse')
        exec('def formula(value, pulse):\n' + '\n'.join('    '+line for line in body.splitlines()), namespace)
        for pulse in pulses:
            for value in ('', None, 'unknown', 'Custom experiment', '15N_HSQC', 'r1'):
                with self.subTest(pulse=pulse, value=value):
                    self.assertEqual(namespace['formula'](value, pulse), category(value, pulse))
        self.assertIsNone(classify_pulse('C_SOFAST'))
