"""Acquisition axis widths and carriers, independent of processing metadata."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from nmr_register.bruker import MAX_METADATA_BYTES, _spectral_metadata, observe


def status(**values):
    return ''.join('##$'+key+'= '+str(value)+'\n' for key, value in values.items())


class SpectralMetadataTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.path = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def write(self, dimension=1, **values):
        name = 'acqus' if dimension == 1 else 'acqu'+str(dimension)+'s'
        text = status(**values)
        (self.path/name).write_text(text)
        return text

    def read(self, pulse='unknown'):
        return _spectral_metadata(self.path, (self.path/'acqus').read_text(), pulse)[0]

    def test_1d_recorded_ppm_needs_no_frequency_and_keeps_zero(self):
        self.write(PARMODE=0, NUC1='<1H>', SW=12.033, O1P=0, O1=100, BF1=600)
        self.assertEqual(self.read(), {'ObservedSpectralWidths_ppm':'Direct (1H): 12.033',
                                      'ObservedSpectralOffsets_ppm':'Direct (1H): 0'})
        self.write(PARMODE=0, NUC1='<1H>', SW=12.033, O1P=-1.25, O1=100, BF1=600)
        self.assertEqual(self.read()['ObservedSpectralOffsets_ppm'], 'Direct (1H): -1.25')

    def test_2d_uses_each_status_axis_and_frequency_formula(self):
        self.write(PARMODE=1, NUC1='<1H>', NUC2='<13C>', SW=12.033, O1=2820, BF1=600)
        self.write(2, NUC1='<15N>', SW_h=1380.12, SFO1=60, O1=-120, BF1=60)
        fields = self.read()
        self.assertEqual(fields['ObservedSpectralWidths_ppm'], 'Direct (1H): 12.033\nIndirect 1 (15N): 23.002')
        self.assertEqual(fields['ObservedSpectralOffsets_ppm'], 'Direct (1H): 4.7\nIndirect 1 (15N): -2')

    def test_3d_assignment_retains_single_point_third_axis(self):
        self.write(PARMODE=2, NUC1='<1H>', SW=12)
        self.write(2, NUC1='<15N>', SW=22, FnMODE=6)
        self.write(3, NUC1='<13C>', SW=70, FnMODE=5, TD=1)
        fields = self.read('BT_HNCACB')
        self.assertEqual(fields['ObservedSpectralWidths_ppm'],
                         'Direct (1H): 12\nIndirect 1 (15N): 22\nIndirect 2 (13C): 70')

    def test_reviewed_relaxation_qf_third_axis_is_array(self):
        for pulse, nucleus in [('15N_T1_HSQC','13C'), ('hsqctretf3gpsitc3d.3.ts3','off')]:
            with self.subTest(pulse=pulse):
                text = self.write(PARMODE=2, NUC1='<1H>', SW=12)
                self.write(2, NUC1='<15N>', SW=22, FnMODE=6)
                self.write(3, NUC1='<'+nucleus+'>', SW=70, O1P=30, FnMODE=1, TD=7)
                fields, evidence = _spectral_metadata(self.path, text, pulse)
                self.assertNotIn('Indirect 2', fields['ObservedSpectralWidths_ppm'])
                self.assertNotIn('Indirect 2', fields['ObservedSpectralOffsets_ppm'])
                self.assertIn('QF relaxation array excluded', '; '.join(evidence))
        # QF is not universally a relaxation array.
        self.write(3, NUC1='<13C>', SW=70, FnMODE=1)
        self.assertIn('Indirect 2 (13C)', self.read('unknown_qf')['ObservedSpectralWidths_ppm'])

    def test_missing_frequencies_do_not_use_nominal_or_other_dimensions(self):
        self.write(PARMODE=1, NUC1='<1H>', SW_h=7200, BF1=600, SFO2=60, O1=2800, SFO1=0)
        self.write(2, NUC1='<15N>', SW_h=1200, O1=100, SFO1=60)
        fields = self.read()
        self.assertEqual(fields['ObservedSpectralWidths_ppm'], '')
        self.assertEqual(fields['ObservedSpectralOffsets_ppm'], '')

    def test_invalid_values_and_unusable_nuclei_remain_unknown(self):
        for bad in ['nan','inf','-inf','nonsense','(0..0)\n1']:
            with self.subTest(bad=bad):
                self.write(PARMODE=0, NUC1='<1H>', SW=bad, O1P=bad)
                self.assertEqual(set(self.read().values()), {''})
        for nucleus in ['off', '', '1H\nforged', '(0..0)\n1H']:
            self.write(PARMODE=0, NUC1=nucleus, SW=12, O1P=4.7)
            self.assertEqual(set(self.read().values()), {''})

    def test_status_files_only_and_parmode_bounds(self):
        self.write(PARMODE=0, NUC1='<1H>', SW=12)
        self.write(2, NUC1='<15N>', SW=20)
        (self.path/'acqu').write_text(status(NUC1='<1H>', SW=99))
        (self.path/'acqu2').write_text(status(NUC1='<13C>', SW=99))
        (self.path/'procs').write_text(status(SW_p=99, OFFSET=999))
        self.assertEqual(self.read()['ObservedSpectralWidths_ppm'], 'Direct (1H): 12')
        for mode in ['-1','8','1.5','nan','unknown']:
            self.write(PARMODE=mode, NUC1='<1H>', SW=12)
            self.assertEqual(self.read()['ObservedSpectralWidths_ppm'], '')
        self.write(NUC1='<1H>', SW=12)
        self.assertIn('Indirect 1 (15N)', self.read()['ObservedSpectralWidths_ppm'])

    def test_four_and_eight_dimensions_use_bounded_explicit_names(self):
        for count in [4,8]:
            with self.subTest(count=count):
                text = self.write(PARMODE=count-1, NUC1='<1H>', SW=12)
                for dimension in range(2,10):
                    self.write(dimension, NUC1='<13C>', SW=dimension)
                import nmr_register.bruker as bruker
                with patch.object(bruker, '_read_metadata', wraps=bruker._read_metadata) as read:
                    fields, _ = _spectral_metadata(self.path, text, 'unknown')
                self.assertEqual(len(fields['ObservedSpectralWidths_ppm'].splitlines()), count)
                self.assertEqual([call.args[0].name for call in read.call_args_list],
                                 ['acqu'+str(i)+'s' for i in range(2,count+1)])

    def test_unsafe_oversized_duplicate_and_unreadable_indirect_are_skipped(self):
        self.write(PARMODE=1, NUC1='<1H>', SW=12)
        target = self.path/'elsewhere'
        target.write_text(status(NUC1='<15N>', SW=20))
        (self.path/'acqu2s').symlink_to(target)
        self.assertEqual(self.read()['ObservedSpectralWidths_ppm'], '')
        (self.path/'acqu2s').unlink()
        (self.path/'acqu2s').write_text('x'*(MAX_METADATA_BYTES+1))
        self.assertEqual(self.read()['ObservedSpectralWidths_ppm'], '')
        (self.path/'acqu2s').write_text(status(NUC1='<15N>', SW=20)+'##$SW= 30\n')
        self.assertEqual(self.read()['ObservedSpectralWidths_ppm'], '')

    def test_declared_missing_axis_keeps_previous_complete_summary(self):
        self.write(PARMODE=1, NUC1='<1H>', SW=12, O1P=4.7)
        self.write(2, NUC1='<15N>', SW=20, O1P=118)
        previous = self.read()
        self.assertEqual(len(previous['ObservedSpectralWidths_ppm'].splitlines()), 2)
        (self.path/'acqu2s').unlink()
        fields = self.read()
        self.assertEqual(set(fields.values()), {''})
        # The collector retains old complete facts when an observation is unknown.
        from nmr_register.sync import _canonical_payload
        _, payload = _canonical_payload(fields)
        self.assertNotIn('SpectralWidths_ppm', payload)
        self.assertNotIn('SpectralOffsets_ppm', payload)

    def test_metric_completeness_is_independent(self):
        self.write(PARMODE=1, NUC1='<1H>', SW=12, O1P=4.7)
        self.write(2, NUC1='<15N>', SW='nan', O1P=118)
        fields = self.read()
        self.assertEqual(fields['ObservedSpectralWidths_ppm'], '')
        self.assertEqual(fields['ObservedSpectralOffsets_ppm'],
                         'Direct (1H): 4.7\nIndirect 1 (15N): 118')
        self.write(2, NUC1='<15N>', SW=20, O1P='nan')
        fields = self.read()
        self.assertEqual(fields['ObservedSpectralWidths_ppm'],
                         'Direct (1H): 12\nIndirect 1 (15N): 20')
        self.assertEqual(fields['ObservedSpectralOffsets_ppm'], '')

    def test_declared_unknown_nucleus_blocks_both_but_off_is_explicitly_skipped(self):
        self.write(PARMODE=1, NUC1='<1H>', SW=12, O1P=4.7)
        for nucleus in [None, '(0..0)\n15N', '<invalid>']:
            values = {'SW':20, 'O1P':118}
            if nucleus is not None:
                values['NUC1'] = nucleus
            self.write(2, **values)
            self.assertEqual(set(self.read().values()), {''})
        self.write(2, NUC1='<off>', SW=20, O1P=118)
        self.assertEqual(self.read(), {'ObservedSpectralWidths_ppm':'Direct (1H): 12',
                                      'ObservedSpectralOffsets_ppm':'Direct (1H): 4.7'})

    def test_missing_parmode_optional_files_but_unsafe_file_blocks_partial_summary(self):
        self.write(NUC1='<1H>', SW=12, O1P=4.7)
        self.assertEqual(self.read()['ObservedSpectralWidths_ppm'], 'Direct (1H): 12')
        target = self.path/'elsewhere'
        target.write_text(status(NUC1='<15N>', SW=20))
        (self.path/'acqu2s').symlink_to(target)
        self.assertEqual(set(self.read().values()), {''})

    def test_spectral_changes_do_not_change_identity_or_protocol_settings(self):
        self.write(DATE_START=1734600000, PARMODE=1, PULPROG='<hsqc>', NUC1='<1H>',
                   SW=12, O1P=4.7, BF1=600, SFO1=600.00282, TD=2048)
        self.write(2, NUC1='<15N>', SW=20, O1P=118)
        before = observe(self.path)
        self.write(2, NUC1='<15N>', SW=30, O1P=120)
        direct = self.path/'acqus'
        direct.write_text(direct.read_text().replace('##$SW= 12\n', '##$SW= 14\n')
                          .replace('##$O1P= 4.7\n', '##$O1P= -2\n'))
        after = observe(self.path)
        self.assertNotEqual(before['ObservedSpectralWidths_ppm'], after['ObservedSpectralWidths_ppm'])
        for key in ['ObservedFingerprint','_AcquisitionIdentity','ObservedD1_s',
                    'ObservedDelayList_s','ObservedSpinLock_Hz','ObservedField_MHz']:
            self.assertEqual(before[key], after[key], key)
        self.assertIn('acqu2s: SW, O1P ppm', after['ObservedProtocolEvidence'])


if __name__ == '__main__':
    unittest.main()
