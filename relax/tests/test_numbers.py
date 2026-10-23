"""Numbers-table reconstruction and comment handling, without external commands."""
import copy
import csv
from pathlib import Path
import tempfile
import unittest

from test_reporting import assessment, bundle
from test_relax_tools import tools


def read_numbers(path):
    """Parse CSV records before separating footers, preserving multiline text."""
    rows, metadata = [], {}
    with Path(path).open(newline='') as handle:
        reader = csv.reader(handle)
        fields = next(reader)
        for record in reader:
            if record[0].startswith('# '):
                key, value = ','.join(record)[2:].split(': ', 1)
                metadata[key] = value
            else:
                if len(record) != len(fields):
                    raise AssertionError('Incomplete peak record')
                rows.append(dict(zip(fields, record)))
    return fields, rows, metadata


class NumbersTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name)/'numbers.csv'

    def test_shared_values_and_measurements_reconstruct_exactly(self):
        b = bundle(count=3, overlap=(1,))
        warning = '# Warning: "quoted, value"\n# Warning: uncertainties scaled with 1.1702705.2'
        b['fit_text'] = warning+'\n'+b['fit_text']
        report = assessment(b)
        tools.write_numbers_csv(self.path, report)
        fields, rows, meta = read_numbers(self.path)
        self.assertEqual(len(rows), 3)
        self.assertEqual(fields[:6], ['assignment','residue','rate_s-1','rate_error_s-1','overlap','valid_fit'])
        self.assertEqual(len({r['peak_id'] for r in rows}), 3)
        self.assertEqual(rows[0]['assignment'], rows[1]['assignment'])
        self.assertNotEqual(rows[0]['w2'], rows[1]['w2'])
        self.assertEqual(meta['source_folder'], 'collection')
        self.assertEqual(meta['expno'], '3')
        self.assertEqual(meta['experiment_type'], 'r1rho')
        self.assertEqual(meta['peak_list'], 'peaks.list')
        self.assertEqual(meta['w1'], '120.1')
        self.assertNotIn('w1', fields)
        self.assertNotIn('w3', fields)
        self.assertNotIn('all_points_rate_s-1', fields)
        for forbidden in ['analysis_id','report_version','source_collection','source_experiment',
                          'peak_index','plane_01_noise_sd','plane_01_delay_s','plane_01_used_in_fit']:
            self.assertNotIn(forbidden, self.path.read_text())
        self.assertEqual(list(map(float, meta['plane_delays_s'].split(','))), b['times'])
        self.assertEqual(list(map(float, meta['plane_noise_sd'].split(','))), b['noises'])
        self.assertEqual(meta['plane_used_in_fit'], ', '.join(['True']*len(b['times'])))
        for row, original in zip(rows, report['rows']):
            self.assertEqual(row['warnings'], warning)
            self.assertEqual(row['overlap'], str(original['overlap']))
            self.assertEqual(row['flags'], ';'.join(original['flags']))
            self.assertEqual(float(row['rate_s-1']), original['values'][0])
            self.assertEqual(float(row['rate_error_s-1']), original['values'][1])
            self.assertEqual(int(meta['dof']), original['dof'])
            self.assertEqual([float(row[f'plane_{i:02d}_intensity'])
                              for i in range(1,len(b['times'])+1)], b['extracted'][row['peak_id']])

    def test_attempted_shortest_comparison_remains_per_peak(self):
        report = assessment(bundle(count=3))
        report['selection'] = {'audit': dict(candidate_attempted=True, fit_basis='exclude_shortest',
            reason='selected_shortest_exclusion', excluded_delay_s=0., excluded_planes_1based=[1])}
        report['used'][0] = False
        for peak in report['rows']:
            peak['all_points_fit'] = copy.deepcopy(peak)
            peak['all_points_fit']['values'][0] = 2.1
            peak['candidate_fit'] = copy.deepcopy(peak)
            peak['flags'].append('shortest_delay_excluded')
        tools.write_numbers_csv(self.path, report)
        fields, rows, meta = read_numbers(self.path)
        self.assertEqual(meta['fit_basis'], 'exclude_shortest')
        self.assertEqual(float(meta['excluded_shortest_delay_s']), 0.)
        self.assertEqual(meta['plane_used_in_fit'].split(', ')[0], 'False')
        self.assertNotIn('fit_basis', fields)
        self.assertNotIn('shortest_delay_excluded', rows[0]['flags'])
        for row in rows:
            self.assertEqual(float(row['rate_s-1']), 2.)
            self.assertEqual(float(row['all_points_rate_s-1']), 2.1)
            self.assertEqual(float(row['without_shortest_rate_s-1']), 2.)
            self.assertAlmostEqual(float(row['rate_change_percent']), 100*(2/2.1-1))

    def test_noe_has_ratios_and_codes_without_empty_decay_columns(self):
        report = assessment(bundle(count=2, mode='noe', times=[0.,1.]))
        tools.write_numbers_csv(self.path, report)
        fields, rows, meta = read_numbers(self.path)
        self.assertEqual(float(rows[0]['noe_ratio']), .8)
        self.assertEqual(float(rows[0]['noe_ratio_error']), .04)
        self.assertEqual(meta['plane_codes'], '0.0, 1.0')
        for key in ['i0','i0_error','chi2_input','p_holm','all_points_rate_s-1']:
            self.assertNotIn(key, fields)
        self.assertNotIn('plane_delays_s', meta)
        self.assertNotIn('dof', meta)

    def test_unusual_provenance_stays_on_one_metadata_line(self):
        report = assessment(bundle(count=1))
        report['context']['source_collection'] = 'a, b\n# not_a_key: text'
        tools.write_numbers_csv(self.path, report)
        _, rows, meta = read_numbers(self.path)
        self.assertEqual(len(rows), 1)
        self.assertEqual(meta['source_folder'], 'a, b\\n# not_a_key: text')
        self.assertNotIn('not_a_key', meta)


if __name__ == '__main__':
    unittest.main()
