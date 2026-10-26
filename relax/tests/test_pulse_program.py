"""Pulse sequence provenance without processing spectra or refitting data."""
import argparse
import csv
import io
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
import zipfile

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import hetnoe_analyze as hn
from test_reporting import assessment, bundle, prepare, write_saved
from test_relax_tools import tools as rt, working, header
from test_numbers import read_numbers


class PulseProgramTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
    def test_saved_pulprog_is_not_inferred_from_experiment_type(self):
        for text,name in [('##$PULPROG= <15N_HETNOE>\n','15N_HETNOE'),
                          ('##$PULPROG= <N15HetNoe_rk.nl> $$ comment\n','N15HetNoe_rk.nl'),
                          ('##$PULPROG= hsqctretf3gpsitc3d.3.ts3\n','hsqctretf3gpsitc3d.3.ts3'),
                          ('##$TD= 128\n','unknown')]:
            self.assertEqual(rt.pulse_program_from_text(text),name)
        for text in ['##$PULPROG= <>\n','##$PULPROG= <one\n','##$PULPROG= <one>\n##$PULPROG= <two>\n']:
            with self.assertRaises(ValueError):rt.pulse_program_from_text(text)
    def test_sequence_is_one_numbers_footer_and_in_pdf_and_support_csv(self):
        b=bundle(count=2);b['context']['pulse_program']='my_real_sequence'
        report=assessment(b);rt.write_numbers_csv(self.root/'numbers.csv',report)
        fields,rows,meta=read_numbers(self.root/'numbers.csv')
        self.assertNotIn('pulse_program',fields)
        self.assertEqual(meta['pulse_program'],'my_real_sequence')
        self.assertEqual((self.root/'numbers.csv').read_text().count('# pulse_program:'),1)
        self.assertEqual(len(rows),2)
        prepare(self.root,b)
        self.assertIn('Pulse sequence: my_real_sequence',(self.root/'plot.gnu').read_text())
        loaded=rt.load_analysis_csv(self.root/'sample_analysis.csv')
        self.assertEqual(loaded['context']['pulse_program'],'my_real_sequence')
        self.assertEqual(loaded['extracted'],b['extracted'])
    def test_old_portable_schema_without_pulse_still_loads(self):
        b=bundle(count=2);prepare(self.root,b)
        path=self.root/'sample_analysis.csv'
        with path.open(newline='') as f:rows=list(csv.DictReader(f))
        with path.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=rt.PRE_PULSE_CSV_FIELDS,extrasaction='ignore');writer.writeheader();writer.writerows(rows)
        loaded=rt.load_analysis_csv(path)
        self.assertEqual(loaded['context']['pulse_program'],'unknown')
        self.assertEqual(loaded['extracted'],b['extracted'])
    def test_inconsistent_repeated_pulse_is_rejected(self):
        b=bundle(count=2);prepare(self.root,b);path=self.root/'sample_analysis.csv'
        with path.open(newline='') as f:reader=csv.DictReader(f);fields=reader.fieldnames;rows=list(reader)
        rows[1]['pulse_program']='other'
        with path.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
        with self.assertRaisesRegex(ValueError,'provenance'):rt.load_analysis_csv(path)
    def test_csv_rebuild_retains_sequence_without_acquisition_or_fitting(self):
        b=bundle(count=2);b['context']['pulse_program']='original_sequence';prepare(self.root,b)
        (self.root/'acqus').write_text('##$PULPROG= <unrelated>\n')
        destination=self.root/'portable';destination.mkdir()
        args=argparse.Namespace(csv=str(self.root/'sample_analysis.csv'),saved_prefix=None,
            experiment_type=None,peak_list_name=None,output_dir=str(destination))
        def render(command,cwd,stdout=None):
            self.assertEqual(command[0],'gnuplot')
            self.assertIn('Pulse sequence: original_sequence',(Path(cwd)/'plot.gnu').read_text())
            (Path(cwd)/b['context']['pdf_name']).write_bytes(b'%PDF-test')
        with mock.patch.object(rt,'checked_run',side_effect=render):rt.refresh_report(args)
        self.assertEqual(read_numbers(destination/'sample_numbers.csv')[2]['pulse_program'],'original_sequence')
    def test_acquisition_enrichment_requires_matching_source_and_detects_conflict(self):
        source=self.root/'collection'/'3';source.mkdir(parents=True)
        (source/'acqus').write_text('##$PULPROG= <actual>\n')
        c=rt.report_context('r1','sample','sample.pdf','collection/3')
        rt.add_acquisition_pulse(c,source);self.assertEqual(c['pulse_program'],'actual')
        c['pulse_program']='different'
        with self.assertRaisesRegex(ValueError,'disagrees'):rt.add_acquisition_pulse(c,source)
        c['source_experiment']='4';rt.add_acquisition_pulse(c,source)
        self.assertEqual(c['pulse_program'],'different')
    def test_old_saved_diagnostics_without_pulse_are_accepted(self):
        b=bundle(count=2);b['context']['peak_list']='sample_T1rho.list';b['context']['pdf_name']='sample_T1rho.list_intensity.pdf'
        write_saved(self.root,b)
        path=self.root/'sample_fit_diagnostics.tsv'
        with path.open(newline='') as f:reader=csv.DictReader(f,delimiter='\t');fields=[k for k in reader.fieldnames if k!='pulse_program'];rows=list(reader)
        with path.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
        loaded=rt.load_saved_analysis(self.root,'sample','r1rho')
        self.assertEqual(loaded['context']['pulse_program'],'unknown')
        self.assertEqual(loaded['extracted'],b['extracted'])
    def test_normal_noe_analysis_reads_acquisition_and_keeps_three_outputs(self):
        peaks=self.root/'sample_HetNOE.list';peaks.write_text('assignment w1 w2 intensity overlap\n?-A10N-H 120 8 0 False\n')
        (self.root/'acqus').write_text('##$PULPROG= <15N_HETNOE>\n')
        (self.root/'m.ft3').write_bytes(header(planes=2))
        (self.root/'header.txt').write_text('2 0 1\n');(self.root/'noise.txt').write_text('2 1 1\n')
        args=argparse.Namespace(peaks=str(peaks),mode=None,data=str(self.root/'m.ft3'),header=str(self.root/'header.txt'),
            noise=str(self.root/'noise.txt'),format='nmrpipe',peak_format='sparky2d',prefix='sample',early_point_policy='off')
        def command(cmd,cwd,stdout=None):
            if cmd[0]=='data_xtr':stdout.write('FT000001-H 100 80\n')
            elif cmd[0]=='rlx_dfit':(Path(cwd)/'rlx_dfit.out').write_text('FT000001-H 0.8 0.04\n')
            elif cmd[0]=='gnuplot':(Path(cwd)/(peaks.name+'_intensity.pdf')).write_bytes(b'%PDF-test')
            else:raise AssertionError(cmd)
        with working(self.root),mock.patch.object(rt,'checked_run',side_effect=command):rt.analyze(args)
        self.assertEqual(read_numbers(self.root/'sample_numbers.csv')[2]['pulse_program'],'15N_HETNOE')
        with zipfile.ZipFile(self.root/'sample_support.zip') as archive:
            records=list(csv.DictReader(io.StringIO(archive.read('sample_analysis.csv').decode())))
            self.assertEqual(records[0]['pulse_program'],'15N_HETNOE')
        self.assertTrue((self.root/(peaks.name+'_intensity.pdf')).is_file())
    def test_paired_reader_uses_metadata_acquisition_or_explicit_fallback(self):
        path=self.root/'reference.list';path.write_text('assignment w1 w2 intensity overlap\n# pulse_program: seqA\n')
        self.assertEqual(hn.peak_pulse_program(path),'seqA')
        self.assertEqual(hn.peak_pulse_program(path,'seqA'),'seqA')
        with self.assertRaises(ValueError):hn.peak_pulse_program(path,'seqB')
        (self.root/'acqus').write_text('##$PULPROG= <seqA>\n')
        path.write_text('assignment w1 w2 intensity overlap\n')
        self.assertEqual(hn.peak_pulse_program(path),'seqA')
        (self.root/'acqus').unlink()
        self.assertEqual(hn.peak_pulse_program(path),'unknown')
        self.assertEqual(hn.peak_pulse_program(path,'seqB'),'seqB')
    def test_paired_reader_records_different_sequences_and_preserves_noise_metadata(self):
        for name,pulse,intensity in [('reference.list','seq_ref',100),('saturated.list','seq_noe',80)]:
            (self.root/name).write_text('assignment w1 w2 intensity overlap\n?-A10N-H 120 8 '+str(intensity)+' False\n'
                '# noise_sd: 1\n# noise_method: MAD\n# spectrum_scale: 1\n# pulse_program: '+pulse+'\n')
        with working(self.root),mock.patch.object(sys,'argv',['hetnoe_analyze.py','-ref','reference.list','-noe','saturated.list']):hn.main()
        for path in self.root.glob('*results*.out'):
            text=path.read_text()
            self.assertIn('# reference_pulse_program: seq_ref',text)
            self.assertIn('# saturated_pulse_program: seq_noe',text)
            self.assertTrue(text.startswith('10 0.80000'))
        self.assertEqual(len(list(self.root.glob('*results*.out'))),2)

    def test_paired_output_preserves_numeric_lines_and_adds_one_comment(self):
        rows=[dict(resid=10,ratio=-.25,error=.02)]
        path=self.root/'numbers.out';hn.write_out(path,rows,metadata={'pulse_program':'sequence'})
        self.assertEqual(path.read_text(),'10 -0.25000 0.02000\n# pulse_program: sequence\n')


if __name__=='__main__':unittest.main()
