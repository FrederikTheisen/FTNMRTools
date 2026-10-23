"""Saved-results reporting tests: no real processing, extraction or fitting."""
import argparse
import copy
import csv
import math
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
import uuid
from unittest import mock

from test_relax_tools import tools, RELAX, working


def bundle(count=12, times=None, mode='r1rho', perturb=None, overlap=()):
    times = [0., .2, .025, .1, .025] if times is None else times
    peaks=[];extracted={};lines=[]
    for i in range(count):
        tag=f'FT{i+1:06d}-N-H'
        peaks.append(dict(peak_id=tag,assignment=f'?-E{174 if i<2 else 175+i}N-H',
                          residue=174 if i<2 else 175+i,coordinates=['120.1',str(8+i*.01)],
                          overlap=i in overlap,duplicate_assignment=i<2 and count>1))
        lines.append(tag+(' 0.8 0.04' if mode=='noe' else ' 2 0.1 100 3 5 5'))
        extracted[tag]=[(100*math.exp(-2*t) if mode!='noe' else 80.)+
                        (perturb(i,j) if perturb else 0.) for j,t in enumerate(times)]
    return dict(peaks=peaks,extracted=extracted,fit_text='\n'.join(lines)+'\n',times=times,
                noises=[1.+.1*j for j in range(len(times))],
                context=tools.report_context(mode,'sample','peaks.list_intensity.pdf','collection/3',peak_list='peaks.list'))


def assessment(b):
    return tools.build_assessment(tools.fitted_rows(b['peaks'],b['fit_text'],b['context']['experiment_type']),
                                  b['extracted'],b['times'],b['noises'],b['context'])


def prepare(root,b):
    c=b['context']
    return tools.prepare_reports(root,b['peaks'],b['extracted'],b['fit_text'],b['times'],b['noises'],
                               c['experiment_type'],c['output_prefix'],c['pdf_name'],'collection/3',
                               context=c,warning_map=b.get('warning_map'))


def write_saved(root,b):
    prepare(root,b)
    (root/'rlx_dfit.inp').write_text(str(len(b['times']))+' '+' '.join(map(str,b['times']))+'\n'+
        ''.join(tag+' '+' '.join(map(str,y))+'\n' for tag,y in b['extracted'].items()))
    (root/'rlx_dfit.out').write_text(b['fit_text'])
    for name,values in [('relax_header.txt',b['times']),('vEstNoise.txt',b['noises'])]:
        (root/name).write_text(str(len(values))+' '+' '.join(map(str,values))+'\n')


class ProbabilityTests(unittest.TestCase):
    def test_nist_published_chi_square_critical_values(self):
        # https://www.itl.nist.gov/div898/handbook/eda/section3/eda3674.htm
        # NIST table rounds the quantiles to three decimal places.
        for df,x,p in [(1,3.841,.05),(2,5.991,.05),(5,11.070,.05),
                       (10,18.307,.05),(30,59.703,.001),(100,135.807,.01)]:
            with self.subTest(df=df):
                self.assertAlmostEqual(tools.chi_square_sf(x,df),p,delta=3e-5)
    def test_closed_form_probabilities_and_tails(self):
        for x in [0.,.01,1.,4.,30.,500.]:
            self.assertAlmostEqual(tools.chi_square_sf(x,2),math.exp(-x/2),delta=1e-14)
            self.assertAlmostEqual(tools.chi_square_sf(x,1),math.erfc(math.sqrt(x/2)),delta=1e-14)
        for x,df in [(-1,2),(math.nan,2),(1,0)]:
            with self.assertRaises(ValueError):tools.chi_square_sf(x,df)
    def test_holm_stepdown_reference_unsorted_and_ties(self):
        # R stats documented Holm step-down: cummax((m-i+1)*p_ordered), capped at one.
        given={'c':.04,'a':.001,'e':.8,'b':.01,'d':.04}
        self.assertEqual(tools.holm_adjust(given),{'a':.005,'b':.04,'c':.12,'d':.12,'e':.8})
        self.assertEqual(tools.holm_adjust({}),{})
        self.assertEqual(tools.holm_adjust({'a':.5,'b':.7}),{'a':1.,'b':1.})


class AssessmentTests(unittest.TestCase):
    def test_clean_repeated_unsorted_decays(self):
        b=bundle();r=assessment(b)
        self.assertEqual(r['conclusion'],'no_evidence_against');self.assertEqual(r['primary_count'],12)
        self.assertEqual(r['times'],b['times'])
        for row in r['rows']:
            self.assertEqual(row['dof'],3);self.assertEqual(row['p_holm'],1.)
            self.assertEqual(row['repeat_dof'],1);self.assertEqual(row['repeat_chi2'],0.)
    def test_shared_plane_bias_and_individual_outlier(self):
        r=assessment(bundle(perturb=lambda i,j:4 if j==1 else 0))
        self.assertEqual(r['conclusion'],'departures_detected');self.assertEqual(r['coherent'],[2])
        r=assessment(bundle(perturb=lambda i,j:20 if i==0 and j==2 else 0))
        self.assertEqual(r['rejected'],1);self.assertEqual(r['coherent'],[])
        self.assertEqual(r['rows'][0]['values'][4],5.) # fitter statistic remains unrelated
    def test_quality_flags_do_not_silently_exclude_otherwise_valid_fits(self):
        b=bundle(count=1);b['fit_text']=b['fit_text'].replace(' 2 0.1',' 0 0.1')
        r=assessment(b)
        self.assertEqual(r['primary_count'],1)
        self.assertIn('nonpositive_rate',r['rows'][0]['flags'])
        self.assertTrue(r['rows'][0]['valid'])
    def test_original_noise_sum_uses_accurate_reproducible_accumulation(self):
        b=bundle(count=1,times=[0.,.1,.2,.3]);b['noises']=[1.]*4
        tag=b['peaks'][0]['peak_id']
        b['extracted'][tag]=[100*math.exp(-2*t)+z for t,z in zip(b['times'],[1e8,1.,1.,1.])]
        r=assessment(b);row=r['rows'][0]
        self.assertEqual(row['chi2_input'],math.fsum(p['z']**2 for p in row['points']))
        self.assertEqual(row['chi2_input'],1e16+4.)
    def test_bias_rule_boundaries(self):
        self.assertFalse(tools.residual_summary([3]*9)['bias'])
        self.assertFalse(tools.residual_summary([2]*10)['bias'])
        self.assertTrue(tools.residual_summary([3]*8+[-1]*2)['bias'])
        self.assertFalse(tools.residual_summary([3]*7+[-1]*3)['bias'])
    def test_overlap_and_rescaling_are_separate(self):
        b=bundle(perturb=lambda i,j:20 if i==0 and j==2 else 0,overlap=(0,))
        tag=b['peaks'][1]['peak_id']
        warning='# Warning: Data uncertainties scaled with 1.1702705.2 for site '+tag
        b['fit_text']=warning+'\n'+b['fit_text'];r=assessment(b)
        self.assertEqual(r['primary_count'],11);self.assertEqual(r['rejected'],0)
        self.assertEqual(r['rows'][0]['exclusion_reason'],'overlap')
        self.assertIn('uncertainty_rescaled',r['rows'][1]['flags'])
        self.assertEqual(r['rows'][1]['warnings'],warning)
        self.assertEqual(r['rows'][1]['values'][1],.1)
    def test_inadequate_degrees_of_freedom_and_noe(self):
        r=assessment(bundle(times=[0.,.1]))
        self.assertEqual(r['conclusion'],'inconclusive');self.assertEqual(r['primary_count'],0)
        self.assertIsNone(r['rows'][0]['p_input'])
        r=assessment(bundle(mode='noe'))
        self.assertEqual(r['conclusion'],'not_applicable');self.assertIsNone(r['rows'][0]['chi2_input'])
    def test_weighted_repeat_scatter_and_all_pair_groups(self):
        b=bundle(count=1,times=[.1,.2,.1,.2],perturb=lambda i,j:j)
        r=assessment(b);row=r['rows'][0]
        self.assertEqual(row['repeat_dof'],2)
        for t,indices in [(.1,[0,2]),(.2,[1,3])]:
            i,j=indices;expected=(b['extracted'][row['peak_id']][i]-b['extracted'][row['peak_id']][j])/math.hypot(b['noises'][i],b['noises'][j])
            self.assertAlmostEqual(row['repeat_groups'][t]['pair_z'],expected)
            self.assertAlmostEqual(row['repeat_groups'][t]['chi2'],expected**2)
        self.assertIsNone(row['repeat_pair_z'])
    def test_empty_primary_coverage_and_invalid_i0(self):
        b=bundle(count=1);b['fit_text']=b['fit_text'].replace('100 3','0 3')
        r=assessment(b);self.assertEqual(r['conclusion'],'inconclusive')
        self.assertIsNone(r['rows'][0]['points'][0]['prediction'])


class PortableTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
        self.b=bundle(count=3);prepare(self.root,self.b);self.csv=self.root/'sample_analysis.csv'
    def rewrite(self,mutate):
        with self.csv.open(newline='') as f:rows=list(csv.DictReader(f))
        mutate(rows)
        with self.csv.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=tools.CSV_FIELDS);w.writeheader();w.writerows(rows)
    def test_roundtrip_reordered_rows_and_duplicate_assignments(self):
        before=self.csv.read_bytes();r=tools.load_analysis_csv(self.csv)
        prepare(self.root,r);self.assertEqual(self.csv.read_bytes(),before)
        self.rewrite(lambda rows:random.Random(42).shuffle(rows))
        r=tools.load_analysis_csv(self.csv)
        self.assertEqual(r['extracted'],self.b['extracted']);self.assertEqual(r['times'],self.b['times'])
        self.assertEqual(r['peaks'][0]['assignment'],r['peaks'][1]['assignment'])
        self.assertNotEqual(r['peaks'][0]['coordinates'],r['peaks'][1]['coordinates'])
        prepare(self.root,r);self.assertEqual(self.csv.read_bytes(),before)
    def test_compact_peak_summary_has_one_row_per_peak_and_provenance(self):
        path=self.root/'sample_peak_summary.csv'
        with path.open(newline='') as handle:rows=list(csv.DictReader(handle))
        self.assertEqual(len(rows),len(self.b['peaks']))
        self.assertEqual(len({row['peak_id'] for row in rows}),len(rows))
        row=rows[0]
        self.assertEqual(row['source_collection'],'collection')
        self.assertEqual(row['source_experiment'],'3')
        self.assertEqual(row['peak_list'],'peaks.list')
        self.assertEqual(row['experiment_type'],'r1rho')
        self.assertEqual(row['selected_fit_basis'],'all_points')
        self.assertTrue(row['rate_s-1'])
        self.assertTrue(row['rate_error_s-1'])
        self.assertIn('overlap',rows[0])
        self.assertEqual(row['all_points_rate_s-1'],'')
        self.assertEqual(row['without_shortest_rate_s-1'],'')

    def test_warning_text_csv_quoting_and_global_warning(self):
        warning='# Warning: "quoted, value"\n# Warning: global uncertainties scaled with 1.1702705.2'
        self.b['fit_text']=warning+'\n'+self.b['fit_text'];prepare(self.root,self.b)
        loaded=tools.load_analysis_csv(self.csv);prepare(self.root,loaded)
        with self.csv.open(newline='') as f:rows=list(csv.DictReader(f))
        self.assertTrue(all(r['warnings']==warning for r in rows))
    def test_reject_missing_rows_inconsistent_metadata_noise_and_tampering(self):
        before=self.csv.read_bytes()
        for mutation in [lambda rows:rows.pop(2),lambda rows:rows[0].update(value='9'),
                         lambda rows:rows[0].update(noise_sd='4'),lambda rows:rows[0].update(report_version='different'),
                         lambda rows:rows[0].update(intensity='99'),lambda rows:rows[0].update(plane_index='20'),
                         lambda rows:rows[0].update(plane_bias='True')]:
            with self.subTest(mutation=mutation):
                self.csv.write_bytes(before);self.rewrite(mutation)
                with self.assertRaises(ValueError):tools.load_analysis_csv(self.csv)
    def test_stable_uuid_separate_report_version(self):
        with self.csv.open() as f:original=next(csv.DictReader(f))
        with mock.patch.object(tools,'REPORT_VERSION','next'):
            b=tools.load_analysis_csv(self.csv);prepare(self.root,b)
        with self.csv.open() as f:changed=next(csv.DictReader(f))
        self.assertEqual(original['analysis_id'],changed['analysis_id'])
        self.assertEqual(changed['report_version'],'next')
        self.assertEqual(uuid.UUID(changed['analysis_id']).version,4)
    def test_saved_legacy_consistency_and_native_bytes(self):
        write_saved(self.root,self.b)
        names=['rlx_dfit.inp','rlx_dfit.out','vEstNoise.txt','relax_header.txt']
        before={n:(self.root/n).read_bytes() for n in names}
        loaded=tools.load_saved_analysis(self.root,'sample','r1rho');prepare(self.root,loaded)
        self.assertEqual(before,{n:(self.root/n).read_bytes() for n in names})
        self.assertEqual(loaded['extracted'],self.b['extracted'])
        for name in ['rlx_dfit.out','vEstNoise.txt','intensity.out','sample_fit_diagnostics.tsv']:
            write_saved(self.root,self.b);path=self.root/name;old=path.read_text()
            if name=='rlx_dfit.out':path.write_text(old.replace(' 2 0.1',' 3 0.1',1))
            elif name=='vEstNoise.txt':path.write_text(old.replace('1.0','2.0',1))
            elif name=='intensity.out':path.write_text(old.replace('0 1 0.01','0 0.7 0.01',1))
            else:path.write_text(old.replace('120.1','nan',1))
            with self.subTest(name=name),self.assertRaises(ValueError):tools.load_saved_analysis(self.root,'sample','r1rho')
    def test_saved_refresh_retains_uuid_and_fit_origin(self):
        self.b['context'].update(fit_origin='analysis')
        write_saved(self.root,self.b)
        loaded=tools.load_saved_analysis(self.root,'sample','r1rho')
        self.assertEqual(loaded['context']['analysis_id'],self.b['context']['analysis_id'])
        self.assertEqual(loaded['context']['fit_origin'],'analysis')
    def test_mapping_and_identifier_mismatches_fail_before_render(self):
        write_saved(self.root,self.b)
        mapping=self.root/'sample_fit_diagnostics.tsv';before=mapping.read_text()
        for altered in [before.replace('FT000002-N-H','FT000001-N-H'),
                        before.replace('False\tTrue','False\tFalse',1)]:
            mapping.write_text(altered)
            with self.assertRaises(ValueError):tools.load_saved_analysis(self.root,'sample','r1rho')
            mapping.write_text(before)
        fit=self.root/'rlx_dfit.out';fit.write_text(fit.read_text().replace('FT000002-N-H','UNEXPECTED-N-H'))
        with self.assertRaises(ValueError):tools.load_saved_analysis(self.root,'sample','r1rho')
    def test_standalone_csv_render_without_processing_or_fitting(self):
        dest=self.root/'portable';dest.mkdir();local=dest/'only.csv';shutil.copyfile(self.csv,local)
        binary=self.root/'bin';binary.mkdir();calls=self.root/'calls'
        # gnuplot is the sole allowed external executable, and has no fitting capability here.
        plot=binary/'gnuplot';plot.write_text('#!/usr/bin/python3\nfrom pathlib import Path\nimport shlex\n'
            f'Path({str(calls)!r}).write_text("gnuplot only")\n'
            'line=next(l for l in Path("plot.gnu").read_text().splitlines() if l.startswith("set output "))\n'
            'Path(shlex.split(line)[2]).write_bytes(b"%PDF-STUB")\n');plot.chmod(0o755)
        for name in ['data_xtr','rlx_dfit','nmrPipe','bruk2pipe']:
            script=binary/name;script.write_text('#!/bin/sh\nexit 99\n');script.chmod(0o755)
        command=[sys.executable,str(RELAX/'relax_tools.py'),'report','--csv',str(local),'--output-dir',str(dest)]
        result=subprocess.run(command,env=dict(os.environ,PATH=str(binary)),text=True,capture_output=True)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        self.assertEqual(calls.read_text(),'gnuplot only')
        self.assertTrue((dest/'peaks.list_intensity.pdf').exists());self.assertFalse((dest/'rlx_dfit.out').exists())
        self.assertEqual(tools.load_analysis_csv(dest/'sample_analysis.csv')['extracted'],self.b['extracted'])
    def test_render_and_publication_failures_preserve_reports(self):
        pdf=self.root/'peaks.list_intensity.pdf';pdf.write_bytes(b'old-pdf')
        before={p.name:p.read_bytes() for p in self.root.iterdir() if p.is_file()}
        args=argparse.Namespace(csv=str(self.csv),saved_prefix=None,experiment_type=None,output_dir=None,peak_list_name=None)
        with mock.patch.object(tools,'checked_run',side_effect=OSError('render failed')):
            with self.assertRaises(OSError):tools.refresh_report(args)
        self.assertEqual(before,{n:(self.root/n).read_bytes() for n in before})
        def render(command,stage): (stage/'peaks.list_intensity.pdf').write_bytes(b'new-pdf')
        original=Path.rename
        def broken(path,target):
            if path.name=='sample_results_filtered.out' and path.parent.name.startswith('.relax-report-'):
                raise OSError('publication failed')
            return original(path,target)
        with mock.patch.object(tools,'checked_run',side_effect=render),mock.patch.object(Path,'rename',broken):
            with self.assertRaises(OSError):tools.refresh_report(args)
        self.assertEqual(before,{n:(self.root/n).read_bytes() for n in before})
    def test_successful_refresh_backups_and_native_bytes(self):
        write_saved(self.root,self.b)
        original=(self.root/'sample_results.out').read_bytes()
        native={n:(self.root/n).read_bytes() for n in ['rlx_dfit.inp','rlx_dfit.out','vEstNoise.txt','relax_header.txt']}
        def render(command,stage): (stage/'peaks.list_intensity.pdf').write_bytes(b'new-pdf')
        args=argparse.Namespace(csv=None,saved_prefix='sample',experiment_type='r1rho',output_dir=None,peak_list_name=None)
        with working(self.root),mock.patch.object(tools,'checked_run',side_effect=render):tools.refresh_report(args)
        backups=list((self.root/'.report-backups').glob('*/*/sample_results.out'))
        self.assertEqual(len(backups),1);self.assertEqual(backups[0].read_bytes(),original)
        with self.csv.open() as f:row=next(csv.DictReader(f))
        self.assertEqual(native,{n:(self.root/n).read_bytes() for n in native})
        self.assertEqual(row['analysis_id'],self.b['context']['analysis_id'])
        self.assertFalse(any('sha256' in key for key in row))
    def legacy_rows(self):
        with self.csv.open(newline='') as handle: rows=list(csv.DictReader(handle))
        for row in rows:
            row.update(schema_version='1',analysis_id='a1-obsolete-content-id',
                       fit_revision='sha256:old-fit-code',report_revision='sha256:old-report-code',
                       fit_input_sha256='unused input digest',fit_output_sha256='unused output digest',
                       noise_sha256='unused noise digest',assessment_method='old detailed methods paragraph')
        return rows
    def write_legacy_csv(self, rows):
        with self.csv.open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=tools.LEGACY_CSV_FIELDS,extrasaction='ignore')
            writer.writeheader();writer.writerows(rows)
    def test_schema1_csv_migrates_once_with_unchanged_numbers(self):
        legacy=self.legacy_rows();self.write_legacy_csv(legacy)
        loaded=tools.load_analysis_csv(self.csv)
        self.assertEqual(loaded['context']['analysis_id'],'')
        prepare(self.root,loaded)
        first=self.csv.read_bytes()
        with self.csv.open(newline='') as handle: migrated=list(csv.DictReader(handle))
        self.assertEqual(uuid.UUID(migrated[0]['analysis_id']).version,4)
        self.assertEqual(migrated[0]['schema_version'],'2')
        for old,new in zip(legacy,migrated):
            for key in tools.CSV_FIELDS:
                if key not in ('schema_version','analysis_id','report_version'):
                    self.assertEqual(old[key],new[key],key)
        prepare(self.root,tools.load_analysis_csv(self.csv))
        self.assertEqual(self.csv.read_bytes(),first)
        for removed in ('sha256','assessment_method','fit_revision','report_revision','a1-obsolete'):
            self.assertNotIn(removed,self.csv.read_text())
    def test_schema1_csv_rejects_inconsistent_measurements_and_metadata(self):
        rows=self.legacy_rows()
        for field,value in [('intensity','99'),('value','99'),('noise_sd','99'),
                            ('prediction','99'),('source_collection','other')]:
            changed=copy.deepcopy(rows);changed[0][field]=value;self.write_legacy_csv(changed)
            with self.subTest(field=field),self.assertRaises(ValueError):tools.load_analysis_csv(self.csv)
    def test_schema1_saved_diagnostics_migrate_and_preserve_uuid_after_export(self):
        write_saved(self.root,self.b)
        path=self.root/'sample_fit_diagnostics.tsv'
        with path.open(newline='') as handle:
            reader=csv.DictReader(handle,delimiter='\t');fields=reader.fieldnames;rows=list(reader)
        old_fields=fields[:-len(tools.TRACE_FIELDS)]+tools.LEGACY_TRACE_FIELDS
        for row in rows:
            row.update(schema_version='1',analysis_id='a1-obsolete',fit_revision='sha256:old',
                       report_revision='sha256:old',fit_input_sha256='unused',fit_output_sha256='unused',noise_sha256='unused')
        with path.open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=old_fields,delimiter='\t',extrasaction='ignore')
            writer.writeheader();writer.writerows(rows)
        loaded=tools.load_saved_analysis(self.root,'sample','r1rho');prepare(self.root,loaded)
        again=tools.load_saved_analysis(self.root,'sample','r1rho')
        self.assertEqual(loaded['context']['analysis_id'],again['context']['analysis_id'])
        self.assertEqual(again['extracted'],self.b['extracted'])
    def test_removed_pages_and_hash_fields_absent_from_all_exports(self):
        plot=(self.root/'plot.gnu').read_text()
        for removed in ['Peak results','Methods and assessment coverage','Repeated-delay scatter',
                        'NOE ratio methods','Report code:', 'sha256']:
            self.assertNotIn(removed,plot)
        self.assertIn('Shared residual summaries by acquisition plane',plot)
        self.assertIn('Shared residual summaries by delay',plot)
        self.assertIn('report v2',plot);self.assertIn('I0 =',plot)
        self.assertIn('>=10 eligible peaks, >=80% same sign',plot)
        self.assertIn('Gaussian measurement noise',plot)
        for name in ['sample_analysis.csv','sample_fit_diagnostics.tsv','sample_results.out','intensity.out']:
            text=(self.root/name).read_text()
            self.assertNotIn('sha256',text);self.assertNotIn('assessment_method',text)
            self.assertNotIn('report_revision',text);self.assertNotIn('fit_revision',text)
        source=(RELAX/'relax_tools.py').read_text()
        self.assertNotIn('import hashlib',source);self.assertNotIn('def file_sha256',source)
    def test_uuid_validation_and_input_changes_during_render(self):
        before=self.csv.read_bytes()
        self.rewrite(lambda rows:[r.update(analysis_id='a1-obsolete') for r in rows])
        with self.assertRaises(ValueError):tools.load_analysis_csv(self.csv)
        self.csv.write_bytes(before)
        pdf=self.root/'peaks.list_intensity.pdf';pdf.write_bytes(b'old-pdf')
        def rendering_changes_input(command,stage):
            (stage/pdf.name).write_bytes(b'new-pdf')
            self.csv.write_bytes(before+b'\n')
        args=argparse.Namespace(csv=str(self.csv),saved_prefix=None,experiment_type=None,output_dir=None,peak_list_name=None)
        with mock.patch.object(tools,'checked_run',side_effect=rendering_changes_input):
            with self.assertRaisesRegex(ValueError,'changed during rendering'):tools.refresh_report(args)
        self.assertEqual(pdf.read_bytes(),b'old-pdf')
    def test_noe_and_invalid_amplitude_roundtrips(self):
        for mode in ['noe','r1rho']:
            b=bundle(count=3,mode=mode)
            if mode!='noe':b['fit_text']=b['fit_text'].replace('100 3','nan 3',1)
            prepare(self.root,b);loaded=tools.load_analysis_csv(self.csv);prepare(self.root,loaded)
            self.assertEqual(loaded['extracted'],b['extracted'])
    def test_plot_area_ratio_and_reported_prediction(self):
        self.assertEqual(tools.decay_panel_layout(),['set lmargin at screen 0.14','set rmargin at screen 0.96',
            'set tmargin at screen 0.9','set bmargin at screen 0.3','set format x ""','unset xlabel'])
        def margins(resid):
            lines=tools.decay_panel_layout(resid,top=.82,bottom=.11)
            return float(lines[2].split()[-1])-float(lines[3].split()[-1])
        self.assertAlmostEqual(margins(False)/margins(True),4.)
        with self.csv.open() as f:rows=list(csv.DictReader(f))
        for r in rows:
            predicted=float(r['i0'])*math.exp(-float(r['value'])*float(r['delay_s']))
            self.assertEqual(float(r['prediction']),predicted)
            self.assertEqual(float(r['residual']),float(r['intensity'])-predicted)


if __name__=='__main__': unittest.main()
