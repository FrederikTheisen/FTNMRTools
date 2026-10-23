"""Synthetic selection and publication tests; no experimental data are changed."""
import argparse
import copy
import csv
import json
import math
from pathlib import Path
import random
import shutil
import subprocess
import tempfile
import unittest
from unittest import mock

from test_relax_tools import tools, header, working


def synthetic(count=24, times=None, perturb=None):
    times = [.075,.001,.2,.025,.15,.025,.1] if times is None else times
    peaks=[];data={}
    for i in range(count):
        tag=f'FT{i+1:06d}-N-H'
        peaks.append(dict(peak_id=tag,assignment=f'?-E{160+i}N-E{160+i}H',residue=160+i,
                          coordinates=['10',str(120+i*.01),str(8+i*.001)],
                          overlap=False,duplicate_assignment=False))
        a=100+i;rate=2+i*.02
        data[tag]=[a*math.exp(-rate*t)+(perturb(i,t)*a if perturb else 0) for t in times]
    return peaks,data,times,[.05]*len(times)


def profile_fit(times, y, noises):
    """Independent scalar profiled least squares, solely a native-command test stub."""
    def evaluate(rate):
        x=[math.exp(-rate*t) for t in times]
        a=math.fsum(v*z/s**2 for v,z,s in zip(y,x,noises))/math.fsum(z*z/s**2 for z,s in zip(x,noises))
        return math.fsum(((v-a*z)/s)**2 for v,z,s in zip(y,x,noises)),a
    lo,hi=0.,50.
    ratio=(math.sqrt(5)-1)/2
    c=hi-ratio*(hi-lo);d=lo+ratio*(hi-lo)
    for _ in range(100):
        if evaluate(c)[0]<evaluate(d)[0]:hi=d;d=c;c=hi-ratio*(hi-lo)
        else:lo=c;c=d;d=lo+ratio*(hi-lo)
    rate=(lo+hi)/2;ss,a=evaluate(rate)
    return [rate,.01,a,.01,ss,float(len(times)-2)]


def fit_text(peaks,data,times,noises):
    return ''.join(p['peak_id']+' '+' '.join(f'{v:.17g}' for v in profile_fit(times,data[p['peak_id']],noises))+'\n'
                   for p in reversed(peaks))


def assess(bundle, policy='auto', candidate_mutation=None):
    peaks,data,times,noises=bundle
    full=tools.fitted_rows(peaks,fit_text(*bundle),'r1rho')
    first=tools.evaluate_early_point(peaks,data,full,times,noises,policy)
    if first['reason']!='candidate_required':return first
    used=[i for i,t in enumerate(times) if t!=min(times)]
    cand=tools.fitted_rows(peaks,fit_text(peaks,{k:[v[i] for i in used] for k,v in data.items()},
                                       [times[i] for i in used],[noises[i] for i in used]),'r1rho')
    if candidate_mutation:candidate_mutation(cand)
    return tools.evaluate_early_point(peaks,data,full,times,noises,policy,cand)


class DecisionTests(unittest.TestCase):
    def test_pure_exponential_retains_all_points(self):
        result=assess(synthetic())
        self.assertEqual(result['reason'],'shortest_residual_within_limit')
        self.assertEqual(result['fit_basis'],'all_points');self.assertFalse(result['candidate_attempted'])

    def test_noise_only_datasets_do_not_select_exclusion(self):
        for seed in range(10):
            rng=random.Random(seed)
            b=synthetic(perturb=lambda i,t:rng.gauss(0,.0004))
            self.assertEqual(assess(b)['fit_basis'],'all_points')

    def test_shared_excess_selects_shortest_not_first_plane(self):
        result=assess(synthetic(perturb=lambda i,t:.03 if t==.001 else 0))
        self.assertEqual(result['fit_basis'],'exclude_shortest')
        self.assertEqual(result['excluded_planes_1based'],[2])
        self.assertEqual(result['fit_count'],6)
        self.assertGreater(result['score_reduction_fraction'],.99)

    def test_repeated_shortest_delay_excludes_all_replicates(self):
        b=synthetic(times=[.075,.001,.2,.025,.15,.025,.1,.001],perturb=lambda i,t:.03 if t==.001 else 0)
        result=assess(b)
        self.assertEqual(result['excluded_planes_1based'],[2,8])
        self.assertEqual(result['remaining_distinct_delays'],5)

    def test_weighted_repeats_use_one_mean_and_fixed_amplitude(self):
        row=dict(peak_id='p',values=[2.,.1,100.,.1,1.,1.])
        times=[.1,.2,.1];noises=[1.,1.,2.]
        values=[100*math.exp(-.2)+1,100*math.exp(-.4),100*math.exp(-.2)+5]
        groups=tools.early_groups(times,noises,{'p':values},['p'],{'p':row},{'p':200.})
        self.assertAlmostEqual(groups[0]['median'],(1+5/4)/1.25/200)
        self.assertAlmostEqual(groups[0]['normalized_noise'],1/math.sqrt(1.25)/200)

    def test_few_outliers_and_overlap_only_bias_do_not_trigger(self):
        self.assertFalse(assess(synthetic(perturb=lambda i,t:.03 if i<3 and t==.001 else 0))['candidate_attempted'])
        b=synthetic(count=30,perturb=lambda i,t:.03 if i<5 and t==.001 else 0)
        for p in b[0][:5]:p['overlap']=True
        self.assertFalse(assess(b)['candidate_attempted'])

    def test_cohort_support_unique_residue_and_validity(self):
        self.assertEqual(assess(synthetic(count=19))['reason'],'insufficient_cohort')
        b=synthetic(count=21)
        b[0][0]['residue']=b[0][1]['residue']
        self.assertEqual(assess(b)['reason'],'insufficient_cohort')
        b=synthetic(times=[.001,.025,.075,.1,.2])
        self.assertEqual(assess(b)['reason'],'insufficient_remaining_delays')
        b=synthetic(count=20);full=tools.fitted_rows(b[0],fit_text(*b),'r1rho')
        full[0]['valid']=False
        self.assertEqual(tools.evaluate_early_point(b[0],b[1],full,b[2],b[3])['reason'],'insufficient_cohort')

    def test_sign_agreement_and_negative_early_bias(self):
        b=synthetic(perturb=lambda i,t:(.03 if i<16 else -.03) if t==.001 else 0)
        self.assertEqual(assess(b)['reason'],'insufficient_sign_agreement')
        self.assertEqual(assess(synthetic(perturb=lambda i,t:-.03 if t==.001 else 0))['fit_basis'],'exclude_shortest')

    def test_remaining_curvature_and_large_rate_change_rejected(self):
        b=synthetic(perturb=lambda i,t:.04 if t==.001 else .06 if t==.075 else 0)
        result=assess(b)
        self.assertEqual(result['fit_basis'],'all_points')
        self.assertIn(result['reason'],['insufficient_structure_improvement','remaining_systematic_residual'])
        result=assess(synthetic(perturb=lambda i,t:.4 if t==.001 else 0))
        self.assertEqual(result['reason'],'rate_change_exceeds_limit')

    def test_invalid_candidate_keeps_fixed_cohort(self):
        b=synthetic(perturb=lambda i,t:.03 if t==.001 else 0)
        result=assess(b,candidate_mutation=lambda rows:rows[0].update(valid=False))
        self.assertEqual(result['reason'],'invalid_candidate_cohort_fit')
        self.assertEqual(result['cohort_size'],24)

    def test_report_and_off_policies(self):
        b=synthetic(perturb=lambda i,t:.03 if t==.001 else 0)
        report=assess(b,'report');off=assess(b,'off')
        self.assertTrue(report['candidate_passed']);self.assertEqual(report['fit_basis'],'all_points')
        self.assertFalse(off['candidate_attempted']);self.assertEqual(off['reason'],'disabled')


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.b=synthetic(perturb=lambda i,t:.03 if t==.001 else 0)
        self.calls=[];self.failure=None
    def command(self,command,stage,stdout=None):
        self.calls.append(list(command))
        if command[0]=='data_xtr':
            for p in reversed(self.b[0]):stdout.write(p['peak_id']+' '+' '.join(map(str,self.b[1][p['peak_id']]))+'\n')
        elif command[0]=='rlx_dfit':
            inp=Path(stage)/command[command.index('-data')+1]
            lines=inp.read_text().splitlines();t=list(map(float,lines[0].split()[1:]))
            noises=tools.counted_vector(Path(stage)/command[command.index('-d_esd')+1],positive=True)
            data={l.split()[0]:list(map(float,l.split()[1:])) for l in lines[1:]}
            text=fit_text(self.b[0],data,t,noises)
            if 'without_shortest' in inp.name:
                if self.failure=='command':raise subprocess.CalledProcessError(7,command)
                if self.failure=='missing':text='\n'.join(text.splitlines()[1:])+'\n'
                if self.failure=='invalid_overlap':
                    tag=self.b[0][0]['peak_id']
                    text='\n'.join(tag+' 2 .1 0 .1 5 4' if l.startswith(tag+' ') else l for l in text.splitlines())+'\n'
            (Path(stage)/command[command.index('-out')+1]).write_text(text)
        elif command[0]=='gnuplot':(Path(stage)/'peaks.list_intensity.pdf').write_bytes(b'%PDF synthetic')
        else:raise AssertionError(command)
    def run_analysis(self,policy='auto'):
        peaks,data,times,noises=self.b
        (self.root/'peaks.list').write_text('assignment w1 w2 w3 intensity overlap comment\n'+''.join(
            p['assignment']+' '+' '.join(p['coordinates'])+' 1 '+str(p['overlap'])+'\n' for p in peaks))
        (self.root/'m.ft3').write_bytes(header(planes=len(times)))
        for name,vector in [('relax_header.txt',times),('vEstNoise.txt',noises)]:
            (self.root/name).write_text(str(len(vector))+' '+' '.join(map(str,vector))+'\n')
        args=argparse.Namespace(peaks='peaks.list',data='m.ft3',header='relax_header.txt',noise='vEstNoise.txt',
                                format='nmrpipe',mode='r1rho',peak_format='sparky3d',prefix='sample',early_point_policy=policy)
        with working(self.root),mock.patch.object(tools,'checked_run',side_effect=self.command):tools.analyze(args)
    def test_selected_outputs_preserve_observations_noise_and_provenance(self):
        self.run_analysis()
        audit=json.loads((self.root/'sample_early_point_summary.json').read_text())['selection']
        self.assertEqual(audit['fit_basis'],'exclude_shortest')
        self.assertEqual(len(tools.counted_vector(self.root/'vEstNoise.txt')),7)
        self.assertEqual(len(tools.counted_vector(self.root/'rlx_dfit_noise.txt')),6)
        self.assertEqual((self.root/'rlx_dfit.out').read_bytes(),(self.root/'rlx_dfit_without_shortest.out').read_bytes())
        self.assertEqual(len((self.root/'data_xtr.out').read_text().splitlines()),24)
        values=[l.split() for l in (self.root/'intensity.out').read_text().splitlines() if l and not l.startswith('#')]
        self.assertEqual(len(values),24*7);self.assertEqual(sum(v[6]=='0' for v in values),24)
        for name in ['sample_results.out','sample_results_filtered.out','sample_results_peak_overlap.out']:
            self.assertIn('# fit_basis: exclude_shortest',(self.root/name).read_text())
            self.assertIn('# excluded_planes_1based: [2]',(self.root/name).read_text())
        plot=(self.root/'plot.gnu').read_text()
        self.assertIn('($7==0?',plot);self.assertIn('same retained delays',plot)
        with (self.root/'sample_fit_diagnostics.tsv').open() as fp:rows=list(csv.DictReader(fp,delimiter='\t'))
        self.assertTrue(all(r['fit_basis']=='exclude_shortest' for r in rows))
        self.assertTrue(all(r['all_points_value'] and r['candidate_value'] for r in rows))
        loaded=tools.load_analysis_csv(self.root/'sample_analysis.csv')
        self.assertEqual(loaded['extracted'],self.b[1])
        self.assertEqual(loaded['context']['_early_point_selection']['audit'],audit)
        self.assertEqual(len([c for c in self.calls if c[0]=='rlx_dfit']),2)
    def test_csv_and_saved_refresh_retain_mask_without_refitting(self):
        self.run_analysis();before=(self.root/'sample_analysis.csv').read_bytes()
        bundle=tools.load_saved_analysis(self.root,'sample','r1rho')
        c=bundle['context'];dest=self.root/'portable';dest.mkdir()
        with mock.patch.object(tools,'checked_run',side_effect=AssertionError('No refitting')):
            tools.prepare_reports(dest,bundle['peaks'],bundle['extracted'],bundle['fit_text'],bundle['times'],
                 bundle['noises'],'r1rho','sample',c['pdf_name'],'source/3',context=c,warning_map=bundle['warning_map'])
        self.assertEqual((dest/'sample_analysis.csv').read_bytes(),before)
        loaded=tools.load_analysis_csv(dest/'sample_analysis.csv')
        self.assertEqual(loaded['context']['_early_point_selection']['audit']['excluded_planes_1based'],[2])
        tampered=self.root/'sample_analysis.csv'
        with tampered.open() as fp:rows=list(csv.DictReader(fp))
        rows[1]['used_in_fit']='True'
        with tampered.open('w',newline='') as fp:
            w=csv.DictWriter(fp,fieldnames=tools.EARLY_CSV_FIELDS);w.writeheader();w.writerows(rows)
        with self.assertRaises(ValueError):tools.load_analysis_csv(tampered)
    def test_report_only_and_off_keep_native_full_fit(self):
        for policy in ['report','off']:
            self.calls=[];self.run_analysis(policy)
            self.assertEqual((self.root/'rlx_dfit.out').read_bytes(),(self.root/'rlx_dfit_all_points.out').read_bytes())
            self.assertEqual(len(tools.counted_vector(self.root/'rlx_dfit_noise.txt')),7)
            audit=json.loads((self.root/'sample_early_point_summary.json').read_text())['selection']
            self.assertEqual(audit['policy'],policy)
            if policy=='off':self.assertIn('not attempted',(self.root/'rlx_dfit_without_shortest.out').read_text())
            self.assertEqual(len([c for c in self.calls if c[0]=='rlx_dfit']),2 if policy=='report' else 1)
            tools.load_analysis_csv(self.root/'sample_analysis.csv')
    def test_standalone_selected_csv_refresh_runs_renderer_only(self):
        self.run_analysis()
        dest=self.root/'standalone';dest.mkdir()
        path=dest/'sample_analysis.csv';shutil.copy2(self.root/'sample_analysis.csv',path)
        args=argparse.Namespace(csv=str(path),saved_prefix=None,experiment_type=None,
                                output_dir=str(dest),peak_list_name=None)
        self.calls=[]
        with mock.patch.object(tools,'checked_run',side_effect=self.command):tools.refresh_report(args)
        self.assertEqual([c[0] for c in self.calls],['gnuplot'])
        audit=tools.load_analysis_csv(path)['context']['_early_point_selection']['audit']
        self.assertEqual(audit['excluded_planes_1based'],[2])
        self.assertFalse((dest/'rlx_dfit.out').exists())
    def test_saved_archive_tampering_rejected(self):
        self.run_analysis()
        path=self.root/'rlx_dfit_all_points.out'
        lines=path.read_text().splitlines();fields=lines[0].split();fields[3]=str(float(fields[3])+1)
        lines[0]=' '.join(fields);path.write_text('\n'.join(lines)+'\n')
        with self.assertRaisesRegex(ValueError,'Archived native fit disagrees'):
            tools.load_saved_analysis(self.root,'sample','r1rho')
    def test_legacy_selected_schema_migrates_without_hashes_or_refitting(self):
        self.run_analysis()
        csv_path=self.root/'sample_analysis.csv';diagnostic=self.root/'sample_fit_diagnostics.tsv'
        with csv_path.open(newline='') as fp:rows=list(csv.DictReader(fp))
        legacy=dict(schema_version='2',analysis_id='a2-obsolete-id',fit_revision='sha256:old',
                    report_revision='sha256:old',fit_input_sha256='old',fit_output_sha256='old',
                    noise_sha256='old',assessment_method='old long method paragraph')
        audit=json.loads(rows[0]['early_point_audit_json'])
        audit['archive_sha256']={name:'obsolete digest' for name in tools.early_archive_names()}
        legacy['early_point_audit_json']=json.dumps(audit,sort_keys=True,separators=(',',':'))
        for row in rows:row.update(legacy)
        with csv_path.open('w',newline='') as fp:
            writer=csv.DictWriter(fp,fieldnames=tools.LEGACY_EARLY_CSV_FIELDS,extrasaction='ignore')
            writer.writeheader();writer.writerows(rows)
        with diagnostic.open(newline='') as fp:
            reader=csv.DictReader(fp,delimiter='\t');fields=reader.fieldnames;mapping=list(reader)
        before=fields.index('schema_version');after=fields.index('fit_basis')
        fields=fields[:before]+tools.LEGACY_TRACE_FIELDS+['early_point_audit_json']+fields[after:]
        for row in mapping:row.update(legacy)
        with diagnostic.open('w',newline='') as fp:
            writer=csv.DictWriter(fp,fieldnames=fields,delimiter='\t',extrasaction='ignore')
            writer.writeheader();writer.writerows(mapping)
        with mock.patch.object(tools,'checked_run',side_effect=AssertionError('No refitting')):
            bundle=tools.load_saved_analysis(self.root,'sample','r1rho')
            c=bundle['context'];tools.prepare_reports(self.root,bundle['peaks'],bundle['extracted'],bundle['fit_text'],
                bundle['times'],bundle['noises'],'r1rho','sample',c['pdf_name'],
                c['source_collection']+'/'+c['source_experiment'],context=c,warning_map=bundle['warning_map'])
        new=tools.load_analysis_csv(csv_path)
        self.assertEqual(new['extracted'],self.b[1])
        self.assertEqual(new['context']['_early_point_selection']['audit']['excluded_planes_1based'],[2])
        self.assertNotIn('sha256',csv_path.read_text())
        self.assertNotIn('sha256',diagnostic.read_text())
        self.assertNotIn('sha256',(self.root/'sample_early_point_summary.json').read_text())
        repeated=tools.load_saved_analysis(self.root,'sample','r1rho')
        self.assertEqual(repeated['context']['analysis_id'],new['context']['analysis_id'])
    def test_selection_audit_and_full_fit_tampering_rejected(self):
        self.run_analysis();path=self.root/'sample_analysis.csv';original=path.read_bytes()
        for kind in ['audit','full_fit']:
            path.write_bytes(original)
            with path.open() as fp:rows=list(csv.DictReader(fp))
            for row in rows:
                if kind=='audit':
                    audit=json.loads(row['early_point_audit_json']);audit['after_score']=999
                    row['early_point_audit_json']=json.dumps(audit,sort_keys=True,separators=(',',':'))
                elif row['peak_id']==rows[0]['peak_id']:row['all_points_i0']='99'
            with path.open('w',newline='') as fp:
                w=csv.DictWriter(fp,fieldnames=tools.EARLY_CSV_FIELDS);w.writeheader();w.writerows(rows)
            with self.assertRaises(ValueError):tools.load_analysis_csv(path)
    def test_candidate_command_or_malformed_output_preserves_old_results(self):
        old=self.root/'sample_results.out';old.write_bytes(b'previous result')
        for failure in ['command','missing']:
            self.failure=failure
            with self.assertRaises((ValueError,subprocess.CalledProcessError)):self.run_analysis()
            self.assertEqual(old.read_bytes(),b'previous result')
    def test_invalid_noncohort_fit_never_falls_back_to_full_fit(self):
        self.b[0][0]['overlap']=True;self.failure='invalid_overlap'
        self.run_analysis()
        with (self.root/'sample_fit_diagnostics.tsv').open() as fp:rows=list(csv.DictReader(fp,delimiter='\t'))
        self.assertEqual(rows[0]['fit_basis'],'exclude_shortest');self.assertEqual(rows[0]['valid_fit'],'False')
        self.assertEqual(float(rows[0]['i0']),0.)
    def test_candidate_artifact_publication_failure_rolls_back_every_output(self):
        self.run_analysis()
        before={p.name:p.read_bytes() for p in self.root.iterdir() if p.is_file()}
        rename=Path.rename
        def broken(path,target):
            if path.name=='rlx_dfit_without_shortest.out' and path.parent.name.startswith('.relax-analyze-'):
                raise OSError('candidate archive publish failed')
            return rename(path,target)
        with mock.patch.object(Path,'rename',broken),self.assertRaises(OSError):self.run_analysis('report')
        self.assertEqual(before,{name:(self.root/name).read_bytes() for name in before})


if __name__=='__main__':unittest.main()
