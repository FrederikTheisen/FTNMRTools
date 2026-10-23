"""Synthetic-only regression tests. External processing/fitting commands are stubs."""
from contextlib import contextmanager
import csv
import importlib.util
import json
import math
import os
from pathlib import Path
import signal
import struct
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

RELAX = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('relax_tools', RELAX / 'relax_tools.py')
tools = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tools)


def header(dim=3, planes=3, vectors=6, endian='<'):
    values = [0.0] * 512
    values[2], values[9], values[15], values[26], values[99], values[219] = 2.345, dim, planes, 3, 4, vectors
    return struct.pack(endian + '512f', *values) + struct.pack(endian + '4f', 0, 1, 2, 3)


@contextmanager
def working(path):
    old = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)


PEAKS = '''assignment w1 w2 w3 intensity overlap comment
?-D151N-D151H 10 120.1 8.0 100 False
?-E174N-E174H 10 122.1 8.4 80 True 174
?-E174N-E174H 10 122.2 8.5 90 False
# noise_sd: 2
# noise_method: MAD
# spectrum_scale: 1.0
# source_folder: 20250101_sample
'''

STUB = r'''#!/usr/bin/python3
import json, math, os
from pathlib import Path
import signal, struct, sys
name = Path(sys.argv[0]).name
args = sys.argv[1:]
fn = args[args.index('-fn')+1] if '-fn' in args else ''
with open(os.environ['FT_TEST_LOG'], 'a') as f: f.write(json.dumps([name, args])+'\n')
fail = os.environ.get('FT_TEST_FAIL')
if fail and fail in (name, fn): sys.exit(7)
def header(dim=2, count=3):
    h=[0.0]*512
    h[2]=2.345;h[9]=dim;h[15]=count;h[26]=3;h[99]=4;h[219]=6
    return struct.pack('<512f',*h)+struct.pack('<4f',0,1,2,3)
if name == 'bruk2pipe':
    Path(args[args.index('-out')+1]).write_bytes(header())
elif name == 'nmrPipe':
    data=Path(args[args.index('-in')+1]).read_bytes() if '-in' in args else sys.stdin.buffer.read()
    if '-out' in args:
        output=Path(args[args.index('-out')+1])
        if fail=='grid_count' and str(output)=='ft/test002.ft2':
            h=list(struct.unpack('<512f',data[:2048]));h[99]=8
            data=struct.pack('<512f',*h)+data[2048:]
        output.write_bytes(data)
    else:sys.stdout.buffer.write(data)
elif name == 'xyz2pipe':
    n=len(list(Path('ft').glob('test*.ft2')))
    if fail == 'plane_count':n+=1
    sys.stdout.buffer.write(header(3,n))
elif name == 'pipe2ucsf':Path(args[-1]).write_bytes(b'SYNTHETIC UCSF')
elif name == 'specStat.com':print('vEstNoise', 'nan' if fail == 'bad_noise' else '2.0')
elif name == 'sethdr' or name == 'setfdata':
    path=Path(args[0]);data=path.read_bytes();h=list(struct.unpack('<512f',data[:2048]))
    if fail=='series_header' and '-ndim' in args and args[args.index('-ndim')+1]=='3':sys.exit(7)
    for flag,index in (('-ndim',9),('-zN',15),('-fileCount',442)):
        if flag in args:h[index]=float(args[args.index(flag)+1])
    path.write_bytes(struct.pack('<512f',*h)+data[2048:])
elif name in ('series.com','ucsfdata'):pass
elif name == 'data_xtr':
    lines=Path(args[args.index('-peak')+1]).read_text().splitlines()[1:]
    for line in reversed(lines):
        tag=line.split()[0]
        print(tag, '100 80 90' if fail != 'intensity_count' else '100 80')
elif name == 'rlx_dfit':
    if fail == 'fit_interrupted':os.kill(os.getpid(),signal.SIGTERM)
    if fail == 'empty_fit':sys.exit(0)
    mode=args[args.index('-mode')+1]
    lines=Path(args[args.index('-data')+1]).read_text().splitlines()[1:]
    if fail == 'missing_fit':lines=lines[:-1]
    if fail == 'unexpected_fit':lines.append('FT999999-N-H 1 2 3')
    with open(args[args.index('-out')+1],'w') as f:
        f.write('# Synthetic fitter output\n')
        for line in reversed(lines):
            tag=line.split()[0]
            f.write('# Warning: Data uncertainties scaled with 1.1702705.2 for site '+tag+'\n')
            if mode=='noe':f.write(tag+' -0.3 0.04\n')
            elif fail == 'zero_i0' and tag.startswith('FT000001'):
                f.write(tag+' 2 0.1 0 3 5 7\n')
            else:f.write(tag+' 2 0.1 100 3 5 7\n')
elif name == 'gnuplot':
    text=Path('plot.gnu').read_text()
    assert '\nfit ' not in text
    import re
    output=re.search(r'^set output "([^"]+)"',text,re.M).group(1)
    Path(output).write_bytes(b'%PDF synthetic plot\n')
else:raise RuntimeError(name)
'''


class DelayTests(unittest.TestCase):
    def test_units_and_scientific_notation(self):
        for value, expected in [('0.25', .25), ('1',1), ('1s',1), ('100m',.1), ('100ms',.1),
                                ('100u',.0001), ('100us',.0001), ('100µs',.0001), ('100μs',.0001),
                                ('2.5e-1s',.25), ('1e3us',.001), ('0m',0)]:
            with self.subTest(value=value):self.assertAlmostEqual(tools.parse_delay(value),expected)
    def test_rejects_bad_delays(self):
        for value in ('nan','inf','-1m','1e999s','1hour','','1 2','2ss'):
            with self.subTest(value=value), self.assertRaises(ValueError):tools.parse_delay(value)
    def test_order_repeats_and_cycle_configuration(self):
        with tempfile.TemporaryDirectory() as d, working(d):
            p=Path(d);(p/'vdlist').write_text('100m\n0m\n100m\n25m\n')
            self.assertEqual(tools.create_delays(p,'vclist','vdlist','ncyc_list'),[.1,0,.1,.025])
            (p/'vclist').write_text('0\n2\n1\n')
            with self.assertRaises(ValueError):tools.create_delays(p,'vclist','vdlist','ncyc_list')
            actual=tools.create_delays(p,'vclist','vdlist','ncyc_list','vc',.01,.01)
            self.assertEqual(actual,[.01,.03,.02])
            (p/'ncyc_list').write_text('0\n2\n1\n')
            with self.assertRaises(ValueError):tools.create_delays(p,'vclist','vdlist','ncyc_list','ncyc')
            actual=tools.create_delays(p,'vclist','vdlist','ncyc_list','ncyc',cycle_seconds='.005')
            self.assertEqual(actual,[0,.01,.005])
    def test_missing_empty_and_multiple_sources(self):
        with tempfile.TemporaryDirectory() as d, working(d):
            p=Path(d)
            with self.assertRaises(ValueError):tools.create_delays(p,'vc','vd','nc')
            (p/'vd1').write_text('')
            with self.assertRaises(ValueError):tools.create_delays(p,'vc','vd1','nc')
            (p/'vd2').write_text('1s\n')
            with self.assertRaises(ValueError):tools.create_delays(p,'vc','vd*','nc')
    def test_vector_and_noise_validation(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'noise'
            for text in ('3 1 2','2 1 -1','2 nan 1','2 0 1'):
                p.write_text(text)
                with self.assertRaises(ValueError):tools.counted_vector(p,positive=True)
            p.write_text('vEstNoise 1.25e+09\n')
            self.assertEqual(tools.noise_value(p),1.25e9)
            for text in ('vEstNoise -2','vEstNoise 0','vEstNoise NaN','unrelated 2'):
                p.write_text(text)
                with self.assertRaises(ValueError):tools.noise_value(p)
    def test_nmrpipe_header_endian_and_plane_count(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'series'
            for endian in ('<','>'):
                p.write_bytes(header(endian=endian));self.assertEqual(tools.plane_count(p),3)
            p.write_bytes(header(dim=2))
            with self.assertRaises(ValueError):tools.plane_count(p)
            p.write_bytes(b'bad header')
            with self.assertRaises(ValueError):tools.plane_count(p)


class PeakAndFitTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.path=Path(self.temp.name);self.peakfile=self.path/'peaks.list';self.peakfile.write_text(PEAKS)
        self.peaks=tools.read_peaks(self.peakfile,'sparky3d')
    def test_duplicate_labels_preserved_and_exact_duplicates_rejected(self):
        self.assertEqual(len(self.peaks),3)
        self.assertEqual(len({p['peak_id'] for p in self.peaks}),3)
        self.assertTrue(self.peaks[1]['duplicate_assignment'])
        self.assertFalse(self.peaks[2]['overlap'])
        self.peakfile.write_text(PEAKS+'?-D151N-D151H 10.0 120.10 8.00 100 False\n')
        with self.assertRaises(ValueError):tools.read_peaks(self.peakfile,'sparky3d')
    def test_reordered_extraction_and_missing_or_unexpected_ids(self):
        text=''.join(p['peak_id']+' 100 80 90\n' for p in reversed(self.peaks))
        result=tools.parse_extracted(text,self.peaks,3)
        self.assertEqual(result[self.peaks[0]['peak_id']],[100,80,90])
        for bad in (text.splitlines()[0], text+'WRONG-N-H 1 2 3\n',text+text.splitlines()[0],text.replace('100 80 90','100 80')):
            with self.assertRaises(ValueError):tools.parse_extracted(bad,self.peaks,3)
    def test_native_warning_is_preserved_verbatim(self):
        text=(Path(__file__).parent/'fixtures/native_exp2.out').read_text()
        values,warnings=tools.parse_fit(text,'r1')
        self.assertAlmostEqual(values['?-A155N-A155H'][0],1.19894)
        self.assertEqual(warnings['?-A155N-A155H'],['# Warning: Data uncertainties scaled with 1.1702705.2 for site ?-A155N-A155H'])
    def fit_text(self,zero=False):
        return ''.join('# Warning: Data uncertainties scaled with 1.1702705.2 for site '+p['peak_id']+'\n'+
                       p['peak_id']+f' 2 0.1 {0 if zero and i==0 else 100} 3 5 7\n'
                       for i,p in reversed(list(enumerate(self.peaks))))
    def report(self,mode='r1',zero=False):
        extracted={p['peak_id']:[100,80,90] for p in self.peaks}
        text=self.fit_text(zero) if mode!='noe' else ''.join(p['peak_id']+' -0.3 0.04\n' for p in reversed(self.peaks))
        return tools.prepare_reports(self.path,self.peaks,extracted,text,[0,.2,.1],[2,3,4],mode,'test','test.pdf','sample/5')
    def test_reported_parameters_normalization_and_residuals(self):
        self.report()
        numeric=[list(map(float,l.split())) for l in (self.path/'intensity.out').read_text().splitlines() if l and not l.startswith('#')]
        self.assertAlmostEqual(numeric[1][1],.8)
        self.assertAlmostEqual(numeric[1][2],.03)
        prediction=100*math.exp(-2*.2)
        self.assertAlmostEqual(numeric[1][3],prediction/100)
        self.assertAlmostEqual(numeric[1][5],(80-prediction)/3)
        plot=(self.path/'plot.gnu').read_text()
        self.assertIn('exp(-(2)*x)',plot);self.assertNotIn('\nfit ',plot)
        filtered=(self.path/'test_results_filtered.out').read_text().splitlines()[1:]
        self.assertEqual([l.split()[0] for l in filtered],['151','174'])
        with (self.path/'test_fit_diagnostics.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        self.assertEqual(rows[1]['overlap'],'True');self.assertEqual(rows[2]['overlap'],'False')
        self.assertIn('uncertainty_rescaled',rows[0]['flags'])
    def test_invalid_i0_raw_panel_and_filter(self):
        self.report(zero=True)
        self.assertIn('Raw intensity (invalid fit/I0)',(self.path/'plot.gnu').read_text())
        filtered=(self.path/'test_results_filtered.out').read_text()
        self.assertNotIn('\n151 ',filtered)
        plot=(self.path/'plot.gnu').read_text()
        overview=plot.split('set xlabel "Residue"')[1]
        self.assertNotIn('\n151 ',overview)
        self.assertIn('Residual unavailable',(self.path/'plot.gnu').read_text())
    def test_nonfinite_i0_keeps_raw_panel(self):
        extracted={p['peak_id']:[100,80,90] for p in self.peaks}
        text=self.fit_text().replace(self.peaks[0]['peak_id']+' 2 0.1 100',self.peaks[0]['peak_id']+' 2 0.1 nan')
        tools.prepare_reports(self.path,self.peaks,extracted,text,[0,.2,.1],[2,3,4],'r1','x','x.pdf','sample')
        self.assertIn('Raw intensity (invalid fit/I0)',(self.path/'plot.gnu').read_text())

    def test_noe_never_uses_exponential_or_decay_residual(self):
        self.report(mode='noe')
        plot=(self.path/'plot.gnu').read_text()
        self.assertNotIn('exp(',plot);self.assertNotIn('Residual',plot);self.assertNotIn('\nfit ',plot)
        self.assertIn('hetNOE',plot)
    def test_bad_fit_records_rejected(self):
        for text in ('','BAD 1 2','x 1 2 3 4 5 6\nx 1 2 3 4 5 6\n'):
            with self.assertRaises(ValueError):
                tools.prepare_reports(self.path,self.peaks,{},text,[0,1,2],[1,1,1],'r1','x','x.pdf','sample')


class PublicationTests(unittest.TestCase):
    def test_rollback_on_publication_failure(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);stage=root/'stage';dest=root/'dest';stage.mkdir();dest.mkdir()
            for n in ('a','b'):(stage/n).write_text('new');(dest/n).write_text('old')
            original=Path.rename
            def fail_one(path,target):
                if path==stage/'b' and Path(target)==dest/'b':raise OSError('simulated publication failure')
                return original(path,target)
            with mock.patch.object(Path,'rename',fail_one),self.assertRaises(OSError):tools.publish(stage,dest,['a','b'])
            self.assertEqual((dest/'a').read_text(),'old');self.assertEqual((dest/'b').read_text(),'old')
    def test_concurrent_run_rejected(self):
        with tempfile.TemporaryDirectory() as d,tools.run_directory(d,'proc'):
            with self.assertRaises(ValueError):
                with tools.run_directory(d,'analyze'):pass
    def test_symlink_outputs_rejected_without_modifying_target(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);stage=p/'stage';dest=p/'dest';stage.mkdir();dest.mkdir()
            (stage/'a').write_text('new');(p/'target').write_text('old');(dest/'a').symlink_to(p/'target')
            with self.assertRaises(ValueError):tools.publish(stage,dest,['a'])
            self.assertEqual((p/'target').read_text(),'old')


class ScriptIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='ftnmr-stub-tests-');self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.source=self.root/'sample data';self.source.mkdir()
        self.bin=self.root/'bin';self.bin.mkdir();self.log=self.root/'calls.jsonl'
        for name in ('bruk2pipe','nmrPipe','series.com','sethdr','setfdata','xyz2pipe','pipe2ucsf','ucsfdata','specStat.com','data_xtr','rlx_dfit','gnuplot'):
            p=self.bin/name;p.write_text(STUB);p.chmod(0o755)
        self.env=dict(os.environ,PATH=str(self.bin)+':/usr/bin:/bin',NMRTXT=str(self.root),
                      FT_TEST_LOG=str(self.log),PYTHONDONTWRITEBYTECODE='1')
        (self.source/'ser').write_text('synthetic input only')
        (self.source/'vdlist').write_text('0s\n200ms\n100m\n')
        (self.source/'m.ft3').write_bytes(header())
        (self.source/'relax_header.txt').write_text('3 0 .2 .1\n')
        (self.source/'vEstNoise.txt').write_text('3 2 3 4\n')
        (self.source/'sample.list').write_text(PEAKS)
    def run_script(self,name,fail=None,script=None):
        env=dict(self.env)
        if fail:env['FT_TEST_FAIL']=fail
        return subprocess.run(['/bin/tcsh','-f',str(script or RELAX/name)],cwd=self.source,env=env,text=True,capture_output=True)
    def calls(self):
        return [json.loads(l) for l in self.log.read_text().splitlines()] if self.log.exists() else []
    def test_worker_cannot_run_directly_in_an_experiment(self):
        old=(self.source/'m.ft3').read_bytes()
        result=subprocess.run(['/bin/tcsh','-f',str(RELAX/'relax_proc.com'),'--worker',str(self.source)],
                              cwd=self.source,env=self.env,text=True,capture_output=True)
        self.assertNotEqual(result.returncode,0)
        self.assertEqual((self.source/'m.ft3').read_bytes(),old)
        self.assertEqual(self.calls(),[])

    def test_successful_processing_and_repeat_run(self):
        result=self.run_script('relax_proc.com');self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        self.assertEqual(tools.plane_count(self.source/'m.ft3'),3)
        self.assertEqual(tools.counted_vector(self.source/'vEstNoise.txt'),[2,2,2])
        self.assertEqual(len(list((self.source/'ft').glob('test*.ft2'))),3)
        result=self.run_script('relax_proc.com');self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        self.assertEqual((self.source/'ser').read_text(),'synthetic input only')
    def test_processing_failures_preserve_existing_outputs(self):
        old=(self.source/'m.ft3').read_bytes()
        (self.source/'ft').mkdir();(self.source/'ft/old.ft2').write_text('old spectrum')
        for failure in ('bruk2pipe','COADD','SOL','series_header','sethdr','setfdata','xyz2pipe','pipe2ucsf','ucsfdata','specStat.com','bad_noise','plane_count','grid_count'):
            with self.subTest(failure=failure):
                result=self.run_script('relax_proc.com',failure)
                self.assertNotEqual(result.returncode,0,result.stdout+result.stderr)
                self.assertEqual((self.source/'m.ft3').read_bytes(),old)
                self.assertEqual((self.source/'ft/old.ft2').read_text(),'old spectrum')
                self.assertFalse(list(self.source.glob('.relax-proc-*')))
    def test_analysis_reorders_ids_and_ignores_leftover_temporary_rows(self):
        (self.source/'pre_results.tmp').write_text('OLD 9 9\n')
        interrupted=self.source/'.relax-analyze-interrupted';interrupted.mkdir()
        (interrupted/'rlx_dfit.out').write_text('old interrupted fit')
        result=self.run_script('relax_analyze.com');self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        results=(self.source/'sample_results.out').read_text().splitlines()[1:]
        self.assertEqual([l.split()[0] for l in results],['151','174','174'])
        self.assertEqual((self.source/'pre_results.tmp').read_text(),'OLD 9 9\n')
        self.assertTrue((self.source/'sample_fit_diagnostics.tsv').is_file())
        before=(self.source/'sample_results.out').read_bytes()
        result=self.run_script('relax_analyze.com');self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        self.assertEqual((self.source/'sample_results.out').read_bytes(),before)
    def test_analysis_failure_and_interruption_preserve_old_results(self):
        (self.source/'sample_results.out').write_text('old results\n')
        (self.source/'rlx_dfit.out').write_text('old raw fit\n')
        for failure in ('data_xtr','rlx_dfit','fit_interrupted','empty_fit','missing_fit','unexpected_fit','intensity_count','gnuplot'):
            with self.subTest(failure=failure):
                result=self.run_script('relax_analyze.com',failure)
                self.assertNotEqual(result.returncode,0,result.stdout+result.stderr)
                self.assertEqual((self.source/'sample_results.out').read_text(),'old results\n')
                self.assertEqual((self.source/'rlx_dfit.out').read_text(),'old raw fit\n')
                self.assertFalse(list(self.source.glob('.relax-analyze-*')))
    def test_count_mismatch_rejected_before_external_commands(self):
        for header_text in ('2 0 1\n','3 0 nan 1\n','3 0 0 0\n'):
            (self.source/'relax_header.txt').write_text(header_text)
            result=self.run_script('relax_analyze.com');self.assertNotEqual(result.returncode,0)
        self.assertEqual(self.calls(),[])
        (self.source/'relax_header.txt').write_text('3 0 .2 .1\n')
        (self.source/'m.ft3').write_bytes(header(planes=2))
        result=self.run_script('relax_analyze.com');self.assertNotEqual(result.returncode,0)
        self.assertEqual(self.calls(),[])
    def test_ambiguous_peak_lists_rejected(self):
        (self.source/'extra.list').write_text(PEAKS)
        result=self.run_script('relax_analyze.com');self.assertNotEqual(result.returncode,0)
        self.assertEqual(self.calls(),[])
    def test_r2_and_noe_shell_modes(self):
        original=(RELAX/'relax_analyze.com').read_text()
        for mode in ('r2','noe'):
            script=self.root/f'{mode}.com';script.write_text(original.replace('set DATA_TYPE = r1',f'set DATA_TYPE = {mode}').replace('set HELPER = "$0:h/relax_tools.py"',f'set HELPER = "{RELAX}/relax_tools.py"'))
            if mode=='noe':(self.source/'relax_header.txt').write_text('3 1 -1 1\n')
            result=self.run_script('relax_analyze.com',script=script)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            if mode=='noe':self.assertNotIn('exp(', (self.source/'plot.gnu').read_text())
        fit_calls=[args for name,args in self.calls() if name=='rlx_dfit']
        self.assertEqual(fit_calls[-1][fit_calls[-1].index('-mode')+1],'noe')
        self.assertIn('10.0',fit_calls[0])


if __name__=='__main__':unittest.main()
