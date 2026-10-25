"""Synthetic etaXY mathematics, acquisition layout, extraction and publication."""
import argparse
import csv
import math
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
import zipfile

RELAX = Path(__file__).resolve().parents[1]
if (Path(__file__).parent/'eta_xy.py').exists():
    RELAX = Path(__file__).parent  # also allow testing the prepared source in /tmp
sys.path.insert(0,str(RELAX))
import eta_xy as eta


PULSE = '''"l5=1"
"d15=5u"
"d16=300u"
"d24 = d20*0.5-d15-p16-d16-p2-p2*0.637"
"d27 = 0.5*(d24+d15+p16+d16-p14)"
if "l5 %4 == 1" goto 20
if "l5 %4 == 2" goto 20
30 d24
 d24
if "l5 %4 == 1" goto 40
if "l5 %4 == 3" goto 40
 MCWRK do:f3 wr #0 if #0 zd iu5
 lo to LBLF1I1 times 4
ph10 = {0}*4 {2}*4
ph11 = {0}*1 {2}*1
ph14 = {0}*2 {2}*2
'''


def acquisition(root):
    root.mkdir(parents=True,exist_ok=True)
    pulses=[1.]*17;pulses[2]=27.;pulses[16]=1000.
    delays=[0.]*21;delays[20]=.05
    (root/'acqus').write_text('''##$PULPROG= <nico_nnh_tr_noSE>
##$TD= 16
##$PARMODE= 1
##$AQ_mod= 3
##$DTYPA= 0
##$BYTORDA= 0
##$DECIM= 2760
##$DSPFVS= 20
##$GRPDLY= 67.98
##$SFO1= 600.48
##$O1= 2820
##$SW_h= 7211.5
##$P= (0..16)
'''+ ' '.join(map(str,pulses))+'\n##$D= (0..20)\n'+' '.join(map(str,delays))+'\n##END=\n')
    (root/'acqu2s').write_text('''##$TD= 16
##$FnMODE= 5
##$SFO1= 60.85
##$O1= 7200
##$SW_h= 1200
##END=
''')
    (root/'pulseprogram').write_text(PULSE)
    (root/'ser').write_bytes(b''.join(bytes([i])*1024 for i in range(16)))


def pipe_bytes(dim=2,planes=1):
    h=[0.]*512;h[2]=2.345;h[9]=dim;h[15]=planes;h[26]=3;h[99]=4;h[219]=4
    return struct.pack('<512f',*h)+bytes(64)


def records(path):
    with Path(path).open(newline='') as handle:
        reader=csv.reader(handle);fields=next(reader)
        return [dict(zip(fields,row)) for row in reader if not row[0].startswith('# ')]


class EtaMathTests(unittest.TestCase):
    def test_known_rate_and_transfer_efficiency_cancellation(self):
        for rate in [.01,2.,15.]:
            q=math.tanh(rate*.05)
            # Excitation/detection gains factorize by initial/final state.
            a,b,c,d=2.,7.,3.,5.
            signals=[a*c,a*d*q,b*c*q,b*d]
            r=eta.calculate(signals,[.0001]*4,.05)
            self.assertTrue(r['valid'])
            self.assertAlmostEqual(r['eta_xy_abs_s_1'],rate,places=10)
            negative=eta.calculate([-v for v in signals],[.0001]*4,.05)
            self.assertEqual(r['eta_xy_abs_s_1'],negative['eta_xy_abs_s_1'])
            # A sign change in one preparation/detection channel also cancels.
            signed=eta.calculate([signals[0],signals[1],-signals[2],-signals[3]], [.0001]*4,.05)
            self.assertTrue(signed['valid'])
            self.assertEqual(r['eta_xy_abs_s_1'],signed['eta_xy_abs_s_1'])
    def test_uncertainty_matches_independent_numerical_derivatives(self):
        y=[1000.,120.,80.,900.];s=[2.,3.,4.,5.];time=.05
        derivatives=[]
        for i in range(4):
            lo=y.copy();hi=y.copy();step=abs(y[i])*1e-5
            lo[i]-=step;hi[i]+=step
            derivatives.append((eta.calculate(hi,s,time)['eta_xy_abs_s_1']-
                                eta.calculate(lo,s,time)['eta_xy_abs_s_1'])/(2*step))
        expected=math.sqrt(sum((d*sd)**2 for d,sd in zip(derivatives,s)))
        self.assertAlmostEqual(eta.calculate(y,s,time)['eta_xy_abs_error_s_1'],expected,places=9)
    def test_signed_products_domain_and_weak_signals(self):
        cases=[([100.,0.,1.,100.],'zero_intensity'),
               ([100.,-10.,10.,100.],'inconsistent_product_sign'),
               ([100.,10.,10.,-100.],'inconsistent_product_sign'),
               ([100.,110.,110.,100.],'ratio_outside_model_domain')]
        for y,flag in cases:
            with self.subTest(flag=flag):
                r=eta.calculate(y,[1.]*4,.05)
                self.assertFalse(r['valid']);self.assertIsNone(r['eta_xy_abs_s_1'])
                self.assertIn(flag,r['flags'])
        r=eta.calculate([100.,-10.,-10.,100.],[1.]*4,.05)
        self.assertTrue(r['valid'])  # both negative cross peaks are not repaired or discarded
        r=eta.calculate([100.,.1,.1,100.],[1.]*4,.05)
        self.assertTrue(r['valid']);self.assertIn('low_cross_snr',r['flags'])
        self.assertIn('nonlinear_uncertainty',r['flags'])
    def test_invalid_inputs_and_single_delay_factor(self):
        for noise,time in [([0.]*4,.05),([math.nan]*4,.05),([1.]*4,0.)]:
            with self.assertRaises(ValueError):eta.calculate([100.,10.,10.,100.],noise,time)
        r=eta.calculate([100.,10.,10.,100.],[1.]*4,.05)
        self.assertAlmostEqual(r['eta_xy_abs_s_1'],math.atanh(.1)/.05)


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)/'collection'/'6';acquisition(self.root)
        self.peaks=self.root/'etaXY.list'
        self.peaks.write_text('assignment w1 w2 intensity overlap\n'
            '?-E174N-E174H 120.1 8.1 999 False\n'
            '?-E174N-E174H 120.2 8.2 999 True\n'
            '# noise_sd: 9999\n# noise_method: MAD\n# spectrum_scale: 1\n')
        self.noise=self.root/'noise.txt';self.noise.write_text('4 1 2 3 4\n')
        self.data=self.root/'m.ft3';self.data.write_bytes(pipe_bytes(3,4))
        self.args=argparse.Namespace(experiment=str(self.root),peaks=str(self.peaks),data=str(self.data),
            extracted=None,noise=str(self.noise),output_dir=str(self.root),prefix='eta')
        self.calls=[]
    def command(self,command,cwd,stdout=None):
        self.calls.append(list(map(str,command)))
        if command[0]=='data_xtr':
            lines=(Path(cwd)/'peaks.list').read_text().splitlines()[1:]
            for i,line in reversed(list(enumerate(lines))):
                stdout.write(line.split()[0]+(' 1000 100 100 1000\n' if i==0 else ' 900 90 90 900\n'))
        elif command[0]=='gnuplot':
            (Path(cwd)/'eta.pdf').write_bytes(b'%PDF-synthetic')
        else:raise AssertionError('Unexpected external command '+str(command))
    def test_acquisition_profile_and_raw_deinterleaving(self):
        p=eta.experiment_profile(self.root)
        self.assertEqual(p['time'],.05);self.assertEqual(p['rows']//8,2)
        raw=(self.root/'ser').read_bytes()
        folder=self.root/'stage';folder.mkdir()
        outputs=eta.split_raw(p,folder)
        for i,path in enumerate(outputs):
            self.assertEqual(path.read_bytes(),b''.join(bytes([j])*1024 for j in range(i,16,4)))
        self.assertEqual((self.root/'ser').read_bytes(),raw)
        (self.root/'ser').write_bytes(raw[:-1])
        with self.assertRaisesRegex(ValueError,'ser size'):eta.split_raw(p,folder)
    def test_modified_profile_split_and_missing_records_are_rejected(self):
        path=self.root/'pulseprogram';original=path.read_text()
        path.write_text(original.replace('d20*0.5','d20'))
        with self.assertRaises(ValueError):eta.experiment_profile(self.root)
        path.write_text(original)
        (self.root/'audita.txt').write_text('created by split 4: dataset 2 of 4 from EXPNO 6')
        with self.assertRaisesRegex(ValueError,'split plane 2'):eta.experiment_profile(self.root)
        (self.root/'audita.txt').unlink();path.unlink()
        with self.assertRaises(FileNotFoundError):eta.experiment_profile(self.root)
    def test_one_2d_list_reordered_extraction_duplicates_and_three_outputs(self):
        with mock.patch.object(eta.rt,'checked_run',side_effect=self.command):eta.analyze(self.args)
        rows=records(self.root/'eta_numbers.csv')
        self.assertEqual(len(rows),2)
        self.assertEqual(rows[0]['assignment'],rows[1]['assignment'])
        self.assertNotEqual(rows[0]['w1'],rows[1]['w1'])
        self.assertEqual(rows[0]['overlap'],'False');self.assertEqual(rows[1]['overlap'],'True')
        self.assertAlmostEqual(float(rows[0]['eta_xy_abs_s-1']),math.atanh(.1)/.05)
        self.assertEqual(rows[0]['plane_01_intensity'],'1000.0')
        self.assertEqual([c[0] for c in self.calls],['data_xtr','gnuplot'])
        self.assertEqual(self.calls[0][self.calls[0].index('-pfmt')+1],'sparky2d')
        text=(self.root/'eta_numbers.csv').read_text()
        self.assertIn('# plane_noise_sd: 1.0, 2.0, 3.0, 4.0',text)
        self.assertNotIn('analysis_id',text)
        with zipfile.ZipFile(self.root/'eta_support.zip') as archive:
            self.assertEqual(archive.read('original_peak_list.list'),self.peaks.read_bytes())
            self.assertEqual(archive.read('pulseprogram'),(self.root/'pulseprogram').read_bytes())
        self.assertEqual({p.name for p in self.root.glob('eta*') if p.name != 'etaXY.list'},
                         {'eta.pdf','eta_numbers.csv','eta_support.zip'})
    def test_saved_intensities_need_no_extractor_or_fitter(self):
        with mock.patch.object(eta.rt,'checked_run',side_effect=self.command):eta.analyze(self.args)
        portable=self.root/'portable';portable.mkdir()
        with zipfile.ZipFile(self.root/'eta_support.zip') as archive:archive.extractall(portable)
        self.calls=[]
        self.args.extracted=str(portable/'extracted_ids.out')
        self.args.experiment=str(portable)
        self.args.peaks=str(portable/'original_peak_list.list')
        self.args.noise=str(portable/'vEstNoise.txt')
        with mock.patch.object(eta.rt,'checked_run',side_effect=self.command):eta.analyze(self.args)
        self.assertEqual([c[0] for c in self.calls],['gnuplot'])
        text=(self.root/'eta_numbers.csv').read_text()
        self.assertIn('# source_folder: collection',text)
        self.assertIn('# expno: 6',text)
        self.assertIn('# peak_list: etaXY.list',text)
        mapping=portable/'peaks.list';mapping.write_text(mapping.read_text().replace('120.1','121.1'))
        with mock.patch.object(eta.rt,'checked_run') as command,self.assertRaises(ValueError):eta.analyze(self.args)
        command.assert_not_called()
    def test_three_coordinate_list_ignores_plane_coordinate(self):
        self.peaks.write_text('assignment w1 w2 w3 intensity overlap\n'
            '?-E174N-E174H 10 120.1 8.1 999 False\n'
            '?-E174N-E174H 10 120.2 8.2 999 True\n')
        with mock.patch.object(eta.rt,'checked_run',side_effect=self.command):eta.analyze(self.args)
        self.assertEqual(self.calls[0][self.calls[0].index('-pfmt')+1],'sparky3d')
        rows=records(self.root/'eta_numbers.csv')
        self.assertEqual(rows[0]['w2'],'120.1')
        self.assertEqual(rows[0]['w3'],'8.1')
        self.assertIn('# w1: 10',(self.root/'eta_numbers.csv').read_text())
    def test_missing_extracted_identifier_aborts_before_rendering(self):
        (self.root/'eta_numbers.csv').write_bytes(b'previous')
        def missing(command,cwd,stdout=None):
            self.assertEqual(command[0],'data_xtr')
            stdout.write('FT000001-H 1000 100 100 1000\n')
        with mock.patch.object(eta.rt,'checked_run',side_effect=missing),self.assertRaises(ValueError):
            eta.analyze(self.args)
        self.assertEqual((self.root/'eta_numbers.csv').read_bytes(),b'previous')
    def test_render_failure_and_bad_counts_keep_existing_reports(self):
        for name in ('eta.pdf','eta_numbers.csv','eta_support.zip'):(self.root/name).write_bytes(b'previous')
        def broken(command,cwd,stdout=None):
            if command[0]=='gnuplot':raise subprocess.CalledProcessError(9,command)
            return self.command(command,cwd,stdout)
        with mock.patch.object(eta.rt,'checked_run',side_effect=broken),self.assertRaises(subprocess.CalledProcessError):
            eta.analyze(self.args)
        for name in ('eta.pdf','eta_numbers.csv','eta_support.zip'):
            self.assertEqual((self.root/name).read_bytes(),b'previous')
        self.noise.write_text('3 1 2 3\n')
        with mock.patch.object(eta.rt,'checked_run') as command,self.assertRaises(ValueError):eta.analyze(self.args)
        command.assert_not_called()
    def test_intermediate_pipeline_failure_is_not_hidden(self):
        with self.assertRaises(subprocess.CalledProcessError):
            eta.pipeline([[sys.executable,'-c','raise SystemExit(7)'],
                          [sys.executable,'-c','import sys; sys.stdin.read()']],self.root)
    def test_failed_processing_leaves_old_outputs_and_raw_untouched(self):
        out=self.root/'eta_xy_processed';out.mkdir();(out/'m.ft3').write_bytes(b'old stack')
        raw=(self.root/'ser').read_bytes()
        args=argparse.Namespace(experiment=str(self.root),output_dir=str(out),h_zf=4096,n_zf=1024,
                               h_p0=0.,h_p1=0.,n_p0=0.,n_p1=0.)
        def conversion(command,cwd,stdout=None):
            self.assertEqual(command[0],'bruk2pipe')
            self.assertEqual(command[command.index('-yN')+1],'4')
            self.assertEqual(command[command.index('-yT')+1],'2')
            self.assertEqual(command[command.index('-yMODE')+1],'States-TPPI')
            Path(command[command.index('-out')+1]).write_bytes(pipe_bytes())
        with mock.patch.object(eta.rt,'checked_run',side_effect=conversion),\
                mock.patch.object(eta,'pipeline',side_effect=subprocess.CalledProcessError(7,'nmrPipe')),\
                self.assertRaises(subprocess.CalledProcessError):eta.process(args)
        self.assertEqual((out/'m.ft3').read_bytes(),b'old stack')
        self.assertEqual((self.root/'ser').read_bytes(),raw)
    def test_successful_processing_uses_common_recipe_and_states_tppi_once(self):
        out=self.root/'eta_xy_processed'
        args=argparse.Namespace(experiment=str(self.root),output_dir=str(out),h_zf=4096,n_zf=1024,
                               h_p0=10.,h_p1=20.,n_p0=30.,n_p1=40.)
        pipelines=[]
        def command(command,cwd,stdout=None):
            name=command[0]
            if name=='bruk2pipe':
                self.assertEqual(command[command.index('-yMODE')+1],'States-TPPI')
                self.assertEqual(command[command.index('-yN')+1],'4')
                Path(command[command.index('-out')+1]).write_bytes(pipe_bytes())
            elif name=='sethdr':
                p=Path(command[1]);data=p.read_bytes();h=list(struct.unpack('<512f',data[:2048]))
                for flag,index in [('-ndim',9),('-zN',15)]:
                    if flag in command:h[index]=float(command[command.index(flag)+1])
                p.write_bytes(struct.pack('<512f',*h)+data[2048:])
            elif name=='setfdata':pass
            elif name=='xyz2pipe':stdout.write(pipe_bytes(3,4))
            elif name=='pipe2ucsf':(Path(cwd)/'proc.ucsf').write_bytes(b'synthetic ucsf')
            elif name=='specStat.com':stdout.write('vEstNoise 10.0\n')
            else:raise AssertionError('Unexpected external command '+str(command))
        def pipe(commands,cwd):
            pipelines.append(commands)
            Path(commands[-1][commands[-1].index('-out')+1]).write_bytes(pipe_bytes())
        with mock.patch.object(eta.rt,'checked_run',side_effect=command),\
                mock.patch.object(eta,'pipeline',side_effect=pipe),\
                mock.patch.object(eta.rt.subprocess,'run',side_effect=lambda cmd,**kwargs:command(cmd,self.root)):
            eta.process(args)
        self.assertEqual(len(pipelines),4)
        for pipeline in pipelines:
            self.assertEqual(sum('-alt' in c for c in pipeline),1)
            phases=[c for c in pipeline if 'PS' in c]
            self.assertEqual(phases[0][phases[0].index('-p0')+1],10.)
            self.assertEqual(phases[1][phases[1].index('-p0')+1],30.)
            self.assertFalse(any('MAC' in c for c in pipeline))
        self.assertEqual(len(list((out/'planes').glob('*.ft2'))),4)
        self.assertEqual(eta.rt.plane_count(out/'m.ft3'),4)
        self.assertEqual(eta.rt.counted_vector(out/'vEstNoise.txt'),[10.]*4)
        self.assertFalse((out/'raw').exists())
        self.assertFalse((out/'fid').exists())


if __name__=='__main__':unittest.main()
