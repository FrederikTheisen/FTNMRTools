"""HetNOE demultiplexing, explicit N baseline, and isolated publication."""
import argparse
import array
import os
import math
import random
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

RELAX = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RELAX))
import hetnoe_proc as hn

PULSE = '''"l0=1"
"ST1CNT=trunc(td1/(2*2))"
if "l0 %2 == 1" { d1 50u } else { d2
4 (p0 ph1) 5m lo to 4 times l4 }
MCWRK do:f3 wr #0 if #0 zd iu0
lo to LBLF1I1 times 2
"ph1loop+=1"
lo to LBLF1 times 2
"t1loop+=1"
lo to LBLST1 times ST1CNT
'''
LEGACY = '''"l2=1"
if "l2 == 1" goto 9
8 11m (p1*2 ph0):f1
lo to 8 times l8
goto 10
9 11m d8*2 11m
lo to 9 times l8
10 go=2 ph31
wr #0 if #0 zd
1m iu2
lo to 3 times 2
1m ip8*2
1m igrad EA
1m ru2
lo to 4 times 2
1m id0
lo to 5 times l3
'''


def acquisition(root, datatype=0, byte_order=0, program='15N_HETNOE'):
    root.mkdir(parents=True, exist_ok=True)
    (root/'acqus').write_text(f'''##$PULPROG= <{program}>
##$TD= 18
##$PARMODE= 1
##$AQ_mod= 3
##$FnTYPE= 0
##$DTYPA= {datatype}
##$BYTORDA= {byte_order}
##$DECIM= 2760
##$DSPFVS= 20
##$GRPDLY= 67.98
##$SFO1= 600.48
##$O1= 2820
##$SW_h= 7211.5
##$WAVEMAKER= (0..1)
<not needed> <a quoted name with spaces>
##END=
''')
    (root/'acqu2s').write_text('''##$TD= 16
##$FnMODE= 6
##$SFO1= 60.85
##$O1= 7200
##$SW_h= 1200
##END=
''')
    (root/'pulseprogram').write_text(PULSE if program == '15N_HETNOE' else LEGACY)
    endian = '<' if byte_order == 0 else '>'
    raw = []
    for i in range(16):
        row = struct.pack(endian+('18i' if datatype == 0 else '18d'), *([i+0.25]*18 if datatype == 2 else [i]*18))
        raw.append(row+bytes((-len(row)) % 1024))
    (root/'ser').write_bytes(b''.join(raw))


def pipe_bytes(dim=2, planes=1):
    h = [0.]*512; h[2]=2.345; h[9]=dim; h[15]=planes; h[26]=3; h[99]=4; h[219]=4
    return struct.pack('<512f', *h)+bytes(64)


class HetNOEProcessingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)/'collection'/'5'; acquisition(self.root)
        self.out = self.root/'hetnoe_processed'
        self.args = argparse.Namespace(experiment=str(self.root), raw_rows=None, output_dir=str(self.out),
            h_zf=4096, n_zf=1024, h_p0=10., h_p1=20., n_p0=30., n_p1=40.,
            n_baseline_order=4, h_baseline_order=3)
    def test_reference_saturation_deinterleave_preserves_echo_pairs_and_padding(self):
        p=hn.profile(self.root); raw=(self.root/'ser').read_bytes()
        folder=self.root/'stage'; folder.mkdir()
        paths=hn.split_raw(p, folder)
        for j,path in enumerate(paths):
            self.assertEqual(path.read_bytes(), b''.join(raw[i*1024:(i+1)*1024] for i in range(j,16,2)))
        self.assertEqual(p['rows']//4,4)
        self.assertEqual((self.root/'ser').read_bytes(),raw)
    def test_float64_in_both_byte_orders_becomes_unscaled_native_float32(self):
        for endian in (0,1):
            acquisition(self.root,2,endian)
            folder=self.root/f'stage{endian}'; folder.mkdir()
            paths=hn.split_raw(hn.profile(self.root),folder)
            for j,path in enumerate(paths):
                data=path.read_bytes(); self.assertEqual(len(data),8*1024)
                for k in range(8):
                    values=array.array('f'); values.frombytes(data[k*1024:k*1024+18*4])
                    self.assertEqual(list(values),[j+k*2+.25]*18)
    def test_nonfinite_float_raw_rejected(self):
        acquisition(self.root,2)
        with (self.root/'ser').open('r+b') as f:f.write(struct.pack('<d',float('nan')))
        folder=self.root/'stage';folder.mkdir()
        with self.assertRaisesRegex(ValueError,'Nonfinite'):hn.split_raw(hn.profile(self.root),folder)
    def test_legacy_loop_and_changed_counter_or_loop_rejection(self):
        acquisition(self.root,program='N15HetNoe_rk.nl')
        self.assertEqual(hn.profile(self.root)['program'],'N15HetNoe_rk.nl')
        for old,new in [('"l2=1"','"l2=2"'),('lo to 3 times 2','lo to 3 times 4'),('ru2','ru0')]:
            (self.root/'pulseprogram').write_text(LEGACY.replace(old,new))
            with self.assertRaises(ValueError):hn.profile(self.root)
    def test_split_copy_unknown_program_and_wrong_modes_rejected(self):
        (self.root/'audita.txt').write_text('created by split 2: dataset 1 of 2 from EXPNO 5')
        with self.assertRaisesRegex(ValueError,'split copy'):hn.profile(self.root)
        (self.root/'audita.txt').unlink()
        for name,old,new in [('acqus','15N_HETNOE','unknown'),('acqu2s','FnMODE= 6','FnMODE= 5'),('acqus','FnTYPE= 0','FnTYPE= 2')]:
            path=self.root/name;original=path.read_text();path.write_text(original.replace(old,new))
            with self.assertRaises(ValueError):hn.profile(self.root)
            path.write_text(original)
    def test_counts_require_explicit_full_file_override_never_truncation(self):
        path=self.root/'acqu2s';path.write_text(path.read_text().replace('TD= 16','TD= 12'))
        with self.assertRaisesRegex(ValueError,'--raw-rows'):hn.profile(self.root)
        self.assertEqual(hn.profile(self.root,16)['rows'],16)
        for n in (12,14,20,0):
            with self.assertRaises(ValueError):hn.profile(self.root,n)
        with (self.root/'ser').open('ab') as f:f.write(b'x')
        with self.assertRaisesRegex(ValueError,'trailing bytes'):hn.profile(self.root,16)
    def test_inspect_does_not_invoke_commands_or_write(self):
        before={p.name:p.read_bytes() for p in self.root.iterdir()}
        with mock.patch.object(hn.rt,'checked_run') as run,mock.patch.object(hn,'pipeline') as pipe:
            self.assertEqual(hn.main(['inspect','--experiment',str(self.root)]),0)
        run.assert_not_called();pipe.assert_not_called()
        self.assertEqual({p.name:p.read_bytes() for p in self.root.iterdir()},before)
    def test_baseline_runs_in_N_then_H_with_signed_real_data_and_no_extra_axis_flip(self):
        commands=hn.processing_commands('input','output',self.args)
        active='H'; corrected=[]; phases=[]
        for cmd in commands:
            if 'TP' in cmd:active='N' if active=='H' else 'H'
            if 'PS' in cmd:
                self.assertIn('-di',cmd); phases.append(cmd[cmd.index('-p0')+1])
            if 'POLY' in cmd:corrected.append((active,cmd[cmd.index('-ord')+1],'-auto' in cmd))
            self.assertFalse(any(flag in cmd for flag in ('MC','-alt','-neg','MAC')))
        self.assertEqual(corrected,[('N',4,True),('H',3,True)])
        self.assertEqual(phases,[10.,30.])
    def fake_commands(self,command,cwd,stdout=None):
        cmd=list(map(str,command));name=cmd[0]
        if name=='bruk2pipe':
            self.conversions.append(cmd)
            Path(cmd[cmd.index('-out')+1]).write_bytes(pipe_bytes())
        elif name=='xyz2pipe':stdout.write(pipe_bytes(3,2))
        elif name=='pipe2ucsf':(Path(cwd)/'proc.ucsf').write_bytes(b'synthetic ucsf')
        elif name=='specStat.com':stdout.write('vEstNoise 12.5\n')
        elif name in ('sethdr','setfdata'):pass
        else:raise AssertionError(command)
    def fake_pipeline(self,commands,cwd):
        self.pipelines.append(commands)
        cmd=commands[-1];Path(cmd[cmd.index('-out')+1]).write_bytes(pipe_bytes())
    def test_successful_process_two_identical_recipes_noise_and_roles(self):
        self.conversions=[];self.pipelines=[]
        before={name:(self.root/name).read_bytes() for name in ('ser','acqus','acqu2s','pulseprogram')}
        with mock.patch.object(hn.rt,'checked_run',side_effect=self.fake_commands),\
             mock.patch.object(hn.rt,'set_series_metadata'),\
             mock.patch.object(hn,'pipeline',side_effect=self.fake_pipeline):hn.process(self.args)
        self.assertEqual(len(self.pipelines),2)
        self.assertEqual(self.pipelines[0][1:-1],self.pipelines[1][1:-1])
        for cmd in self.conversions:
            self.assertEqual(cmd[cmd.index('-yMODE')+1],'Echo-AntiEcho')
            self.assertEqual(cmd[cmd.index('-yN')+1],'8')
            self.assertEqual(cmd[cmd.index('-yT')+1],'4')
            self.assertNotIn('-noi2f',cmd)
        self.assertEqual((self.out/'relax_header.txt').read_text(),'2 0 1\n')
        self.assertEqual(hn.rt.counted_vector(self.out/'vEstNoise.txt'),[12.5,12.5])
        text=(self.out/'processing.txt').read_text()
        self.assertIn('15N_baseline: POLY -auto -ord 4',text)
        self.assertIn('H_baseline: POLY -auto -ord 3',text)
        self.assertIn('plane_roles: reference, saturated',text)
        self.assertFalse((self.out/'raw').exists())
        for name,data in before.items():self.assertEqual((self.root/name).read_bytes(),data)
    def test_float_conversion_uses_native_noi2f(self):
        acquisition(self.root,2,1);self.conversions=[];self.pipelines=[]
        with mock.patch.object(hn.rt,'checked_run',side_effect=self.fake_commands),\
             mock.patch.object(hn.rt,'set_series_metadata'),\
             mock.patch.object(hn,'pipeline',side_effect=self.fake_pipeline):hn.process(self.args)
        for cmd in self.conversions:self.assertIn('-noi2f',cmd);self.assertIn('-noswap',cmd)
    def test_failure_after_reference_preserves_previous_outputs(self):
        self.out.mkdir()
        for name in hn.OUTPUTS:
            if name=='planes':(self.out/name).mkdir();(self.out/name/'old').write_bytes(b'old')
            else:(self.out/name).write_bytes(b'old')
        self.conversions=[];self.pipelines=[]
        def fail_second(commands,cwd):
            if self.pipelines:raise subprocess.CalledProcessError(7,'N baseline')
            self.fake_pipeline(commands,cwd)
        with mock.patch.object(hn.rt,'checked_run',side_effect=self.fake_commands),\
             mock.patch.object(hn,'pipeline',side_effect=fail_second),self.assertRaises(subprocess.CalledProcessError):hn.process(self.args)
        for name in hn.OUTPUTS:
            path=self.out/name/('old' if name=='planes' else '')
            self.assertEqual(path.read_bytes(),b'old')
    def test_mixed_input_during_run_and_invalid_noise_do_not_publish(self):
        self.conversions=[];self.pipelines=[]
        def changed(commands,cwd):
            self.fake_pipeline(commands,cwd)
            with (self.root/'acqus').open('a') as f:f.write('$$ changed\n')
        with mock.patch.object(hn.rt,'checked_run',side_effect=self.fake_commands),\
             mock.patch.object(hn.rt,'set_series_metadata'),\
             mock.patch.object(hn,'pipeline',side_effect=changed),self.assertRaisesRegex(ValueError,'inputs changed'):hn.process(self.args)
        self.assertFalse((self.out/'m.ft3').exists())
    def test_intermediate_pipeline_failure_is_checked(self):
        with self.assertRaises(subprocess.CalledProcessError):
            hn.pipeline([[sys.executable,'-c','raise SystemExit(9)'],[sys.executable,'-c','import sys; sys.stdin.read()']],self.root)
    @unittest.skipUnless(os.environ.get('FTNMR_RUN_NMRPIPE_TESTS') == '1', 'optional installed-NMRPipe synthetic smoke test')
    def test_installed_tools_signed_ratio_and_float64_match(self):
        # Synthetic only: echo/antiecho encoding and exactly negative NOE.
        path=self.root/'acqus';path.write_text(path.read_text().replace('TD= 18','TD= 512').replace('2760','1').replace('67.98','0'))
        path=self.root/'acqu2s';path.write_text(path.read_text().replace('TD= 16','TD= 128'))
        rng=random.Random(31)
        with (self.root/'ser').open('wb') as handle:
            for t in range(32):
                for ea in range(2):
                    values=[]
                    for x in range(256):
                        phase=2*math.pi*(50*x/256+5*t/32*(1 if ea==0 else -1))
                        amp=1e6*math.exp(-x/70-t/15)*(1 if ea==0 else -1)
                        values.extend([4*round((amp*math.cos(phase)+rng.gauss(0,120))/4),
                                       4*round((amp*math.sin(phase)+rng.gauss(0,120))/4)])
                    handle.write(struct.pack('<512i',*values))
                    handle.write(struct.pack('<512i',*[-v//4 for v in values]))
        self.args.h_zf=1024;self.args.n_zf=256
        self.args.h_p0=self.args.h_p1=self.args.n_p0=self.args.n_p1=0.
        self.args.h_baseline_order=4
        hn.process(self.args)
        payloads=[(self.out/f'planes/plane_{i:02d}.ft2').read_bytes()[2048:] for i in (1,2)]
        vectors=[]
        for payload in payloads:
            a=array.array('f');a.frombytes(payload);vectors.append(a)
        a,b=vectors;index=max(range(len(a)),key=lambda i:abs(a[i]))
        self.assertAlmostEqual(b[index]/a[index],-.25,places=5)
        self.assertLess(math.sqrt(sum((v+.25*u)**2 for u,v in zip(a,b))/sum((.25*u)**2 for u in a)),1e-5)
        # Same inputs stored in big-endian float64 must decode identically.
        raw=(self.root/'ser').read_bytes();integers=struct.unpack('<'+'i'*(len(raw)//4),raw)
        (self.root/'ser').write_bytes(struct.pack('>'+'d'*len(integers),*integers))
        path=self.root/'acqus';path.write_text(path.read_text().replace('DTYPA= 0','DTYPA= 2').replace('BYTORDA= 0','BYTORDA= 1'))
        hn.process(self.args)
        for i,expected in enumerate(payloads,1):self.assertEqual((self.out/f'planes/plane_{i:02d}.ft2').read_bytes()[2048:],expected)

    def test_entrypoint_absolute_bare_path_and_symlink_help(self):
        script=RELAX/'hetnoe_proc.com'
        tcsh='/bin/tcsh'
        for path in [str(script),'hetnoe_proc.com']:
            command=[tcsh,'-f','-c',f'set path = ( {RELAX} $path ); {path} --help']
            result=subprocess.run(command,cwd=self.root,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            self.assertIn('--n-baseline-order',result.stdout)
        link=self.root/'linked.com';link.symlink_to(script)
        result=subprocess.run([tcsh,'-f',str(link),'--help'],cwd=self.root,capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        self.assertIn('--n-baseline-order',result.stdout)


if __name__=='__main__':unittest.main()
