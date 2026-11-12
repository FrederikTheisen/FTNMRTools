"""Instrument outcome from exact saved TopSpin audit records, never raw size."""
from datetime import datetime, timezone, timedelta
import hashlib
import json
from pathlib import Path
import re
import tempfile
import unittest
from unittest.mock import patch
from nmr_register.bruker import acquisition_status, observe, MAX_METADATA_BYTES
from nmr_register.state import State

START = 1748448031
END = 1748545685
UUID = 'e1c20164-3bdc-11f0-a2ab-0030d625a2a2'


def timestamp(value, iso=False, offset=0):
    dt=datetime.fromtimestamp(value,timezone(timedelta(hours=offset)))
    return dt.strftime('%Y-%m-%d'+('T' if iso else ' ')+'%H:%M:%S.%f'+('' if iso else ' ')+'%z')


def record(number, when, body, version='4.4.1', process='go4', offset=0):
    return (f'( {number:3d},<{timestamp(when,iso=version=="4.4.1",offset=offset)}>,<user>,<instrument>,<{process}>,<TopSpin {version}>,\n'
            f'      <{body}>)\n')


def terminal(start=START+.962, end=END+1.238, status='Completed', uuid=UUID, version='4.4.1', offset=0):
    action=('completed at '+timestamp(end,offset=offset) if status=='Completed' else
            "terminated by command 'stop' at "+timestamp(end,offset=offset))
    intro=("created by zg on 'Avance' (instrument) by user 'user'\n       on data set with UUID '" if version=='4.4.1'
           else "created by zg on data set with UUID '")
    return (intro+uuid+"'\n       started at "+timestamp(start,offset=offset)+', POWCHK enabled, PULCHK disabled,\n'
            +'       '+action+'\n       configuration hash MD5:\n       AA BB CC\n'
            +'       data hash MD5: 2K * 256 * 7\n       DD EE FF')


def audit(*records, version='4.4.1'):
    return ('##TITLE= Audit trail, TopSpin '+version+'\n##JCAMPDX= 5.01\n'
            +'##AUDIT TRAIL= $$ (NUMBER, WHEN, WHO, WHERE, PROCESS, VERSION, WHAT)\n'
            +''.join(records)+'##END=\n\n$$ hash MD5\n$$ AA BB CC\n')


def acqus(start=START,end=END):
    return f'##$DATE_START= {start}\n##$DATE= {end}\n##$PULPROG= <15N_T1_HSQC>\n##$BF1= 950.23\n##$NUC1= <1H>\n'


class AcquisitionCompletionTests(unittest.TestCase):
    def status(self,text,metadata=None):
        result=acquisition_status(acqus() if metadata is None else metadata,text)
        self.assertIn(result['ObservedAcquisitionStatus'],('Completed','Aborted','Unknown'))
        self.assertLess(len(result['ObservedAcquisitionEvidence']),250)
        return result

    def test_saved_topspin3_and4_formats_finish_with_explicit_offsets(self):
        for version,offset in [('3.7.0',2),('4.4.1',0)]:
            with self.subTest(version=version):
                text=audit(record(1,END+1.239,terminal(version=version,offset=offset),version=version,offset=offset),version=version)
                result=self.status(text)
                self.assertEqual(result['ObservedAcquisitionStatus'],'Completed')
                self.assertIn('DATE_START matched',result['ObservedAcquisitionEvidence'])

    def test_reviewed_additional_versions_keep_completion_and_terminal_guards(self):
        for version in ('3.6.5','5.0.0'):
            with self.subTest(version=version):
                text=audit(record(1,END+1.239,terminal(version=version),version=version),version=version)
                self.assertEqual(self.status(text)['ObservedAcquisitionStatus'],'Completed')
                self.assertEqual(self.status(text,acqus(start=START+1))['ObservedAcquisitionStatus'],'Unknown')
                self.assertEqual(self.status(text,acqus(end=END-10))['ObservedAcquisitionStatus'],'Unknown')
                self.assertEqual(self.status(audit(record(1,END,'acquisition in progress',version=version),version=version))['ObservedAcquisitionStatus'],'Unknown')
        for version in ('3.5','3.2','4.4.2','5.0.1'):
            text=audit(record(1,END+1.239,terminal(),version=version),version=version)
            self.assertEqual(self.status(text)['ObservedAcquisitionStatus'],'Unknown')

    def test_jcamp_source_comments_between_records_are_inert(self):
        comments='$$ /redacted/source/experiment/audita.txt\n  $$ ordinary source comment\n'
        text=audit(comments,record(1,END+1.239,terminal()),comments,
                   record(2,END+2,'changed parameter F1 O1 to 3290.60',process='cprserver'),comments)
        self.assertEqual(self.status(text)['ObservedAcquisitionStatus'],'Completed')
        bad=audit(comments+'unexpected noncomment text\n',record(1,END+1.239,terminal()))
        self.assertEqual(self.status(bad)['ObservedAcquisitionStatus'],'Unknown')

    def test_fake_record_in_comment_never_proves_completion(self):
        fake=''.join('$$ '+line+'\n' for line in record(1,END+1.239,terminal()).splitlines())
        self.assertEqual(self.status(audit(fake))['ObservedAcquisitionStatus'],'Unknown')
        copied=audit(record(1,END,'wser 10 2',process='mvrow'),fake)
        self.assertEqual(self.status(copied)['ObservedAcquisitionStatus'],'Unknown')
        real=audit(fake,record(2,END,'acquisition in progress'))
        self.assertEqual(self.status(real)['ObservedAcquisitionStatus'],'Unknown')

    def test_stop_is_distinct_from_completion_and_requires_matching_metadata(self):
        text=audit(record(1,END+1.239,terminal(status='Aborted')))
        self.assertEqual(self.status(text)['ObservedAcquisitionStatus'],'Aborted')
        self.assertEqual(self.status(text,acqus(start=START+1))['ObservedAcquisitionStatus'],'Unknown')
        self.assertEqual(self.status(text,acqus(end=END-10))['ObservedAcquisitionStatus'],'Unknown')

    def test_old_stop_deleted_restart_then_latest_completion(self):
        text=audit(
            record(53,START-1800,'acquisition in progress'),
            record(55,START-1700,terminal(start=START-1800,end=START-1700,status='Aborted',uuid='da412dba-3bd8-11f0-ad47-0030d625a2a2')),
            record(56,START-100,'changed parameter F3 NS to 16',process='cprserver',offset=2),
            record(57,START+.959,'raw data deleted by zg'),
            record(58,START+.962,'acquisition in progress'),
            record(59,END+1.239,terminal()))
        self.assertEqual(self.status(text)['ObservedAcquisitionStatus'],'Completed')
        self.assertIn('entry 59',self.status(text)['ObservedAcquisitionEvidence'])

    def test_later_restart_or_deletion_cannot_reuse_previous_finished(self):
        completed=record(1,END+1.239,terminal())
        for message in ['acquisition in progress','raw data deleted by zg','acquisition aborted by user']:
            with self.subTest(message=message):
                result=self.status(audit(completed,record(2,END+60,message)))
                self.assertEqual(result['ObservedAcquisitionStatus'],'Unknown')
        # A different latest completed generation must not fall back to the old
        # matching completed record either.
        newer=record(2,END+200,terminal(start=END+100,end=END+200,uuid='da412dba-3bd8-11f0-ad47-0030d625a2a2'))
        self.assertEqual(self.status(audit(completed,newer))['ObservedAcquisitionStatus'],'Unknown')

    def test_in_progress_is_never_live_running_even_with_present_raw(self):
        text=audit(record(1,START+.962,'acquisition in progress'))
        result=self.status(text)
        self.assertEqual(result['ObservedAcquisitionStatus'],'Unknown')
        self.assertIn('cannot establish live running',result['ObservedAcquisitionEvidence'])

    def test_processing_completion_substrings_and_nonmonotonic_processing_times(self):
        misleading=record(1,END+1.239,'processing completed at '+timestamp(END+1.238),process='proc')
        self.assertEqual(self.status(audit(misleading))['ObservedAcquisitionStatus'],'Unknown')
        # Processing records can have timezone-valid earlier clocks. They do not
        # supersede the latest instrument acquisition outcome.
        processing=record(2,START-60,'processing completed; previous acquisition in progress',process='proc',offset=2)
        self.assertEqual(self.status(audit(record(1,END+1.239,terminal()),processing))['ObservedAcquisitionStatus'],'Completed')
        quoted=record(1,END+1.239,'user note: '+terminal())
        self.assertEqual(self.status(audit(quoted))['ObservedAcquisitionStatus'],'Unknown')

    def test_missing_stale_invalid_and_conflicting_metadata_or_log(self):
        good=audit(record(1,END+1.239,terminal()))
        for metadata in ['',acqus().replace('DATE_START','OLD_DATE_START'),acqus().replace(str(START),'NaN'),
                         acqus(start=START+1),acqus(end=END+10),acqus()+'##$DATE_START= '+str(START)+'\n']:
            with self.subTest(metadata=metadata):
                self.assertEqual(self.status(good,metadata)['ObservedAcquisitionStatus'],'Unknown')
        for text in ['',good.replace('##END=','##UNFINISHED='),good[:-40],
                     good.replace('+0000',''),good.replace('TopSpin 4.4.1','TopSpin 5.0'),
                     good.replace('completed at','aborted at'),good.replace('completed at','not completed at'),
                     good+'( 2,<broken later acquisition',
                     good.replace('##END=','( 2,<broken later acquisition\n##END='),
                     audit(record(1,END+1.239,terminal()),record(1,END+1.239,terminal()))]:
            with self.subTest(text=text[:100]):
                self.assertEqual(self.status(text)['ObservedAcquisitionStatus'],'Unknown')
        conflicting=terminal()+'\n       terminated by command \'stop\' at '+timestamp(END+1.238)
        self.assertEqual(self.status(audit(record(1,END+1.239,conflicting)))['ObservedAcquisitionStatus'],'Unknown')

    def test_unclosed_terminal_cannot_swallow_a_later_restart_record(self):
        broken=record(1,END+1.239,terminal()).replace('>)\n','\n',1)
        text=audit(broken,record(2,END+60,'acquisition in progress'))
        result=self.status(text)
        self.assertEqual(result['ObservedAcquisitionStatus'],'Unknown')
        self.assertIn('nested',result['ObservedAcquisitionEvidence'])

    def test_invalid_timezone_offsets_are_not_normalized_into_valid_offsets(self):
        good=audit(record(1,END+1.239,terminal()))
        for offset in ['+0060','+1500','+1401','-0060']:
            with self.subTest(offset=offset):
                self.assertEqual(self.status(good.replace('+0000',offset))['ObservedAcquisitionStatus'],'Unknown')

    def test_same_second_ambiguous_generations_and_split_are_unknown(self):
        first=record(1,END+1.239,terminal(start=START+.1))
        second=record(2,END+1.239,terminal(start=START+.9,uuid='da412dba-3bd8-11f0-ad47-0030d625a2a2'))
        self.assertEqual(self.status(audit(first,second))['ObservedAcquisitionStatus'],'Unknown')
        split=record(2,END+60,'created by split 4: dataset 1 of 4 from EXPNO 3',process='cprserver')
        self.assertEqual(self.status(audit(first,split))['ObservedAcquisitionStatus'],'Unknown')

    def test_terminal_date_tolerance_is_bounded_to_two_seconds(self):
        for delta,expected in [(0,'Completed'),(1.238,'Completed'),(2,'Completed'),(2.001,'Unknown')]:
            with self.subTest(delta=delta):
                text=audit(record(1,END+delta,terminal(end=END+delta)))
                self.assertEqual(self.status(text)['ObservedAcquisitionStatus'],expected)
        text=audit(record(1,END+10,terminal()))
        self.assertEqual(self.status(text)['ObservedAcquisitionStatus'],'Unknown')

    def test_observe_keeps_raw_transfer_assessment_and_fingerprint_separate(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'1';path.mkdir()
            metadata=acqus();(path/'acqus').write_text(metadata)
            (path/'audita.txt').write_text(audit(record(1,END+1.239,terminal())))
            legacy={m[1]:m[2].strip() for m in re.finditer(r'^##\$([A-Za-z0-9_]+)=\s*([^\n]*)',metadata,re.M)}
            identity={k:legacy.get(k,'') for k in ('DATE_START','DATE','PULPROG','BF1','SFO1','NUC1','TD')}
            expected=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()
            original_open=Path.open
            def metadata_only(file,*args,**kwargs):
                if file.name in ('fid','ser'):raise AssertionError('raw bytes must not be read')
                return original_open(file,*args,**kwargs)
            with patch.object(Path,'open',metadata_only):
                fields=observe(path)
            self.assertEqual(fields['ObservedAcquisitionStatus'],'Completed')
            self.assertFalse(fields['RawPresent'])
            self.assertEqual(fields['ObservedFingerprint'],expected)
            self.assertFalse({'Assessment','AssessmentReason','AcquisitionStatus'} & fields.keys())
            (path/'ser').write_bytes(b'partial transfer')
            state=State(Path(folder)/'state')
            try:
                fields['AcquisitionKey']='source:sample/1'
                state.observe([fields]);state.observe([fields]);state.acknowledge([key for key,_ in state.pending()])
                (path/'audita.txt').write_text(audit(record(1,END+1.239,terminal()),record(2,END+60,'acquisition in progress')))
                changed=observe(path);changed['AcquisitionKey']='source:sample/1'
                with patch.object(Path,'open',metadata_only):
                    self.assertTrue(observe(path)['RawPresent'])
                state.observe([changed]);payload=state.pending()[0][1]
                self.assertEqual(payload['ObservedAcquisitionStatus'],'Unknown')
                self.assertEqual(payload['ObservedFingerprint'],expected)
                self.assertEqual(payload['ReconciliationReason'],'')
                self.assertFalse(payload['TransferStable'])
            finally:state.close()

    def test_timestamped_halt_is_aborted_but_untimed_stop_is_unproven(self):
        body=terminal(status='Aborted').replace("command 'stop'", "command 'halt'")
        result=self.status(audit(record(1,END+1.239,body)))
        self.assertEqual(result['ObservedAcquisitionStatus'],'Aborted')
        self.assertIn("command 'halt'",result['ObservedAcquisitionEvidence'])
        untimed=re.sub(r"terminated by command 'halt' at [^\n]+", "terminated by command 'stop'",body)
        self.assertEqual(self.status(audit(record(1,END+1.239,untimed)))['ObservedAcquisitionStatus'],'Unknown')

    def test_setup_is_clean_local_copy_without_raw_not_live_running(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'1';path.mkdir()
            (path/'acqus').write_text(acqus())
            initial=observe(path)
            self.assertEqual(initial['ObservedAcquisitionStatus'],'Setup')
            self.assertIn('local copy',initial['ObservedAcquisitionEvidence'])
            self.assertFalse(initial['RawPresent'])
            (path/'fid').touch()
            self.assertEqual(observe(path)['ObservedAcquisitionStatus'],'Setup')
            (path/'fid').write_bytes(b'partial')
            self.assertEqual(observe(path)['ObservedAcquisitionStatus'],'Unknown')
            (path/'audita.txt').write_text(audit(record(1,START,'acquisition in progress')))
            result=observe(path)
            self.assertEqual(result['ObservedAcquisitionStatus'],'Unknown')
            self.assertTrue(result['ObservedAcquisitionEvidence'].startswith('audita.txt: latest attempt unresolved: '))
            self.assertEqual(result['ObservedFingerprint'],initial['ObservedFingerprint'])
            (path/'fid').unlink()
            self.assertEqual(observe(path)['ObservedAcquisitionStatus'],'Setup')
            (path/'acqus').unlink()
            self.assertEqual(observe(path)['ObservedAcquisitionStatus'],'Unknown')

    def test_unsafe_raw_or_invalid_setup_parameters_cannot_establish_setup(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'1';path.mkdir()
            (path/'acqus').write_text(acqus())
            (path/'fid').symlink_to(Path(folder)/'elsewhere')
            self.assertEqual(observe(path)['ObservedAcquisitionStatus'],'Unknown')
            (path/'fid').unlink()
            (path/'acqus').write_text('##$PULPROG= <zg>\n##$D= (0..2)\n1 2\n')
            self.assertEqual(observe(path)['ObservedAcquisitionStatus'],'Unknown')

    def test_bounded_reader_rejects_oversized_and_symlink_audit(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);path=root/'1';path.mkdir()
            (path/'acqus').write_text(acqus())
            external=root/'other-audit';external.write_text(audit(record(1,END+1.239,terminal())))
            (path/'audita.txt').symlink_to(external)
            self.assertEqual(observe(path)['ObservedAcquisitionStatus'],'Unknown')
            (path/'audita.txt').unlink()
            (path/'audita.txt').write_bytes(b' '*(MAX_METADATA_BYTES+1))
            self.assertEqual(observe(path)['ObservedAcquisitionStatus'],'Unknown')


if __name__=='__main__':unittest.main()
