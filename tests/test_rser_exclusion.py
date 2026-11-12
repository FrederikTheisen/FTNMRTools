"""Explicit structured rser history excludes extracted FIDs, never originals."""
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from nmr_register.bruker import observe
from nmr_register.discovery import scan_source
from nmr_register.split import audit_is_derived_split,is_derived_split
from nmr_register.state import State,serialized,MATCHING_METADATA_PREFIX

# Redacted exact record structure from the reviewed TopSpin 3.5 pl 7 audit.
RSER='''##TITLE= Audit trail, TopSpin 3.5 pl 7
##JCAMPDX= 5.01
##OWNER= redacted
$$ /redacted/nmr/folder/8501/audita.txt
##AUDIT TRAIL=  $$ (NUMBER, WHEN, WHO, WHERE, PROCESS, VERSION, WHAT)
(   1,<2024-12-28 16:19:01.012 +0100>,<user>,<instrument>,<mvrow>,<TopSpin 3.5.7>,
      <rser 1 from
       85 "folder" "/redacted/nmr"
       data hash MD5: 1280
       F9 B2 FF 93 ED 04 0A 4F 77 33 60 0F 2C 5E 12 F0>)
##END=

$$ hash MD5
$$ 59 38 A1 A1 E7 22 C7 F2 A3 B9 F4 92 4C 46 34 65
'''

class RserExclusionTests(unittest.TestCase):
    def test_only_complete_structured_rser_provenance_counts(self):
        self.assertTrue(audit_is_derived_split(RSER))
        for text in ('rser 1 from 85 "folder" "/redacted/nmr"',
                     ''.join('$$ '+line+'\n' for line in RSER.splitlines()),
                     RSER.replace('<mvrow>','<cprserver>'),
                     RSER.replace('##END=',''),
                     RSER.replace('       85 "folder"','       unknown "folder"'),
                     RSER.replace('##END=','malformed record\n##END=')):
            self.assertFalse(audit_is_derived_split(text))
        for name in ('Title','Notes','AcquisitionStatus','SourcePath'):
            self.assertFalse(is_derived_split({name:RSER}))

    def test_acquisition_process_record_prevents_exclusion(self):
        for process in ('go1','go4','go64'):
            actual='( 2,<2024-12-28 16:20:01.012 +0100>,<user>,<instrument>,<'+process+'>,<TopSpin 3.5.7>,<acquisition in progress>)\n'
            self.assertFalse(audit_is_derived_split(RSER.replace('##END=',actual+'##END=')))
        self.assertFalse(audit_is_derived_split(RSER.replace('<mvrow>','<go4>')))

    def test_discovery_and_cache_exclusion_keep_original_files_and_hard_warnings(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);folder=root/'folder'; original=folder/'85';derived=folder/'8501'
            for path in (original,derived):
                path.mkdir(parents=True)
                (path/'acqus').write_text('##$DATE_START= 1735380000\n##$DATE= 1735380100\n##$PULPROG= <t1>\n')
            (derived/'audita.txt').write_text(RSER)
            source=SimpleNamespace(id='src',path=root,targets=None,exclude=(),instrument='')
            rows,errors=scan_source(source)
            self.assertFalse(errors);self.assertEqual([r['ObservedExpno'] for r in rows],['85'])
            copied=observe(derived)
            copied.update(AcquisitionKey='src:folder/8501',ObservedSourceID='src',ObservedFolder='folder',ObservedExpno='8501')
            self.assertTrue(is_derived_split(copied))
            state=State(root/'cache');self.addCleanup(state.close)
            survivor=dict(rows[0]);survivor['ReconciliationReason']=MATCHING_METADATA_PREFIX+copied['AcquisitionKey']
            old_cached_copy=dict(copied,ObservedAcquisitionEvidence='audita.txt: unsupported or missing audit format/version',ScanError='')
            with state.db:
                for record in (survivor,old_cached_copy):
                    state.db.execute('INSERT INTO observations(key,fingerprint,signature,payload,first_seen,last_seen) VALUES (?,?,?,?,?,?)',
                                     (record['AcquisitionKey'],'fingerprint','signature',serialized(record),1,1))
                    state.db.execute('INSERT INTO outbox VALUES (?,?,?)',(record['AcquisitionKey'],serialized(record),1))
                state.db.execute('DELETE FROM outbox WHERE key=?',(survivor['AcquisitionKey'],))
            state.exclude_derived_splits([copied],now=2)
            self.assertEqual([key for key,_ in state.cached_observations()],[survivor['AcquisitionKey']])
            self.assertEqual([key for key,_ in state.pending()],[survivor['AcquisitionKey']])
            self.assertEqual(state.pending()[0][1]['ReconciliationReason'],'')
            hard='raw acquisition file replaced; '+MATCHING_METADATA_PREFIX+copied['AcquisitionKey']
            survivor['ReconciliationReason']=hard
            with state.db:state.db.execute('UPDATE observations SET payload=? WHERE key=?',(serialized(survivor),survivor['AcquisitionKey']))
            state.exclude_derived_splits([copied],now=3)
            self.assertEqual(state.cached_observations()[0][1]['ReconciliationReason'],hard)
            self.assertEqual((derived/'audita.txt').read_text(),RSER)
            self.assertTrue((original/'acqus').exists());self.assertTrue((derived/'acqus').exists())
