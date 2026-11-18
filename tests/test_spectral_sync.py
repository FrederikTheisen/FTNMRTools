"""Spectral metadata follows acquisition identity guards and preserves human reviews."""
import tempfile
import unittest
from pathlib import Path
from nmr_register.state import State
from nmr_register.sync import _canonical_payload,flush
from test_sync import Fake,fields

class SpectralSyncTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.state=State(Path(self.temp.name));self.addCleanup(self.state.close)
        self.incoming={**fields(),'ObservedSpectralWidths_ppm':'Direct (1H): 12.03','ObservedSpectralOffsets_ppm':'Direct (1H): 0'}
    def test_canonical_mapping_and_unknown_retention(self):
        machine,facts=_canonical_payload(self.incoming)
        self.assertEqual(facts['SpectralOffsets_ppm'],'Direct (1H): 0')
        self.assertEqual(facts['SpectralWidths_ppm'],'Direct (1H): 12.03')
        self.assertNotIn('ObservedSpectralOffsets_ppm',machine)
        _,unknown=_canonical_payload({**self.incoming,'ObservedSpectralOffsets_ppm':'','ObservedSpectralWidths_ppm':''})
        self.assertNotIn('SpectralWidths_ppm',unknown);self.assertNotIn('SpectralOffsets_ppm',unknown)
    def test_sync_preserves_humans_and_rescan_is_unchanged(self):
        self.state.observe([self.incoming]);client=Fake();flush(self.state,client)
        saved=client.rows[0]['fields'];saved.update(Sample=8,Assessment='Use',QCSNR='Adequate',Notes='keep')
        self.assertEqual(saved['SpectralWidths_ppm'],'Direct (1H): 12.03')
        self.assertEqual(saved['SpectralOffsets_ppm'],'Direct (1H): 0')
        count=len(client.patches);flush(self.state,client)
        self.assertEqual(len(client.patches),count)
        changed={**self.incoming,'ObservedSpectralOffsets_ppm':'Direct (1H): -1'}
        self.state.observe([changed]);flush(self.state,client)
        self.assertEqual(saved['SpectralOffsets_ppm'],'Direct (1H): -1')
        self.assertEqual((saved['Sample'],saved['Assessment'],saved['QCSNR'],saved['Notes']),(8,'Use','Adequate','keep'))
        self.assertEqual(saved['ReconciliationReason'],'')
    def test_ambiguous_or_mismatched_identity_does_not_replace_spectral_data(self):
        remote={**self.incoming,'SourceID':'local','Folder':'different','Expno':'1','SpectralWidths_ppm':'preserve','SpectralOffsets_ppm':'preserve','Assessment':'Use'}
        self.state.observe([self.incoming]);client=Fake([{'id':1,'fields':remote}]);result=flush(self.state,client)
        self.assertEqual(result['conflicts'],1)
        self.assertEqual(client.rows[0]['fields']['SpectralWidths_ppm'],'preserve')
        self.assertEqual(client.rows[0]['fields']['SpectralOffsets_ppm'],'preserve')
        self.assertEqual(client.rows[0]['fields']['Assessment'],'Use')
