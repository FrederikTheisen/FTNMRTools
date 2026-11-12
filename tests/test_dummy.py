import re
import unittest
from nmr_register.schema import columns, ATTENTION_FORMULA
from nmr_register.sync import _status_facts_for_remote
from nmr_register.registration_widget import experiment_actions, registration_html


def evaluate(formula, values):
    source=re.sub(r'\$(\w+)',lambda m:'values.get('+repr(m[1])+')',formula)
    namespace={}
    exec('def run(values):\n'+''.join('    '+line+'\n' for line in source.splitlines()),namespace)
    return namespace['run'](values)


class DummyTests(unittest.TestCase):
    def test_empty_dummy_suppresses_all_attention(self):
        self.assertEqual(evaluate(ATTENTION_FORMULA,{'AcquisitionStatus':'Dummy'}),'')
        cols={c['id']:c for c in columns()['Experiments']}
        self.assertEqual(eval(re.sub(r'\$(\w+)',lambda m:'values.get('+repr(m[1])+')',cols['EffectiveAcquisitionStatus']['formula']),{'values':{'AcquisitionStatus':'Dummy'}}),'Dummy')

    def test_sync_repeatedly_retains_dummy_and_metadata(self):
        remote={'AcquisitionStatus':'Dummy'}
        for observed in ('Unknown','Completed','Aborted','Setup'):
            result=_status_facts_for_remote({'AcquisitionStatus':observed,'Instrument':'new'}, {},remote)
            self.assertEqual(result,{'AcquisitionStatus':'Dummy','Instrument':'new'})
            remote.update(result)
        facts={'AcquisitionStatus':'Aborted'}
        self.assertEqual(_status_facts_for_remote(facts,{},{}),facts)

    def test_link_is_distinct_from_raw_or_typed_folder(self):
        formula=next(c['formula'] for c in columns()['Experiments'] if c['id']=='DesktopLink')
        def run(values):return eval(re.sub(r'\$(\w+)',lambda m:'values.get('+repr(m[1])+')',formula),{'values':values})
        self.assertEqual(run({'Folder':'typed','RawPresent':True}),'Not linked')
        self.assertEqual(run({'AcquisitionKey':'key','SourcePath':'folder','RawPresent':False}),'')
        self.assertEqual(run({'AcquisitionStatus':'Dummy'}),'')

    def test_dialog_defaults_and_limited_initial_statuses(self):
        tables={'Experiments':[],'Samples':[],'Series':[]}
        for value,wanted in [({},'Planned'),({'AcquisitionStatus':'Dummy'},'Dummy')]:
            self.assertEqual(experiment_actions(value,tables)['actions'][0][3]['AcquisitionStatus'],wanted)
        with self.assertRaises(ValueError):experiment_actions({'AcquisitionStatus':'Completed'},tables)
        self.assertIn("['Dummy','Dummy entry']",registration_html('experiment'))
