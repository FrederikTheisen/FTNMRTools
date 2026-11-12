import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from nmr_register.cli import install_electrostatics
from nmr_register.grist import GristError

class InstallTests(unittest.TestCase):
    def setUp(self):
        self.events=[]
        self.client=SimpleNamespace(records=lambda table: self.events.append('read') or [])
        self.config=SimpleNamespace(home=Path('/tmp/example'),api_base='https://example/api/docs/doc')
        self.plan={'issues':[],'records':[{'id':1,'fields':{'Study':'DNA'}}],'catalog_sha256':'hash'}
    def test_dry_run_is_read_only(self):
        with patch('nmr_register.migration.reviewed_pilot_patch',return_value=self.plan):
            result=install_electrostatics(self.client,self.config,'catalog')
        self.assertTrue(result['dry_run']); self.assertEqual(self.events,['read'])
    def test_apply_backs_up_before_writes_and_recomputes_plan(self):
        def event(name,result=None):
            def run(*args,**kwargs):
                self.events.append(name)
                return result
            return run
        with patch('nmr_register.migration.reviewed_pilot_patch',side_effect=event('plan',self.plan)), \
             patch('nmr_register.migration.apply_reviewed_pilot_patch',side_effect=event('seed',{'records':[]})), \
             patch('nmr_register.cli.create_backup',side_effect=event('backup',{'path':'backup'})), \
             patch('nmr_register.setup.upgrade_acquisition_schema',side_effect=event('acquisition')), \
             patch('nmr_register.setup.upgrade_labelling_choices',side_effect=event('labels')), \
             patch('nmr_register.setup.upgrade_electrostatics_schema',side_effect=event('schema')), \
             patch('nmr_register.setup.build_layout',side_effect=event('page')), \
             patch('nmr_register.setup.apply_electrostatics_ui',side_effect=event('layout',{})), \
             patch('nmr_register.setup.apply_experiment_registration_widget',side_effect=event('form',1)), \
             patch('nmr_register.setup.apply_experiments_ui',side_effect=event('experiment_layout',{})), \
             patch('nmr_register.setup.apply_sample_registration_widget',side_effect=event('sample_form',2)), \
             patch('nmr_register.samples_ui.apply_samples_ui',side_effect=event('sample_layout',{})):
            result=install_electrostatics(self.client,self.config,'catalog',True)
        self.assertEqual(self.events,['read','plan','backup','acquisition','labels','schema','page','layout','form','experiment_layout','sample_form','sample_layout','read','plan','seed'])
        self.assertEqual(result['backup'],'backup')
    def test_quota_halts_before_seed_with_resume_message(self):
        with patch('nmr_register.migration.reviewed_pilot_patch',return_value=self.plan), \
             patch('nmr_register.migration.apply_reviewed_pilot_patch') as seed, \
             patch('nmr_register.cli.create_backup',return_value={'path':'backup'}), \
             patch('nmr_register.setup.upgrade_acquisition_schema',side_effect=GristError('Grist HTTP 429')):
            with self.assertRaisesRegex(GristError,'rerun electrostatics --apply'):
                install_electrostatics(self.client,self.config,'catalog',True)
            seed.assert_not_called()
    def test_pilot_changes_during_install_stop_the_seed(self):
        changed={**self.plan,'records':[]}
        with patch('nmr_register.migration.reviewed_pilot_patch',side_effect=[self.plan,changed]), \
             patch('nmr_register.migration.apply_reviewed_pilot_patch') as seed, \
             patch('nmr_register.cli.create_backup',return_value={'path':'backup'}), \
             patch('nmr_register.setup.upgrade_acquisition_schema'), \
             patch('nmr_register.setup.upgrade_labelling_choices'), \
             patch('nmr_register.setup.upgrade_electrostatics_schema'), \
             patch('nmr_register.setup.build_layout'), \
             patch('nmr_register.setup.apply_electrostatics_ui',return_value={}), \
             patch('nmr_register.setup.apply_experiment_registration_widget',return_value=1), \
             patch('nmr_register.setup.apply_experiments_ui',return_value={}), \
             patch('nmr_register.setup.apply_sample_registration_widget',return_value=2), \
             patch('nmr_register.samples_ui.apply_samples_ui',return_value={}):
            with self.assertRaisesRegex(ValueError,'changed during setup'):
                install_electrostatics(self.client,self.config,'catalog',True)
            seed.assert_not_called()
