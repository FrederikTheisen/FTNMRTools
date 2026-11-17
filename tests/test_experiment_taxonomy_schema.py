import re
import unittest

from nmr_register.schema import columns, PRE_CATEGORY_ATTENTION_FORMULA
from nmr_register.setup import upgrade_experiment_taxonomy, show_experiment_category
from nmr_register.experiment_taxonomy import category


class Fake:
    def __init__(self):
        self.cols = [{'id': i + 1, 'fields': {'parentId': 1, 'colId': c['id'], **{k: v for k, v in c.items() if k != 'id'}}}
                     for i, c in enumerate(columns()['Experiments']) if c['id'] != 'ExperimentCategory']
        for row in self.cols:
            if row['fields']['colId'] == 'AttentionReason': row['fields']['formula'] = PRE_CATEGORY_ATTENTION_FORMULA
            if row['fields']['colId'] == 'ExperimentType': row['fields']['label'] = 'ExperimentType'
            if row['fields']['colId'] == 'EffectiveType': row['fields']['label'] = 'Type'
        self.fields = []
    def records(self, name):
        return {'_grist_Tables': [{'id': 1, 'fields': {'tableId': 'Experiments'}}],
                '_grist_Tables_column': self.cols,
                '_grist_Views': [{'id': 10, 'fields': {'name': 'Experiments'}}],
                '_grist_Views_section': [{'id': 20, 'fields': {'parentId': 10, 'tableRef': 1, 'title': 'Experiments'}}],
                '_grist_Views_section_field': self.fields}[name]
    def apply(self, actions):
        for action in actions:
            if action[0] == 'AddColumn':
                self.cols.append({'id': 999, 'fields': {'parentId': 1, 'colId': action[2], **action[3]}})
            elif action[1] == '_grist_Tables_column':
                next(r for r in self.cols if r['id'] == action[2])['fields'].update(action[3])
            else:
                self.fields.append({'id': 1000 + len(self.fields), 'fields': action[3]})


class TaxonomySchemaTests(unittest.TestCase):
    def test_category_formula_agrees_with_taxonomy(self):
        col = next(c for c in columns()['Experiments'] if c['id'] == 'ExperimentCategory')
        source = re.sub(r'\$(\w+)', lambda m: 'values.get(' + repr(m[1]) + ')', col['formula'])
        ns = {}
        exec('def run(values):\n' + ''.join('    ' + line + '\n' for line in source.splitlines()), ns)
        for value, pulse in [('T1', ''), ('HSQC', ''), ('13C HSQC', ''), ('HNCA', ''), ('zg30', ''), ('1D', ''), ('custom', 'zg'), ('', 'zg')]:
            self.assertEqual(ns['run']({'ExperimentType': value, 'PulseProgram': pulse}), category(value, pulse))

    def test_upgrade_is_guarded_and_idempotent(self):
        client = Fake()
        self.assertTrue(upgrade_experiment_taxonomy(client, apply=True))
        self.assertEqual(upgrade_experiment_taxonomy(client, apply=True), [])
        experiment_id = next(r['id'] for r in client.cols if r['fields']['colId'] == 'ExperimentType')
        client.fields.extend([{'id': 1, 'fields': {'parentId': 20, 'colRef': 123, 'parentPos': 2}},
                              {'id': 2, 'fields': {'parentId': 20, 'colRef': experiment_id, 'parentPos': 4}}])
        added = show_experiment_category(client, apply=True)
        self.assertEqual(added[0][3]['parentPos'], 3)
        self.assertEqual(added[0][3]['width'], 120)
        self.assertEqual(show_experiment_category(client, apply=True), [])
        next(r for r in client.cols if r['fields']['colId'] == 'AttentionReason')['fields']['formula'] = 'custom'
        with self.assertRaises(ValueError): upgrade_experiment_taxonomy(client, apply=True)
