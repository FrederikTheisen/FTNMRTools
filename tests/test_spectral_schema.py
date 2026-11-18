import copy
import json
import unittest

from nmr_register.schema import CANONICAL_ACQUISITION_FIELDS, columns
from nmr_register.setup import upgrade_spectral_schema


class Client:
    def __init__(self):
        self.data = {
            '_grist_Tables': [{'id': 1, 'fields': {'tableId': 'Experiments'}}],
            '_grist_Tables_column': [{'id': 10, 'fields': {'parentId': 1, 'colId': 'D1_s', 'type': 'Numeric'}}],
            '_grist_Views': [{'id': i, 'fields': {'name': n}} for i, n in [(2, 'Experiments'), (3, 'Needs attention'), (4, 'Electrostatics')]],
            '_grist_Views_section': [
                {'id': 23, 'fields': {'parentId': 2, 'tableRef': 1, 'parentKey': 'single', 'title': 'Experiment details', 'layoutSpec': json.dumps({'children': [{'leaf': 100}], 'custom': 'keep'})}},
                {'id': 30, 'fields': {'parentId': 3, 'tableRef': 1, 'parentKey': 'single', 'title': 'Attention details', 'layoutSpec': ''}},
                {'id': 41, 'fields': {'parentId': 4, 'tableRef': 1, 'parentKey': 'single', 'title': 'Acquisition details', 'layoutSpec': ''}},
                {'id': 22, 'fields': {'parentId': 2, 'tableRef': 1, 'parentKey': 'record', 'title': 'Experiments'}},
            ],
            '_grist_Views_section_field': [{'id': 100, 'fields': {'parentId': 23, 'colRef': 10, 'parentPos': 1}}],
            'Experiments': [{'id': 1, 'fields': {'D1_s': 2}}],
        }
        self.actions = []

    def records(self, table): return copy.deepcopy(self.data[table])

    def apply(self, actions):
        self.actions.extend(copy.deepcopy(actions))
        for action in actions:
            kind, table, ident, values = action
            if kind == 'AddColumn':
                ident = max(r['id'] for r in self.data['_grist_Tables_column']) + 1
                self.data['_grist_Tables_column'].append({'id': ident, 'fields': dict(values, parentId=1, colId=action[2])})
            elif kind == 'AddRecord':
                ident = max(r['id'] for r in self.data[table]) + 1
                self.data[table].append({'id': ident, 'fields': values})
            elif kind == 'UpdateRecord':
                next(r for r in self.data[table] if r['id'] == ident)['fields'].update(values)
            else: raise AssertionError(action)


class SpectralSchemaTests(unittest.TestCase):
    def test_migration_preserves_data_and_layout_and_is_idempotent(self):
        c = Client(); before = copy.deepcopy(c.data['Experiments']); list_section = copy.deepcopy(c.data['_grist_Views_section'][-1])
        self.assertTrue(upgrade_spectral_schema(c))
        self.assertEqual(c.actions, [])
        upgrade_spectral_schema(c, apply=True)
        self.assertEqual(c.data['Experiments'], before)
        self.assertEqual(c.data['_grist_Views_section'][-1], list_section)
        layout = json.loads(c.data['_grist_Views_section'][0]['fields']['layoutSpec'])
        self.assertEqual(layout['custom'], 'keep')
        self.assertEqual(layout['children'][0], {'leaf': 100})
        for section in (23, 30, 41):
            rows = [r for r in c.data['_grist_Views_section_field'] if r['fields']['parentId'] == section and r['fields']['colRef'] in (11, 12)]
            self.assertEqual(len(rows), 2)
            if section == 23: self.assertEqual({r['id'] for r in rows}, {r['leaf'] for r in layout['children'][1]['children']})
        self.assertEqual(c.data['_grist_Views_section'][1]['fields']['layoutSpec'], '')
        self.assertEqual(upgrade_spectral_schema(c, apply=True), [])

    def test_custom_column_rejected_before_changes(self):
        c = Client(); c.data['_grist_Tables_column'].append({'id': 11, 'fields': {'parentId': 1, 'colId': 'SpectralWidths_ppm', 'type': 'Numeric'}})
        with self.assertRaises(ValueError): upgrade_spectral_schema(c, apply=True)
        self.assertEqual(c.actions, [])

    def test_fresh_columns_are_canonical_stored_text(self):
        names = {'SpectralWidths_ppm', 'SpectralOffsets_ppm'}
        self.assertTrue(names <= CANONICAL_ACQUISITION_FIELDS)
        self.assertEqual({c['id'] for c in columns()['Experiments'] if c['id'] in names and c['type'] == 'Text' and not c['isFormula']}, names)


if __name__ == '__main__': unittest.main()
