import copy
import json
import unittest

from nmr_register.experiments_ui import configure_experiment_type_filter


class Client:
    def __init__(self, condition=None, pinned=False):
        self.data = {
            '_grist_Views': [{'id': 8, 'fields': {'name': 'Experiments'}}],
            '_grist_Tables': [{'id': 1, 'fields': {'tableId': 'Experiments'}}],
            '_grist_Tables_column': [{'id': 58, 'fields': {'parentId': 1, 'colId': 'ExperimentType'}}],
            '_grist_Views_section': [{'id': 22, 'fields': {'parentId': 8, 'parentKey': 'record', 'title': 'Experiments', 'tableRef': 1}}],
            '_grist_Filters': [] if condition is None else [
                {'id': 37, 'fields': {'viewSectionRef': 22, 'colRef': 58,
                                      'filter': json.dumps(condition), 'pinned': pinned, 'other': 'keep'}}],
        }
        self.calls = []

    def records(self, name):
        return copy.deepcopy(self.data[name])

    def apply(self, actions):
        self.calls.extend(copy.deepcopy(actions))
        for item in actions:
            action, table, row_id = item[:3]
            if action == 'RemoveRecord':
                self.data[table] = [row for row in self.data[table] if row['id'] != row_id]
                continue
            fields = item[3]
            if action == 'AddRecord':
                self.data[table].append({'id': 100, 'fields': copy.deepcopy(fields)})
            else:
                next(row for row in self.data[table] if row['id'] == row_id)['fields'].update(fields)


class TypeFilterTests(unittest.TestCase):
    def test_neutral_legacy_filter_removed_only_from_main_grid(self):
        client = Client({'excluded': ['Calibration', '1D', 'Other']}, pinned=True)
        client.data['_grist_Tables_column'].append(
            {'id': 96, 'fields': {'parentId': 1, 'colId': 'EffectiveType',
                                  'isFormula': True, 'formula': '$ExperimentType'}})
        client.data['_grist_Filters'].extend([
            {'id': 2, 'fields': {'viewSectionRef': 22, 'colRef': 96, 'filter': '{}', 'pinned': True}},
            {'id': 3, 'fields': {'viewSectionRef': 23, 'colRef': 96, 'filter': '{}', 'pinned': True}},
            {'id': 4, 'fields': {'viewSectionRef': 22, 'colRef': 99, 'filter': '{}', 'pinned': True}},
        ])
        before = copy.deepcopy(client.data['_grist_Filters'])
        self.assertEqual(configure_experiment_type_filter(client)['actions'],
                         [['RemoveRecord', '_grist_Filters', 2]])
        self.assertEqual(client.data['_grist_Filters'], before)
        configure_experiment_type_filter(client, apply=True)
        self.assertEqual(client.data['_grist_Filters'], [r for r in before if r['id'] != 2])
        self.assertEqual(configure_experiment_type_filter(client, apply=True)['actions'], [])

    def test_active_legacy_filter_requires_review_before_writes(self):
        client = Client()
        client.data['_grist_Tables_column'].append(
            {'id': 96, 'fields': {'parentId': 1, 'colId': 'EffectiveType',
                                  'isFormula': True, 'formula': '$ExperimentType'}})
        client.data['_grist_Filters'].append(
            {'id': 2, 'fields': {'viewSectionRef': 22, 'colRef': 96,
                                  'filter': '{"included":["T1"]}', 'pinned': True}})
        with self.assertRaises(ValueError):
            configure_experiment_type_filter(client, apply=True)
        self.assertEqual(client.calls, [])

    def test_existing_exclusions_and_metadata_are_preserved(self):
        client = Client({'excluded': ['Other', 'Calibration'], 'flag': 'keep'})
        result = configure_experiment_type_filter(client)
        self.assertEqual(result['filter_value'], {'excluded': ['Other', 'Calibration', '1D'], 'flag': 'keep'})
        self.assertEqual(client.calls, [])
        configure_experiment_type_filter(client, apply=True)
        fields = client.data['_grist_Filters'][0]['fields']
        self.assertTrue(fields['pinned'])
        self.assertEqual(fields['other'], 'keep')
        self.assertEqual(configure_experiment_type_filter(client, apply=True)['actions'], [])

    def test_included_filter_preserves_other_selections(self):
        client = Client({'included': ['1D', 'T1', 'Calibration'], 'extra': 3}, pinned=True)
        configure_experiment_type_filter(client, apply=True)
        self.assertEqual(json.loads(client.data['_grist_Filters'][0]['fields']['filter']),
                         {'included': ['T1'], 'extra': 3})
        self.assertEqual(configure_experiment_type_filter(client, apply=True)['actions'], [])

    def test_missing_filter_is_created_and_idempotent(self):
        client = Client()
        configure_experiment_type_filter(client, apply=True)
        fields = client.data['_grist_Filters'][0]['fields']
        self.assertEqual(json.loads(fields['filter']), {'excluded': ['Calibration', '1D']})
        self.assertEqual((fields['viewSectionRef'], fields['colRef'], fields['pinned']), (22, 58, True))
        self.assertEqual(configure_experiment_type_filter(client, apply=True)['actions'], [])

    def test_rejects_ambiguous_or_linked_configuration(self):
        client = Client({'included': ['T1'], 'excluded': ['Other']})
        with self.assertRaises(ValueError):
            configure_experiment_type_filter(client, apply=True)
        client = Client({})
        client.data['_grist_Views_section'][0]['fields']['linkSrcSectionRef'] = 47
        with self.assertRaises(ValueError):
            configure_experiment_type_filter(client, apply=True)
        self.assertEqual(client.calls, [])


if __name__ == '__main__':
    unittest.main()
