"""Only visible human experiment fields are substituted in place."""
import unittest
from nmr_register.experiments_ui import display_partial_assessment


class Client:
    def __init__(self):
        self.applied = []
        self.data = {
            '_grist_Tables': [{'id': 1, 'fields': {'tableId': 'Experiments'}}, {'id': 2, 'fields': {'tableId': 'Samples'}}],
            '_grist_Tables_column': [{'id': 10, 'fields': {'parentId': 1, 'colId': 'Assessment'}}, {'id': 11, 'fields': {'parentId': 1, 'colId': 'AssessmentDisplay'}}],
            '_grist_Views_section': [
                {'id': 20, 'fields': {'tableRef': 1, 'parentId': 3, 'parentKey': 'record'}},
                {'id': 21, 'fields': {'tableRef': 1, 'parentId': 4, 'parentKey': 'single'}},
                {'id': 22, 'fields': {'tableRef': 1, 'parentId': 0, 'parentKey': 'record'}},
                {'id': 23, 'fields': {'tableRef': 1, 'parentId': 3, 'parentKey': 'custom'}},
                {'id': 24, 'fields': {'tableRef': 2, 'parentId': 3, 'parentKey': 'record'}},
            ],
            '_grist_Views_section_field': [
                {'id': 30, 'fields': {'parentId': 20, 'colRef': 10, 'width': 100, 'parentPos': 2, 'widgetOptions': 'style'}},
                {'id': 31, 'fields': {'parentId': 21, 'colRef': 10, 'width': 300, 'parentPos': 1}},
                *[{'id': 32 + i, 'fields': {'parentId': section, 'colRef': 10, 'width': 100}} for i, section in enumerate((22, 23, 24))],
                {'id': 35, 'fields': {'parentId': 20, 'colRef': 99, 'width': 100}},
            ],
        }
    def records(self, name): return self.data[name]
    def apply(self, actions): self.applied.extend(actions)


class PartialAssessmentUITest(unittest.TestCase):
    def test_targeted_in_place_updates_preserve_metadata(self):
        client = Client()
        expected = [
            ['UpdateRecord', '_grist_Views_section_field', 30, {'colRef': 11, 'width': 220}],
            ['UpdateRecord', '_grist_Views_section_field', 31, {'colRef': 11}],
        ]
        self.assertEqual(display_partial_assessment(client), expected)
        self.assertEqual(client.applied, [])
        self.assertEqual(display_partial_assessment(client, apply=True), expected)
        self.assertEqual(client.applied, expected)
        for action in expected:
            row = next(r for r in client.data['_grist_Views_section_field'] if r['id'] == action[2])
            row['fields'].update(action[3])
        self.assertEqual(display_partial_assessment(client), [])
        self.assertEqual(client.data['_grist_Views_section_field'][0]['fields']['parentPos'], 2)
        self.assertEqual(client.data['_grist_Views_section_field'][0]['fields']['widgetOptions'], 'style')


if __name__ == '__main__': unittest.main()
