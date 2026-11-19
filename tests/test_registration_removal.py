"""Removing the retired field preserves layouts and rejects unknown dependencies."""
import copy
import json
import unittest
from nmr_register.experiments_ui import remove_registration_timestamp
from nmr_register.schema import ATTENTION_FORMULA, PRE_REGISTRATION_REMOVAL_ATTENTION_FORMULA

class Metadata:
    def __init__(self):
        self.data = {
            '_grist_Filters': [], '_grist_ACLResources': [], '_grist_ACLRules': [],
            '_grist_Tables': [{'id': 4, 'fields': {'tableId': 'Experiments'}}],
            '_grist_Tables_column': [
                {'id': 74, 'fields': {'parentId': 4, 'colId': 'RegisteredAt', 'type': 'DateTime', 'isFormula': False, 'formula': ''}},
                {'id': 99, 'fields': {'parentId': 4, 'colId': 'AttentionReason', 'type': 'Text', 'isFormula': True, 'formula': PRE_REGISTRATION_REMOVAL_ATTENTION_FORMULA}},
            ],
            '_grist_Views_section_field': [{'id': 450, 'fields': {'parentId': 23, 'colRef': 74}}],
            '_grist_Views_section': [
                {'id': 23, 'fields': {'layoutSpec': json.dumps({'children': [{'children': [{'leaf': 440}, {'leaf': 450}, {'leaf': 441}]}, {'children': [{'leaf': 450}]}], 'custom': 'keep'}), 'sortColRefs': '[74, -88]'}},
                {'id': 30, 'fields': {'layoutSpec': '', 'sortColRefs': '', 'linkSrcSectionRef': 29}},
            ],
        }
        self.actions = []
    def records(self, table): return copy.deepcopy(self.data[table])
    def apply(self, actions): self.actions.append(copy.deepcopy(actions))

class RegistrationRemovalTests(unittest.TestCase):
    def test_plan_preserves_other_fields_layout_and_sorting(self):
        c = Metadata(); before = copy.deepcopy(c.data)
        actions = remove_registration_timestamp(c)
        self.assertEqual(c.actions, []); self.assertEqual(c.data, before)
        self.assertEqual(actions[0], ['UpdateRecord', '_grist_Tables_column', 99, {'formula': ATTENTION_FORMULA}])
        self.assertEqual(actions[-1], ['RemoveColumn', 'Experiments', 'RegisteredAt'])
        card = next(a[3] for a in actions if a[1] == '_grist_Views_section')
        self.assertEqual(json.loads(card['layoutSpec']), {'children': [{'children': [{'leaf': 440}, {'leaf': 441}]}], 'custom': 'keep'})
        self.assertEqual(json.loads(card['sortColRefs']), [-88])
        self.assertFalse(any(a[1] == 'Experiments' and a[0] != 'RemoveColumn' for a in actions))
        remove_registration_timestamp(c, apply=True)
        self.assertEqual(c.actions, [actions])
    def test_already_removed_is_unchanged(self):
        c = Metadata(); c.data['_grist_Tables_column'].pop(0)
        c.data['_grist_Tables_column'][0]['fields']['formula'] = ATTENTION_FORMULA
        self.assertEqual(remove_registration_timestamp(c, apply=True), [])
        self.assertEqual(c.actions, [])
    def test_custom_dependency_rejected_before_any_write(self):
        for dependency in ('custom_formula', 'summary', 'linked_widget', 'filter', 'acl_resource', 'acl_formula'):
            c = Metadata()
            if dependency == 'custom_formula':
                c.data['_grist_Tables_column'].append({'id': 120, 'fields': {'parentId': 4, 'colId': 'Custom', 'formula': '$RegisteredAt', 'isFormula': True}})
            elif dependency == 'summary': c.data['_grist_Tables_column'].append({'id': 120, 'fields': {'summarySourceCol': 74}})
            elif dependency == 'linked_widget': c.data['_grist_Views_section'][0]['fields']['linkSrcColRef'] = 74
            elif dependency == 'filter': c.data['_grist_Filters'].append({'id': 1, 'fields': {'colRef': 74}})
            elif dependency == 'acl_resource': c.data['_grist_ACLResources'].append({'id': 1, 'fields': {'tableId': 'Experiments', 'colIds': 'RegisteredAt'}})
            else: c.data['_grist_ACLRules'].append({'id': 1, 'fields': {'aclFormula': 'rec.RegisteredAt > 0'}})
            with self.subTest(dependency=dependency):
                with self.assertRaises(ValueError): remove_registration_timestamp(c, apply=True)
                self.assertEqual(c.actions, [])

if __name__ == '__main__': unittest.main()
