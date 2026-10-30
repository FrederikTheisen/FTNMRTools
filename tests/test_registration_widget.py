import unittest
from nmr_register.registration_widget import (sample_actions,experiment_actions,registration_html,
    NUMERIC_RULES,number,prepared,JS_LOGIC)


class WidgetTests(unittest.TestCase):
    def setUp(self):
        self.tables={'Series':[],'Samples':[],'Experiments':[],'StorageGrid':[]}
        self.values={'Prepared':'2026-10','Protein':'HeTDR1','Isotope':'15N','DistinguishingCondition':'DNA_titration','SampleNumber':'1'}

    def test_single_initializes_unknowns_null_keeps_literal_zero(self):
        values={**self.values,'NaCl_mM':'0','DNA_uM':'0','pH':''}
        plan=sample_actions(values,self.tables);fields=plan['actions'][0][3]
        self.assertEqual(fields['NaCl_mM'],0);self.assertEqual(fields['DNA_uM'],0)
        self.assertTrue(all(name in fields for name in NUMERIC_RULES))
        self.assertIsNone(fields['pH']);self.assertIsNone(fields['TFE_pct']);self.assertIsNone(fields['StorageColumn'])
        self.assertEqual(fields['PreparedDateText'],'2026-10');self.assertIsNone(fields['PreparedAt'])
        self.assertEqual(fields['SampleID'],'202610_15N_HeTDR1_DNA_titration_1')

    def test_bulk_series_has_complete_shared_values_and_distinct_tubes(self):
        values={**self.values,'mode':'bulk','Series':'new','SeriesID':'DNA series','NaCl_mM':'50','VaryField':'DNA_uM','VaryValues':'0, 100\n200;500'}
        plan=sample_actions(values,self.tables);series=plan['actions'][0]
        self.assertEqual(series[:3],['AddRecord','Series',1]);self.assertEqual(series[3]['NaCl_mM'],50)
        self.assertIsNone(series[3]['DNA_uM'])
        tubes=[a[3] for a in plan['actions'] if a[1]=='Samples']
        self.assertEqual([r['DNA_uM'] for r in tubes],[0,100,200,500])
        self.assertEqual(len(set(r['SampleID'] for r in tubes)),4)
        self.assertTrue(all(r['Series']==1 for r in tubes))
        self.assertEqual(tubes[0]['SampleID'],'202610_15N_HeTDR1_DNA_titration_0uMDNA_1')
        self.assertEqual([a[2] for a in plan['actions'] if a[1]=='Samples'],[1,2,3,4])

    def test_duplicate_sample_and_series_ids_fail(self):
        first=sample_actions(self.values,self.tables)['actions'][0][3]
        self.tables['Samples']=[{'id':3,**first}]
        with self.assertRaisesRegex(ValueError,'Sample ID already exists'):sample_actions(self.values,self.tables)
        self.tables['Series']=[{'id':4,'SeriesID':'existing'}]
        with self.assertRaisesRegex(ValueError,'Series ID already exists'):
            sample_actions({**self.values,'Series':'new','SeriesID':'existing'},self.tables)

    def test_invalid_values_and_repeated_bulk_values_fail(self):
        for value in ['nan','Infinity','-1','0x10','1_000']:
            with self.subTest(value=value),self.assertRaises(ValueError):sample_actions({**self.values,'DNA_uM':value},self.tables)
        for key,value in [('pH','15'),('TFE_pct','101'),('Prepared','2026-02-30'),('SampleNumber','0')]:
            with self.subTest(key=key),self.assertRaises(ValueError):sample_actions({**self.values,key:value},self.tables)
        with self.assertRaisesRegex(ValueError,'Repeated varying'):
            sample_actions({**self.values,'mode':'bulk','VaryValues':'0,0','VaryField':'DNA_uM'},self.tables)

    def test_storage_grid_created_idempotently_and_slots_checked(self):
        self.tables['StorageGrid']=[{'id':4,'StorageBox':'Freezer1','Row':'A'}]
        values={**self.values,'StorageBox':'Freezer1','StorageRow':'A','StorageColumn':'12','mode':'bulk','VaryValues':'0,100','VaryField':'DNA_uM'}
        plan=sample_actions(values,self.tables)
        grid=[a for a in plan['actions'] if a[1]=='StorageGrid'];self.assertEqual(len(grid),7)
        self.assertEqual(grid[0][2],5);self.assertEqual(grid[0][3]['Row'],'B')
        tubes=[a[3] for a in plan['actions'] if a[1]=='Samples']
        self.assertEqual([(s['StorageRow'],s['StorageColumn']) for s in tubes],[('A',12),('B',1)])
        self.tables['Samples']=[{'id':2,'StorageBox':'Freezer1','StorageRow':'A','StorageColumn':12}]
        with self.assertRaisesRegex(ValueError,'occupied'):sample_actions(values,self.tables)
        with self.assertRaises(ValueError):sample_actions({**self.values,'StorageBox':'box','StorageRow':'AB','StorageColumn':'1'},self.tables)

    def test_experiment_all_unknown_numeric_null_and_temperature_zero(self):
        plan=experiment_actions({'Temperature_C':'0','AcquisitionStatus':'Planned'},self.tables,now=123)
        fields=plan['actions'][0][3]
        self.assertEqual(fields['Temperature_C'],0);self.assertEqual(fields['RegisteredAt'],123)
        for name in ['Field_MHz','ObservedField_MHz','ObservedTemperature_C','RawBytes']+['Snapshot'+k for k in NUMERIC_RULES]:
            self.assertIn(name,fields);self.assertIsNone(fields[name])

    def test_linked_experiment_omits_snapshot_so_native_trigger_copies(self):
        self.tables['Series']=[{'id':2,'SeriesID':'S'}]
        self.tables['Samples']=[{'id':5,'SampleID':'tube','Series':2,'DNA_uM':0}]
        fields=experiment_actions({'Series':'2','Sample':'5','AcquisitionStatus':'Acquired'},self.tables)['actions'][0][3]
        self.assertEqual(fields['Sample'],5);self.assertFalse(any(k.startswith('Snapshot') for k in fields))
        with self.assertRaisesRegex(ValueError,'not in the selected series'):
            experiment_actions({'Series':'2','Sample':'6','AcquisitionStatus':'Acquired'},self.tables)

    def test_exact_duplicate_acquisition_or_protected_target_fail(self):
        values={'SourceID':'nmrdata','Folder':'folder','Expno':'12','AcquisitionStatus':'Acquired'}
        self.tables['Experiments']=[{'id':7,'ObservedSourceID':'nmrdata','ObservedFolder':'folder','ObservedExpno':'12'}]
        with self.assertRaisesRegex(ValueError,'already exists'):experiment_actions(values,self.tables)
        with self.assertRaisesRegex(ValueError,'protected'):experiment_actions({'Folder':'211@950','AcquisitionStatus':'Acquired'},self.tables)

    def test_experiment_sample_filter_unknown_available_and_historical_opt_in(self):
        self.tables['Samples']=[{'id':5,'SampleID':'tube','Protein':'TDR1','Available':'Unknown'}]
        values={'Sample':'5','ProteinFilter':'TDR1','AcquisitionStatus':'Acquired'}
        self.assertEqual(experiment_actions(values,self.tables)['actions'][0][3]['Sample'],5)
        self.tables['Samples'][0]['Available']='Unavailable'
        with self.assertRaisesRegex(ValueError,'unavailable'):experiment_actions(values,self.tables)
        self.assertEqual(experiment_actions({**values,'IncludeUnavailable':True},self.tables)['actions'][0][3]['Sample'],5)
        with self.assertRaisesRegex(ValueError,'selected protein'):
            experiment_actions({**values,'IncludeUnavailable':True,'ProteinFilter':'other'},self.tables)
        with self.assertRaisesRegex(ValueError,'Planned or Acquired'):experiment_actions({},self.tables)

    def test_embed_has_only_supported_sdk_and_safe_dom_operations(self):
        html=registration_html('experiment')
        self.assertIn('<script src="grist-plugin-api.js"></script>',html)
        self.assertIn("grist.ready({requiredAccess:'full'})",html)
        self.assertIn('grist.docApi.applyUserActions(plan.actions)',html)
        self.assertIn('grist.docApi.fetchTable(name)',html)
        self.assertNotIn('innerHTML',html);self.assertNotIn('localStorage',html)
        self.assertNotIn('/*LOGIC*/',html);self.assertIn('.value="experiment";',html)
        with self.assertRaises(ValueError):registration_html('</script>')
