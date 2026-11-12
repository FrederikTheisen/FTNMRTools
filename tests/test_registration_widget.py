import unittest
from nmr_register.registration_widget import (sample_actions,experiment_actions,registration_html,
    NUMERIC_RULES,number,prepared,JS_LOGIC,labelling_state,labelling)


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

    def test_study_is_explicit_and_never_guessed(self):
        action=experiment_actions({'AcquisitionStatus':'Planned','Study':'DNA','AnalysisStatus':'Analyzed'},self.tables,now=123)['actions'][0]
        self.assertEqual(action[3]['Study'],'DNA')
        self.assertEqual(action[3]['AnalysisStatus'],'Analyzed')
        self.assertEqual(experiment_actions({'AcquisitionStatus':'Planned'},self.tables,now=123)['actions'][0][3]['Study'],'')
        with self.assertRaisesRegex(ValueError,'Invalid study'):
            experiment_actions({'AcquisitionStatus':'Planned','Study':'unexpected'},self.tables,now=123)
        with self.assertRaisesRegex(ValueError,'Invalid analysis status'):
            experiment_actions({'AcquisitionStatus':'Planned','AnalysisStatus':'unexpected'},self.tables,now=123)
        self.assertIn("['Electrostatics','Electrostatics']",registration_html('experiment'))

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

    def test_counter_auto_allocates_without_changing_existing_ids(self):
        first=sample_actions(self.values,self.tables)['actions'][0][3]
        self.tables['Samples']=[{'id':3,**first}]
        next_plan=sample_actions(self.values,self.tables)
        self.assertEqual(next_plan['sample_ids'],['202610_15N_HeTDR1_DNA_titration_2'])
        self.assertEqual(next_plan['actions'][0][3]['SampleNumber'],'2')
        self.assertEqual(self.tables['Samples'][0]['SampleID'],first['SampleID'])
        self.assertTrue(all(action[0]=='AddRecord' for action in next_plan['actions']))
        self.tables['Series']=[{'id':4,'SeriesID':'existing'}]
        with self.assertRaisesRegex(ValueError,'Series ID already exists'):
            sample_actions({**self.values,'Series':'new','SeriesID':'existing'},self.tables)

    def test_invalid_values_and_repeated_bulk_values_fail(self):
        for value in ['nan','Infinity','-1','0x10','1_000']:
            with self.subTest(value=value),self.assertRaises(ValueError):sample_actions({**self.values,'DNA_uM':value},self.tables)
        for key,value in [('pH','15'),('TFE_pct','101'),('Prepared','2026-02-30')]:
            with self.subTest(key=key),self.assertRaises(ValueError):sample_actions({**self.values,key:value},self.tables)
        with self.assertRaisesRegex(ValueError,'Repeated varying'):
            sample_actions({**self.values,'mode':'bulk','VaryValues':'0,0','VaryField':'DNA_uM'},self.tables)

    def test_storage_grid_created_idempotently_and_slots_checked(self):
        self.tables['StorageGrid']=[{'id':4,'StorageBox':'legacy','Row':'A'}]
        values={**self.values,'StorageRow':'A','StorageColumn':'12','mode':'bulk','VaryValues':'0,100','VaryField':'DNA_uM'}
        plan=sample_actions(values,self.tables)
        grid=[a for a in plan['actions'] if a[1]=='StorageGrid'];self.assertEqual(grid,[])
        tubes=[a[3] for a in plan['actions'] if a[1]=='Samples']
        self.assertEqual([(s['StorageRow'],s['StorageColumn']) for s in tubes],[('A',12),('B',1)])
        self.tables['Samples']=[{'id':2,'StorageBox':'old box','StorageRow':'a','StorageColumn':'12'}]
        with self.assertRaisesRegex(ValueError,'occupied'):sample_actions(values,self.tables)
        with self.assertRaises(ValueError):sample_actions({**self.values,'StorageRow':'AB','StorageColumn':'1'},self.tables)

    def test_visual_rack_position_keeps_existing_storage_encoding(self):
        plan=sample_actions({**self.values,'StorageSlot':'A6'},self.tables)
        fields=plan['actions'][0][3]
        self.assertEqual((fields['StorageRow'],fields['StorageColumn']),('A',6))
        self.assertEqual(plan['storage_slots'],['A6'])
        bulk=sample_actions({**self.values,'mode':'bulk','VaryValues':'0,100,200','VaryField':'DNA_uM','StorageSlot':'A11'},self.tables)
        self.assertEqual(bulk['storage_slots'],['A11','A12','B1'])
        self.assertEqual([(a[3]['StorageRow'],a[3]['StorageColumn']) for a in bulk['actions']], [('A',11),('A',12),('B',1)])

    def test_bulk_rack_preview_does_not_skip_occupied_slot(self):
        self.tables['Samples']=[{'id':2,'StorageRow':'B','StorageColumn':1,'Available':'Dead'}]
        values={**self.values,'mode':'bulk','VaryValues':'0,100','VaryField':'DNA_uM','StorageSlot':'A12'}
        with self.assertRaisesRegex(ValueError,'occupied'):
            sample_actions(values,self.tables)
        unplaced=sample_actions({**self.values,'mode':'bulk','VaryValues':'0,100','VaryField':'DNA_uM'},self.tables)
        self.assertEqual(unplaced['storage_slots'],['',''])
        self.assertTrue(all(a[3]['StorageRow']=='' and a[3]['StorageColumn'] is None for a in unplaced['actions']))

    def test_experiment_all_unknown_numeric_null_and_temperature_zero(self):
        plan=experiment_actions({'Temperature_C':'0','AcquisitionStatus':'Planned'},self.tables,now=123)
        fields=plan['actions'][0][3]
        self.assertEqual(fields['Temperature_C'],0);self.assertEqual(fields['RegisteredAt'],123)
        self.assertEqual(fields['AcquisitionStatus'],'Planned')
        for name in ['Field_MHz','D1_s','SpinLock_Hz','SpinLockOffset_Hz']+['Snapshot'+k for k in NUMERIC_RULES]:
            self.assertIn(name,fields);self.assertIsNone(fields[name])
        self.assertEqual(fields['DelayList_s'],'')
        self.assertFalse(any(name.startswith('Observed') for name in fields))

    def test_known_manual_metadata_uses_canonical_fact_fields(self):
        values={'Instrument':'Bruker 600','ExperimentType':'T1','PulseProgram':'t1','Acquired':'2026-11',
            'Field_MHz':'600.13','Temperature_C':'25','D1_s':'0','DelayList_s':'0.1, 0.5, 1.0',
            'SpinLock_Hz':'1000','SpinLockOffset_Hz':'-25'}
        fields=experiment_actions(values,self.tables,now=123)['actions'][0][3]
        self.assertEqual({k:fields[k] for k in ('Instrument','ExperimentType','PulseProgram','AcquisitionStatus','AcquisitionDateText','Field_MHz','Temperature_C','D1_s','DelayList_s','SpinLock_Hz','SpinLockOffset_Hz')},
            {'Instrument':'Bruker 600','ExperimentType':'T1','PulseProgram':'t1','AcquisitionStatus':'Planned','AcquisitionDateText':'2026-11','Field_MHz':600.13,'Temperature_C':25,'D1_s':0,'DelayList_s':'0.1, 0.5, 1.0','SpinLock_Hz':1000,'SpinLockOffset_Hz':-25})
        self.assertFalse(any(name.startswith('Observed') for name in fields))

    def test_experiment_form_has_plain_canonical_facts_without_status_override(self):
        html=registration_html('experiment')
        self.assertIn("field('experiment','AcquisitionStatus'",html)
        self.assertNotIn('Acquisition status override',html)
        for name in ('Instrument','ExperimentType','Acquired','AcquisitionStatus','Field_MHz','Temperature_C','PulseProgram','D1_s','DelayList_s','SpinLock_Hz','SpinLockOffset_Hz'):
            self.assertIn(name,html)
        self.assertIn("tables.Experiments.map(r=>r[name])",html)
        self.assertNotIn("tables.Experiments.flatMap(r=>[r[name],r[observed]])",html)

    def test_linked_experiment_omits_snapshot_so_native_trigger_copies(self):
        self.tables['Series']=[{'id':2,'SeriesID':'S'}]
        self.tables['Samples']=[{'id':5,'SampleID':'tube','Series':2,'DNA_uM':0}]
        fields=experiment_actions({'Series':'2','Sample':'5','AcquisitionStatus':'Planned'},self.tables)['actions'][0][3]
        self.assertEqual(fields['Sample'],5);self.assertFalse(any(k.startswith('Snapshot') for k in fields))
        self.assertEqual(fields['AcquisitionStatus'],'Planned')
        with self.assertRaisesRegex(ValueError,'not in the selected series'):
            experiment_actions({'Series':'2','Sample':'6','AcquisitionStatus':'Planned'},self.tables)

    def test_exact_duplicate_acquisition_or_protected_target_fail(self):
        values={'SourceID':'nmrdata','Folder':'folder','Expno':'12','AcquisitionStatus':'Planned'}
        self.tables['Experiments']=[{'id':7,'ObservedSourceID':'nmrdata','ObservedFolder':'folder','ObservedExpno':'12'}]
        with self.assertRaisesRegex(ValueError,'already exists'):experiment_actions(values,self.tables)
        with self.assertRaisesRegex(ValueError,'protected'):experiment_actions({'Folder':'211@950','AcquisitionStatus':'Planned'},self.tables)

    def test_experiment_sample_filter_unknown_available_and_historical_opt_in(self):
        self.tables['Samples']=[{'id':5,'SampleID':'tube','Protein':'TDR1','Available':'Unknown'}]
        values={'Sample':'5','ProteinFilter':'TDR1','AcquisitionStatus':'Planned'}
        self.assertEqual(experiment_actions(values,self.tables)['actions'][0][3]['Sample'],5)
        self.tables['Samples'][0]['Available']='Unavailable'
        with self.assertRaisesRegex(ValueError,'dead'):experiment_actions(values,self.tables)
        self.assertEqual(experiment_actions({**values,'IncludeUnavailable':True},self.tables)['actions'][0][3]['Sample'],5)
        with self.assertRaisesRegex(ValueError,'selected protein'):
            experiment_actions({**values,'IncludeUnavailable':True,'ProteinFilter':'other'},self.tables)
        self.assertEqual(experiment_actions({},self.tables)['actions'][0][3]['AcquisitionStatus'],'Planned')
        with self.assertRaisesRegex(ValueError,'initial acquisition status'):experiment_actions({'AcquisitionStatus':'Running'},self.tables)

    def test_standard_labelling_combinations_and_custom_series_text_are_preserved(self):
        self.assertEqual(labelling({'Label15N':True,'Label13C':True,'Label2H':True,'IsotopeCustom':''}),'15N+13C+2H')
        self.assertEqual(labelling_state('15N / 13C'),{'labels':['15N','13C'],'custom':''})
        custom='U-15N, 70% 2H'
        self.assertEqual(labelling_state(custom),{'labels':[],'custom':custom})
        values={**self.values,'Series':'4'};values.pop('Isotope')
        for isotope in ('15N+13C+2H',custom):
            self.tables['Series']=[{'id':4,'Isotope':isotope,'NaCl_mM':0}]
            fields=sample_actions(values,self.tables)['actions'][0][3]
            self.assertEqual(fields['Isotope'],isotope)
        fields=sample_actions({**self.values,'Label15N':False,'IsotopeCustom':custom},self.tables)['actions'][0][3]
        self.assertEqual(fields['Isotope'],custom)

    def test_optional_sample_label_derives_composition_or_basic_identity(self):
        values={**self.values,'DistinguishingCondition':'','NaCl_mM':'0','DNA_uM':'100'}
        fields=sample_actions(values,self.tables)['actions'][0][3]
        self.assertEqual(fields['SampleID'],'202610_15N_HeTDR1_0mMNaCl_100uMDNA_1')
        self.assertEqual(fields['DistinguishingCondition'],'0mMNaCl_100uMDNA')
        values.pop('NaCl_mM');values.pop('DNA_uM')
        fields=sample_actions(values,self.tables)['actions'][0][3]
        self.assertEqual(fields['SampleID'],'202610_15N_HeTDR1_1')
        self.assertEqual(fields['Available'],'Usable')
        self.assertTrue(all(fields[name] is None for name in NUMERIC_RULES))

    def test_full_date_and_month_keep_their_precision(self):
        full=sample_actions({**self.values,'Prepared':'2026-10-30'},self.tables)['actions'][0][3]
        month=sample_actions({**self.values,'Prepared':'2026-10'},self.tables)['actions'][0][3]
        self.assertIsNotNone(full['PreparedAt']);self.assertEqual(full['PreparedDateText'],'')
        self.assertIsNone(month['PreparedAt']);self.assertEqual(month['PreparedDateText'],'2026-10')
        self.assertTrue(full['SampleID'].startswith('20261030_'));self.assertTrue(month['SampleID'].startswith('202610_'))

    def test_embed_has_only_supported_sdk_and_safe_dom_operations(self):
        html=registration_html('experiment')
        self.assertIn('<script src="grist-plugin-api.js"></script>',html)
        self.assertIn("grist.ready({requiredAccess:'full'})",html)
        self.assertIn('grist.docApi.applyUserActions(plan.actions)',html)
        self.assertIn('grist.docApi.fetchTable(name)',html)
        self.assertNotIn('innerHTML',html);self.assertNotIn('localStorage',html)
        self.assertNotIn("field('shared','SampleNumber'",html)
        self.assertNotIn("field('identity','SampleNumber'",html)
        self.assertIn('Additional composition (optional)',html)
        self.assertIn("'Tube status'",html);self.assertIn('Dead means permanently unavailable',html)
        self.assertNotIn("['Retired','Retired']",html)
        self.assertIn('Suspected degradation',html)
        self.assertIn("'NaCl (mM)'",html)
        self.assertNotIn('/*LOGIC*/',html);self.assertIn('.value="experiment";',html)
        self.assertIn('const Rack=(()=>{',html)
        self.assertIn("Rack.mountRack(byId('rackPicker')",html)
        self.assertIn('columns A to H, rows 1 to 12',html)
        self.assertNotIn('First row A–H',html);self.assertNotIn('First column 1–12',html)
        self.assertIn('Bulk positions:',html);self.assertIn('will not skip occupied slots',html)
        with self.assertRaises(ValueError):registration_html('</script>')
