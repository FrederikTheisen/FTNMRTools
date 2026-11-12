"""The reviewed pilot only: existing rows, exact identities, current human fields."""
import copy
import csv
from pathlib import Path
import tempfile
import unittest
from nmr_register.migration import (REVIEWED_PILOT_IDENTITIES,reviewed_pilot_patch,
    apply_reviewed_pilot_patch,reconcile_sources,build_pilot)


class Client:
    def __init__(self,rows):self.rows=copy.deepcopy(rows);self.writes=[]
    def records(self,table):
        assert table=='Experiments'
        return copy.deepcopy(self.rows)
    def update_records(self,table,records):
        self.writes.append(copy.deepcopy(records))
        for record in records:next(row for row in self.rows if row['id']==record['id'])['fields'].update(record['fields'])


class AnalysisEvidenceMigrationTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.catalog=self.root/'acquisition-catalog.csv'
        self.rows=[]
        for rowid,(source,folder,expno) in REVIEWED_PILOT_IDENTITIES.items():
            self.rows.append({'id':rowid,'fields':{'SourceID':source,'Folder':folder,'Expno':expno,
                'SourceEvidence':'workbook row '+str(rowid),'Assessment':'Do not use' if rowid==1 else '',
                'AssessmentReason':'keep reason','Notes':'keep notes','Sample':10+rowid,
                'MigrationIssue':'Sample link unresolved' if rowid in (2,3,4,5) else '',
                'MetadataDifference':'observed Temperature_C differs from confirmed value' if rowid in (6,7,8) else '',
                'AcquisitionKey':source+':'+folder+'/'+expno if rowid<=8 else ''}})
        self.catalog_rows=[]
        for rowid in (1,6,7,8):
            _,folder,expno=REVIEWED_PILOT_IDENTITIES[rowid]
            self.catalog_rows.append({'relative_path':folder+'/'+expno,'status':'analysed' if rowid%2 else 'ANALYZED',
                'primary_record_status':'primary raw acquisition','related_primary_path':''})
        # Actual catalog spelling intentionally differs from reviewed workbook.
        for expno in ('5','6','7'):
            self.catalog_rows.append({'relative_path':'20251218_15N_heTDR1_50mMNaCl/'+expno,'status':'analysed',
                'primary_record_status':'primary raw acquisition','related_primary_path':''})
        self.write_catalog()
    def tearDown(self):self.tmp.cleanup()
    def write_catalog(self):
        names=('relative_path','status','primary_record_status','related_primary_path','SourceID')
        with self.catalog.open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=names);writer.writeheader();writer.writerows(self.catalog_rows)

    def test_reviewed_patch_membership_current_analysis_and_exact_typo_mismatch(self):
        patch=reviewed_pilot_patch(self.rows,self.catalog)
        self.assertEqual(patch['issues'],[])
        updates={record['id']:record['fields'] for record in patch['records']}
        self.assertEqual(set(updates),set(range(1,12)))
        self.assertEqual([rowid for rowid,fields in updates.items() if fields.get('AnalysisStatus')=='Analyzed'],[1,6,7,8])
        self.assertTrue(all(updates[rowid]['Study']=='DNA' for rowid in (2,3,4,5)))
        self.assertTrue(all(updates[rowid]['Study']=='Electrostatics' for rowid in (1,6,7,8,9,10,11)))
        self.assertNotIn('AnalysisStatus',updates[9])
        self.assertIn(str(self.catalog.resolve()),updates[6]['SourceEvidence'])
        self.assertIn('nmrdata:20251105_15N_heTDR1_NaCl/405',updates[6]['SourceEvidence'])
        self.assertTrue(all(set(fields)<= {'Study','AnalysisStatus','SourceEvidence'} for fields in updates.values()))

    def test_preserves_nonblank_human_status_study_and_all_assessments(self):
        self.rows[0]['fields'].update(Study='Other',AnalysisStatus='Not started')
        self.rows[5]['fields'].update(AnalysisStatus='In progress')
        self.rows[6]['fields'].update(AnalysisStatus='Analyzed')
        before=copy.deepcopy(self.rows)
        patch=reviewed_pilot_patch(self.rows,self.catalog)
        client=Client(self.rows)
        apply_reviewed_pilot_patch(patch,client,dry_run=False)
        self.assertEqual(client.rows[0]['fields']['Study'],'Other')
        self.assertEqual(client.rows[0]['fields']['AnalysisStatus'],'Not started')
        self.assertEqual(client.rows[5]['fields']['AnalysisStatus'],'In progress')
        for original,current in zip(before,client.rows):
            for name in ('Assessment','AssessmentReason','Notes','Sample'):
                self.assertEqual(original['fields'][name],current['fields'][name])
        self.assertEqual(len(client.rows),11)
        self.assertEqual(reviewed_pilot_patch(client.rows,self.catalog)['records'],[])

    def test_zero_dna_control_is_curated_dna_without_guessing_other_folders(self):
        try:import openpyxl
        except ImportError:self.skipTest('openpyxl optional')
        book=openpyxl.Workbook();book.active.title='Samples'
        book['Samples'].append(['position','id','protein','comment','unused','buffer','comment'])
        for name in ('Database','ANP32','N-Tail','TDR1','CtIP'):
            sheet=book.create_sheet(name);sheet.append(['date','instrument','user','folder','expno','sample','condition','temp','type','pulse','options','status','comment'])
        book['TDR1'].append(['2025-10','600_1300','','20251017_15N_HeTDR1_DNA','104','','0 uM DNA',25,'15N_HSQC','','','Done',''])
        book['TDR1'].append(['2025-10','600_1300','','unreviewed_DNA_name','1','','0 uM DNA',25,'15N_HSQC','','','Done',''])
        workbook=self.root/'pilot.xlsx';book.save(workbook)
        reviews=self.root/'reviews.tsv';reviews.write_text('target\tassessment\tcomment\n')
        plan=build_pilot(workbook,reviews,[],['20251017_15N_HeTDR1_DNA/104','unreviewed_DNA_name/1'],'nmrdata')
        values={entry['fields']['Folder']:entry['fields'] for entry in plan['experiments']}
        self.assertEqual(values['20251017_15N_HeTDR1_DNA']['Study'],'DNA')
        self.assertNotIn('Study',values['unreviewed_DNA_name'])
        other=build_pilot(workbook,reviews,[],['20251017_15N_HeTDR1_DNA/104'],'other-source')
        self.assertNotIn('Study',other['experiments'][0]['fields'])

    def test_source_ids_inconsistent_keys_and_duplicate_identities_fail_closed(self):
        self.rows[0]['fields']['SourceID']='other'
        self.rows[5]['fields']['AcquisitionKey']='nmrdata:other/405'
        self.rows.append(copy.deepcopy(self.rows[6]));self.rows[-1]['id']=20
        patch=reviewed_pilot_patch(self.rows,self.catalog)
        self.assertEqual(len(patch['issues']),3)
        self.assertFalse(set((1,6,7))&{record['id'] for record in patch['records']})
        with self.assertRaises(ValueError):apply_reviewed_pilot_patch(patch,Client(self.rows),False)
        with self.assertRaises(ValueError):reviewed_pilot_patch(self.rows,self.catalog,source_id='other')

    def test_catalog_copy_derived_duplicate_or_source_mismatch_never_promote(self):
        self.catalog_rows[0]['related_primary_path']='some/parent'
        self.catalog_rows[1]['primary_record_status']='processed only; raw absent'
        self.catalog_rows[2]['SourceID']='other-source'
        self.catalog_rows.append(copy.deepcopy(self.catalog_rows[3]))
        self.write_catalog()
        patch=reviewed_pilot_patch(self.rows,self.catalog)
        self.assertTrue(patch['issues'])
        self.assertFalse(any('AnalysisStatus' in record['fields'] for record in patch['records']))

    def test_existing_identity_ambiguity_blocks_current_analysis_but_metadata_warning_does_not(self):
        self.rows[0]['fields']['ReconciliationReason']='raw acquisition replaced'
        patch=reviewed_pilot_patch(self.rows,self.catalog)
        updates={record['id']:record['fields'] for record in patch['records']}
        self.assertNotIn('AnalysisStatus',updates[1])
        self.assertEqual(updates[6]['AnalysisStatus'],'Analyzed')

    def test_apply_rechecks_human_edits_and_catalog_fingerprint_before_any_write(self):
        patch=reviewed_pilot_patch(self.rows,self.catalog)
        client=Client(self.rows)
        self.assertTrue(apply_reviewed_pilot_patch(patch,client)['dry_run']);self.assertFalse(client.writes)
        client.rows[5]['fields']['AnalysisStatus']='In progress'
        with self.assertRaisesRegex(ValueError,'changed after patch review'):
            apply_reviewed_pilot_patch(patch,client,False)
        self.assertFalse(client.writes)
        self.catalog.write_text(self.catalog.read_text()+'\n')
        with self.assertRaisesRegex(ValueError,'Catalog changed'):
            apply_reviewed_pilot_patch(patch,Client(self.rows),False)

    def test_reconcile_sources_only_fills_blank_exact_current_status(self):
        plan={'scope':'pilot','samples':[],'experiments':[]}
        for record in self.rows:
            plan['experiments'].append({'fields':copy.deepcopy(record['fields']),'sample_id':None})
        plan['experiments'][0]['fields']['AnalysisStatus']='Not started'
        reconcile_sources(plan,catalog=self.catalog)
        fields={entry['fields']['Folder']+'/'+entry['fields']['Expno']:entry['fields'] for entry in plan['experiments']}
        self.assertEqual(fields['20250528_15N_hANP32a_50mMNaCl/5']['AnalysisStatus'],'Not started')
        self.assertEqual(fields['20251105_15N_heTDR1_NaCl/405']['AnalysisStatus'],'Analyzed')
        self.assertNotIn('AnalysisStatus',fields['20251218_15N_heDR1_50mMNaCl/5'])
        self.assertEqual(fields['20250528_15N_hANP32a_50mMNaCl/5']['Assessment'],'Do not use')
        other=copy.deepcopy(plan);other['experiments'][5]['fields'].update(SourceID='other',AnalysisStatus='')
        reconcile_sources(other,catalog=self.catalog)
        self.assertEqual(other['experiments'][5]['fields']['AnalysisStatus'],'')
