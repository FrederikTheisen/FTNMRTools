"""Saved-result overviews, metadata-only refresh and rollback tests."""
import argparse
import copy
import csv
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
import zipfile

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import condition_overview as co
import relax_tools as rt


def decay_bundle():
    peaks=[dict(peak_id='FT000001-H',assignment='?-A10N-H',residue=10,coordinates=['120','8'],overlap=False,duplicate_assignment=False)]
    times=[0.,.1,.2,.4]
    import math
    y=[100*math.exp(-2*t) for t in times]
    return dict(peaks=peaks,extracted={'FT000001-H':y},fit_text='FT000001-H 2 0.1 100 3 0 5\n',times=times,noises=[1.]*4)


class OverviewTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.project=self.root/'project';self.project.mkdir()
        self.data=self.root/'data';self.data.mkdir()
        self.catalog=self.root/'catalog.csv'
    def catalog_rows(self, rows):
        fields=list(dict.fromkeys(k for r in rows for k in r))
        with self.catalog.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
    def fixture(self,mode='r1',name='collection/3',date='2025-01-01',salt=50,temp=25,condition_extra=None):
        folder=self.project/'ANP32aIDR'/f'{temp}C'/f'{mode}_{salt}mM_600MHz'/'run'/name;folder.mkdir(parents=True)
        acq=self.data/name;acq.mkdir(parents=True,exist_ok=True);(acq/'acqus').write_text('##$PULPROG= <actual_sequence>\n')
        prefix='sample';numbers=folder/'sample_numbers.csv';archive=folder/'sample_support.zip'
        if mode=='noe':
            fields=['residue','hetNOE','hetNOE_sd','I_ref','I_NOE','I_ref_sd','I_NOE_sd','ref_overlap','NOE_overlap','review_flags','comparison_eligible']
            with numbers.open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerow(dict(zip(fields,[10,-.25,.02,100,-25,1,2,False,True,'overlap',False])))
            (folder/'sample_intensity.pdf').write_bytes(b'old PDF')
            provenance={'definition':'signed ratio','analysis_driver_sha256':'original-fit-code','reference_metadata':{'noise_sd':1}}
            (folder/'provenance.json').write_text(json.dumps(provenance))
            with zipfile.ZipFile(archive,'w') as z:
                z.writestr('inputs/original.list',b'assigned peaks unchanged')
                z.writestr('provenance.json',json.dumps(provenance))
        else:
            b=decay_bundle();c=rt.report_context(mode,'sample','sample_T1.list_intensity.pdf',name,peak_list='sample_T1.list')
            names=rt.prepare_reports(folder,**b,mode=mode,prefix='sample',pdf_name=c['pdf_name'],run_tag=name,context=c)
            (folder/c['pdf_name']).write_bytes(b'old PDF')
            with zipfile.ZipFile(archive,'w') as z:
                for file in names:
                    if file!=numbers.name:z.write(folder/file,file)
                z.writestr('rlx_dfit.inp',b'native input unchanged');z.writestr('rlx_dfit.out',b'native fit unchanged')
                z.writestr('rlx_dfit_noise.txt',b'native noise unchanged')
                z.writestr('data_xtr.out',b'native measurements unchanged')
        spec=dict(relative_path=name)
        job=dict(spec=spec,frozen_sha256='original-fitter-code',driver_sha256='original-driver',output_hashes={str(numbers):'old'})
        (folder/'job.json').write_text(json.dumps(job))
        c={'nominal_temperature_C':temp};c.update(condition_extra or {})
        entry=dict(relative_path=name,protein='ANP32A',construct='IDR',experiment_type=mode,NaCl_mM=str(salt),temperature_C=str(temp),temperature_basis='reviewed nominal',field_MHz='600',conditions=json.dumps(c),status='analysed',output_directory=str(folder),investigator_exclusion='',acquisition_date=date)
        return entry
    def records(self,*entries,**kw):
        self.catalog_rows(entries);return co.discover(self.catalog,self.project,self.data,**kw)
    def review_fixture(self,mode='etaXY',name='extra/6',kind='primary raw acquisition'):
        acq=self.data/name;acq.mkdir(parents=True,exist_ok=True)
        pulse='nico_nnh_tr_noSE' if mode=='etaXY' else 'nico_t1sat_v2s'
        (acq/'acqus').write_text('##$PULPROG= <'+pulse+'>\n')
        return dict(relative_path=name,protein='SeV Ntail',construct='Ntail',experiment_type=mode,
                    pulse_program=pulse,primary_record_status=kind,related_primary_path='',
                    status='blocked_method_inputs',output_directory='',NaCl_mM='50',temperature_C='25',
                    temperature_basis='reviewed nominal',field_MHz='600.48',acquisition_date='2025-06-01',
                    conditions=json.dumps({'nominal_temperature_C':'25'}),investigator_exclusion='')
    def fake_render(self,cmd,cwd,stdout=None):
        self.assertEqual(cmd[0],'gnuplot','No processing/extraction/fitting executable may run')
        script=(Path(cwd)/'plot.gnu').read_text();import re
        name=re.search(r'set output "([^"]+)"',script)[1]
        self.assertIn('actual_sequence',script)
        (Path(cwd)/name).write_bytes(b'%PDF new')
    def test_numbers_metadata_quoting_and_exact_byte_preservation(self):
        data=b'"residue","warnings"\r\n10,"a, quoted\n# warning: embedded"\r\n# source_folder: collection\r\n'
        fields,rows,meta=co.read_numbers_bytes(data)
        self.assertEqual(rows[0]['warnings'],'a, quoted\n# warning: embedded')
        new=co.add_pulse_footer(data,'sequence');self.assertTrue(new.startswith(data))
        self.assertEqual(co.read_numbers_bytes(new)[:2],(fields,rows))
        self.assertEqual(co.add_pulse_footer(new,'sequence'),new)
        with self.assertRaises(ValueError):co.add_pulse_footer(new,'other')
    def test_condition_grouping_construct_additives_and_provisional_temperature(self):
        a=self.fixture();b=copy.deepcopy(a);b['NaCl_mM']='500'
        self.assertNotEqual(co.confirmed_condition(a)[0],co.confirmed_condition(b)[0])
        b=copy.deepcopy(a);b['construct']='FL';self.assertNotEqual(co.confirmed_condition(a)[0],co.confirmed_condition(b)[0])
        b=copy.deepcopy(a);b['conditions']=json.dumps({'nominal_temperature_C':25,'additives':'30bpDNA'})
        self.assertNotEqual(co.confirmed_condition(a)[0],co.confirmed_condition(b)[0])
        a['temperature_basis']='provisional controller group';self.assertIsNone(co.confirmed_condition(a)[0])
        a['temperature_basis']='reviewed';a['NaCl_mM']='unknown';self.assertIsNone(co.confirmed_condition(a)[0])
    def test_groups_across_proteins_fields_observables_and_six_header_values(self):
        original=self.fixture();base_key=co.confirmed_condition(original)[0]
        for key,value in [('protein','SeV Ntail'),('field_MHz','700'),('experiment_type','noe'),('temperature_C','5')]:
            changed=dict(original);changed[key]=value
            self.assertNotEqual(base_key,co.confirmed_condition(changed)[0])
        changed=dict(original,protein='heTDR1',construct='heTDR1')
        row=self.records(changed)[0]
        path=self.root/'hetdr-overview.tsv';co.write_overview(path,[row],[],self.project)
        text=path.read_text();self.assertIn('HeTDR1',text);self.assertIn('actual_sequence',text)
        row['mode']='r1rho';row['protocol']={'spinlock_nominal_Hz':'1500','spinlock_offset_from_carrier_Hz':'0'}
        info=co.header_cells(row,self.project)
        self.assertIn('spinlock_Hz=1500',info[1]);self.assertIn('offset_Hz=0',info[1])

    def test_duplicate_catalog_and_archive_paths_rejected(self):
        a=self.fixture()
        with self.assertRaisesRegex(ValueError,'Duplicate current acquisition'):self.records(a,a)
        a['output_directory']=str(self.project/'workers')
        with self.assertRaisesRegex(ValueError,'archived/diagnostic'):self.records(a)
        with self.assertRaisesRegex(ValueError,'protected'):co.safe_source('211@950/1')
    def test_qc_keeps_review_unassessed_and_excludes_only_reviewed_decisions(self):
        entries=[self.fixture(name=f'collection/{i}') for i in (3,4,5)]
        path=self.root/'fid.csv'
        with path.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=['relative_path','FID_screen_status','notes']);w.writeheader();w.writerows([
                dict(relative_path='collection/3',FID_screen_status='major_signature_or_logged_event',notes='review only'),
                dict(relative_path='collection/4',FID_screen_status='startup_pattern_only',notes='startup review')])
        decisions=self.root/'decisions.json';decisions.write_text(json.dumps({'known_strong_growth_experiments':['collection/4']}))
        records=self.records(*entries,fid_inventory=path,decisions=decisions)
        self.assertEqual([r['qc'][0] for r in records],['review','review','not_assessed'])
        self.assertFalse(records[0]['exclusion']);self.assertTrue(records[1]['exclusion'])
    def test_cells_missing_duplicates_overlap_invalid_and_negative_noe(self):
        p=dict(residue='10',rate_s=2,**{'rate_s-1':'2','rate_error_s-1':'.1','reduced_chi2_input':'1.2','overlap':'True','valid_fit':'True'})
        self.assertEqual(co.peak_cells([p],'r1')[:4],['2','.1','True','1.2'])
        self.assertEqual(co.peak_cells([],'r1'),['','','','','missing'])
        self.assertEqual(co.peak_cells([p,p],'r1')[-1],'multiple_peaks')
        p['duplicate_assignment']='True';self.assertEqual(co.peak_cells([p],'r1')[0],'')
        p.pop('duplicate_assignment');p['valid_fit']='False';self.assertEqual(co.peak_cells([p],'r1')[0],'')
        self.assertIn('invalid_fit',co.peak_cells([p],'r1')[-1])
        self.assertEqual(co.peak_cells([{'hetNOE':'-.25','hetNOE_sd':'.02'}],'noe')[:2],['-.25','.02'])
    def test_overview_has_six_headers_and_sorted_union_with_source_paths(self):
        a,b=self.fixture(name='collection/3'),self.fixture(name='collection/4',date='2024-12-01')
        records=self.records(a,b);records[0]['rows'].append(dict(records[0]['rows'][0],residue='9'))
        dest=self.root/'overview.tsv';co.write_overview(dest,records,[],self.project)
        with dest.open() as f:rows=list(csv.reader(f,delimiter='\t'))
        self.assertEqual([r[0] for r in rows[3:9]],['source','pulse / protocol','delays used','noise / median first-plane height','acquisition QC','fit assessment'])
        self.assertEqual('collection/4',rows[3][1]);self.assertEqual(rows[9][0],'residue')
        self.assertEqual([r[0] for r in rows[10:12]],['9','10'])
        self.assertEqual(len(rows[9]),11);self.assertIn('not_assessed',rows[7][1]);self.assertEqual(rows[9][1:3],['R1','error'])
    def test_decay_refresh_preserves_numbers_native_fits_original_provenance_and_backups(self):
        row=self.records(self.fixture())[0];old=row['numbers'].read_bytes()
        support=row['directory']/'sample_support.zip'
        with zipfile.ZipFile(support) as z:native=co.native_members(z)
        with mock.patch.object(rt,'checked_run',side_effect=self.fake_render):self.assertEqual(co.refresh(row),'refreshed')
        self.assertEqual(co.read_numbers_bytes(row['numbers'].read_bytes())[:2],co.read_numbers_bytes(old)[:2])
        self.assertEqual(row['metadata']['pulse_program'],'actual_sequence')
        with zipfile.ZipFile(support) as z:
            self.assertEqual(co.native_members(z),native)
            self.assertEqual(json.loads(z.read('report_refresh.json'))['pulse_program'],'actual_sequence')
        job=co.load_json(row['directory']/'job.json');self.assertEqual(job['frozen_sha256'],'original-fitter-code')
        self.assertFalse(job['report_refresh']['refitted'])
        self.assertTrue(list((row['directory']/'.report-backups').rglob('sample_numbers.csv')))
        again=self.records(self.fixture_for_existing(row))[0]
        with mock.patch.object(rt,'checked_run') as cmd:self.assertEqual(co.refresh(again),'already_current');cmd.assert_not_called()
    def fixture_for_existing(self,row):
        return {k:v for k,v in row.items() if k in ['relative_path','protein','construct','experiment_type','NaCl_mM','temperature_C','temperature_basis','field_MHz','conditions','status','output_directory','investigator_exclusion','acquisition_date']}
    def test_custom_noe_refresh_preserves_numeric_bytes_inputs_and_fit_origin(self):
        row=self.records(self.fixture(mode='noe'))[0];old=row['numbers'].read_bytes()
        with mock.patch.object(rt,'checked_run',side_effect=self.fake_render):co.refresh(row)
        self.assertTrue(row['numbers'].read_bytes().startswith(old))
        with zipfile.ZipFile(row['directory']/'sample_support.zip') as z:
            self.assertEqual(z.read('inputs/original.list'),b'assigned peaks unchanged')
            self.assertEqual(json.loads(z.read('provenance.json'))['analysis_driver_sha256'],'original-fit-code')
        self.assertIn('pulse_program',co.load_json(row['directory']/'provenance.json'))
    def test_render_failure_leaves_every_existing_report_unchanged(self):
        row=self.records(self.fixture())[0];before={p.name:p.read_bytes() for p in row['directory'].iterdir() if p.is_file()}
        with mock.patch.object(rt,'checked_run',side_effect=subprocess.CalledProcessError(9,'gnuplot')):
            with self.assertRaises(subprocess.CalledProcessError):co.refresh(row)
        for name,data in before.items():self.assertEqual((row['directory']/name).read_bytes(),data)
    def test_publication_failure_rolls_back_complete_acquisition(self):
        row=self.records(self.fixture())[0];folder=row['directory']
        before={p.name:p.read_bytes() for p in folder.iterdir() if p.is_file()}
        rename=Path.rename
        def broken(path,target):
            if path.name=='job.json' and path.parent!=folder and path.parent.name!='.previous' and Path(target)==folder/'job.json':raise OSError('simulated publication failure')
            return rename(path,target)
        with mock.patch.object(rt,'checked_run',side_effect=self.fake_render),mock.patch.object(Path,'rename',new=broken):
            with self.assertRaises(OSError):co.refresh(row)
        for name,data in before.items():self.assertEqual((folder/name).read_bytes(),data)
    def test_source_acquisition_change_aborts_before_publication(self):
        row=self.records(self.fixture())[0];before=row['numbers'].read_bytes()
        def changed(cmd,cwd,stdout=None):
            self.fake_render(cmd,cwd,stdout);row['acqus_path'].write_text('##$PULPROG= <other>\n')
        with mock.patch.object(rt,'checked_run',side_effect=changed),self.assertRaisesRegex(ValueError,'input changed'):co.refresh(row)
        self.assertEqual(row['numbers'].read_bytes(),before)
    def test_plain_parameter_headers_contiguous_residues_and_short_source(self):
        row=self.records(self.fixture())[0]
        row['rows'].append(dict(row['rows'][0],residue='13'))
        path=self.root/'overview.tsv';co.write_overview(path,[row],[],self.project)
        with path.open() as f:rows=list(csv.reader(f,delimiter='\t'))
        self.assertEqual(rows[0][0],'# generated');self.assertIn('DO NOT EDIT',rows[0][1])
        source=next(r for r in rows if r[0]=='source');self.assertEqual(source[1],'collection/3')
        start=next(i for i,r in enumerate(rows) if r[0]=='residue')
        self.assertEqual(rows[start],['residue','R1','error','overlap','reduced_chi2','flags'])
        body=[r for r in rows[start+1:] if r[0].isdigit()]
        self.assertEqual([r[0] for r in body],['10','11','12','13'])
        self.assertEqual(body[1],['11','','','','','missing'])
        self.assertTrue(any(r[0]=='# numbers_file' and r[1]=='collection/3' for r in rows))
    def test_each_source_labels_its_own_parameter_block_after_date_ordering(self):
        a=self.fixture(name='collection/3',date='2025-02-01')
        b=self.fixture(name='collection/4',date='2025-01-01')
        records=self.records(a,b)
        for row in records:
            if row['relative_path']=='collection/4':row['rows'][0]['rate_s-1']='9.123456789012345'
        path=self.root/'overview.tsv';co.write_overview(path,list(reversed(records)),[],self.project)
        with path.open() as f:lines=list(csv.reader(f,delimiter='\t'))
        source=next(r for r in lines if r[0]=='source')
        values=next(r for r in lines if r[0]=='10')
        self.assertEqual((source[1],values[1]),('collection/4','9.123456789012345'))
        self.assertEqual((source[6],values[6]),('collection/3','2.0'))
        self.assertEqual(len([r for r in lines if r[0]=='# numbers_file']),2)

    def test_delay_labels_only_state_what_was_used(self):
        row=self.records(self.fixture())[0]
        row['metadata'].update(fit_basis='all_points',selection_reason='shortest_residual_within_limit')
        self.assertEqual(co.header_cells(row,self.project)[2],'all delays used')
        row['metadata'].update(fit_basis='exclude_shortest',excluded_shortest_delay_s='0.001')
        self.assertEqual(co.header_cells(row,self.project)[2],'1 ms delay omitted')
    def test_relative_noise_uses_unique_valid_nonoverlap_absolute_first_heights(self):
        row=self.records(self.fixture())[0]
        row['rows']=[dict(residue='10',plane_01_intensity='100',valid_fit='True',overlap='False'),
                     dict(residue='11',plane_01_intensity='-300',valid_fit='True',overlap='False'),
                     dict(residue='12',plane_01_intensity='100000',valid_fit='True',overlap='True'),
                     dict(residue='13',plane_01_intensity='1',valid_fit='False',overlap='False')]
        row['metadata'].update(plane_noise_sd='2, 4',fit_basis='exclude_shortest')
        text=co.relative_noise_header(row)
        self.assertEqual(text,'1%, 2%; median |first-plane height|=200 (n=2)')
        row['rows']=[dict(residue='10',plane_01_intensity='0',valid_fit='True')]
        self.assertIn('unavailable',co.relative_noise_header(row))
        noe=self.records(self.fixture(mode='noe',name='collection/4'))[0]
        noe['rows'][0]['NOE_overlap']='False'
        self.assertIn('ref=1%, sat=2%',co.relative_noise_header(noe))
    def test_error_scaling_label_is_explained_without_changing_numbers(self):
        peak={'uncertainty_rescaled':'True','flags':'original_noise_lack_of_fit;uncertainty_rescaled'}
        self.assertEqual(co.row_flags(peak),['error_adjusted_by_fitter','lack_of_fit'])
    def test_user_reviews_are_persistent_and_apply_condition_then_acquisition(self):
        a,b=self.fixture(name='collection/3'),self.fixture(name='collection/4')
        records=self.records(a,b);path=self.project/'condition-reviews.tsv'
        co.initialise_reviews(path,records)
        self.assertEqual(len(co.read_tsv(path)),2)
        condition=str(Path(*records[0]['directory'].relative_to(self.project).parts[:3]))
        with path.open('w',newline='') as f:
            w=csv.writer(f,delimiter='\t');w.writerow(['target','assessment','comment'])
            w.writerow([condition,'review','condition assessment'])
            w.writerow(['collection/3','usable','my inspection: acceptable'])
            w.writerow(['collection/4','exclude','damaged FID'])
        before=path.read_bytes();co.read_user_reviews(path,records,self.project);co.initialise_reviews(path,records)
        self.assertEqual(path.read_bytes(),before)
        self.assertEqual(records[0]['user_assessment'],'usable')
        self.assertIn('my inspection',records[0]['qc'][1]);self.assertIn('condition assessment',records[0]['qc'][1])
        self.assertTrue(records[1]['exclusion']);self.assertEqual(records[1]['qc'][0],'user_exclude')
        path.write_text('target\tassessment\tcomment\ncollection/3\tBAD\tMy rating\n')
        co.read_user_reviews(path,records,self.project)
        self.assertEqual(records[0]['user_assessment'],'BAD')
        self.assertFalse(records[0]['exclusion'])
    def test_new_review_sheet_prefills_machine_context_and_keeps_user_cells_blank(self):
        records=self.records(self.fixture())
        path=self.project/'condition-reviews.tsv';co.initialise_reviews(path,records)
        row=co.read_tsv(path)[0]
        self.assertEqual(row['experiment_type'],'R1')
        self.assertEqual(row['pulse_program'],'actual_sequence')
        self.assertEqual(row['NaCl_mM'],'50')
        self.assertEqual(row['nominal_temperature_C'],'25')
        self.assertEqual(row['field'],'600 MHz')
        self.assertEqual(row['machine_disposition'],'included')
        self.assertEqual(row['machine_assessment'],'not_assessed')
        self.assertNotIn('delays_used',row)
        self.assertNotIn('fit_assessment',row)
        self.assertEqual((row['assessment'],row['comment']),('',''))
        self.assertEqual(row['machine_reason'],'')
    def test_brief_machine_reasons_keep_counts_and_unassessed_rows_blank(self):
        row=self.records(self.fixture())[0]
        row['qc']=('review', "retained-range trace flags: {'direct_startup_discontinuity_review': 8, 'water_frequency_change_review': 3}")
        self.assertEqual(co.brief_machine_reason(row),'water-region frequency outlier (3 traces)')
        row['pending']='nominal temperature provisional'
        row['exclusion']='reviewed strong late-FID growth: findings.json'
        row['qc']=('review','Strong water-dominated late-FID growth at 150,200 ms; 14 direct-startup review traces; ADC overflow recorded')
        self.assertEqual(co.brief_machine_reason(row),'late-FID growth (150,200 ms); condition pending; ADC overflow')
        row['qc']=('not_assessed','no saved acquisition-level FID screen')
        self.assertEqual(co.brief_machine_reason(row),'')
        self.assertEqual(co.machine_review_context(row)['machine_disposition'],'excluded')
    def test_field_display_rounding_and_legacy_column_migration(self):
        records=self.records(self.fixture());row=records[0]
        for value,expected in [('599.8','600 MHz'),('950.23','950 MHz'),('unknown','unknown')]:
            row['field_MHz']=value
            self.assertEqual(co.machine_review_context(row)['field'],expected)
        row['field_MHz']='700.13'
        path=self.project/'review.tsv'
        path.write_text('target\tfield_MHz\tassessment\tcomment\ncollection/3\t599.8\treview\tmy note\n')
        co.initialise_reviews(path,records,update_context=True)
        result=co.read_tsv(path)[0]
        self.assertEqual(result['field'],'700 MHz')
        self.assertNotIn('field_MHz',result)
        self.assertEqual((result['assessment'],result['comment']),('review','my note'))
    def test_lookup_rows_precede_header_and_quoted_hash_comments_survive(self):
        records=self.records(self.fixture());path=self.project/'review.tsv'
        co.initialise_reviews(path,records)
        with path.open(newline='') as handle:
            physical=list(csv.reader(handle,delimiter='\t'))
        preamble=physical[:physical.index(co.REVIEW_FIELDS)]
        definitions=[r for r in preamble if r]
        self.assertEqual(len(definitions),len(co.REVIEW_LOOKUP))
        self.assertTrue(all(len(r)==2 and r[0].startswith('# ') and '\n' not in r[1] for r in definitions))
        self.assertIn(['# later FID discontinuity',dict(co.REVIEW_LOOKUP)['later FID discontinuity']],definitions)
        rows=co.read_tsv(path);rows[0]['comment']='one line\n# this is my comment, not a lookup row'
        with path.open('w',newline='') as handle:
            w=csv.writer(handle,delimiter='\t');w.writerows(preamble)
            w.writerow(co.REVIEW_FIELDS);w.writerows([[row.get(k,'') for k in co.REVIEW_FIELDS] for row in rows])
        before=path.read_bytes()
        co.read_user_reviews(path,records,self.project)
        self.assertIn('# this is my comment',records[0]['user_comment'])
        co.initialise_reviews(path,records)
        self.assertEqual(path.read_bytes(),before)
        co.initialise_reviews(path,records,update_context=True)
        self.assertEqual(co.read_tsv(path)[0]['comment'],rows[0]['comment'])
    def test_unknown_qc_notes_are_reviewed_and_no_major_flags_need_no_reason(self):
        row=self.records(self.fixture())[0]
        row['qc']=('review','saved screening error needing inspection')
        self.assertEqual(co.brief_machine_reason(row),'QC review')
        row['qc']=('screened_no_major_flags','retained-range trace flags: none')
        self.assertEqual(co.brief_machine_reason(row),'')
    def test_review_context_needs_no_peak_results_or_fit_headers(self):
        row=self.records(self.fixture())[0]
        for key in ('rows','metadata','directory'):
            row.pop(key)
        with mock.patch.object(co,'header_cells',side_effect=AssertionError('Review must not prepare fit summaries')):
            context=co.machine_review_context(row)
        self.assertNotIn('fit_assessment',context)
        self.assertNotIn('delays_used',context)
        self.assertEqual(context['experiment_type'],'R1')
    def test_review_refresh_drops_old_analysis_columns_preserving_user_notes(self):
        records=self.records(self.fixture());path=self.project/'review.tsv'
        path.write_text('target\tdelays_used\tfit_assessment\tassessment\tcomment\tcustom\ncollection/3\tall delays used\t144 valid peaks; Holm=10/20\treview\tkeep my assessment\tkeep custom note\n')
        co.initialise_reviews(path,records,update_context=True)
        row=co.read_tsv(path)[0]
        self.assertNotIn('delays_used',row);self.assertNotIn('fit_assessment',row)
        self.assertEqual((row['assessment'],row['comment'],row['custom']),('review','keep my assessment','keep custom note'))
        text=path.read_text()
        self.assertNotIn('Holm',text);self.assertNotIn('valid peaks',text)
        self.assertNotIn('departures_detected',text)
    def test_screen_reason_labels_distinguish_fid_time_and_increment_features(self):
        row=self.records(self.fixture())[0]
        row['qc']=('review',"retained-range trace flags: {'isolated_direct_discontinuity_review': 2, 'isolated_water_amplitude_change_review': 1}")
        self.assertEqual(co.brief_machine_reason(row),'later FID discontinuity (2 traces); water-region amplitude outlier (1 trace)')
        lookup=dict(co.REVIEW_LOOKUP)
        self.assertNotIn('FID start discontinuity',lookup)
        self.assertIn('after 0.5 ms',lookup['later FID discontinuity'])
        self.assertIn('across 15N increments',lookup['water-region frequency outlier'])
    def test_startup_only_is_omitted_without_residual_review_or_clean_claim(self):
        row=self.records(self.fixture())[0]
        for note in ("retained-range trace flags: {'direct_startup_discontinuity_review': 20}",
                     'startup_pattern_only; 20 direct-startup review traces (maximum jump <=0.5 ms)'):
            row['qc']=('review',note)
            before=row['qc']
            context=co.machine_review_context(row)
            self.assertEqual(context['machine_assessment'],'screened_no_retained_flags')
            self.assertEqual(context['machine_reason'],'')
            self.assertEqual(row['qc'],before)
            self.assertEqual(context['machine_disposition'],'included')
    def test_omitting_startup_preserves_mixed_flags_and_reviewed_exclusions(self):
        row=self.records(self.fixture())[0]
        row['qc']=('review',"retained-range trace flags: {'direct_startup_discontinuity_review': 20, 'water_frequency_change_review': 2}")
        self.assertEqual(co.review_machine_assessment(row),'review')
        self.assertEqual(co.brief_machine_reason(row),'water-region frequency outlier (2 traces)')
        row['qc']=('review','startup_pattern_only; 20 direct-startup review traces')
        row['exclusion']='reviewed strong late-FID growth: findings.json'
        context=co.machine_review_context(row)
        self.assertEqual(context['machine_disposition'],'excluded')
        self.assertEqual(context['machine_assessment'],'review')
        self.assertEqual(context['machine_reason'],'late-FID growth')
    def test_unknown_or_malformed_flags_remain_review(self):
        row=self.records(self.fixture())[0]
        for note in ("retained-range trace flags: {'direct_startup_discontinuity_review': 20, 'unrecognised_flag': 2}",
                     'retained-range trace flags: malformed',
                     "retained-range trace flags: {'direct_startup_discontinuity_review': 'unknown'}"):
            row['qc']=('review',note)
            self.assertEqual(co.review_machine_assessment(row),'review')
            self.assertEqual(co.brief_machine_reason(row),'QC review')
    def test_explicit_context_update_preserves_user_values_and_custom_condition_notes(self):
        records=self.records(self.fixture());path=self.project/'condition-reviews.tsv'
        condition=str(Path(*records[0]['directory'].relative_to(self.project).parts[:3]))
        with path.open('w',newline='') as f:
            w=csv.writer(f,delimiter='\t');w.writerow(['target','assessment','comment','custom_note'])
            w.writerow(['collection/3','review','my "quoted" comment, with\nsecond line','keep this'])
            w.writerow([condition,'','condition-wide note','also keep'])
        before=path.read_bytes();co.read_user_reviews(path,records,self.project)
        co.initialise_reviews(path,records,update_context=True)
        rows=co.read_tsv(path)
        self.assertEqual(rows[0]['assessment'],'review')
        self.assertEqual(rows[0]['comment'],'my "quoted" comment, with\nsecond line')
        self.assertEqual(rows[0]['custom_note'],'keep this')
        self.assertEqual(rows[1]['comment'],'condition-wide note')
        self.assertEqual(rows[1]['custom_note'],'also keep')
        self.assertEqual(rows[0]['machine_assessment'],'not_assessed')
        backups=list((self.project/'.report-backups').rglob('condition-reviews.tsv'))
        self.assertEqual(len(backups),1);self.assertEqual(backups[0].read_bytes(),before)
        current=path.read_bytes();co.initialise_reviews(path,records)
        self.assertEqual(path.read_bytes(),current)
    def test_machine_context_never_becomes_a_user_override(self):
        records=self.records(self.fixture());records[0]['exclusion']='reviewed damage'
        path=self.project/'condition-reviews.tsv';co.initialise_reviews(path,records)
        rows=co.read_tsv(path);rows[0]['machine_disposition']='included'
        with path.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=co.REVIEW_FIELDS,delimiter='\t');w.writeheader();w.writerows(rows)
        co.read_user_reviews(path,records,self.project)
        self.assertEqual(records[0]['exclusion'],'reviewed damage')
        self.assertEqual(records[0]['user_assessment'],'')
    def test_context_snapshot_is_automated_even_after_user_acceptance(self):
        records=self.records(self.fixture());records[0]['exclusion']='reviewed damage'
        path=self.project/'condition-reviews.tsv'
        path.write_text('target\tassessment\tcomment\ncollection/3\tusable\tI reviewed it.\n')
        co.read_user_reviews(path,records,self.project)
        self.assertEqual(records[0]['exclusion'],'')
        co.initialise_reviews(path,records,update_context=True)
        row=co.read_tsv(path)[0]
        self.assertEqual(row['machine_disposition'],'excluded')
        self.assertEqual(row['machine_assessment'],'not_assessed')
        self.assertEqual(row['assessment'],'usable')
        self.assertEqual(row['comment'],'I reviewed it.')
    def test_review_reader_rejects_duplicate_headers_but_accepts_old_and_enriched_tables(self):
        path=self.project/'review.tsv'
        path.write_text('target\tassessment\tcomment\ncollection/3\t\told note\n')
        self.assertEqual(co.read_tsv(path)[0]['comment'],'old note')
        path.write_text('target\tassessment\tcomment\tcomment\ncollection/3\t\tnote\tnote\n')
        with self.assertRaises(ValueError):co.read_tsv(path)

    def test_generator_never_changes_existing_review_comments(self):
        a=self.fixture();self.catalog_rows([a]);path=self.project/'condition-reviews.tsv'
        path.write_text('target\tassessment\tcomment\ncollection/3\treview\tPreserve my comment, exactly.\n')
        before=path.read_bytes()
        args=argparse.Namespace(project_root=str(self.project),data_root=str(self.data),catalog=str(self.catalog),fid_inventory=None,r1_qc=None,qc_decisions=None,user_reviews=str(path),dry_run=False,refresh_reports=False)
        co.run(args);self.assertEqual(path.read_bytes(),before)
        overview=co.destination_for(self.records(a),self.project)
        self.assertIn('Preserve my comment, exactly.',overview.read_text())

    def test_dry_run_has_no_writes_and_full_run_publishes_index(self):
        a=self.fixture();self.catalog_rows([a]);args=argparse.Namespace(project_root=str(self.project),data_root=str(self.data),catalog=str(self.catalog),fid_inventory=None,r1_qc=None,qc_decisions=None,dry_run=True,refresh_reports=False)
        with mock.patch.object(co,'refresh') as refresh:counts=co.run(args);refresh.assert_not_called()
        self.assertEqual(counts['included'],1);self.assertFalse((self.project/'.relax.lock').exists())
        args.dry_run=False;co.run(args)
        path=co.destination_for(self.records(a),self.project)
        self.assertTrue(path.is_file());self.assertTrue((self.project/'condition-overviews.tsv').is_file())

    def test_extra_observables_need_no_numbers_and_show_pulse_nominal_conditions_and_qc(self):
        eta=self.review_fixture();water=self.review_fixture('water 1H T1','extra/7')
        water['protein']='heTDR1';water['NaCl_mM']=''
        eta['temperature_C']='15';eta['temperature_basis']='provisional controller group'
        self.catalog_rows([eta,water])
        records=co.discover_reviews(self.catalog,self.project,self.data)
        self.assertEqual(len(records),2)
        contexts={r['relative_path']:co.machine_review_context(r) for r in records}
        self.assertEqual(contexts['extra/6']['experiment_type'],'etaXY')
        self.assertEqual(contexts['extra/6']['pulse_program'],'nico_nnh_tr_noSE')
        self.assertEqual(contexts['extra/6']['nominal_temperature_C'],'25')
        self.assertEqual(contexts['extra/6']['machine_disposition'],'condition_review')
        self.assertEqual(contexts['extra/7']['experiment_type'],'water T1sat')
        self.assertEqual(contexts['extra/7']['pulse_program'],'nico_t1sat_v2s')
        self.assertTrue(all(c['machine_assessment']=='not_assessed' and c['machine_reason']=='' for c in contexts.values()))
        self.assertFalse(any(self.project.rglob('*_numbers.csv')))
    def test_extra_discovery_excludes_components_setups_and_duplicates_and_marks_processed_only(self):
        entries=[self.review_fixture(name='extra/6'),
                 self.review_fixture(name='extra/6001',kind='inferred split/derived product'),
                 self.review_fixture(name='extra/9',kind='parameters only'),
                 self.review_fixture(name='extra/17',kind='processed only; raw absent'),
                 self.review_fixture(name='extra/18')]
        entries[1]['related_primary_path']='extra/6';entries[-1]['status']='excluded_duplicate_raw'
        self.catalog_rows(entries)
        records=co.discover_reviews(self.catalog,self.project,self.data)
        self.assertEqual({r['relative_path'] for r in records},{'extra/6','extra/17'})
        only=next(r for r in records if r['relative_path']=='extra/17')
        self.assertEqual(co.machine_review_context(only)['machine_disposition'],'processed_only_review')
    def test_reviews_only_refresh_preserves_notes_and_never_reads_results_or_runs_commands(self):
        eta=self.review_fixture();self.catalog_rows([eta])
        path=self.project/'condition-reviews.tsv'
        path.write_text('target\tassessment\tcomment\nextra/6\treview\tMy etaXY assessment\n')
        args=argparse.Namespace(project_root=str(self.project),data_root=str(self.data),catalog=str(self.catalog),
                               fid_inventory=None,r1_qc=None,qc_decisions=None,reviews_only=True,
                               dry_run=True,refresh_reports=False)
        before=path.read_bytes()
        with mock.patch.object(co,'discover',side_effect=AssertionError('Must not read numerical reports')), \
             mock.patch.object(co,'refresh',side_effect=AssertionError('Must not refresh results')), \
             mock.patch.object(co.subprocess,'run',side_effect=AssertionError('No external commands')):
            co.run(args);self.assertEqual(path.read_bytes(),before)
            self.assertFalse((self.project/'.relax.lock').exists())
            args.dry_run=False;co.run(args)
        result=co.read_tsv(path)[0]
        self.assertEqual((result['assessment'],result['comment']),('review','My etaXY assessment'))
        self.assertEqual(path.read_bytes(),before)
        self.assertFalse((self.project/'condition-overviews.tsv').exists())
    def test_standard_overview_accepts_extra_review_targets_without_creating_extra_overviews(self):
        a=self.fixture();eta=self.review_fixture();water=self.review_fixture('water 1H T1','extra/7')
        self.catalog_rows([a,eta,water])
        path=self.project/'condition-reviews.tsv'
        path.write_text('target\tassessment\tcomment\nextra/6\treview\tKeep this note\n')
        args=argparse.Namespace(project_root=str(self.project),data_root=str(self.data),catalog=str(self.catalog),
                               fid_inventory=None,r1_qc=None,qc_decisions=None,update_review_context=True,
                               dry_run=False,refresh_reports=False)
        counts=co.run(args)
        self.assertEqual(counts['acquisitions'],1);self.assertEqual(counts['overviews'],1)
        records=co.read_tsv(path)
        self.assertEqual(len(records),3)
        self.assertEqual(next(r for r in records if r['target']=='extra/6')['comment'],'Keep this note')
        self.assertEqual(len(list(self.project.rglob('*_overview.tsv'))),1)
    def test_extra_duplicate_conflicting_pulse_and_protected_path_fail_before_update(self):
        eta=self.review_fixture();self.catalog_rows([eta,eta])
        with self.assertRaisesRegex(ValueError,'Duplicate assessment'):
            co.discover_reviews(self.catalog,self.project,self.data)
        eta['pulse_program']='wrong';self.catalog_rows([eta])
        with self.assertRaisesRegex(ValueError,'Pulse-program conflict'):
            co.discover_reviews(self.catalog,self.project,self.data)
        eta['relative_path']='211@950/6';self.catalog_rows([eta])
        with self.assertRaisesRegex(ValueError,'protected'):
            co.discover_reviews(self.catalog,self.project,self.data)

    def test_all_primary_relaxation_types_are_inventoried_without_completed_analysis(self):
        entries=[]
        for i,(mode,pulse,status) in enumerate([
            ('15N T1','15N_T1_HSQC','blocked_inputs'),
            ('15N T1rho','15N_T1rho','blocked_inputs'),
            ('1H-15N hetNOE','15N_HETNOE','blocked_noe_inputs'),
            ('15N CPMG','hsqcrexetf3gpsitc3d','blocked_method_inputs'),
            ('15N T2','15N_T2_HSQC','pending_readiness'),
            ('unresolved relaxation candidate','zg30','blocked_processed_only_review')],3):
            entry=self.review_fixture(name=f'all/{i}')
            entry.update(experiment_type=mode,pulse_program=pulse,status=status)
            if mode=='unresolved relaxation candidate':entry['primary_record_status']='processed only; raw absent'
            (self.data/entry['relative_path']/'acqus').write_text('##$PULPROG= <'+pulse+'>\n')
            entries.append(entry)
        self.catalog_rows(entries)
        records=co.discover_reviews(self.catalog,self.project,self.data)
        self.assertEqual(len(records),6)
        self.assertEqual({co.machine_review_context(r)['experiment_type'] for r in records},
                         {'R1','R1rho','hetNOE','15N CPMG','R2','unresolved relaxation'})
        self.assertTrue(all(r['qc'][0]=='not_assessed' for r in records))
    def test_append_preserves_edited_document_bytes_custom_columns_and_free_ratings(self):
        a=self.review_fixture(name='extra/6');b=self.review_fixture('water 1H T1','extra/7')
        self.catalog_rows([a,b]);records=co.discover_reviews(self.catalog,self.project,self.data)
        path=self.project/'review.tsv'
        before=(b'# My edited explanation\tkeep this\t\t\r\n\t\t\t\r\n'
                b'target\tassessment\tcomment\tcustom\r\n'
                b'extra/6\tOK\t"My comment\r\n# inside the comment"\tmy custom value\r\n'
                b'\t\t\t\r\n')
        path.write_bytes(before)
        co.initialise_reviews(path,records,append_missing=True)
        after=path.read_bytes();self.assertTrue(after.startswith(before))
        rows=co.read_tsv(path)
        self.assertEqual(len(rows),2)
        self.assertEqual(rows[0],{'target':'extra/6','assessment':'OK','comment':'My comment\r\n# inside the comment','custom':'my custom value'})
        self.assertEqual(rows[1],{'target':'extra/7','assessment':'','comment':'','custom':''})
        co.initialise_reviews(path,records,append_missing=True)
        self.assertEqual(path.read_bytes(),after)
        backups=list((self.project/'.report-backups').rglob('review.tsv'))
        self.assertEqual(len(backups),1);self.assertEqual(backups[0].read_bytes(),before)
    def test_append_failures_and_concurrent_edit_leave_the_edited_document_untouched(self):
        entry=self.review_fixture();self.catalog_rows([entry])
        records=co.discover_reviews(self.catalog,self.project,self.data)
        path=self.project/'review.tsv';before=b'target\tassessment\tcomment\n'
        path.write_bytes(before)
        with mock.patch.object(co.rt,'publish',side_effect=OSError('publication failed')):
            with self.assertRaises(OSError):co.initialise_reviews(path,records,append_missing=True)
        self.assertEqual(path.read_bytes(),before)
        original=co.read_review_table
        def edit_after_read(p):
            result=original(p);p.write_bytes(before+b'extra/6\tBAD\tFresh edit\n');return result
        with mock.patch.object(co,'read_review_table',side_effect=edit_after_read):
            with self.assertRaisesRegex(ValueError,'changed during append'):
                co.initialise_reviews(path,records,append_missing=True)
        self.assertEqual(path.read_bytes(),before+b'extra/6\tBAD\tFresh edit\n')
    def test_live_scan_adds_new_r1_and_r2_but_omits_split_and_known_copy_records(self):
        existing=self.review_fixture();self.catalog_rows([existing])
        def acquired(name,pulse,audit=''):
            folder=self.data/name;folder.mkdir(parents=True)
            (folder/'acqus').write_text('##$PULPROG= <'+pulse+'>\n##$BF1= 850.13\n')
            (folder/'ser').write_bytes(b'raw file must never be processed')
            if audit:(folder/'audita.txt').write_text(audit)
        acquired('new/5','hsqct1etf3gpsitc3d.3')
        acquired('new/6','15N_T2_HSQC')
        acquired('new/10005','hsqct1etf3gpsitc3d.3','created by split 7: dataset 1 of 7 from EXPNO 5')
        acquired('copy/5','hsqct1etf3gpsitc3d.3')
        with self.catalog.with_name('source_record_dispositions.csv').open('w',newline='') as h:
            w=csv.DictWriter(h,fieldnames=['relative_path','experiment_type','primary_record_status','status','related_primary_path'])
            w.writeheader();w.writerow(dict(relative_path='copy/5',experiment_type='15N T1',primary_record_status='byte-identical raw copy',status='excluded_nonprimary',related_primary_path='new/5'))
        with mock.patch.object(co.subprocess,'run',side_effect=AssertionError('No processing/extraction/fitting')):
            records=co.discover_reviews(self.catalog,self.project,self.data)
        targets={r['relative_path']:r for r in records}
        self.assertEqual(set(targets),{'extra/6','new/5','new/6'})
        self.assertEqual(targets['new/5']['mode'],'r1');self.assertEqual(targets['new/6']['mode'],'r2')
        self.assertEqual(co.machine_review_context(targets['new/5'])['field'],'850 MHz')
        self.assertEqual(targets['new/5']['qc'][0],'not_assessed')
        self.assertTrue(targets['new/5']['pending'])


if __name__=='__main__':unittest.main()
