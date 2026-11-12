import copy,itertools,re,unittest
from nmr_register.schema import (QC_CRITERIA,QC_DEFAULT,assessment_from_values,
    ASSESSMENT_FORMULA,ASSESSMENT_EXPLANATION_FORMULA,PRE_TOLERANCE_ASSESSMENT_FORMULA,PRE_TOLERANCE_ASSESSMENT_EXPLANATION_FORMULA,ATTENTION_FORMULA,LEGACY_QC_ATTENTION_FORMULA,columns)
from nmr_register.setup import upgrade_assessment_schema

def formula(source,values):
    import re
    source=re.sub(r'\$([A-Za-z_][A-Za-z_0-9]*)',r"values.get('\1')",source)
    ns={};exec('def calc(values):\n'+''.join('    '+line+'\n' for line in source.splitlines()),ns)
    return ns['calc'](values)
class FormulaTests(unittest.TestCase):
    def test_all_valid_combinations(self):
        names=list(QC_CRITERIA)
        for choices in itertools.product(*[(QC_DEFAULT,*QC_CRITERIA[n][1]) for n in names]):
            values=dict(zip(names,choices));expected=assessment_from_values(values)
            self.assertEqual(formula(ASSESSMENT_FORMULA,values),expected,values)
            self.assertTrue(formula(ASSESSMENT_EXPLANATION_FORMULA,values))
    def test_caution_threshold_and_every_worst_veto(self):
        good={k:v[1][0] for k,v in QC_CRITERIA.items()}
        good['QCAcquisitionCoverage']='Sufficient'
        cautions=[k for k in QC_CRITERIA if k!='QCAcquisitionCoverage']
        for count in range(6):
            values=dict(good)
            for key in cautions[:count]:values[key]=QC_CRITERIA[key][1][1]
            expected='Limited use' if count>=2 else 'Use'
            self.assertEqual(assessment_from_values(values),expected)
            self.assertEqual(formula(ASSESSMENT_FORMULA,values),expected)
            explanation=formula(ASSESSMENT_EXPLANATION_FORMULA,values)
            if count==1:self.assertIn('Single caution accepted:',explanation)
            for key in cautions[:count]:self.assertIn(QC_CRITERIA[key][0],explanation)
        for key,(_,choices) in QC_CRITERIA.items():
            values={key:choices[-1]}
            self.assertEqual(assessment_from_values(values),'Do not use')
            self.assertEqual(formula(ASSESSMENT_FORMULA,values),'Do not use')
    def test_imported_fallback_unchanged(self):
        for decision in ('Use','Limited use','Do not use'):
            for blank in (None,'',QC_DEFAULT):
                values={k:blank for k in QC_CRITERIA};values['ImportedAssessment']=decision
                self.assertEqual(assessment_from_values(values),decision)
                self.assertEqual(formula(ASSESSMENT_FORMULA,values),decision)
                self.assertEqual(formula(ASSESSMENT_EXPLANATION_FORMULA,values),'Imported assessment: '+decision+'; checklist not yet completed')
            values['QCSettings']='Unknown'
            self.assertEqual(assessment_from_values(values),'')
            self.assertEqual(formula(ASSESSMENT_FORMULA,values),'')
    def test_failure_precedes_missing_and_unknown_does_not_pass(self):
        self.assertEqual(formula(ASSESSMENT_FORMULA,{'QCFID':'Unusable'}),'Do not use')
        self.assertIn('FIDs and spectrum: Unusable',formula(ASSESSMENT_EXPLANATION_FORMULA,{'QCFID':'Unusable'}))
        good={k:v[1][0] for k,v in QC_CRITERIA.items()};good['QCFID']='Marginal'
        self.assertEqual(formula(ASSESSMENT_FORMULA,good),'')
        good['QCFID']='Usable with exclusions';self.assertEqual(formula(ASSESSMENT_FORMULA,good),'Use')
        good['QCAcquisitionCoverage']='Sufficient';good['QCFID']='Clean';self.assertEqual(formula(ASSESSMENT_FORMULA,good),'Use')
        self.assertIn('$AssessmentExplanation',ATTENTION_FORMULA)

class Fake:
    def __init__(self,records):
        self.rows=records;self.actions=[]
        excluded=set(QC_CRITERIA)|{'AssessmentExplanation'}
        self.cols=[{'id':i+1,'fields':{'parentId':4,'colId':c['id'],**{k:v for k,v in c.items() if k!='id'}}} for i,c in enumerate(columns()['Experiments']) if c['id'] not in excluded]
        for c in self.cols:
            if c['fields']['colId']=='Assessment':c['fields'].update(isFormula=False,formula='')
            if c['fields']['colId']=='AttentionReason':c['fields']['formula']=LEGACY_QC_ATTENTION_FORMULA
    def records(self,name):return {'_grist_Tables':[{'id':4,'fields':{'tableId':'Experiments'}}], '_grist_Tables_column':self.cols,'Experiments':self.rows}[name]
    def apply(self,actions):
        self.actions.append(actions)
        for a in actions:
            if a[0]=='AddColumn':self.cols.append({'id':1000+len(self.cols),'fields':{'parentId':4,'colId':a[2],**a[3]}})
            elif a[0]=='ModifyColumn':next(c for c in self.cols if c['fields']['colId']==a[2])['fields'].update(a[3])
            elif a[0]=='UpdateRecord':next(c for c in self.cols if c['id']==a[2])['fields'].update(a[3])
            elif a[0]=='BulkUpdateRecord':
                for j,rid in enumerate(a[2]):next(r for r in self.rows if r['id']==rid)['fields'].update({k:v[j] for k,v in a[3].items()})
class UpgradeTests(unittest.TestCase):
    def test_fails_before_writes_if_unreviewed(self):
        c=Fake([{'id':1,'fields':{'Assessment':'Do not use','AssessmentReason':'No signal'}}])
        with self.assertRaisesRegex(ValueError,'needs reviewed mapping'):upgrade_assessment_schema(c)
        self.assertEqual(c.actions,[])
    def test_reviewed_decision_preserved_atomically_and_repeat(self):
        c=Fake([{'id':1,'fields':{'Assessment':'Do not use','AssessmentReason':'No signal'}},{'id':2,'fields':{'Assessment':'','AssessmentReason':''}}])
        mapping={1:{'expected_assessment':'Do not use','expected_reason':'No signal','checks':{'QCFID':'Unusable'}}}
        a=upgrade_assessment_schema(c,mapping);self.assertEqual(len(c.actions),1);self.assertEqual(len([x for x in a if x[0]=='AddColumn']),7)
        defaults=[x for x in a if x[0]=='AddColumn' and x[2] in QC_CRITERIA]
        self.assertTrue(all(x[3]['recalcDeps']==[['L']] for x in defaults))
        self.assertEqual(c.rows[0]['fields']['QCFID'],'Unusable');self.assertEqual(c.rows[0]['fields']['AssessmentReason'],'No signal')
        self.assertEqual(c.rows[1]['fields']['QCAnalysis'],'Not assessed');self.assertEqual(upgrade_assessment_schema(c),[])
    def test_checklist_labels_cannot_rename_formula_dependencies(self):
        c=Fake([{'id':1,'fields':{'Assessment':'','AssessmentReason':''}}])
        upgrade_assessment_schema(c)
        criteria=[r for r in c.cols if r['fields']['colId'] in QC_CRITERIA]
        self.assertTrue(all(r['fields'].get('untieColIdFromLabel') for r in criteria))
        fid=next(r for r in criteria if r['fields']['colId']=='QCFID')
        fid['fields'].update(label='FID',untieColIdFromLabel=False)
        actions=upgrade_assessment_schema(c)
        self.assertEqual(actions,[['UpdateRecord','_grist_Tables_column',fid['id'],
            {'untieColIdFromLabel':True,'label':'FIDs and spectrum'}]])
        self.assertEqual(fid['fields']['colId'],'QCFID')
        self.assertEqual(upgrade_assessment_schema(c),[])

    def test_pre_tolerance_upgrade_preserves_data_and_is_guarded(self):
        c=Fake([{'id':1,'fields':{'Assessment':'Use','AssessmentReason':'manual note','ImportedAssessment':'Do not use','QCSNR':'Marginal'}}])
        c.cols=[{'id':i+1,'fields':{'parentId':4,'colId':spec['id'],**{k:v for k,v in spec.items() if k!='id'}}} for i,spec in enumerate(columns()['Experiments'])]
        for col in c.cols:
            if col['fields']['colId']=='Assessment':col['fields']['formula']=PRE_TOLERANCE_ASSESSMENT_FORMULA
            if col['fields']['colId']=='AssessmentExplanation':col['fields']['formula']=PRE_TOLERANCE_ASSESSMENT_EXPLANATION_FORMULA
        rows=copy.deepcopy(c.rows)
        actions=upgrade_assessment_schema(c)
        self.assertEqual(len(actions),2)
        self.assertTrue(all(a[0]=='UpdateRecord' and set(a[3])=={'formula'} for a in actions))
        self.assertEqual(c.rows,rows)
        self.assertEqual(upgrade_assessment_schema(c),[])
        next(col for col in c.cols if col['fields']['colId']=='Assessment')['fields']['formula']='unknown formula'
        with self.assertRaisesRegex(ValueError,'requires review'):upgrade_assessment_schema(c)
        self.assertEqual(len(c.actions),1)
        c.cols=[col for col in c.cols if col['fields']['colId']!='QCFID']
        next(col for col in c.cols if col['fields']['colId']=='Assessment')['fields']['formula']=PRE_TOLERANCE_ASSESSMENT_FORMULA
        with self.assertRaisesRegex(ValueError,'incomplete QC schema'):upgrade_assessment_schema(c)
        self.assertEqual(len(c.actions),1)

    def test_changed_or_insufficient_mapping_fails(self):
        c=Fake([{'id':1,'fields':{'Assessment':'Do not use','AssessmentReason':'Bad'}}])
        for mapping in ({1:{'expected_assessment':'Do not use','expected_reason':'Wrong','checks':{'QCFID':'Unusable'}}}, {1:{'expected_assessment':'Do not use','expected_reason':'Bad','checks':{'QCFID':'Clean'}}}):
            with self.assertRaises(ValueError):upgrade_assessment_schema(c,mapping)
        self.assertEqual(c.actions,[])
if __name__=='__main__':unittest.main()
