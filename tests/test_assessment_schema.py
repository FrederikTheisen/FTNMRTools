import copy,itertools,re,unittest
from functools import lru_cache
from nmr_register.schema import (QC_CRITERIA,QC_DEFAULT,QC_CAUTION,QC_FAILING,assessment_from_values,
    ASSESSMENT_FORMULA,ASSESSMENT_EXPLANATION_FORMULA,PRE_TOLERANCE_ASSESSMENT_FORMULA,PRE_TOLERANCE_ASSESSMENT_EXPLANATION_FORMULA,ATTENTION_FORMULA,LEGACY_QC_ATTENTION_FORMULA,columns)
from nmr_register.setup import upgrade_assessment_schema

@lru_cache(maxsize=None)
def compiled_formula(source):
    import re
    source=re.sub(r'\$([A-Za-z_][A-Za-z_0-9]*)',r"values.get('\1')",source)
    ns={};exec('def calc(values):\n'+''.join('    '+line+'\n' for line in source.splitlines()),ns)
    return ns['calc']

def formula(source, values):
    return compiled_formula(source)(values)
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
            for key in cautions[:count]:values[key]=QC_CAUTION[key][0]
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
class AnalysisClassificationTests(unittest.TestCase):
    def test_historical_formulas_are_exact(self):
        import hashlib
        from nmr_register import schema
        expected = {'LEGACY_ASSESSMENT_FORMULA': '9a22f5b28edba61eaaac8def53bf21f21516502ddd6da8a5e296cd3df426b0fe', 'LEGACY_ASSESSMENT_EXPLANATION_FORMULA': '07e372676f8c79f39b8ebee1a7388d1d743a19f9648ce12ea22d58ef8d12bfaa', 'PRE_TOLERANCE_ASSESSMENT_FORMULA': 'f6d4ceeabafa28cfff9dee6d40b9d31075023fcceee17e6d256688ee4177f45d', 'PRE_TOLERANCE_ASSESSMENT_EXPLANATION_FORMULA': 'e1b8f718bd7d6b1363cc36ca78ce6a13f29bef2ea3939cd73b7212fdc657e231', 'PRE_PARTIAL_ASSESSMENT_EXPLANATION_FORMULA': '9c2cad1c3543572fc99e7bb54c5362c2b810eacd5d10eb2a7911bedd86d59cfe', 'PRE_COMPACT_ASSESSMENT_DISPLAY_FORMULA': '0af9b33c737eedc7ad5871e93c3a1d6137f9a5048d550b8ce9707225079dd8c7', 'PRE_ANALYSIS_ASSESSMENT_FORMULA': '1082501a0bbc3522ad36fd55dfc71e5d3caf3c1a97dc24d22a41b3b7e9fabc81', 'PRE_ANALYSIS_ASSESSMENT_EXPLANATION_FORMULA': '7b1f3ffa7b6f58c43bb9ddbf435eb6c04789bfbd1411bfd15e551a352c78f93a', 'PRE_ANALYSIS_ASSESSMENT_DISPLAY_FORMULA': 'a51270df4bb2b0d58f5b46f9d3cd3da84b29d9a6d94d4fdb54f87680e45f6537'}
        for name, digest in expected.items():
            self.assertEqual(hashlib.sha256(getattr(schema, name).encode()).hexdigest(), digest, name)

    def test_every_analysis_choice_and_failure_precedence(self):
        from nmr_register.schema import ASSESSMENT_DISPLAY_FORMULA, assessment_display_from_values
        good = {k:v[1][0] for k,v in QC_CRITERIA.items()}
        for choice in QC_CRITERIA['QCAnalysis'][1]:
            for extra_caution in (False, True):
                values = dict(good, QCAnalysis=choice, ExperimentType='Other')
                if extra_caution: values['QCSNR'] = 'Marginal'
                expected = 'Do not use' if choice in QC_FAILING['QCAnalysis'] else 'Limited use' if choice == 'Limited' and extra_caution else 'Use'
                self.assertEqual(assessment_from_values(values), expected)
                self.assertEqual(formula(ASSESSMENT_FORMULA, values), expected)
                self.assertEqual(assessment_display_from_values(values), expected)
                self.assertEqual(formula(ASSESSMENT_DISPLAY_FORMULA, values), expected)
                values['QCSettings'] = QC_DEFAULT
                if choice in QC_FAILING['QCAnalysis']:
                    self.assertEqual(formula(ASSESSMENT_FORMULA, values), 'Do not use')
                    self.assertIn(choice, formula(ASSESSMENT_EXPLANATION_FORMULA, values))

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

    def test_analysis_upgrade_is_metadata_only_preserves_options_and_idempotent(self):
        import json
        from nmr_register import schema
        c = Fake([])
        c.cols = [{'id':i+1,'fields':{'parentId':4,'colId':spec['id'],**{k:v for k,v in spec.items() if k!='id'}}} for i,spec in enumerate(columns()['Experiments'])]
        c.rows = [{'id':i+1,'fields':{'QCAnalysis':choice,'AssessmentReason':'keep','ImportedAssessment':'Limited use'}} for i,choice in enumerate(('Reliable','Limited','Unreliable'))]
        before = copy.deepcopy(c.rows)
        for name in ('Assessment','AssessmentExplanation','AssessmentDisplay'):
            next(r for r in c.cols if r['fields']['colId']==name)['fields']['formula'] = getattr(schema, 'PRE_ANALYSIS_' + {'Assessment':'ASSESSMENT_FORMULA','AssessmentExplanation':'ASSESSMENT_EXPLANATION_FORMULA','AssessmentDisplay':'ASSESSMENT_DISPLAY_FORMULA'}[name])
        analysis = next(r for r in c.cols if r['fields']['colId']=='QCAnalysis')
        analysis['fields']['widgetOptions'] = json.dumps({'choices':['Not assessed','Reliable','Limited','Unreliable'],'choiceOptions':{'Reliable':{'fillColor':'green'}},'alignment':'left'})
        actions = upgrade_assessment_schema(c)
        self.assertEqual(len(actions), 4)
        self.assertTrue(all(a[0:2]==['UpdateRecord','_grist_Tables_column'] for a in actions))
        self.assertEqual(c.rows, before)
        options = json.loads(analysis['fields']['widgetOptions'])
        self.assertEqual(options['choices'], [QC_DEFAULT,*QC_CRITERIA['QCAnalysis'][1]])
        self.assertEqual(options['choiceOptions'], {'Reliable':{'fillColor':'green'}})
        self.assertEqual(options['alignment'], 'left')
        self.assertEqual(upgrade_assessment_schema(c), [])
        analysis['fields']['widgetOptions'] = '{bad'
        with self.assertRaisesRegex(ValueError,'choices require review'): upgrade_assessment_schema(c)
        self.assertEqual(len(c.actions), 1)

    def test_changed_or_insufficient_mapping_fails(self):
        c=Fake([{'id':1,'fields':{'Assessment':'Do not use','AssessmentReason':'Bad'}}])
        for mapping in ({1:{'expected_assessment':'Do not use','expected_reason':'Wrong','checks':{'QCFID':'Unusable'}}}, {1:{'expected_assessment':'Do not use','expected_reason':'Bad','checks':{'QCFID':'Clean'}}}):
            with self.assertRaises(ValueError):upgrade_assessment_schema(c,mapping)
        self.assertEqual(c.actions,[])
if __name__=='__main__':unittest.main()
