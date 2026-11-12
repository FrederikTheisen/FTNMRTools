import unittest
from types import SimpleNamespace as Row
from nmr_register.setup import electrostatics_formulas

class CoverageTests(unittest.TestCase):
    def test_uncertainty_duplicates_and_tubes(self):
        def experiment(n, sample=1, **kwargs):
            data=dict(SourceID='source',Folder='folder',Expno=str(n),AcquisitionKey='key',
                ObservedSourceID='source',ObservedFolder='folder',ObservedExpno=str(n),
                RawPresent=True,AcquisitionStatus='Done',EffectiveAcquisitionStatus='Completed',ReconciliationReason='',MigrationIssue='',
                Sample=Row(id=sample) if sample else None,ResultCSV=[],AnalysisStatus='',Assessment='')
            data.update(kwargs)
            return Row(**data)
        rows=[experiment(1,1,ResultCSV=[1]),experiment(2,1),experiment(3,None),
              experiment(4,2,ReconciliationReason='duplicate'),experiment(5,2),experiment(5,3),
              experiment(6,2,ObservedExpno='different')]
        class Experiments:
            @staticmethod
            def lookupRecords(**k):
                return [r for r in rows if all(getattr(r,key)==value for key,value in k.items())]
        def evaluate(name):
            return eval(electrostatics_formulas()[name].replace('$group','group'),
                        {'Experiments':Experiments,'group':rows})
        self.assertEqual(evaluate('VerifiedAcquisitions'),3)
        self.assertEqual(evaluate('KnownTubes'),3)
        self.assertEqual(evaluate('UnknownSamples'),1)
        self.assertEqual(evaluate('SelectedCSVs'),1)
        self.assertEqual(evaluate('AnalysisSummary'),'0 / 1')
        rows[0].AnalysisStatus='Analyzed'
        self.assertEqual(evaluate('AnalysisSummary'),'1 / 1')
        rows[0].AcquisitionStatus='Planned';rows[0].EffectiveAcquisitionStatus='Planned'
        self.assertEqual(evaluate('RecordsDone'),'7 / 6')
        rows[1].EffectiveAcquisitionStatus='Aborted'
        rows[2].EffectiveAcquisitionStatus='Unknown'
        self.assertEqual(evaluate('RecordsDone'),'7 / 4')
        self.assertEqual(evaluate('VerifiedAcquisitions'),3)

class ConditionTests(unittest.TestCase):
    def test_same_tube_and_missing_axes(self):
        import re
        from nmr_register.schema import columns
        formula=next(c['formula'] for c in columns()['Experiments'] if c['id']=='CoverageUncertaintyKey')
        source=re.sub(r'\$([A-Za-z_][A-Za-z_0-9]*)',r"r['\1']",formula)
        namespace={}
        exec('def evaluate(r):\n'+''.join('    '+line+'\n' for line in source.splitlines()),namespace)
        evaluate=namespace['evaluate']
        fields={key:1 for key in re.findall(r'\$([A-Za-z_][A-Za-z_0-9]*)',formula)}
        fields.update(id=1,Sample=Row(id=10),SnapshotNaCl_mM=0,SnapshotpH=0,EffectiveTemperature_C=0)
        self.assertEqual(evaluate(fields),'')
        fields['SnapshotBuffer']=None
        self.assertEqual(evaluate(fields),'sample:10')
        fields['id']=2
        self.assertEqual(evaluate(fields),'sample:10')
        fields['Sample']=Row(id=11)
        self.assertEqual(evaluate(fields),'sample:11')
        fields['Sample']=None
        self.assertEqual(evaluate(fields),'row:2')
        fields.update(Sample=Row(id=10),EffectiveType='')
        self.assertEqual(evaluate(fields),'row:2')
        fields['id']=3
        self.assertEqual(evaluate(fields),'row:3')
