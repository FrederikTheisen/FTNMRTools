import re,unittest
from types import SimpleNamespace as Row
from nmr_register.schema import SAMPLE_TAG_TRIGGER,SAMPLE_TAG_ISSUE_FORMULA,ATTENTION_FORMULA,MACHINE_FIELDS
from nmr_register.setup import upgrade_sample_tag_schema

def evaluate(formula, row, samples, experiments, value=0):
    source=re.sub(r'\$([A-Za-z_][A-Za-z_0-9]*)',r"row['\1']",formula)
    class Table:
        def __init__(self,rows):self.rows=rows
        def lookupRecords(self,**terms):return [r for r in self.rows if all(getattr(r,k,None)==v for k,v in terms.items())]
    ns={'Samples':Table(samples),'Experiments':Table(experiments)}
    exec('def calc(row,value):\n'+''.join('    '+line+'\n' for line in source.splitlines()),ns)
    return ns['calc'](row,value)

class TagFormulaTests(unittest.TestCase):
    def setUp(self):
        self.row=dict(id=1,ObservedSampleID='tube-A',ObservedSampleTagIssue='',ObservedSourceID='nmr',ObservedFolder='run',ObservedExpno='1',AcquisitionKey='nmr:run/1',SourceID='',Folder='',Expno='',ReconciliationReason='',Sample=None)
        self.exp=[Row(id=1,**{k:v for k,v in self.row.items() if k!='id'})];self.samples=[Row(id=8,SampleID='tube-A')]
    def trig(self,value=0):return evaluate(SAMPLE_TAG_TRIGGER,self.row,self.samples,self.exp,value)
    def issue(self):return evaluate(SAMPLE_TAG_ISSUE_FORMULA,self.row,self.samples,self.exp)
    def test_fill_only_and_issues_track_changes(self):
        self.assertEqual(self.trig(),8);self.assertEqual(self.trig(7),7);self.assertEqual(self.issue(),'')
        self.row['Sample']=Row(id=7);self.assertIn('existing sample differs',self.issue())
        self.samples.append(Row(id=9,SampleID='tube-A'));self.assertIn('duplicate Sample ID',self.issue());self.assertEqual(self.trig(),0)
        self.samples=[];self.assertIn('not found',self.issue());self.assertEqual(self.trig(),0)
    def test_identity_and_absence(self):
        self.row['Folder']='other';self.assertEqual(self.trig(),0);self.assertIn('identity differs',self.issue())
        self.row['Folder']='';self.row['ReconciliationReason']='copy';self.assertEqual(self.trig(),0)
        self.row['ReconciliationReason']='';self.row['ObservedSampleID']='';self.assertEqual(self.trig(),0);self.assertEqual(self.issue(),'')
        self.row['ObservedSampleTagIssue']='conflicting title tags';self.assertEqual(self.issue(),'conflicting title tags')
    def test_duplicate_and_mismatched_identities_fail_closed(self):
        self.exp.append(Row(id=2,AcquisitionKey='nmr:run/1',ObservedSourceID='nmr',ObservedFolder='run',ObservedExpno='1'))
        self.assertEqual(self.trig(),0);self.assertIn('duplicate acquisition key',self.issue())
        self.exp.pop();self.row['AcquisitionKey']='wrong';self.assertEqual(self.trig(),0);self.assertIn('key differs',self.issue())
        self.row['AcquisitionKey']='nmr:run/1';self.exp.append(Row(id=3,AcquisitionKey='other',SourceID='nmr',Folder='run',Expno='1'))
        self.assertEqual(self.trig(),0);self.assertIn('duplicate registered',self.issue())
    def test_root_level_folder_and_invalid_expno(self):
        self.row.update(ObservedFolder='',AcquisitionKey='nmr:1')
        self.exp=[Row(id=1,**{k:v for k,v in self.row.items() if k!='id'})]
        self.assertEqual(self.trig(),8);self.assertEqual(self.issue(),'')
        self.row.update(ObservedExpno='not-numeric',AcquisitionKey='nmr:not-numeric')
        self.assertEqual(self.trig(),0);self.assertIn('identity incomplete',self.issue())
    def test_snapshot_trigger_excludes_manual_sample_column(self):
        from nmr_register.schema import columns,trigger_actions
        ids={(table,c['id']):1000*i+j for i,(table,cs) in enumerate(columns().items(),1) for j,c in enumerate(cs,1)}
        actions=trigger_actions({},ids);col=ids[('Experiments','Sample')]
        deps=next(a[3]['recalcDeps'] for a in actions if a[2]==col and 'recalcDeps' in a[3])
        self.assertNotIn(col,deps);self.assertIn(ids[('Experiments','ObservedSampleID')],deps)
    def test_machine_allowlist_keeps_sample_human_owned(self):
        self.assertNotIn('Sample',MACHINE_FIELDS['Experiments']);self.assertTrue({'ObservedSampleID','ObservedSampleTagEvidence','ObservedSampleTagIssue'} <= MACHINE_FIELDS['Experiments']);self.assertIn('$SampleTagIssue',ATTENTION_FORMULA)

class UpgradeTests(unittest.TestCase):
    def test_restartable_upgrade_preserves_sample_values(self):
        from nmr_register.schema import columns
        class Fake:
            def __init__(self):
                self.cols=[{'id':i+1,'fields':{'parentId':4,'colId':c['id'],**{k:v for k,v in c.items() if k!='id'}}} for i,c in enumerate(columns()['Experiments']) if c['id'] not in ('ObservedSampleID','ObservedSampleTagEvidence','ObservedSampleTagIssue','SampleTagIssue')]
                for c in self.cols:
                    if c['fields']['colId']=='Sample':c['fields']['formula']='';c['fields']['recalcDeps']=None
                    if c['fields']['colId']=='AttentionReason':c['fields']['formula']=ATTENTION_FORMULA.replace('if $SampleTagIssue: issues.append($SampleTagIssue)\n','')
                self.records_data=[{'id':1,'fields':{'Sample':9}},{'id':2,'fields':{'Sample':0}}];self.batches=[]
            def records(self,name):return {'_grist_Tables':[{'id':4,'fields':{'tableId':'Experiments'}}], '_grist_Tables_column':self.cols}[name]
            def apply(self,actions):
                self.batches.append(actions)
                for a in actions:
                    if a[0]=='AddColumn':self.cols.append({'id':1000+len(self.cols),'fields':{'parentId':4,'colId':a[2],**a[3]}})
                    else:next(c for c in self.cols if c['id']==a[2])['fields'].update(a[3])
        client=Fake();actions=upgrade_sample_tag_schema(client);self.assertTrue(actions);self.assertEqual([r['fields']['Sample'] for r in client.records_data],[9,0]);self.assertEqual(upgrade_sample_tag_schema(client),[])
        from nmr_register.schema import LEGACY_QC_ATTENTION_FORMULA
        for old in (LEGACY_QC_ATTENTION_FORMULA, LEGACY_QC_ATTENTION_FORMULA.replace('if $SampleTagIssue: issues.append($SampleTagIssue)\n','')):
            with self.subTest(old=old):
                client=Fake()
                client.cols=[c for c in client.cols if c['fields']['colId']!='AssessmentExplanation']
                attention=next(c for c in client.cols if c['fields']['colId']=='AttentionReason')
                attention['fields']['formula']=old
                upgrade_sample_tag_schema(client)
                self.assertEqual(attention['fields']['formula'], LEGACY_QC_ATTENTION_FORMULA.replace('if not $RegisteredAt: issues.append("Registration incomplete")\n',''))
                self.assertNotIn('RegisteredAt',attention['fields']['formula'])
                self.assertEqual(upgrade_sample_tag_schema(client),[])
if __name__=='__main__':unittest.main()
