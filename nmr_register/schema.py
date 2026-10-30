"""Register schema and ownership contract, expressed as Grist user actions."""
import json

CONDITIONS = {'Protein':'Text', 'Construct':'Text', 'Isotope':'Text', 'Buffer':'Text', 'pH':'Numeric',
              'NaCl_mM':'Numeric', 'DNA_uM':'Numeric', 'TFE_pct':'Numeric',
              'Protein_uM':'Numeric', 'OtherConditions':'Text'}
MACHINE_FIELDS = {'Experiments': frozenset(('AcquisitionKey','SourcePath','ObservedSourceID',
    'ObservedFolder','ObservedExpno','ObservedInstrument','ObservedType','ObservedAcquiredAt',
    'ObservedPulseProgram','ObservedField_MHz','ObservedTemperature_C','RawPresent','RawBytes','TransferStable',
    'ScanError','LastChangedAt','ObservedFingerprint','ReconciliationReason','MetadataDifference','ResultCSV','ResultSHA256','ResultSourcePath')),
    'ScanStatus': frozenset(('Scanner','LastScanAt','LastSyncAt','LastError','ScannedRoot'))}
DECISIONS = ('Use', 'Limited use', 'Do not use')
SAMPLE_ID_FORMULA = '''datepart = $PreparedAt.strftime("%Y%m%d") if $PreparedAt else ""
month = $PreparedDateText or ""
if not datepart and len(month) == 7 and month[4] == "-" and month[:4].isdigit() and month[5:].isdigit() and 1 <= int(month[5:]) <= 12:
  datepart = month[:4] + month[5:]
return value if value else ("_".join(x for x in [datepart, $Isotope, $Protein, $DistinguishingCondition, $SampleNumber] if x) if datepart and $Isotope and $Protein and $DistinguishingCondition else "")'''
ATTENTION_FORMULA = '''issues = []
if $MigrationIssue: issues.append($MigrationIssue)
if $ReconciliationReason: issues.append($ReconciliationReason)
if $MetadataDifference: issues.append($MetadataDifference)
if $ScanError: issues.append($ScanError)
if not $RegisteredAt: issues.append("Registration incomplete")
if not $Sample: issues.append("Sample not assigned")
if $Assessment and $Assessment not in ("Use", "Limited use", "Do not use"): issues.append("Invalid assessment")
if $Assessment in ("Limited use", "Do not use") and not ($AssessmentReason or "").strip(): issues.append("Assessment reason missing")
if not $Assessment and ($RawPresent or $AcquisitionStatus == "Acquired") and $EffectiveAcquiredAt and (NOW() - $EffectiveAcquiredAt).total_seconds() > 86400: issues.append("Acquired over 24 hours ago; unassessed")
return "; ".join(issues)'''

def col(name, typ='Text', formula='', is_formula=False, **extra):
    return dict(id=name, type=typ, isFormula=is_formula, formula=formula, **extra)

def columns():
    common = [col(k,v) for k,v in CONDITIONS.items()]
    snapshot = [col('Snapshot'+k,v, '$Sample.'+k+' if $Sample else None', label={'Protein':'Protein','NaCl_mM':'NaCl (mM)'}.get(k,'Sample '+k)) for k,v in CONDITIONS.items()]
    return {
      'Series':[col('SeriesID'), *common, col('Notes')],
      'Samples':[col('SampleID', 'Text', SAMPLE_ID_FORMULA),
        col('PreparedAt','Date'), col('PreparedDateText'), col('DistinguishingCondition'), col('SampleNumber'), col('Series','Ref:Series'),
        *[col(k,v, 'value if value not in (None, "") else ($Series.'+k+' if $Series else None)') for k,v in CONDITIONS.items()],
        col('StorageBox'),col('StorageRow'),col('StorageColumn','Int'),col('Available','Choice',widgetOptions=json.dumps({'choices':['Available','Unavailable','Unknown']})),col('CurrentLocation'),col('Notes'),col('SourceEvidence')],
      'Experiments':[col('AcquisitionKey'),col('SourcePath'),col('ObservedSourceID'),
        col('ObservedFolder'),col('ObservedExpno'),col('ObservedInstrument'),col('ObservedType'),
        col('ObservedAcquiredAt','DateTime'),col('ObservedPulseProgram'),col('ObservedField_MHz','Numeric'),col('ObservedTemperature_C','Numeric'),
        col('RawPresent','Bool'),col('RawBytes','Numeric'),col('TransferStable','Bool'),col('ScanError'),
        col('LastChangedAt','DateTime'),col('ObservedFingerprint'),col('ReconciliationReason'),col('MetadataDifference'),col('ResultCSV','Attachments'),col('ResultSHA256'),col('ResultSourcePath'),
        col('SourceID'),col('Folder'),col('Expno'),col('Instrument'),col('ExperimentType'),
        col('AcquiredAt','DateTime'),col('AcquisitionDateText'),col('AcquisitionStatus'),col('Field_MHz','Numeric'),
        col('Series','Ref:Series'), col('ProteinSelector'), col('RegisteredAt','DateTime'),
        col('Sample','Ref:Samples',widgetOptions=json.dumps({'dropdownCondition':{'text':'(not $Series or choice.Series == $Series) and (not $ProteinSelector or choice.Protein == $ProteinSelector)'}})),col('SnapshotSampleID','Text','$Sample.SampleID if $Sample else None'),*snapshot,
        col('Temperature_C','Numeric',label='Temperature (°C)'),col('PulseProgram'),col('Assessment','Choice',widgetOptions=json.dumps({'choices':list(DECISIONS)})),
        col('AssessmentReason'),col('Notes'),col('SourceEvidence'),col('MigrationKey'),col('MigrationIssue'),
        col('EffectiveInstrument','Text', '$Instrument or $ObservedInstrument',True,label='Instrument'),
        col('EffectiveType','Text', '$ExperimentType or $ObservedType',True,label='Type'),
        col('EffectiveAcquiredAt','DateTime', '$AcquiredAt or $ObservedAcquiredAt',True,label='Date'),
        col('EffectiveField_MHz','Numeric', '$Field_MHz if $Field_MHz is not None else $ObservedField_MHz',True,label='Field (MHz)'),
        col('EffectiveTemperature_C','Numeric', '$Temperature_C if $Temperature_C is not None else $ObservedTemperature_C',True,label='Temperature (°C)'),
        col('AttentionReason','Text',ATTENTION_FORMULA,True),
        col('NeedsAttention','Bool','bool($AttentionReason)',True)],
      'StorageGrid':[col('StorageBox'),col('Row'), *[col('Col'+str(n), 'Text', '", ".join(s.SampleID for s in Samples.lookupRecords(StorageBox=$StorageBox, StorageRow=$Row, StorageColumn='+str(n)+'))', True, label=str(n)) for n in range(1,13)]],
      'Unclassified':[col('SourceKey'),col('Label'),col('SourceEvidence'),col('Notes'),col('Resolution')],
      'ScanStatus':[col('Scanner'),col('LastScanAt','DateTime'),col('LastSyncAt','DateTime'),col('LastError'),col('ScannedRoot')],
    }

def schema_actions():
    return [['AddTable', name, cols] for name,cols in columns().items()]

def trigger_actions(table_ids, column_ids):
    """Configure defaults on new samples and snapshot only on Sample changes.

    IDs are Grist metadata IDs keyed by table name and (table, colId).
    """
    actions=[]
    sample_ref=column_ids[('Experiments','Sample')]
    for name in ['SnapshotSampleID']+['Snapshot'+k for k in CONDITIONS]:
        actions.append(['UpdateRecord','_grist_Tables_column',column_ids[('Experiments',name)],
                        {'recalcWhen':0,'recalcDeps':['L',sample_ref]}])
    series_ref=column_ids[('Samples','Series')]
    id_deps=[column_ids[('Samples',name)] for name in ('PreparedAt','PreparedDateText','Protein','Isotope','DistinguishingCondition','SampleNumber')]
    actions.append(['UpdateRecord','_grist_Tables_column',column_ids[('Samples','SampleID')],
                    {'recalcWhen':0,'recalcDeps':['L',*id_deps]}])
    for name in CONDITIONS:
        actions.append(['UpdateRecord','_grist_Tables_column',column_ids[('Samples',name)],
                        {'recalcWhen':0,'recalcDeps':['L',series_ref]}])
    for table, name, target, label in (('Samples','Series','Series','SeriesID'),
                                       ('Experiments','Series','Series','SeriesID'),
                                       ('Experiments','Sample','Samples','SampleID')):
        ref=column_ids[(table,name)]
        actions.append(['UpdateRecord','_grist_Tables_column',ref,
                        {'visibleCol':column_ids[(target,label)]}])
        actions.append(['SetDisplayFormula',table,None,ref,'$'+name+'.'+label])
    return actions

def snapshot_values(sample):
    return {'SnapshotSampleID': sample.get('SampleID'), **{'Snapshot'+k:sample.get(k) for k in CONDITIONS}}

def series_defaults(series, sample=None):
    result=dict(sample or {})
    for key in CONDITIONS:
        if key not in result or result[key] is None:
            result[key]=series.get(key)
    return result
