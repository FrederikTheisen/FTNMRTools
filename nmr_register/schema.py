"""Register schema and ownership contract, expressed as Grist user actions."""
import json

CONDITIONS = {'Protein':'Text', 'Construct':'Text', 'Isotope':'Text', 'Buffer':'Text', 'pH':'Numeric',
              'NaCl_mM':'Numeric', 'DNA_uM':'Numeric', 'TFE_pct':'Numeric',
              'Protein_uM':'Numeric', 'OtherConditions':'Text'}
CANONICAL_ACQUISITION_FIELDS = frozenset(('Instrument','ExperimentType','AcquiredAt','AcquisitionStatus','Field_MHz','Temperature_C','PulseProgram','D1_s','DelayList_s','SpinLock_Hz','SpinLockOffset_Hz'))
# Identity remains human-editable; the collector may only fill verified blanks.
FILL_ONLY_IDENTITY_FIELDS = ('SourceID','Folder','Expno')
MACHINE_FIELDS = {'Experiments': CANONICAL_ACQUISITION_FIELDS | frozenset(('AcquisitionKey','SourcePath','ObservedSourceID',
    'ObservedFolder','ObservedExpno','ObservedAcquisitionEvidence','ObservedProtocolEvidence','RawPresent','RawBytes','TransferStable',
    'ObservedSampleID','ObservedSampleTagEvidence','ObservedSampleTagIssue','Title','TitleSource','ScanError','LastChangedAt','ObservedFingerprint','ReconciliationReason','MetadataDifference','ResultCSV','ResultSHA256','ResultSourcePath')),
    'ScanStatus': frozenset(('Scanner','LastScanAt','LastSyncAt','LastError','ScannedRoot'))}
def canonical_experiment_type(value):
    """Human-readable labels for known instrument classifications; unknown stays unset."""
    if not isinstance(value,str):return None
    normalized=''.join(ch for ch in value.lower() if ch not in ' _-\t\r\n')
    return {'1d':'1D','r1':'T1','t1':'T1','15nt1':'T1','r1rho':'T1rho','t1rho':'T1rho','t1ρ':'T1rho',
            'noe':'HetNOE','hetnoe':'HetNOE','heteronuclearnoe':'HetNOE',
            'hsqc':'HSQC','15nhsqc':'HSQC','13chsqc':'13C HSQC','r2':'T2','t2':'T2','etaxy':'etaXY','ηxy':'etaXY',
            'watert1sat':'Water T1 (T1sat)','t1sat':'Water T1 (T1sat)',
            'nicot1sat':'Water T1 (T1sat)','watert1(t1sat)':'Water T1 (T1sat)',
            'watert1':'Water T1 (T1sat)','cpmg':'CPMG','hmbc':'HMBC',
            'cosy':'COSY','noesy':'NOESY','calibration':'Calibration'}.get(normalized)

SOURCE_LABEL_FORMULA = """def text(value): return str(value).strip() if value is not None else ""
hf, he, hs = text($Folder), text($Expno), text($SourceID)
of, oe, os = text($ObservedFolder), text($ObservedExpno), text($ObservedSourceID)
def display(folder, expno, root=False):
  return (expno or "?" ) + " · " + (folder or ("(root)" if root else "(folder not set)"))
human = bool(hf or he)
observed = bool(of or oe)
root = not of and bool(os and oe) and $AcquisitionKey == os + ":" + oe
if human and observed:
  if hf == of and he == oe and not (hs and os and hs != os): return display(of, oe, root)
  source_diff = bool(hs and os and hs != os)
  return "Observed" + (" (different source)" if source_diff else "") + ": " + display(of, oe, root) + " | Recorded" + (" (different source)" if source_diff else "") + ": " + display(hf, he)
if human: return display(hf, he)
if observed: return display(of, oe, root)
return ""
"""

ATTENTION_IDENTITY_PREFIX = """def text(value): return str(value).strip() if value is not None else ""
hf, he, hs = text($Folder), text($Expno), text($SourceID)
of, oe, os = text($ObservedFolder), text($ObservedExpno), text($ObservedSourceID)
human, observed = bool(hf or he), bool(of or oe)
root = not of and bool(os and oe) and $AcquisitionKey == os + ":" + oe
source_diff = bool(hs and os and hs != os)
conflict = human and observed and (hf != of or he != oe or source_diff)
"""

def attention_identity_formula(part):
    if part=='folder':
        values='human_value = hf or "(folder not set)"\nobserved_value = of or ("(root)" if root else "(folder not set)")\n'
    elif part=='expno':
        values='human_value = he or "?"\nobserved_value = oe or "?"\n'
    else:raise ValueError('Unknown identity part')
    return ATTENTION_IDENTITY_PREFIX + values + """if conflict: return "Observed" + (" (different source)" if source_diff else "") + ": " + observed_value + " | Recorded" + (" (different source)" if source_diff else "") + ": " + human_value
if observed: return observed_value
if human: return human_value
return ""
"""

DATE_ONLY_OPTIONS = {'dateFormat': 'YYYY-MM-DD', 'timeFormat': ''}
DATE_ONLY_EXPERIMENT_COLUMNS = frozenset(('AcquiredAt', 'EffectiveAcquiredAt', 'RegisteredAt'))
LABELLING_CHOICES = ('15N', '13C', '2H', '15N+13C', '15N+2H', '13C+2H', '15N+13C+2H')
SAMPLE_STATUSES = ('Usable', 'Suspected degradation', 'Degraded', 'Dead')
DECISIONS = ('Use', 'Limited use', 'Do not use')
ACQUISITION_STATUSES = ('Planned', 'Dummy', 'Setup', 'Running', 'Completed', 'Aborted', 'Unknown')
QC_CRITERIA = {
    'QCAcquisitionCoverage': ('Acquisition coverage', ('Complete','Sufficient','Insufficient')),
    'QCSettings': ('Settings', ('Appropriate','Usable with restrictions','Invalid')),
    'QCFID': ('FIDs and spectrum', ('Clean','Usable with exclusions','Unusable')),
    'QCSNR': ('SNR', ('Adequate','Marginal','Inadequate')),
    'QCSample': ('Sample during acquisition', ('Acceptable','Concern','Unacceptable')),
    'QCAnalysis': ('Analysis quality', ('Reliable','Reliable — early points excluded','Two-phase — unresolved','Limited','Unreliable')),
}
QC_DEFAULT = 'Not assessed'
QC_ACCEPTABLE = {
    'QCAcquisitionCoverage': ('Complete','Sufficient'),
    'QCSettings': ('Appropriate',), 'QCFID': ('Clean',),
    'QCSNR': ('Adequate',), 'QCSample': ('Acceptable',),
    'QCAnalysis': ('Reliable','Reliable — early points excluded'),
}
QC_CAUTION = {
    'QCAcquisitionCoverage': (), 'QCSettings': ('Usable with restrictions',),
    'QCFID': ('Usable with exclusions',), 'QCSNR': ('Marginal',),
    'QCSample': ('Concern',), 'QCAnalysis': ('Limited',),
}
QC_FAILING = {
    'QCAcquisitionCoverage': ('Insufficient',), 'QCSettings': ('Invalid',),
    'QCFID': ('Unusable',), 'QCSNR': ('Inadequate',),
    'QCSample': ('Unacceptable',), 'QCAnalysis': ('Two-phase — unresolved','Unreliable'),
}

def assessment_from_values(values):
    failures=[]; limited=[]; missing=[]
    if all(values.get(name) in (None, '', QC_DEFAULT) for name in QC_CRITERIA):
        imported = values.get('ImportedAssessment')
        if imported in DECISIONS:return imported
    for name,(label,choices) in QC_CRITERIA.items():
        item=values.get(name) or QC_DEFAULT
        if item in QC_FAILING[name]:failures.append(label+': '+item)
        elif item in QC_ACCEPTABLE[name]:pass
        elif item in QC_CAUTION.get(name, ()):limited.append(label+': '+item)
        else:missing.append(label)
    if failures:return 'Do not use'
    if missing:return ''
    if len(limited)>=2:return 'Limited use'
    return 'Use'

QC_MISSING_LABELS = ('acquisitions', 'settings', 'FIDs + spectra', 'SNR', 'sample', 'analysis')

def assessment_display_from_values(values):
    """Show provisional progress without approving an incomplete checklist."""
    decision = assessment_from_values(values)
    missing = []
    rated = False
    for (name, (_, choices)), short in zip(QC_CRITERIA.items(), QC_MISSING_LABELS):
        value = values.get(name)
        if value in choices:
            rated = True
        else:
            missing.append(short)
    if not missing:
        return decision
    if all(values.get(name) in (None, '', QC_DEFAULT) for name in QC_CRITERIA) and decision:
        return decision
    if not rated:
        return 'Not assessed'
    if decision == 'Do not use':
        return decision
    suffix = ' · analysis missing' if missing == ['analysis'] else ''
    cautions = sum(values.get(name) in QC_CAUTION.get(name, ()) for name in QC_CRITERIA)
    return ('Limited use' if cautions >= 2 else 'Use') + ' so far' + suffix

_QC_FORMULA_PREFIX = """criteria = [
  ("Acquisition coverage", $QCAcquisitionCoverage, ("Complete", "Sufficient"), None, "Insufficient"),
  ("Settings", $QCSettings, ("Appropriate",), "Usable with restrictions", "Invalid"),
  ("FIDs and spectrum", $QCFID, ("Clean",), "Usable with exclusions", "Unusable"),
  ("SNR", $QCSNR, ("Adequate",), "Marginal", "Inadequate"),
  ("Sample during acquisition", $QCSample, ("Acceptable",), "Concern", "Unacceptable"),
  ("Analysis quality", $QCAnalysis, ("Reliable",), "Limited", "Unreliable"),
]
failed = [label + ": " + str(value) for label, value, good, middle, bad in criteria if value == bad]
missing = [label for label, value, good, middle, bad in criteria if not value or value == "Not assessed" or (value not in good and value != bad and value != middle)]
limited = [label + ": " + str(value) for label, value, good, middle, bad in criteria if middle and value == middle]
"""
LEGACY_ASSESSMENT_FORMULA = _QC_FORMULA_PREFIX + """if failed: return "Do not use"
if missing: return ""
if limited: return "Limited use"
return 'Use'"""
LEGACY_ASSESSMENT_EXPLANATION_FORMULA = _QC_FORMULA_PREFIX + """if failed: return "; ".join(failed)
if missing: return "Assessment incomplete: " + ", ".join(missing)
if limited: return "; ".join(limited)
return 'All criteria acceptable'"""

_IMPORTED_ASSESSMENT_PREFIX = '''if all(not value or value == "Not assessed" for label, value, good, middle, bad in criteria) and $ImportedAssessment in ("Use", "Limited use", "Do not use"):
  return $ImportedAssessment
'''
PRE_TOLERANCE_ASSESSMENT_FORMULA = _QC_FORMULA_PREFIX + _IMPORTED_ASSESSMENT_PREFIX + LEGACY_ASSESSMENT_FORMULA[len(_QC_FORMULA_PREFIX):]
PRE_TOLERANCE_ASSESSMENT_EXPLANATION_FORMULA = _QC_FORMULA_PREFIX + _IMPORTED_ASSESSMENT_PREFIX.replace('return $ImportedAssessment', 'return "Imported assessment: " + $ImportedAssessment + "; checklist not yet completed"') + LEGACY_ASSESSMENT_EXPLANATION_FORMULA[len(_QC_FORMULA_PREFIX):]
ASSESSMENT_FORMULA = PRE_TOLERANCE_ASSESSMENT_FORMULA.replace('if limited: return "Limited use"', 'if len(limited) >= 2: return "Limited use"')
ASSESSMENT_EXPLANATION_FORMULA = PRE_TOLERANCE_ASSESSMENT_EXPLANATION_FORMULA.replace(
    'if limited: return "; ".join(limited)',
    'if len(limited) >= 2: return "; ".join(limited)\nif limited: return "Single caution accepted: " + limited[0]')
PRE_PARTIAL_ASSESSMENT_EXPLANATION_FORMULA = ASSESSMENT_EXPLANATION_FORMULA
PRE_COMPACT_ASSESSMENT_DISPLAY_FORMULA = _QC_FORMULA_PREFIX + _IMPORTED_ASSESSMENT_PREFIX + """if not missing:
  return "Do not use" if failed else ("Limited use" if len(limited) >= 2 else "Use")
if not any(value in good or value == bad or (middle and value == middle) for label, value, good, middle, bad in criteria):
  return "Not assessed"
short_labels = dict(zip([c[0] for c in criteria], ["acquisitions", "settings", "FIDs + spectra", "SNR", "sample", "analysis"]))
marker = short_labels[missing[0]] + " missing" if len(missing) == 1 else str(len(missing)) + " checks missing"
if failed: return "Do not use · " + marker
return ("Limited use" if len(limited) >= 2 else "Use") + " so far · " + marker"""
ASSESSMENT_DISPLAY_FORMULA = _QC_FORMULA_PREFIX + _IMPORTED_ASSESSMENT_PREFIX + """if failed: return "Do not use"
if not missing: return "Limited use" if len(limited) >= 2 else "Use"
if not any(value in good or (middle and value == middle) for label, value, good, middle, bad in criteria):
  return "Not assessed"
return ("Limited use" if len(limited) >= 2 else "Use") + " so far" + (" · analysis missing" if missing == ["Analysis quality"] else "")"""

ASSESSMENT_EXPLANATION_FORMULA = _QC_FORMULA_PREFIX + _IMPORTED_ASSESSMENT_PREFIX.replace('return $ImportedAssessment', 'return "Imported assessment: " + $ImportedAssessment + "; checklist not yet completed"') + """if missing:
  details = []
  if failed: details.append("Do not use: " + "; ".join(failed))
  elif any(value in good or (middle and value == middle) for label, value, good, middle, bad in criteria):
    details.append(("Limited use" if len(limited) >= 2 else "Use") + " so far")
  else: details.append("Not assessed")
  if limited: details.append("Cautions: " + "; ".join(limited))
  details.append("Missing checks: " + ", ".join(missing))
  return "; ".join(details)
if failed: return "; ".join(failed)
if len(limited) >= 2: return "; ".join(limited)
if limited: return "Single caution accepted: " + limited[0]
return 'All criteria acceptable'"""

# Preserve exact previously shipped formulas for guarded metadata-only upgrades.
PRE_ANALYSIS_ASSESSMENT_FORMULA = ASSESSMENT_FORMULA
PRE_ANALYSIS_ASSESSMENT_EXPLANATION_FORMULA = ASSESSMENT_EXPLANATION_FORMULA
PRE_ANALYSIS_ASSESSMENT_DISPLAY_FORMULA = ASSESSMENT_DISPLAY_FORMULA
_ANALYSIS_QC_FORMULA_PREFIX = "criteria = [\n" + "".join(
    "  (" + repr(label) + ", $" + name + ", " + repr(QC_ACCEPTABLE[name]) + ", "
    + repr(QC_CAUTION.get(name, ())) + ", " + repr(QC_FAILING[name]) + "),\n"
    for name, (label, choices) in QC_CRITERIA.items()) + """ ]
failed = [label + ": " + str(value) for label, value, good, middle, bad in criteria if value in bad]
missing = [label for label, value, good, middle, bad in criteria if not value or value == "Not assessed" or (value not in good and value not in bad and value not in middle)]
limited = [label + ": " + str(value) for label, value, good, middle, bad in criteria if value in middle]
"""
def _analysis_formula(previous):
    return (_ANALYSIS_QC_FORMULA_PREFIX + previous[len(_QC_FORMULA_PREFIX):]).replace(
        'value == bad', 'value in bad').replace('value == middle', 'value in middle')

ASSESSMENT_FORMULA = _analysis_formula(PRE_ANALYSIS_ASSESSMENT_FORMULA)
ASSESSMENT_EXPLANATION_FORMULA = _analysis_formula(PRE_ANALYSIS_ASSESSMENT_EXPLANATION_FORMULA)
ASSESSMENT_DISPLAY_FORMULA = _analysis_formula(PRE_ANALYSIS_ASSESSMENT_DISPLAY_FORMULA)

SAMPLE_ID_FORMULA = '''datepart = $PreparedAt.strftime("%Y%m%d") if $PreparedAt else ""
month = $PreparedDateText or ""
if not datepart and len(month) == 7 and month[4] == "-" and month[:4].isdigit() and month[5:].isdigit() and 1 <= int(month[5:]) <= 12:
  datepart = month[:4] + month[5:]
return value if value else ("_".join(x for x in [datepart, $Isotope, $Protein, $DistinguishingCondition, $SampleNumber] if x) if datepart and $Isotope and $Protein and $DistinguishingCondition else "")'''
ATTENTION_FORMULA = '''if $AcquisitionStatus == "Dummy": return ""
issues = []
if $MigrationIssue: issues.append($MigrationIssue)
if $SampleTagIssue: issues.append($SampleTagIssue)
if $ReconciliationReason: issues.append($ReconciliationReason)
if $MetadataDifference: issues.append($MetadataDifference)
if $ScanError: issues.append($ScanError)
if not $RegisteredAt: issues.append("Registration incomplete")
if not $Sample: issues.append("Sample not assigned")
if $Assessment and $Assessment not in ("Use", "Limited use", "Do not use"): issues.append("Invalid assessment")
if $Assessment in ("Limited use", "Do not use") and not ($AssessmentReason or "").strip() and not ($AssessmentExplanation or "").strip(): issues.append("Assessment reason missing")
if $EffectiveAcquisitionStatus == "Aborted": issues.append("Acquisition aborted")
if $EffectiveAcquisitionStatus == "Unknown": issues.append("Acquisition status unknown")
if not $Assessment and $EffectiveAcquisitionStatus == "Completed" and $EffectiveAcquiredAt and (NOW() - $EffectiveAcquiredAt).total_seconds() > 86400: issues.append("Acquired over 24 hours ago; unassessed")
return "; ".join(issues)'''

PRE_CALIBRATION_ATTENTION_FORMULA = ATTENTION_FORMULA
PRE_DUMMY_ATTENTION_FORMULA = ATTENTION_FORMULA.split("\n", 1)[1]
ATTENTION_FORMULA = ATTENTION_FORMULA.replace('if not $Assessment and $EffectiveAcquisitionStatus == "Completed"', 'if $ExperimentType != "Calibration" and not $Assessment and $EffectiveAcquisitionStatus == "Completed"')

PRE_PARTIAL_ATTENTION_FORMULA = ATTENTION_FORMULA
ATTENTION_FORMULA = ATTENTION_FORMULA.replace(
    'issues.append("Acquired over 24 hours ago; unassessed")',
    'issues.append("Acquired over 24 hours ago; " + $AssessmentDisplay if $AssessmentDisplay and $AssessmentDisplay != "Not assessed" else "Acquired over 24 hours ago; unassessed")')

PRE_STUDY_ASSESSMENT_ATTENTION_FORMULA = ATTENTION_FORMULA
ATTENTION_FORMULA = ATTENTION_FORMULA.replace(
    'if $ExperimentType != "Calibration" and not $Assessment',
    'if ($Study or "").strip() and $ExperimentType != "Calibration" and not $Assessment')

PRE_LIVE_SAMPLE_ATTENTION_FORMULA = ATTENTION_FORMULA
ATTENTION_FORMULA = ATTENTION_FORMULA.replace(
    'if $MigrationIssue: issues.append($MigrationIssue)',
    '''migration_issue = $MigrationIssue or ""
if $Sample and migration_issue:
  kept = []
  removed = False
  skip_condition = False
  for part in migration_issue.replace("\\r\\n", "\\n").replace("\\n", ";").split(";"):
    clause = part.strip()
    if not clause: continue
    if clause.startswith("Sample link unresolved:"):
      removed = True
      skip_condition = True
      continue
    if skip_condition and clause.lower().startswith("workbook condition:"):
      skip_condition = False
      continue
    skip_condition = False
    kept.append(clause)
  if removed: migration_issue = "; ".join(kept)
if migration_issue: issues.append(migration_issue)''')

PRE_ASSIGNED_SAMPLE_ATTENTION_FORMULA = ATTENTION_FORMULA
ATTENTION_FORMULA = ATTENTION_FORMULA.replace('if $Sample and migration_issue:', 'if migration_issue:')

LEGACY_QC_ATTENTION_FORMULA = PRE_DUMMY_ATTENTION_FORMULA.replace(
    ' and not ($AssessmentExplanation or "").strip()', '')

# Both formulas check the same identity invariants. The trigger returns the
# current value first: observed metadata never replaces a human assignment.
TAG_IDENTITY_GUARD = """source, folder, expno = $ObservedSourceID, $ObservedFolder, $ObservedExpno
if not source or folder is None or not expno or not str(expno).isdigit() or not $AcquisitionKey: return "Sample tag: acquisition identity incomplete"
if $AcquisitionKey != source + ":" + (folder + "/" if folder else "") + str(expno): return "Sample tag: acquisition key differs from observed identity"
for confirmed, observed in [($SourceID, source), ($Folder, folder), ($Expno, expno)]:
  if confirmed not in (None, "") and str(confirmed) != str(observed): return "Sample tag: acquisition identity differs from registration"
if $ReconciliationReason: return "Sample tag: acquisition identity requires review"
if len(Experiments.lookupRecords(AcquisitionKey=$AcquisitionKey)) != 1: return "Sample tag: duplicate acquisition key"
observed = Experiments.lookupRecords(ObservedSourceID=source, ObservedFolder=folder, ObservedExpno=expno)
if len(observed) != 1: return "Sample tag: duplicate observed acquisition identity"
manual = Experiments.lookupRecords(SourceID=source, Folder=folder, Expno=expno)
if any(r.id != $id for r in manual): return "Sample tag: duplicate registered acquisition identity"
"""
SAMPLE_TAG_TRIGGER = """if value: return value
if not $ObservedSampleID or $ObservedSampleTagIssue: return value
""" + "\n".join(line.split('return "Sample tag:', 1)[0] + 'return value'
                if 'return "Sample tag:' in line else line
                for line in TAG_IDENTITY_GUARD.splitlines()) + """
matches = Samples.lookupRecords(SampleID=$ObservedSampleID)
return matches[0].id if len(matches) == 1 else value"""
SAMPLE_TAG_ISSUE_FORMULA = """if $ObservedSampleTagIssue: return $ObservedSampleTagIssue
if not $ObservedSampleID: return ""
""" + TAG_IDENTITY_GUARD + """matches = Samples.lookupRecords(SampleID=$ObservedSampleID)
if len(matches) == 0: return "Sample tag: Sample ID not found"
if len(matches) > 1: return "Sample tag: duplicate Sample ID"
if $Sample and $Sample.id != matches[0].id: return "Sample tag: existing sample differs from tag"
return ''"""

def col(name, typ='Text', formula='', is_formula=False, **extra):
    return dict(id=name, type=typ, isFormula=is_formula, formula=formula, **extra)

def acquisition_columns():
    """Canonical settings and hidden compatibility aliases; numeric blanks are null."""
    result=[]
    for name,typ,label in [('D1_s','Numeric','D1 (s)'),('DelayList_s','Text','Delay list (s)'),
                           ('SpinLock_Hz','Numeric','Spin-lock field (Hz)'),('SpinLockOffset_Hz','Numeric','Spin-lock offset (Hz)')]:
        extra={'description':'Ordered seconds, for example [0, 0.1, 0.1].'} if name=='DelayList_s' else {}
        default={'recalcWhen':0,'recalcDeps':['L']} if typ=='Numeric' else {}
        result.extend([col(name,typ,'None' if typ=='Numeric' else '',label=label,**extra,**default),
                       col('Effective'+name,typ,'$'+name,True,label=label,**extra)])
    return result+[col('ObservedProtocolEvidence',label='Acquisition settings source'),
                   col('EffectivePulseProgram','Text','$PulseProgram',True,label='Pulse program')]


def columns():
    common = [col(k,'Choice',widgetOptions=json.dumps({'choices':list(LABELLING_CHOICES) if k=='Isotope' else []})) if k in ('Protein','Isotope') else col(k,v) for k,v in CONDITIONS.items()]
    snapshot = [col('Snapshot'+k,v, '$Sample.'+k+' if $Sample else None', label={'Protein':'Protein','NaCl_mM':'NaCl (mM)'}.get(k,'Sample '+k)) for k,v in CONDITIONS.items()]
    return {
      'Series':[col('SeriesID'), *common, col('Notes')],
      'Samples':[col('SampleID', 'Text', SAMPLE_ID_FORMULA),
        col('PreparedAt','Date'), col('PreparedDateText'), col('DistinguishingCondition'), col('SampleNumber'), col('Series','Ref:Series'),
        *[col(k,'Choice' if k in ('Protein','Isotope') else v, 'value if value not in (None, "") else ($Series.'+k+' if $Series else None)', **({'widgetOptions':json.dumps({'choices':list(LABELLING_CHOICES) if k=='Isotope' else []})} if k in ('Protein','Isotope') else {})) for k,v in CONDITIONS.items()],
        col('StorageBox'),col('StorageRow'),col('StorageColumn','Int'),col('Available','Choice','"Usable"',recalcWhen=0,recalcDeps=['L'],widgetOptions=json.dumps({'choices':list(SAMPLE_STATUSES)})),col('CurrentLocation'),col('Notes'),col('SourceEvidence')],
      'Experiments':[col('AcquisitionKey'),col('SourcePath'),col('ObservedSourceID'),
        col('ObservedFolder'),col('ObservedExpno'),col('Title',label='Title'),col('TitleSource'),
        col('RawPresent','Bool'),col('RawBytes','Numeric'),col('TransferStable','Bool'),col('ScanError'),
        col('ObservedSampleID',label='Observed sample tag'),col('ObservedSampleTagEvidence',label='Sample tag source'),col('ObservedSampleTagIssue'),
        col('LastChangedAt','DateTime'),col('ObservedFingerprint'),col('ReconciliationReason'),col('MetadataDifference'),col('ResultCSV','Attachments'),col('ResultSHA256'),col('ResultSourcePath'),
        col('SourceID'),col('Folder'),col('Expno'),col('Instrument'),col('ExperimentType'),
        col('Study','Choice',widgetOptions=json.dumps({'choices':['Electrostatics','DNA','TFE','Other']})),
        col('AnalysisStatus','Choice',widgetOptions=json.dumps({'choices':['Not started','In progress','Analyzed']})),
        col('AcquiredAt','DateTime',widgetOptions=json.dumps(DATE_ONLY_OPTIONS)),col('AcquisitionDateText'),col('AcquisitionStatus','Choice',label='Acquisition status',widgetOptions=json.dumps({'choices':list(ACQUISITION_STATUSES)})),
        col('ObservedAcquisitionEvidence',label='Acquisition status source'),
        col('EffectiveAcquisitionStatus','Text', '"Completed" if $AcquisitionStatus in ("Completed", "Finished", "Done", "Acquired") else ("Aborted" if $AcquisitionStatus in ("Aborted", "Stopped") else ($AcquisitionStatus if $AcquisitionStatus in ("Planned", "Dummy", "Setup", "Running") else "Unknown"))',True,label='Acquisition status'),col('Field_MHz','Numeric'),
        col('Series','Ref:Series'), col('ProteinSelector'), col('RegisteredAt','DateTime',widgetOptions=json.dumps(DATE_ONLY_OPTIONS)),
        col('Sample','Ref:Samples',SAMPLE_TAG_TRIGGER,widgetOptions=json.dumps({'dropdownCondition':{'text':'not $Series or choice.Series == $Series'}})),col('SnapshotSampleID','Text','$Sample.SampleID if $Sample else None'),*snapshot,
        col('Temperature_C','Numeric',label='Temperature (°C)'),col('PulseProgram',label='Pulse program'),*acquisition_columns(),*[col(name,'Choice', '"Not assessed"', label=label, untieColIdFromLabel=True, recalcWhen=0,recalcDeps=['L'],widgetOptions=json.dumps({'choices':[QC_DEFAULT,*choices]})) for name,(label,choices) in QC_CRITERIA.items()],
        col('ImportedAssessment','Choice',widgetOptions=json.dumps({'choices':list(DECISIONS)})),
        col('Assessment','Choice',ASSESSMENT_FORMULA,True,widgetOptions=json.dumps({'choices':list(DECISIONS)})),
        col('AssessmentDisplay','Text',ASSESSMENT_DISPLAY_FORMULA,True,label='Assessment',untieColIdFromLabel=True),
        col('AssessmentExplanation','Text',ASSESSMENT_EXPLANATION_FORMULA,True,label='Assessment details'),
        col('AssessmentReason'),col('Notes'),col('SourceEvidence'),col('MigrationKey'),col('MigrationIssue'),
        col('SourceLabel','Text', SOURCE_LABEL_FORMULA,True,label='EXPNO · Folder',untieColIdFromLabel=True),
        col('AttentionFolder','Text',attention_identity_formula('folder'),True,label='Folder',untieColIdFromLabel=True),
        col('AttentionExpno','Text',attention_identity_formula('expno'),True,label='EXPNO',untieColIdFromLabel=True),
        col('AnalysisDisplay','Text', '"Review identity" if ($ReconciliationReason or ($AcquisitionKey and ($SourceID != $ObservedSourceID or $Folder != $ObservedFolder or $Expno != $ObservedExpno))) else ("Selected CSV" if $ResultCSV else ($AnalysisStatus or "Not recorded"))',True,label='Analysis'),
        col('CoverageConditions','Text', '''parts = []
for label, val, unit in [("Construct", $SnapshotConstruct, ""), ("Label", $SnapshotIsotope, ""), ("Buffer", $SnapshotBuffer, ""), ("pH", $SnapshotpH, ""), ("Protein", $SnapshotProtein_uM, "µM"), ("Notes", $SnapshotOtherConditions, ""), ("Instrument", $EffectiveInstrument, ""), ("Pulse", $EffectivePulseProgram, ""), ("D1", $EffectiveD1_s, " s"), ("Delays", $EffectiveDelayList_s, " s"), ("Spin lock", $EffectiveSpinLock_Hz, " Hz"), ("Offset", $EffectiveSpinLockOffset_Hz, " Hz")]:
  if val is not None and val != "": parts.append(label + ": " + str(val) + unit)
if $CoverageUncertaintyKey: parts.append("Conditions incomplete")
return "; ".join(parts)''',True,label='Conditions'),
        col('CoverageUncertaintyKey','Text', '''axes = ($EffectiveType, $EffectiveInstrument, $EffectiveField_MHz, $EffectiveTemperature_C)
prep = ($SnapshotProtein, $SnapshotConstruct, $SnapshotIsotope, $SnapshotBuffer, $SnapshotpH, $SnapshotNaCl_mM, $SnapshotDNA_uM, $SnapshotTFE_pct, $SnapshotProtein_uM)
if any(x is None or x == "" for x in axes): return "row:" + str($id)
if any(x is None or x == "" for x in prep): return ("sample:" + str($Sample.id)) if $Sample else ("row:" + str($id))
return ""''',True),
        col('EffectiveInstrument','Text', '$Instrument',True,label='Instrument'),
        col('EffectiveType','Text', '$ExperimentType',True,label='Type'),
        col('EffectiveAcquiredAt','DateTime', '$AcquiredAt',True,label='Date',widgetOptions=json.dumps(DATE_ONLY_OPTIONS)),
        col('EffectiveField_MHz','Numeric', '$Field_MHz',True,label='Field (MHz)'),
        col('EffectiveTemperature_C','Numeric', '$Temperature_C if $Temperature_C is None else (round($Temperature_C) if abs($Temperature_C - round($Temperature_C)) <= 0.200000001 else (round($Temperature_C * 2) / 2 if abs($Temperature_C - round($Temperature_C * 2) / 2) <= 0.020000001 else $Temperature_C))',True,label='Temperature group (°C)',untieColIdFromLabel=True),
        col('SampleTagIssue','Text',SAMPLE_TAG_ISSUE_FORMULA,True),
        col('DesktopLink','Text', '"Not linked" if $AcquisitionStatus != "Dummy" and not ($AcquisitionKey and $SourcePath) else ""',True,label='Desktop data',untieColIdFromLabel=True,description='Whether this record is linked to a collected desktop folder. A link does not imply raw FID data or a finished acquisition.'),
        col('AttentionReason','Text',ATTENTION_FORMULA,True),
        col('NeedsAttention','Bool','bool($AttentionReason)',True)],
      'StorageGrid':[col('StorageBox'),col('Row'), *[col('Col'+str(n), 'Text', '", ".join(s.SampleID for s in Samples.lookupRecords(StorageRow=$Row, StorageColumn='+str(n)+'))', True, label=str(n)) for n in range(1,13)]],
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
    tag_deps=('ObservedSampleID','ObservedSampleTagIssue','AcquisitionKey','ObservedSourceID','ObservedFolder','ObservedExpno','SourceID','Folder','Expno','ReconciliationReason')
    actions.append(['UpdateRecord','_grist_Tables_column',sample_ref,{'recalcWhen':0,'recalcDeps':['L',*[column_ids[('Experiments',name)] for name in tag_deps]]}])
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

PRE_DUMMY_EFFECTIVE_FORMULA = next(c["formula"] for c in columns()["Experiments"] if c["id"] == "EffectiveAcquisitionStatus").replace('("Planned", "Dummy", "Setup", "Running")', '("Planned", "Setup", "Running")')
