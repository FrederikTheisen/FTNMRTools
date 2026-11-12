"""Build an explicit, stable legacy review snapshot without modifying the input."""
import csv
import io
from .migration import read_reviews

TO_LEGACY={'Use':'usable','Limited use':'review','Do not use':'exclude'}

def review_snapshot(source_path, experiments, source_id=None):
    """Return TSV bytes; exact Folder/Expno identity only; unknown rows preserved.

    Never use AcquisitionKey as a guess for a legacy target. Legacy decisions without matching experiments block export. Duplicate identities or incomplete decisions fail closed.
    """
    rows=read_reviews(source_path)
    with open(source_path,encoding='utf-8-sig') as stream:
        source=stream.read()
    lines=source.splitlines(keepends=True)
    index=next(i for i,line in enumerate(lines) if line.startswith('target\t'))
    header=lines[index].rstrip('\r\n').split('\t')
    seen=set()
    for item in experiments:
        fields=item.get('fields',item)
        row_source=(fields.get('SourceID') or '').strip()
        if row_source and not source_id:
            raise ValueError('A source_id is required for source-specific reviews')
        if source_id and row_source != source_id:
            continue
        folder=str(fields.get('Folder') or '').strip(); expno=str(fields.get('Expno') or '').strip()
        if not folder or not expno: continue
        target=folder+'/'+expno
        if target in seen: raise ValueError('Duplicate experiment target: '+target)
        seen.add(target)
        if target not in rows: continue
        decision=fields.get('Assessment') or ''; reason=(fields.get('AssessmentReason') or '').strip()
        if decision in ('Limited use','Do not use') and not reason:
            reason=(fields.get('AssessmentExplanation') or '').strip()
        if fields.get('ReconciliationReason') or fields.get('MigrationIssue') or fields.get('MetadataDifference'):
            raise ValueError('Unresolved experiment cannot be exported: '+target)
        if decision and (decision not in TO_LEGACY or (decision != 'Use' and not reason)):
            raise ValueError('Invalid or incomplete assessment: '+target)
        rows[target]['assessment']=TO_LEGACY.get(decision,'')
        rows[target]['comment']=reason
    unresolved=[target for target,row in rows.items() if target not in seen and ((row.get('assessment') or '').strip() or (row.get('comment') or '').strip())]
    if unresolved:
        raise ValueError('Legacy decisions without a Grist experiment require reconciliation: '+', '.join(unresolved))
    out=io.StringIO(newline=''); out.writelines(lines[:index]); writer=csv.DictWriter(out,fieldnames=header,delimiter='\t',lineterminator='\n',extrasaction='ignore')
    writer.writeheader(); writer.writerows(rows.values())
    return out.getvalue().encode('utf-8')
