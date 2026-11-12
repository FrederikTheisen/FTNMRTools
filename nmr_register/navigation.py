"""Narrow native link from Electrostatics acquisitions to the Experiments grid."""
import json


def apply_experiment_navigation(client, *, dry_run=True):
    """Add a row-specific native hyperlink without changing layouts or filters.

    Grist assigns a column metadata ID when a column is created, so an initial
    install uses two small transactions. Rerunning safely resumes either phase.
    Unknown existing columns or ambiguous widgets fail closed for review.
    """
    tables=client.records('_grist_Tables');views=client.records('_grist_Views')
    sections=client.records('_grist_Views_section');metadata=client.records('_grist_Tables_column')
    ids=[r['id'] for r in tables if r['fields'].get('tableId')=='Experiments']
    if len(ids)!=1:raise ValueError('Experiments table missing or ambiguous')
    table=ids[0]
    def one_view(name):
        rows=[r for r in views if r['fields'].get('name')==name]
        if len(rows)!=1:raise ValueError('Page missing or ambiguous: '+name)
        return rows[0]['id']
    def one_grid(view,title):
        rows=[r for r in sections if r['fields'].get('parentId')==view and r['fields'].get('title')==title and r['fields'].get('parentKey')=='record']
        if len(rows)!=1 or rows[0]['fields'].get('tableRef')!=table:raise ValueError('Unexpected or ambiguous grid: '+title)
        return rows[0]
    page=one_view('Experiments');electro=one_view('Electrostatics')
    grid=one_grid(page,'Experiments'); acquisitions=one_grid(electro,'Individual acquisitions')
    if grid['fields'].get('linkSrcSectionRef'):raise ValueError('Experiments grid is linked; row navigation requires review')
    formula='SELF_HYPERLINK(label="Open", page='+str(page)+') + "#a1.s'+str(grid['id'])+'.r" + str($id)'
    found=[r for r in metadata if r['fields'].get('parentId')==table and r['fields'].get('colId')=='OpenExperiment']
    if len(found)>1:raise ValueError('Duplicate OpenExperiment column')
    actions=[]
    if not found:
        add=['AddColumn','Experiments','OpenExperiment',{'type':'Text','isFormula':True,'formula':formula,'label':'Open experiment','widgetOptions':json.dumps({'widget':'HyperLink'})}]
        actions.append(add)
        if dry_run:return {'actions':actions,'remaining':'Add the link field after column creation','page':page,'grid':grid['id'],'acquisitions':acquisitions['id']}
        client.apply([add])
        found=[r for r in client.records('_grist_Tables_column') if r['fields'].get('parentId')==table and r['fields'].get('colId')=='OpenExperiment']
        if len(found)!=1:raise ValueError('Could not verify link column; rerun after review')
    col=found[0]; f=col['fields']
    try:options=json.loads(f.get('widgetOptions') or '{}')
    except ValueError:raise ValueError('Invalid link column options') from None
    if f.get('type')!='Text' or not f.get('isFormula') or f.get('formula')!=formula or options.get('widget')!='HyperLink':raise ValueError('Existing link column differs; review')
    fields=client.records('_grist_Views_section_field');current=[r for r in fields if r['fields'].get('parentId')==acquisitions['id']]
    matches=[r for r in current if r['fields'].get('colRef')==col['id']]
    if len(matches)>1:raise ValueError('Duplicate link field in acquisitions grid')
    if not matches:
        positions=[r['fields'].get('parentPos') for r in current if isinstance(r['fields'].get('parentPos'),(int,float))]
        field=['AddRecord','_grist_Views_section_field',None,{'parentId':acquisitions['id'],'colRef':col['id'],'parentPos':min(positions)-1 if positions else 1,'width':85}]
        actions.append(field)
        if not dry_run:client.apply([field])
    return {'actions':actions,'page':page,'grid':grid['id'],'acquisitions':acquisitions['id'],'column':col['id']}
