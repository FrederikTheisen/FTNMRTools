"""Install the single-rack creator in the existing Samples page section."""
import copy
import json


def _one(items, predicate, description):
    found=[r for r in items if predicate(r)]
    if len(found)!=1:raise ValueError(description+' missing or ambiguous')
    return found[0]


def plan_storage_widget(client, html):
    """Return narrowly scoped metadata actions without touching sample data."""
    if not isinstance(html,str) or 'data-nmr-storage="1"' not in html:raise ValueError('Storage widget HTML is missing or invalid')
    tables=client.records('_grist_Tables');views=client.records('_grist_Views')
    columns=client.records('_grist_Tables_column');sections=client.records('_grist_Views_section')
    fields=client.records('_grist_Views_section_field')
    sample=_one(tables,lambda r:r['fields'].get('tableId')=='Samples','Samples table')['id']
    old_table=_one(tables,lambda r:r['fields'].get('tableId')=='StorageGrid','StorageGrid table')['id']
    page=_one(views,lambda r:r['fields'].get('name')=='Samples','Samples page')['id']
    grid=_one(sections,lambda r:r['fields'].get('parentId')==page and r['fields'].get('title')=='Samples','Samples grid')
    rack=_one(sections,lambda r:r['fields'].get('parentId')==page and r['fields'].get('title')=='Storage grid','Storage grid')
    card=_one(sections,lambda r:r['fields'].get('parentId')==page and r['fields'].get('title')=='More sample details','More sample details')
    if grid['fields'].get('tableRef')!=sample or grid['fields'].get('parentKey')!='record':raise ValueError('Unexpected Samples grid')
    rf=rack['fields'];native=(rf.get('parentKey')=='record' and rf.get('tableRef')==old_table)
    creator=(rf.get('parentKey') in ('creator','custom') and rf.get('tableRef')==sample)
    if not (native or creator):raise ValueError('Unexpected storage widget type or table')
    if any(s['id']!=rack['id'] and s['fields'].get('linkSrcSectionRef')==rack['id'] for s in sections):raise ValueError('Storage widget controls another widget')
    if native and rf.get('linkSrcSectionRef') not in (None,0):raise ValueError('Legacy storage widget has a custom link')
    if creator and rf.get('linkSrcSectionRef') not in (None,0,grid['id']):raise ValueError('Storage widget has a custom link')
    try:
        options=json.loads(rf.get('options') or '{}')
        if not isinstance(options,dict):raise ValueError()
        custom=json.loads(options.get('customView') or '{}')
        if not isinstance(custom,dict):raise ValueError()
    except (ValueError,TypeError):raise ValueError('Invalid storage widget options') from None
    # Never replace an unrelated creator installed by the user.
    if creator and options.get('html') and (not isinstance(options['html'],str) or 'data-nmr-storage="1"' not in options['html']):
        raise ValueError('Existing storage creator requires review')
    options=copy.deepcopy(options);custom=copy.deepcopy(custom);custom['access']='full'
    options['customView']=json.dumps(custom,separators=(',',':'));options['html']=html
    changes={}
    for key,value in {'parentKey':'creator','tableRef':sample,'linkSrcSectionRef':grid['id'],'linkSrcColRef':0,'linkTargetColRef':0}.items():
        if (rf.get(key) or 0)!=(value or 0):changes[key]=value
    if json.loads(rf.get('options') or '{}')!=options:changes['options']=json.dumps(options,separators=(',',':'))
    actions=[]
    if native:
        old_columns={r['id'] for r in columns if r['fields'].get('parentId')==old_table}
        for f in fields:
            if f['fields'].get('parentId')==rack['id']:
                if f['fields'].get('colRef') not in old_columns:raise ValueError('Storage section contains unfamiliar fields')
                actions.append(['RemoveRecord','_grist_Views_section_field',f['id']])
    if changes:actions.append(['UpdateRecord','_grist_Views_section',rack['id'],changes])
    colrefs={}
    for name,label in (('StorageRow','Storage column'),('StorageColumn','Storage row')):
        col=_one(columns,lambda r:r['fields'].get('parentId')==sample and r['fields'].get('colId')==name,name)
        colrefs[name]=col['id'];f=col['fields']
        if f.get('label')!=label or not f.get('untieColIdFromLabel'):
            actions.append(['UpdateRecord','_grist_Tables_column',col['id'],{'label':label,'untieColIdFromLabel':True}])
    field_ids={}
    for name,col in colrefs.items():
        field=_one(fields,lambda r:r['fields'].get('parentId')==card['id'] and r['fields'].get('colRef')==col,name+' detail field')
        field_ids[name]=field['id']
    try:layout=json.loads(card['fields'].get('layoutSpec') or '{}')
    except (TypeError,ValueError):raise ValueError('Invalid sample details layout') from None
    if not isinstance(layout,dict):raise ValueError('Invalid sample details layout')
    before=copy.deepcopy(layout);pairs=[]
    def reorder(node):
        if isinstance(node,dict):
            children=node.get('children')
            if isinstance(children,list):
                a=[i for i,c in enumerate(children) if isinstance(c,dict) and c.get('leaf')==field_ids['StorageRow']]
                b=[i for i,c in enumerate(children) if isinstance(c,dict) and c.get('leaf')==field_ids['StorageColumn']]
                if a and b:
                    pairs.append((a[0],b[0]))
                    if b[0]>a[0]:children[a[0]],children[b[0]]=children[b[0]],children[a[0]]
                for c in children:reorder(c)
        elif isinstance(node,list):
            for c in node:reorder(c)
    reorder(layout)
    if len(pairs)!=1:raise ValueError('Storage detail fields have an unfamiliar layout')
    if layout!=before:actions.append(['UpdateRecord','_grist_Views_section',card['id'],{'layoutSpec':json.dumps(layout,separators=(',',':'))}])
    return {'page':page,'grid':rack['id'],'samples':grid['id'],'actions':actions}


def apply_storage_widget(client,doc_id,*,dry_run=False):
    """Install or refresh the same creator section; preserve the other widgets."""
    from .storage_widget import storage_html
    sections=client.records('_grist_Views_section');views=client.records('_grist_Views')
    page=_one(views,lambda r:r['fields'].get('name')=='Samples','Samples page')['id']
    grid=_one(sections,lambda r:r['fields'].get('parentId')==page and r['fields'].get('title')=='Storage grid','Storage grid')['id']
    plan=plan_storage_widget(client,storage_html(section_id=grid))
    if not dry_run and plan['actions']:client.apply(plan['actions'])
    return plan
