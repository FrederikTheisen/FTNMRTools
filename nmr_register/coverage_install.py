"""Install the Electrostatics matrix without rebuilding unrelated pages."""
from .setup import _mcp


def install_coverage_matrix(client, doc_id):
    from .coverage_widget import coverage_html
    call=lambda name,**kw:_mcp(client,name,{'doc_id':doc_id,**kw})
    pages=[p for p in call('grist_get_pages') if p.get('name')=='Electrostatics']
    if len(pages)!=1:raise ValueError('Electrostatics page missing or ambiguous')
    page=pages[0]['id']; widgets=call('grist_get_page_widgets',page_id=page)
    def one(title):
        found=[w for w in widgets if w.get('title')==title]
        if len(found)!=1:raise ValueError('Missing or ambiguous widget: '+title)
        return found[0]
    acquisitions=one('Individual acquisitions'); details=one('Acquisition details')
    if acquisitions.get('table_id')!='Experiments' or details.get('table_id')!='Experiments':
        raise ValueError('Unexpected acquisition widget table')
    existing=[w for w in widgets if w.get('title')=='Electrostatics matrix']
    legacy=[w for w in widgets if w.get('title')=='Condition coverage']
    if len(existing)>1 or len(legacy)>1:raise ValueError('Ambiguous coverage widget')
    if existing:
        matrix=existing[0]
        if matrix.get('type')!='creator' or matrix.get('table_id')!='Experiments':raise ValueError('Unexpected matrix widget type or table')
    else:
        matrix=call('grist_add_page_widget',page_id=page,widget_options={'table_id':'Experiments','type':'creator','title':'Electrostatics matrix'})
    wid=matrix.get('id',matrix.get('widget_id'))
    if not isinstance(wid,int):raise ValueError('Missing matrix widget ID')
    html=coverage_html(section_id=wid)
    if call('grist_get_creator_html',widget_id=wid).get('html') != html:
        call('grist_set_creator_html',widget_id=wid,html=html)
    if call('grist_get_creator_html',widget_id=wid).get('html') != html:
        raise ValueError('Could not verify installed matrix HTML')
    settings=call('grist_get_custom_widget_settings',widget_id=wid)
    if settings.get('access')!='full':settings=call('grist_set_custom_widget_settings',widget_id=wid,access='full')
    if settings.get('access')!='full' or settings.get('unmapped_required_columns'):
        raise ValueError('Matrix widget settings incomplete')
    choices=call('grist_get_page_widget_select_by_options',widget_id=acquisitions['id'])
    if isinstance(choices,dict):choices=choices.get('options',[])
    matching=[o for o in choices if o.get('link_from_widget_id')==wid and o.get('link_from_column_id') is None and o.get('link_to_column_id') is None]
    if len(matching)>1:raise ValueError('Ambiguous matrix row-ID link')
    link={'link_from_widget_id':wid,'link_from_column_id':None,'link_to_column_id':None}
    sections=client.records('_grist_Views_section')
    source=[r for r in sections if r['id']==wid]
    target=[r for r in sections if r['id']==acquisitions['id']]
    if len(source)!=1 or len(target)!=1:raise ValueError('Matrix or acquisition section missing or ambiguous')
    sf,tf=source[0]['fields'],target[0]['fields']
    tables=[r['id'] for r in client.records('_grist_Tables') if r['fields'].get('tableId')=='Experiments']
    if len(tables)!=1:raise ValueError('Experiments table missing or ambiguous')
    if (sf.get('parentId')!=page or tf.get('parentId')!=page or
        sf.get('parentKey')!='creator' or tf.get('parentKey')!='record' or
        sf.get('title')!='Electrostatics matrix' or tf.get('title')!='Individual acquisitions' or
        sf.get('tableRef')!=tables[0] or tf.get('tableRef')!=tables[0] or
        sf.get('linkSrcSectionRef')):
        raise ValueError('Unexpected matrix or acquisition section metadata')
    if tf.get('linkSrcSectionRef')!=wid or tf.get('linkSrcColRef') or tf.get('linkTargetColRef'):
        if matching:
            result=call('grist_set_page_widget_select_by',widget_id=acquisitions['id'],widget_select_by=link)
            if result.get('matches_by')!='row id':raise ValueError('Acquisition link is not by row ID')
        else:
            # Custom-source eligibility is an in-memory iframe observable, so the
            # server's options omit it until a browser runs grist.ready(). Both
            # sections are verified above as same-page, same-table, non-cyclic.
            client.apply([['UpdateRecord','_grist_Views_section',acquisitions['id'],
                           {'linkSrcSectionRef':wid,'linkSrcColRef':0,'linkTargetColRef':0}]])
        checked=[r for r in client.records('_grist_Views_section') if r['id']==acquisitions['id']]
        if len(checked)!=1 or any(checked[0]['fields'].get(k)!=v for k,v in
                [('linkSrcSectionRef',wid),('linkSrcColRef',0),('linkTargetColRef',0)]):
            raise ValueError('Could not verify acquisition row-ID link')
    if legacy:
        old=legacy[0]['id']
        dependencies=[r['id'] for r in client.records('_grist_Views_section') if r['fields'].get('linkSrcSectionRef')==old]
        if dependencies:raise ValueError('Old coverage widget still has linked sections: '+str(dependencies))
        call('grist_remove_page_widget',widget_id=old)
    current=call('grist_get_page_widgets',page_id=page)
    all_ids=[w['id'] for w in current]
    if len(set(all_ids))!=len(all_ids) or not {wid,acquisitions['id'],details['id']} <= set(all_ids):
        raise ValueError('Unexpected page widgets after matrix installation')
    collapsed=[{'leaf':i} for i in all_ids if i not in (wid,acquisitions['id'])]
    layout={'children':[{'leaf':wid,'size':7},{'leaf':acquisitions['id'],'size':3}],'collapsed':collapsed}
    call('grist_set_page_layout',page_id=page,layout=layout)
    return {'page':page,'coverage':wid,'acquisitions':acquisitions['id'],'details':details['id']}
