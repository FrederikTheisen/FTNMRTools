"""Produce setup actions; caller must provide authenticated Grist transport."""
from .schema import schema_actions, trigger_actions

def create_schema(client):
    """Create tables in an empty document, then configure snapshot triggers.

    client.apply(actions) and client.records(table) are required; intentionally no
    credentials or network handling. Fails if tables already exist.
    """
    existing={r['fields']['tableId']:r['id'] for r in client.records('_grist_Tables')}
    actions=schema_actions(); wanted={a[1] for a in actions}
    metadata=client.records('_grist_Tables_column')
    for action in actions:
        name=action[1]
        if name not in existing: continue
        actual={r['fields']['colId']:r['fields'] for r in metadata if r['fields']['parentId']==existing[name]}
        for expected in action[2]:
            found=actual.get(expected['id'])
            if not found or found.get('type') != expected['type'] or found.get('formula','') != expected['formula'] or bool(found.get('isFormula')) != expected['isFormula']:
                raise ValueError('Existing register schema differs at '+name+'.'+expected['id'])
    missing=[a for a in actions if a[1] not in existing]
    if missing: client.apply(missing)
    tables={r['fields']['tableId']:r['id'] for r in client.records('_grist_Tables')}
    inverse={v:k for k,v in tables.items()}
    cols={(inverse[r['fields']['parentId']],r['fields']['colId']):r['id']
          for r in client.records('_grist_Tables_column') if r['fields']['parentId'] in inverse}
    client.apply(trigger_actions(tables,cols))
    # AddTable creates a primary page for every table. Hide only the pages just
    # created here; the underlying raw views and tables remain available.
    created=[r for r in client.records('_grist_Tables') if r['fields']['tableId'] in wanted]
    primary={r['fields'].get('primaryViewId') for r in created} - {0,None}
    cleanup=[['UpdateRecord','_grist_Views',r['fields']['primaryViewId'],
              {'name':'Register raw '+r['fields']['tableId']}]
             for r in created if r['fields'].get('primaryViewId')]
    for meta in ('_grist_Pages','_grist_TabBar'):
        cleanup.extend(['RemoveRecord',meta,r['id']] for r in client.records(meta)
                       if r['fields'].get('viewRef') in primary)
    if cleanup: client.apply(cleanup)

def _result(client, action):
    """Return a Grist apply endpoint's first action result."""
    result=client.apply([action])
    if isinstance(result, dict): result=result.get('retValues', result.get('results'))
    if not isinstance(result, list) or len(result) != 1 or not isinstance(result[0], dict):
        raise ValueError('Grist apply transport must return retValues for each user action')
    return result[0]

def build_layout(client):
    """Create four pages and native linked widgets in a freshly-created register.

    Requires records(table), add(table, fields), and apply(actions), where apply
    returns Grist's retValues. Refuses existing named pages to avoid duplication.
    The caller is responsible for authenticating the transport and reviewing the
    resulting view in a disposable document before use.
    """
    from .schema import CONDITIONS
    import json
    tables={r['fields']['tableId']:r['id'] for r in client.records('_grist_Tables')}
    byid={v:k for k,v in tables.items()}
    cols={(byid[r['fields']['parentId']],r['fields']['colId']):r['id']
          for r in client.records('_grist_Tables_column') if r['fields']['parentId'] in byid}
    def page(name, table):
        matches=[r['id'] for r in client.records('_grist_Views') if r['fields'].get('name') == name]
        if len(matches)>1: raise ValueError('Ambiguous page: '+name)
        if matches: return matches[0]
        view=_result(client,['AddView',table,'empty',name])['id']
        return view
    def section(view, table, kind='record', grouping=None, title=''):
        sections=client.records('_grist_Views_section')
        matches=[r['id'] for r in sections if r['fields'].get('parentId') == view and r['fields'].get('title') == title]
        if len(matches)>1: raise ValueError('Ambiguous widget: '+title)
        if matches:
            if grouping is not None:
                current=next(r for r in sections if r['id']==matches[0]); ref=current['fields'].get('tableRef')
                existing={r['fields'].get('summarySourceCol') for r in client.records('_grist_Tables_column')
                          if r['fields'].get('parentId')==ref and r['fields'].get('summarySourceCol')}
                wanted={cols[(table,name)] for name in grouping}
                if existing != wanted:
                    client.apply([['UpdateSummaryViewSection',matches[0],[cols[(table,name)] for name in grouping]]])
            return matches[0]
        # Recover a section created immediately before an interrupted title update.
        candidates=[r for r in sections if r['fields'].get('parentId') == view
                    and not r['fields'].get('title') and r['fields'].get('parentKey') == kind]
        if grouping is None:
            unmatched=[r['id'] for r in candidates if r['fields'].get('tableRef') == tables[table]]
        else:
            source_tables={r['id']:r['fields'].get('summarySourceTable') for r in client.records('_grist_Tables')}
            all_columns=client.records('_grist_Tables_column')
            wanted={cols[(table,c)] for c in grouping}
            unmatched=[]
            for r in candidates:
                ref=r['fields'].get('tableRef')
                refs={c['fields'].get('summarySourceCol') for c in all_columns
                      if c['fields'].get('parentId') == ref and c['fields'].get('summarySourceCol')}
                if source_tables.get(ref) == tables[table] and refs == wanted: unmatched.append(r['id'])
        if len(unmatched)>1: raise ValueError('Ambiguous unnamed widget: '+title)
        if unmatched:
            update('_grist_Views_section',unmatched[0],{'title':title})
            return unmatched[0]
        refs=None if grouping is None else [cols[(table,c)] for c in grouping]
        result=_result(client,['CreateViewSection',tables[table],view,kind,refs,None])['sectionRef']
        update('_grist_Views_section',result,{'title':title})
        return result
    def update(table,row,fields):
        client.apply([['UpdateRecord',table,row,fields]])
    def link(target,source,targetcol=None):
        update('_grist_Views_section',target,{'linkSrcSectionRef':source,'linkTargetColRef':targetcol or 0})
    def arrange(view, sections):
        # Keep already installed registration forms when setup is rerun.
        forms=[r['id'] for r in client.records('_grist_Views_section')
               if r['fields'].get('parentId')==view and r['fields'].get('title') in ('Register samples','Register experiment')]
        all_sections=list(sections)+[ref for ref in forms if ref not in sections]
        update('_grist_Views',view,{'layoutSpec':json.dumps({'children':[{'leaf':s} for s in all_sections]})})
    def fields(section_id, table, names):
        wanted=[cols[(table,name)] for name in names]
        current=[r for r in client.records('_grist_Views_section_field') if r['fields'].get('parentId')==section_id]
        bycol={r['fields'].get('colRef'):r['id'] for r in current}
        actions=[['RemoveRecord','_grist_Views_section_field',r['id']] for r in current if r['fields'].get('colRef') not in wanted]
        for i,colref in enumerate(wanted,1):
            values={'parentPos':i,'width':150}
            if colref in bycol: actions.append(['UpdateRecord','_grist_Views_section_field',bycol[colref],values])
            else: actions.append(['AddRecord','_grist_Views_section_field',None,dict(values,parentId=section_id,colRef=colref)])
        if actions: client.apply(actions)
    def pins(section_id, table, names):
        current={(r['fields'].get('viewSectionRef'),r['fields'].get('colRef')) for r in client.records('_grist_Filters')}
        actions=[]
        for name in names:
            ref=cols[(table,name)]
            if (section_id,ref) not in current:
                actions.append(['AddRecord','_grist_Filters',None,{'viewSectionRef':section_id,'colRef':ref,'filter':'{}','pinned':True}])
        if actions: client.apply(actions)
    experiment_columns=['EffectiveAcquiredAt','SnapshotProtein','Sample','EffectiveType','EffectiveInstrument','EffectiveField_MHz','EffectiveTemperature_C','SnapshotNaCl_mM','Assessment','ResultCSV']
    experiment_detail=['SourceID','Folder','Expno','RegisteredAt','Series','ProteinSelector','Sample','Instrument','ExperimentType','AcquiredAt','AcquisitionDateText','AcquisitionStatus','Field_MHz','Temperature_C','PulseProgram','SnapshotProtein','SnapshotConstruct','SnapshotBuffer','SnapshotpH','SnapshotNaCl_mM','SnapshotDNA_uM','SnapshotTFE_pct','SnapshotProtein_uM','SnapshotOtherConditions','Assessment','AssessmentReason','Notes','AttentionReason','ResultCSV']
    sample_columns=['SampleID','Protein','Series','NaCl_mM','DNA_uM','TFE_pct','StorageBox','StorageRow','StorageColumn','Available']
    sample_detail=['SampleID','PreparedAt','PreparedDateText','Isotope','Protein','Construct','DistinguishingCondition','SampleNumber','Series','Buffer','pH','NaCl_mM','DNA_uM','TFE_pct','Protein_uM','OtherConditions','StorageBox','StorageRow','StorageColumn','Available','CurrentLocation','Notes']
    v=page('Experiments','Experiments'); e=section(v,'Experiments',title='Experiments'); card=section(v,'Experiments','single',title='Experiment details'); link(card,e); arrange(v,[e,card])
    fields(e,'Experiments',experiment_columns); fields(card,'Experiments',experiment_detail)
    pins(e,'Experiments',['Sample','SnapshotProtein','EffectiveType','EffectiveInstrument','EffectiveField_MHz','EffectiveTemperature_C','SnapshotNaCl_mM','SnapshotDNA_uM','SnapshotTFE_pct','Assessment'])
    stale=[r['id'] for r in client.records('_grist_Filters') if r['fields'].get('viewSectionRef')==e and r['fields'].get('colRef')==cols[('Experiments','Temperature_C')]]
    if stale: client.apply([['RemoveRecord','_grist_Filters',i] for i in stale])
    v=page('Samples','Samples'); series=section(v,'Series',title='Series'); sample=section(v,'Samples',title='Samples'); detail=section(v,'Samples','single',title='Sample details'); exp=section(v,'Experiments',title='Sample experiments'); grid=section(v,'StorageGrid',title='Storage grid')
    link(sample,series,cols[('Samples','Series')]); link(detail,sample); link(exp,sample,cols[('Experiments','Sample')]); arrange(v,[series,sample,detail,exp,grid])
    fields(series,'Series',['SeriesID','Protein','Construct','Isotope','Buffer','pH','NaCl_mM','DNA_uM','TFE_pct','Protein_uM','OtherConditions','Notes'])
    fields(sample,'Samples',sample_columns); fields(detail,'Samples',sample_detail); fields(exp,'Experiments',experiment_columns)
    fields(grid,'StorageGrid',['StorageBox','Row']+['Col'+str(i) for i in range(1,13)])
    pins(sample,'Samples',['Protein','Series','NaCl_mM','DNA_uM','TFE_pct','Available'])
    v=page('Needs attention','Experiments'); attention=section(v,'Experiments',title='Needs attention'); card=section(v,'Experiments','single',title='Attention details'); scanner=section(v,'ScanStatus','single',title='Scanner status'); link(card,attention)
    fields(attention,'Experiments',experiment_columns+['AttentionReason']); fields(card,'Experiments',experiment_detail+['MigrationIssue'])
    fields(scanner,'ScanStatus',['Scanner','LastScanAt','LastSyncAt','LastError','ScannedRoot'])
    filters=[r for r in client.records('_grist_Filters') if r['fields'].get('viewSectionRef')==attention and r['fields'].get('colRef')==cols[('Experiments','NeedsAttention')]]
    if not filters: client.apply([['AddRecord','_grist_Filters',None,{'viewSectionRef':attention,'colRef':cols[('Experiments','NeedsAttention')], 'filter':json.dumps({'included':[True]}),'pinned':True}]])
    arrange(v,[attention,card,scanner])
    v=page('Overview','Experiments'); summary=section(v,'Experiments','record',['SnapshotProtein','SnapshotConstruct','SnapshotBuffer','SnapshotpH','SnapshotNaCl_mM','SnapshotDNA_uM','SnapshotTFE_pct','SnapshotProtein_uM','EffectiveTemperature_C','EffectiveField_MHz'],title='Conditions and counts'); arrange(v,[summary])
    summary_record=next(r for r in client.records('_grist_Views_section') if r['id']==summary)
    summary_ref=summary_record['fields']['tableRef']
    summary_table=next(r['fields']['tableId'] for r in client.records('_grist_Tables') if r['id']==summary_ref)
    aggregate_formulas={
        'UseCount':'len([r for r in $group if r.Assessment == "Use"])',
        'LimitedCount':'len([r for r in $group if r.Assessment == "Limited use"])',
        'ExcludedCount':'len([r for r in $group if r.Assessment == "Do not use"])',
        'UnassessedCount':'len([r for r in $group if not r.Assessment])',
        'ResultsCount':'len([r for r in $group if r.ResultCSV])',
    }
    present={r['fields']['colId'] for r in client.records('_grist_Tables_column') if r['fields'].get('parentId')==summary_ref}
    for name,formula in aggregate_formulas.items():
        if name not in present:
            client.apply([['AddVisibleColumn',summary_table,name,{'type':'Int','isFormula':True,'formula':formula}]])
    for r in client.records('_grist_Tables_column'):
        if r['fields'].get('parentId')==summary_ref: cols[(summary_table,r['fields']['colId'])]=r['id']
    fields(summary,summary_table,['SnapshotProtein','SnapshotConstruct','SnapshotBuffer','SnapshotpH','SnapshotNaCl_mM','SnapshotDNA_uM','SnapshotTFE_pct','SnapshotProtein_uM','EffectiveTemperature_C','EffectiveField_MHz','count',*aggregate_formulas])
    # Eight rows form one display grid. The box can be edited on each row or
    # copied to create another box; formula cells resolve the actual samples.
    existing={(r['fields'].get('StorageBox'),r['fields'].get('Row')) for r in client.records('StorageGrid')}
    boxes={'Unspecified box'} | {r['fields'].get('StorageBox') for r in client.records('Samples') if r['fields'].get('StorageBox')}
    for box in sorted(boxes):
        for row in 'ABCDEFGH':
            if (box,row) not in existing: client.add('StorageGrid',{'StorageBox':box,'Row':row})

def _mcp(client, name, arguments):
    """Call a Grist MCP tool and decode its JSON text result, failing closed."""
    import json
    result=client.mcp_call(name,arguments)
    if result.get('isError'):
        raise ValueError('Grist MCP tool failed: '+name)
    blocks=[block.get('text','') for block in result.get('content',[]) if block.get('type')=='text']
    if len(blocks)!=1: raise ValueError('Unexpected Grist MCP response: '+name)
    try: return json.loads(blocks[0])
    except (TypeError,ValueError): raise ValueError('Unexpected Grist MCP result: '+name) from None

def build_registration_widgets(client,doc_id):
    """Create/update the two document-owned forms; safe to rerun by title.

    These forms need full document access to atomically create linked records.
    """
    from .registration_widget import registration_html
    pages=_mcp(client,'grist_get_pages',{'doc_id':doc_id})
    result={}
    for page_name,table,mode,title in (
        ('Samples','Samples','sample','Register samples'),
        ('Experiments','Experiments','experiment','Register experiment')):
        matches=[p for p in pages if p.get('name')==page_name]
        if len(matches)!=1: raise ValueError('Page is missing or ambiguous: '+page_name)
        page_id=matches[0]['id']
        widgets=_mcp(client,'grist_get_page_widgets',{'doc_id':doc_id,'page_id':page_id})
        matches=[w for w in widgets if w.get('title')==title]
        if len(matches)>1: raise ValueError('Registration widget is ambiguous: '+title)
        if matches:
            w=matches[0]
            if w.get('type')!='creator' or w.get('table_id')!=table:
                raise ValueError('Existing registration widget has unexpected type or table: '+title)
        else:
            w=_mcp(client,'grist_add_page_widget',{'doc_id':doc_id,'page_id':page_id,
                'widget_options':{'table_id':table,'type':'creator','title':title}})
        wid=w['id'] if 'id' in w else w['widget_id']
        html=registration_html(mode)
        current=_mcp(client,'grist_get_creator_html',{'doc_id':doc_id,'widget_id':wid})
        if current.get('html') != html:
            _mcp(client,'grist_set_creator_html',{'doc_id':doc_id,'widget_id':wid,'html':html})
        _mcp(client,'grist_set_custom_widget_settings',{'doc_id':doc_id,'widget_id':wid,'access':'full'})
        result[mode]=wid
    return result
