"""Compact Samples page layout; never edits sample records or IDs."""
from .setup import _mcp

DATE_FORMULA = '$PreparedAt.strftime("%Y-%m-%d") if $PreparedAt else ($PreparedDateText or "")'
LEGACY_STORAGE_FORMULA = '''slot = (str($StorageRow) + str($StorageColumn)) if $StorageRow and $StorageColumn else ""
return " / ".join(part for part in [$StorageBox, slot] if part)'''
STORAGE_FORMULA = '(str($StorageRow) + str($StorageColumn)) if $StorageRow and $StorageColumn else ""'

TABLE_FIELDS = ['SampleID','Protein','Isotope','Buffer','NaCl_mM','StorageDisplay','Available']
DETAIL_FIELDS = ['SampleID','PreparationDateDisplay','Protein','Isotope','Protein_uM',
                 'Buffer','pH','NaCl_mM','StorageDisplay','CurrentLocation','Available']
EXTRA_FIELDS = ['SampleID','PreparedAt','Series','Protein','Construct','Isotope',
                'Buffer','pH','NaCl_mM','Protein_uM','DNA_uM','TFE_pct','OtherConditions',
                'StorageColumn','StorageRow','CurrentLocation','Available','Notes','SourceEvidence']
LABELS = {'SampleID':'Sample ID','Isotope':'Labelling','NaCl_mM':'NaCl (mM)',
          'DNA_uM':'DNA (µM)','TFE_pct':'TFE (%)','Protein_uM':'Protein (µM)',
          'Available':'Tube status','CurrentLocation':'Current location','OtherConditions':'Other conditions',
          'PreparedAt':'Preparation date',
          'StorageRow':'Storage column','StorageColumn':'Storage row'}
WIDTHS = {'SampleID':240,'Protein':110,'Isotope':90,'Buffer':150,'NaCl_mM':100,'StorageDisplay':140,'Available':115}


STATUS_ROW_STYLES = (
    {'formula':'$Available == "Suspected degradation"', 'style':{'fillColor':'#FFF1CC','textColor':'#5C4300'}},
    {'formula':'$Available == "Degraded"', 'style':{'fillColor':'#FCE6DF','textColor':'#743125'}},
    {'formula':'$Available == "Dead"', 'style':{'fillColor':'#F8B4B4','textColor':'#7F1D1D'}},
)

def apply_sample_status_styles(client,doc_id):
    """Append native status row styles while retaining existing conditional rules."""
    call=lambda name,**kwargs:_mcp(client,name,{'doc_id':doc_id,**kwargs})
    result=call('grist_get_table_conditional_styles',table_id='Samples')
    styles=result.get('styles',result) if isinstance(result,dict) else None
    if not isinstance(styles,dict):raise ValueError('Unexpected Samples conditional styles')
    rows=styles.get('rowRules',[]);cols=styles.get('columnRules',{})
    if not isinstance(rows,list) or not isinstance(cols,dict):raise ValueError('Unexpected Samples conditional rules')
    for rule in rows:
        if not isinstance(rule,dict) or set(rule)!={'formula','style'} or not isinstance(rule['style'],dict):
            raise ValueError('Existing sample row rule requires review')
    for rules in cols.values():
        if not isinstance(rules,list) or any(not isinstance(rule,dict) or set(rule)!={'formula','style'} or not isinstance(rule['style'],dict) for rule in rules):
            raise ValueError('Existing sample column rule requires review')
    updated=list(rows)
    old_dead={'formula':'$Available == "Dead"','style':{'fillColor':'#EDF0F2','textColor':'#52616B'}}
    dead=[r for r in updated if r['formula']==old_dead['formula']]
    if len(dead)>1 or (dead and dead[0] not in (old_dead,STATUS_ROW_STYLES[-1])):raise ValueError('Existing sample status rule requires review')
    if dead and dead[0]==old_dead:updated[updated.index(dead[0])]=STATUS_ROW_STYLES[-1]
    old={'formula':'$Available == "Retired"','style':{'fillColor':'#EDF0F2','textColor':'#52616B'}}
    previous=[r for r in updated if r['formula']==old['formula']]
    if len(previous)>1 or (previous and previous[0]!=old):raise ValueError('Existing sample status rule requires review')
    if previous:
        current=[r for r in updated if r['formula']==STATUS_ROW_STYLES[-1]['formula']]
        if len(current)>1 or (current and current[0]!=STATUS_ROW_STYLES[-1]):raise ValueError('Existing sample status rule requires review')
        index=updated.index(previous[0])
        if current:updated.pop(index)
        else:updated[index]=STATUS_ROW_STYLES[-1]
    for wanted in STATUS_ROW_STYLES:
        matches=[rule for rule in updated if rule['formula']==wanted['formula']]
        if len(matches)>1 or (matches and matches[0]!=wanted):raise ValueError('Existing sample status rule requires review')
        if not matches:updated.append(wanted)
    if updated!=rows:
        call('grist_set_table_conditional_styles',table_id='Samples',styles={'rowRules':updated,'columnRules':cols})
    return {'rowRules':updated,'columnRules':cols}


def _one(items, key, value):
    matches=[item for item in items if item.get(key)==value]
    if len(matches)!=1: raise ValueError('Missing or ambiguous '+str(value))
    return matches[0]


def unpin_sample_series_filter(client,doc_id):
    """Hide only the pinned Series filter; preserve its filtering condition."""
    call=lambda name,**kwargs:_mcp(client,name,{'doc_id':doc_id,**kwargs})
    page=_one(call('grist_get_pages'),'name','Samples')['id']
    widget=_one(call('grist_get_page_widgets',page_id=page),'title','Samples')
    filters=widget.get('filters',[])
    if not isinstance(filters,list) or any(not isinstance(f,dict) for f in filters):raise ValueError('Unexpected Samples filters')
    updated=[{**f,'pinned':False} if f.get('column_id')=='Series' and f.get('pinned') else f for f in filters]
    if updated!=filters:call('grist_update_page_widget',widget_id=widget['id'],widget_options={'filters':updated})
    return updated


def apply_more_details(client,doc_id):
    """Migrate the former card by ID and arrange its editable fields."""
    call=lambda name,**kwargs:_mcp(client,name,{'doc_id':doc_id,**kwargs})
    page=_one(call('grist_get_pages'),'name','Samples')['id']
    widgets=call('grist_get_page_widgets',page_id=page)
    sample=_one(widgets,'title','Samples')['id']
    matches=[w for w in widgets if w.get('title') in ('Additional conditions','More sample details')]
    if len(matches)>1: raise ValueError('Ambiguous more sample details widget')
    if not matches:
        widget=call('grist_add_page_widget',page_id=page,widget_options={'table_id':'Samples','type':'card','title':'More sample details'})
        widget_id=widget.get('id',widget.get('widget_id'))
    else:
        widget=matches[0];widget_id=widget['id']
        if widget.get('type')!='card' or widget.get('table_id')!='Samples':
            raise ValueError('More sample details widget has an unexpected type or table')
        if widget['title']!='More sample details':
            call('grist_update_page_widget',widget_id=widget_id,widget_options={'title':'More sample details'})
    call('grist_update_table_column',table_id='Samples',column_id='SourceEvidence',column_options={
        'label':'Evidence / source','untie_col_id_from_label':True})
    call('grist_set_page_widget_select_by',widget_id=widget_id,widget_select_by={
        'link_from_widget_id':sample,'link_from_column_id':None,'link_to_column_id':None})
    call('grist_set_widget_fields',widget_id=widget_id,column_ids=EXTRA_FIELDS)
    ids={f['column_id']:f['field_id'] for f in call('grist_get_card_layout',widget_id=widget_id)['fields']}
    groups=[['SampleID'],['PreparedAt','Series'],['Protein','Construct','Isotope'],
            ['Buffer','pH','NaCl_mM'],['Protein_uM','DNA_uM','TFE_pct'],['OtherConditions'],
            ['StorageColumn','StorageRow'],['CurrentLocation','Available'],
            ['Notes'],['SourceEvidence']]
    layout={'children':[{'leaf':ids[row[0]]} if len(row)==1 else {'children':[{'leaf':ids[c]} for c in row]} for row in groups]}
    call('grist_set_card_layout',widget_id=widget_id,layout=layout)
    return widget_id


def apply_samples_ui(client,doc_id):
    """Apply only presentation changes; return IDs for audit and retry.

    Requires an authenticated Grist MCP-capable client. Safe to repeat; existing
    widget IDs and all stored source data are retained. The More sample details
    card is created only if absent.
    """
    call=lambda name,**kwargs:_mcp(client,name,{'doc_id':doc_id,**kwargs})
    page=_one(call('grist_get_pages'),'name','Samples')['id']
    rack=apply_single_rack_ui(client,doc_id)
    metadata=client.records('_grist_Tables')
    sample_table=_one([{'name':r['fields'].get('tableId'),'id':r['id']} for r in metadata],'name','Samples')['id']
    columns={r['fields']['colId']:r for r in client.records('_grist_Tables_column') if r['fields'].get('parentId')==sample_table}
    for name,label,formula in [('PreparationDateDisplay','Preparation date',DATE_FORMULA),('StorageDisplay','Storage',STORAGE_FORMULA)]:
        if name not in columns:
            client.apply([['AddColumn','Samples',name,{'type':'Text','isFormula':True,'formula':formula,'label':label}]])
        elif columns[name]['fields'].get('formula') != formula or not columns[name]['fields'].get('isFormula'):
            raise ValueError('Existing display column has different behavior: '+name)
    for name,label in LABELS.items():
        options={'label':label,'untie_col_id_from_label':True}
        if name=='Available': options['description']='Physical tube condition; Usable is the working default and does not assess an experiment.'
        call('grist_update_table_column',table_id='Samples',column_id=name,column_options=options)
    widgets=call('grist_get_page_widgets',page_id=page)
    samples=_one(widgets,'title','Samples')['id']; detail=_one(widgets,'title','Sample details')['id']
    extra_id=rack['additional_conditions']
    call('grist_set_page_widget_select_by',widget_id=samples,widget_select_by=None)
    call('grist_set_page_widget_select_by',widget_id=detail,widget_select_by={'link_from_widget_id':samples,'link_from_column_id':None,'link_to_column_id':None})
    for widget,fields in ((samples,TABLE_FIELDS),(detail,DETAIL_FIELDS)):
        call('grist_set_widget_fields',widget_id=widget,column_ids=fields)
    sample_widget=_one(widgets,'id',samples)
    filters=sample_widget.get('filters',[])
    if any(f.get('pinned') and f.get('column_id') in ('DNA_uM','TFE_pct','Series') for f in filters):
        updated=[{**f,'pinned':False} if f.get('column_id') in ('DNA_uM','TFE_pct','Series') else f for f in filters]
        call('grist_update_page_widget',widget_id=samples,widget_options={'filters':updated})
    for name,width in WIDTHS.items():
        call('grist_update_widget_field',widget_id=samples,column_id=name,options={'width':width})
    detail_info=call('grist_get_card_layout',widget_id=detail)
    ids={f['column_id']:f['field_id'] for f in detail_info['fields']}
    groups=[['SampleID','PreparationDateDisplay'],['Protein','Isotope','Protein_uM'],
            ['Buffer','pH','NaCl_mM'],['StorageDisplay','CurrentLocation','Available']]
    card_layout={'children':[{'leaf':ids[row[0]]} if len(row)==1 else {'children':[{'leaf':ids[c]} for c in row]} for row in groups]}
    call('grist_set_card_layout',widget_id=detail,layout=card_layout)
    page_info=call('grist_get_page_layout',page_id=page)
    all_ids=[w['id'] for w in page_info['widgets']]
    if samples not in all_ids or detail not in all_ids: raise ValueError('Samples widgets disappeared')
    collapsed=[{'leaf':wid} for wid in all_ids if wid not in (samples,detail)]
    layout={'children':[{'leaf':samples,'size':7},{'leaf':detail,'size':3}],'collapsed':collapsed}
    call('grist_set_page_layout',page_id=page,layout=layout)
    apply_sample_status_styles(client,doc_id)
    return {'page':page,'samples':samples,'detail':detail,'additional_conditions':extra_id,'collapsed':[x['leaf'] for x in collapsed]}


def apply_single_rack_ui(client,doc_id):
    """Refresh the linked rack selector and editable sample details."""
    from .storage_install import apply_storage_widget
    extra=apply_more_details(client,doc_id)
    rack=apply_storage_widget(client,doc_id)
    return {'grid':rack['grid'],'additional_conditions':extra}
