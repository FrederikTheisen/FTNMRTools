"""Guarded Dummy schema upgrade and native experiment presentation."""
import json
from .schema import columns, PRE_ASSIGNED_SAMPLE_ATTENTION_FORMULA, PRE_LIVE_SAMPLE_ATTENTION_FORMULA, PRE_STUDY_ASSESSMENT_ATTENTION_FORMULA, PRE_CALIBRATION_ATTENTION_FORMULA, PRE_PARTIAL_ATTENTION_FORMULA, PRE_DUMMY_ATTENTION_FORMULA, PRE_DUMMY_EFFECTIVE_FORMULA
from .setup import _mcp

ROW_STYLES = (
    {'formula':'$AcquisitionStatus == "Dummy"','style':{'fillColor':'#F0F1F2','textColor':'#697078'}},
    {'formula':'$AcquisitionStatus != "Dummy" and $NeedsAttention','style':{'fillColor':'#FFF6D8','textColor':'#433B22'}},
)
LINK_STYLE = {'formula':'$DesktopLink == "Not linked"','style':{'textColor':'#456A85'}}


def upgrade_dummy_schema(client, apply=False):
    """Accept only the exact shipped formulas; leave user metadata intact."""
    tables=[r for r in client.records('_grist_Tables') if r['fields'].get('tableId')=='Experiments']
    if len(tables)!=1: raise ValueError('Experiments table missing or ambiguous')
    table=tables[0]['id']
    existing={r['fields'].get('colId'):r for r in client.records('_grist_Tables_column') if r['fields'].get('parentId')==table}
    wanted={c['id']:c for c in columns()['Experiments']}
    actions=[]
    for name,previous in [('AttentionReason',PRE_DUMMY_ATTENTION_FORMULA),('EffectiveAcquisitionStatus',PRE_DUMMY_EFFECTIVE_FORMULA)]:
        current=existing[name]; fields=current['fields'];target=wanted[name]['formula']
        if not fields.get('isFormula') or fields.get('formula') not in ((previous,target,PRE_ASSIGNED_SAMPLE_ATTENTION_FORMULA, PRE_LIVE_SAMPLE_ATTENTION_FORMULA,PRE_STUDY_ASSESSMENT_ATTENTION_FORMULA,PRE_CALIBRATION_ATTENTION_FORMULA,PRE_PARTIAL_ATTENTION_FORMULA) if name=='AttentionReason' else (previous,target)): raise ValueError(name+' formula requires review')
        if fields['formula']!=target:actions.append(['ModifyColumn','Experiments',name,{'formula':target}])
    status=existing['AcquisitionStatus']['fields']
    if status.get('isFormula') or status.get('formula') or status.get('type') not in ('Text','Choice'):raise ValueError('AcquisitionStatus requires review')
    options=json.loads(status.get('widgetOptions') or '{}');choices=options.get('choices',[])
    if not isinstance(choices,list) or any(not isinstance(c,str) for c in choices):raise ValueError('Status choices require review')
    if 'Dummy' not in choices:actions.append(['ModifyColumn','Experiments','AcquisitionStatus',{'type':'Choice','widgetOptions':json.dumps({**options,'choices':choices+['Dummy']})}])
    link=wanted['DesktopLink']
    if 'DesktopLink' not in existing:actions.append(['AddColumn','Experiments','DesktopLink',{k:v for k,v in link.items() if k!='id'}])
    elif any(existing['DesktopLink']['fields'].get(k)!=link[k] for k in ('formula','isFormula','type')):raise ValueError('DesktopLink requires review')
    else:
        metadata={}
        if not existing['DesktopLink']['fields'].get('untieColIdFromLabel'):metadata['untieColIdFromLabel']=True
        if not existing['DesktopLink']['fields'].get('description'):metadata['description']=link['description']
        if metadata:actions.append(['UpdateRecord','_grist_Tables_column',existing['DesktopLink']['id'],metadata])
    if apply and actions:
        client.apply(actions)
        if 'DesktopLink' not in existing:
            # AddColumn may discard presentation metadata. Pin it after the new
            # column has a real metadata ID, without guessing that ID in dry runs.
            added=[r for r in client.records('_grist_Tables_column') if r['fields'].get('parentId')==table and r['fields'].get('colId')=='DesktopLink']
            if len(added)!=1:raise ValueError('New DesktopLink column missing or ambiguous')
            metadata={}
            if not added[0]['fields'].get('untieColIdFromLabel'):metadata['untieColIdFromLabel']=True
            if not added[0]['fields'].get('description'):metadata['description']=link['description']
            if metadata:
                followup=['UpdateRecord','_grist_Tables_column',added[0]['id'],metadata]
                client.apply([followup]);actions.append(followup)
    return actions


def apply_experiment_styles(client,doc_id):
    """Append native rules without changing custom rules or column options."""
    call=lambda name,**kw:_mcp(client,name,{'doc_id':doc_id,**kw})
    result=call('grist_get_table_conditional_styles',table_id='Experiments')
    styles=result.get('styles',result)
    rows=styles.get('rowRules',[]);cols=styles.get('columnRules',{})
    if not isinstance(rows,list) or not isinstance(cols,dict):raise ValueError('Unexpected experiment styles')
    updated=list(rows);updated_cols={**cols}
    for rule in ROW_STYLES:
        matching=[r for r in updated if r.get('formula')==rule['formula']]
        if matching and matching!=[rule]:raise ValueError('Experiment rule requires review')
        if not matching:updated.append(rule)
    links=list(cols.get('DesktopLink',[]))
    matching=[r for r in links if r.get('formula')==LINK_STYLE['formula']]
    if matching and matching!=[LINK_STYLE]:raise ValueError('Desktop link style requires review')
    if not matching:links.append(LINK_STYLE)
    updated_cols['DesktopLink']=links
    target={'rowRules':updated,'columnRules':updated_cols}
    if target!={'rowRules':rows,'columnRules':cols}:call('grist_set_table_conditional_styles',table_id='Experiments',styles=target)
    return target


def editable_experiment_status(client,apply=False):
    """Replace the displayed effective status in experiment grids in place."""
    tables={r['id']:r['fields'].get('tableId') for r in client.records('_grist_Tables')}
    cols={r['fields'].get('colId'):r['id'] for r in client.records('_grist_Tables_column') if tables.get(r['fields'].get('parentId'))=='Experiments'}
    sections={r['id'] for r in client.records('_grist_Views_section') if tables.get(r['fields'].get('tableRef'))=='Experiments' and r['fields'].get('parentKey')=='record'}
    fields=client.records('_grist_Views_section_field');actions=[]
    for sid in sections:
        group=[r for r in fields if r['fields'].get('parentId')==sid]
        if any(r['fields'].get('colRef')==cols['AcquisitionStatus'] for r in group):continue
        for r in group:
            if r['fields'].get('colRef')==cols['EffectiveAcquisitionStatus']:actions.append(['UpdateRecord','_grist_Views_section_field',r['id'],{'colRef':cols['AcquisitionStatus']}])
    if apply and actions:client.apply(actions)
    return actions


def display_partial_assessment(client, apply=False):
    """Show the display formula in human experiment views, preserving field metadata."""
    tables = [r for r in client.records('_grist_Tables') if r['fields'].get('tableId') == 'Experiments']
    if len(tables) != 1:
        raise ValueError('Experiments table missing or ambiguous')
    table_id = tables[0]['id']
    columns = [r for r in client.records('_grist_Tables_column') if r['fields'].get('parentId') == table_id]
    ids = {}
    for name in ('Assessment', 'AssessmentDisplay'):
        matches = [r['id'] for r in columns if r['fields'].get('colId') == name]
        if len(matches) != 1:
            raise ValueError(name + ' column missing or ambiguous')
        ids[name] = matches[0]
    sections = {r['id'] for r in client.records('_grist_Views_section')
                if r['fields'].get('tableRef') == table_id
                and r['fields'].get('parentId', 0) > 0
                and r['fields'].get('parentKey') in ('record', 'single')}
    actions = []
    for row in client.records('_grist_Views_section_field'):
        fields = row['fields']
        if fields.get('parentId') in sections and fields.get('colRef') == ids['Assessment']:
            updates = {'colRef': ids['AssessmentDisplay']}
            width = fields.get('width')
            if isinstance(width, (int, float)) and not isinstance(width, bool) and width < 220:
                updates['width'] = 220
            actions.append(['UpdateRecord', '_grist_Views_section_field', row['id'], updates])
    if apply and actions:
        client.apply(actions)
    return actions


def install_experiment_visibility(client, doc_id, *, apply=False):
    """Install the exact selector after validating the original grid and filter."""
    from .experiment_visibility import experiment_visibility_html
    views=client.records('_grist_Views');sections=client.records('_grist_Views_section')
    tables=client.records('_grist_Tables');columns=client.records('_grist_Tables_column');filters=client.records('_grist_Filters')
    def one(items, message):
        if len(items)!=1:raise ValueError(message)
        return items[0]
    view=one([r for r in views if r['fields'].get('name')=='Experiments'],'Experiments page missing or ambiguous')
    table=one([r for r in tables if r['fields'].get('tableId')=='Experiments'],'Experiments table missing or ambiguous')
    page=view['id'];tid=table['id']
    grid=one([r for r in sections if r['fields'].get('parentId')==page and r['fields'].get('title')=='Experiments' and r['fields'].get('parentKey')=='record' and r['fields'].get('tableRef')==tid],'Experiment grid missing or ambiguous')
    col=one([r for r in columns if r['fields'].get('parentId')==tid and r['fields'].get('colId')=='ExperimentType'],'Experiment type column missing or ambiguous')
    flt=one([r for r in filters if r['fields'].get('viewSectionRef')==grid['id'] and r['fields'].get('colRef')==col['id']],'Experiment type filter missing or ambiguous')
    try:layout=json.loads(view['fields'].get('layoutSpec') or '{}');old=json.loads(flt['fields'].get('filter') or '{}')
    except ValueError:raise ValueError('Invalid experiment layout or filter') from None
    def contains(node, leaf):
        return isinstance(node,dict) and (node.get('leaf')==leaf or any(contains(child,leaf) for child in node.get('children',[]) if isinstance(child,dict)))
    if not isinstance(layout,dict) or not isinstance(layout.get('children'),list) or not contains(layout,grid['id']):raise ValueError('Experiment layout requires review')
    if not isinstance(old,dict) or 'included' in old or not isinstance(old.get('excluded',[]),list) or any(not isinstance(x,str) for x in old.get('excluded',[])):raise ValueError('Experiment type filter requires review')
    def prune(node, wid):
        if not isinstance(node,dict):raise ValueError('Experiment layout requires review')
        if node.get('leaf')==wid:return None
        result=dict(node)
        for key in ('children','collapsed'):
            if key in node:
                if not isinstance(node[key],list):raise ValueError('Experiment layout requires review')
                result[key]=[v for child in node[key] if (v:=prune(child,wid)) is not None]
        if 'children' in result and not result['children'] and 'leaf' not in result:return None
        return result
    def placed(current,wid):
        return bool(current['children'] and current['children'][0].get('leaf')==wid and not any(contains(x,wid) for x in current.get('collapsed',[])))
    def attach(current,wid):
        cleaned=prune(current,wid)
        if not cleaned or not contains(cleaned,grid['id']):raise ValueError('Experiment layout requires review')
        return {**cleaned,'children':[{'leaf':wid,'size':0.9},*cleaned['children']]}
    prior=[r for r in sections if r['fields'].get('parentId')==page and r['fields'].get('title')=='Experiment visibility']
    if len(prior)>1:raise ValueError('Visibility controller is ambiguous')
    newfilter={**old,'excluded':[x for x in old.get('excluded',[]) if x not in ('Calibration','1D')]}
    if prior:
        controller=prior[0];wid=controller['id'];f=controller['fields']
        if (f.get('tableRef')!=tid or f.get('parentKey') not in ('custom','creator')
            or f.get('linkSrcSectionRef') or f.get('linkSrcColRef') or f.get('linkTargetColRef')
            or grid['fields'].get('linkSrcColRef') or grid['fields'].get('linkTargetColRef')
            or any(r['id']!=grid['id'] and r['fields'].get('linkSrcSectionRef')==wid for r in sections)):
            raise ValueError('Existing visibility controller requires review')
        try:options=json.loads(f.get('options') or '{}')
        except ValueError:raise ValueError('Invalid controller options') from None
        html=options.get('html'); shipped=experiment_visibility_html(wid,grid['id']);source=grid['fields'].get('linkSrcSectionRef')
        complete=html==shipped and source==wid
        unfinished=html in (None,'') and not source
        if not (complete or unfinished):raise ValueError('Existing controller differs from shipped control')
        actions=[]
        if unfinished:
            options['html']=shipped
            actions.extend([['UpdateRecord','_grist_Views_section',wid,{'options':json.dumps(options,separators=(',',':'))}],['UpdateRecord','_grist_Views_section',grid['id'],{'linkSrcSectionRef':wid}]])
        if newfilter!=old:actions.append(['UpdateRecord','_grist_Filters',flt['id'],{'filter':json.dumps(newfilter,separators=(',',':'))}])
        if not placed(layout,wid):actions.append(['UpdateRecord','_grist_Views',page,{'layoutSpec':json.dumps(attach(layout,wid),separators=(',',':'))}])
        if apply and actions:client.apply(actions)
        return {'page':page,'grid':grid['id'],'widget':wid,'actions':actions}
    if grid['fields'].get('linkSrcSectionRef') or grid['fields'].get('linkSrcColRef') or grid['fields'].get('linkTargetColRef'):raise ValueError('Existing experiment selector requires review')
    if not apply:return {'page':page,'grid':grid['id'],'filter':flt['id'],'filter_value':newfilter,'original_layout':layout}
    w=_mcp(client,'grist_add_page_widget',{'doc_id':doc_id,'page_id':page,'widget_options':{'table_id':'Experiments','type':'creator','title':'Experiment visibility'}})
    wid=w.get('id',w.get('widget_id'))
    if not isinstance(wid,int) or isinstance(wid,bool) or wid<=0:raise ValueError('Could not create visibility control')
    matches=[r for r in client.records('_grist_Views_section') if r['id']==wid]
    if len(matches)!=1 or matches[0]['fields'].get('parentId')!=page or matches[0]['fields'].get('tableRef')!=tid or matches[0]['fields'].get('parentKey') not in ('custom','creator') or matches[0]['fields'].get('linkSrcSectionRef'):raise ValueError('Created visibility control requires review')
    options=json.loads(matches[0]['fields'].get('options') or '{}')
    if options.get('html') not in (None,''):raise ValueError('Created visibility control is not empty')
    options['html']=experiment_visibility_html(wid,grid['id'])
    newlayout=attach(layout,wid)
    actions=[['UpdateRecord','_grist_Views_section',wid,{'options':json.dumps(options,separators=(',',':'))}],['UpdateRecord','_grist_Views_section',grid['id'],{'linkSrcSectionRef':wid}],['UpdateRecord','_grist_Filters',flt['id'],{'filter':json.dumps(newfilter,separators=(',',':'))}],['UpdateRecord','_grist_Views',page,{'layoutSpec':json.dumps(newlayout,separators=(',',':'))}]]
    client.apply(actions)
    return {'page':page,'grid':grid['id'],'widget':wid,'actions':actions}
