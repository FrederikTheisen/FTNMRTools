"""Guarded Dummy schema upgrade and native experiment presentation."""
import json
from .schema import columns, PRE_CATEGORY_ATTENTION_FORMULA, PRE_ASSIGNED_SAMPLE_ATTENTION_FORMULA, PRE_LIVE_SAMPLE_ATTENTION_FORMULA, PRE_STUDY_ASSESSMENT_ATTENTION_FORMULA, PRE_CALIBRATION_ATTENTION_FORMULA, PRE_PARTIAL_ATTENTION_FORMULA, PRE_DUMMY_ATTENTION_FORMULA, PRE_DUMMY_EFFECTIVE_FORMULA, PRE_REGISTRATION_REMOVAL_ATTENTION_FORMULA
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
        if not fields.get('isFormula') or fields.get('formula') not in ((previous,target,PRE_CATEGORY_ATTENTION_FORMULA,PRE_ASSIGNED_SAMPLE_ATTENTION_FORMULA, PRE_LIVE_SAMPLE_ATTENTION_FORMULA,PRE_STUDY_ASSESSMENT_ATTENTION_FORMULA,PRE_CALIBRATION_ATTENTION_FORMULA,PRE_PARTIAL_ATTENTION_FORMULA,PRE_REGISTRATION_REMOVAL_ATTENTION_FORMULA) if name=='AttentionReason' else (previous,target)): raise ValueError(name+' formula requires review')
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


def configure_experiment_type_filter(client, apply=False):
    """Pin the broad experiment type filter and hide calibration / 1D."""
    def one(rows, message):
        if len(rows) != 1:
            raise ValueError(message)
        return rows[0]

    view = one([r for r in client.records('_grist_Views')
                if r['fields'].get('name') == 'Experiments'],
               'Experiments page missing or ambiguous')
    table = one([r for r in client.records('_grist_Tables')
                 if r['fields'].get('tableId') == 'Experiments'],
                'Experiments table missing or ambiguous')
    tid = table['id']
    grid = one([r for r in client.records('_grist_Views_section')
                if r['fields'].get('parentId') == view['id']
                and r['fields'].get('title') == 'Experiments'
                and r['fields'].get('parentKey') == 'record'
                and r['fields'].get('tableRef') == tid],
               'Experiment grid missing or ambiguous')
    if any(grid['fields'].get(k) for k in
           ('linkSrcSectionRef', 'linkSrcColRef', 'linkTargetColRef')):
        raise ValueError('Experiment grid is linked; review before configuring filter')
    columns = client.records('_grist_Tables_column')
    category_columns = [r for r in columns if r['fields'].get('parentId') == tid
                        and r['fields'].get('colId') == 'ExperimentCategory']
    precise = one([r for r in columns if r['fields'].get('parentId') == tid
                   and r['fields'].get('colId') == 'ExperimentType'],
                  'Experiment column missing or ambiguous')
    if category_columns:
        broad = one(category_columns, 'Experiment type column ambiguous')
        targets = [(broad, ('Calibration',)), (precise, ('1D',))]
    else:
        # Pre-taxonomy documents retain their previous native filter behavior.
        targets = [(precise, ('Calibration', '1D'))]
    all_filters = client.records('_grist_Filters')
    aliases = [r for r in columns
               if r['fields'].get('parentId') == tid
               and r['fields'].get('colId') == 'EffectiveType']
    alias_filters = [r for r in all_filters
                     if r['fields'].get('viewSectionRef') == grid['id']
                     and r['fields'].get('colRef') in {c['id'] for c in aliases}]
    if alias_filters:
        if len(aliases) != 1 or not aliases[0]['fields'].get('isFormula') or aliases[0]['fields'].get('formula') != '$ExperimentType':
            raise ValueError('Legacy experiment type column requires review')
        for row in alias_filters:
            try:
                condition = json.loads(row['fields'].get('filter') or '{}')
            except (TypeError, ValueError):
                raise ValueError('Legacy experiment type filter requires review') from None
            if condition not in ({}, {'excluded': []}):
                raise ValueError('Legacy experiment type filter has restrictions; review before removing')
    actions = [['RemoveRecord', '_grist_Filters', r['id']] for r in alias_filters]
    filter_values = {}
    for col, hidden in targets:
        filters = [r for r in all_filters
                   if r['fields'].get('viewSectionRef') == grid['id']
                   and r['fields'].get('colRef') == col['id']]
        if len(filters) > 1:
            raise ValueError('Experiment type filter ambiguous')
        try:
            current = json.loads(filters[0]['fields'].get('filter') or '{}') if filters else {}
        except (TypeError, ValueError):
            raise ValueError('Invalid experiment type filter') from None
        if not isinstance(current, dict) or ('included' in current and 'excluded' in current):
            raise ValueError('Experiment type filter requires review')
        for key in ('included', 'excluded'):
            if key in current and (not isinstance(current[key], list)
                                   or any(not isinstance(x, str) for x in current[key])):
                raise ValueError('Experiment type filter requires review')
        target = dict(current)
        if 'included' in target:
            target['included'] = [x for x in target['included'] if x not in hidden]
        else:
            excluded = list(target.get('excluded', []))
            target['excluded'] = excluded + [x for x in hidden if x not in excluded]
        if filters:
            updates = {}
            if target != current:
                updates['filter'] = json.dumps(target, separators=(',', ':'))
            if not filters[0]['fields'].get('pinned'):
                updates['pinned'] = True
            if updates:
                actions.append(['UpdateRecord', '_grist_Filters', filters[0]['id'], updates])
        else:
            actions.append(['AddRecord', '_grist_Filters', None,
                            {'viewSectionRef': grid['id'], 'colRef': col['id'],
                             'filter': json.dumps(target, separators=(',', ':')), 'pinned': True}])
        filter_values[col['fields']['colId']] = target
    if apply and actions:
        client.apply(actions)
    return {'page': view['id'], 'grid': grid['id'], 'filter_value': filter_values.get('ExperimentCategory', target), 'filter_values': filter_values, 'actions': actions}


def remove_registration_timestamp(client, apply=False):
    """Remove the retired column and its shipped warning without resetting views.

    The caller refreshes embedded forms in the same action batch. Unknown
    formula dependencies are rejected before any document changes.
    """
    import re
    from .schema import ATTENTION_FORMULA
    tables = [r for r in client.records('_grist_Tables') if r['fields'].get('tableId') == 'Experiments']
    if len(tables) != 1: raise ValueError('Experiments table missing or ambiguous')
    columns = client.records('_grist_Tables_column')
    retired = [r for r in columns if r['fields'].get('parentId') == tables[0]['id'] and r['fields'].get('colId') == 'RegisteredAt']
    if len(retired) > 1: raise ValueError('Registration timestamp column is ambiguous')
    actions = []
    for row in columns:
        f = row['fields']
        if re.search(r'\bRegisteredAt\b', f.get('formula') or ''):
            if (f.get('parentId') != tables[0]['id'] or f.get('colId') != 'AttentionReason'
                    or not f.get('isFormula') or f['formula'] != PRE_REGISTRATION_REMOVAL_ATTENTION_FORMULA):
                raise ValueError('Custom registration timestamp dependency requires review')
            actions.append(['UpdateRecord', '_grist_Tables_column', row['id'], {'formula': ATTENTION_FORMULA}])
    if retired:
        column = retired[0]; ref = column['id']; f = column['fields']
        if f.get('isFormula') or f.get('formula') or f.get('type') not in ('DateTime', 'DateTime:UTC'):
            raise ValueError('Registration timestamp schema requires review')
        if any(r['fields'].get(k) == ref for r in columns for k in ('summarySourceCol', 'displayCol', 'visibleCol', 'reverseCol')):
            raise ValueError('Registration timestamp has a linked column dependency')
        if any(r['fields'].get('colRef') == ref for r in client.records('_grist_Filters')):
            raise ValueError('Registration timestamp has a filter dependency')
        for table in ('_grist_ACLResources', '_grist_ACLRules'):
            for row in client.records(table):
                f = row['fields']
                if f.get('aclColumn') == ref or any(re.search(r'\bRegisteredAt\b', str(f.get(k) or ''))
                        for k in ('colIds', 'aclFormula', 'aclFormulaParsed', 'userAttributes')):
                    raise ValueError('Registration timestamp has an access-rule dependency')
        fields = client.records('_grist_Views_section_field')
        removed = {r['id'] for r in fields if r['fields'].get('colRef') == ref}
        def prune(node):
            if not isinstance(node, dict): raise ValueError('Invalid card layout')
            if node.get('leaf') in removed: return None
            result = dict(node)
            if 'children' in node:
                if not isinstance(node['children'], list): raise ValueError('Invalid card layout children')
                result['children'] = [child for child in (prune(c) for c in node['children']) if child is not None]
                if not result['children']: return None
            return result
        for row in client.records('_grist_Views_section'):
            f = row['fields']; update = {}
            if f.get('linkSrcColRef') == ref or f.get('linkTargetColRef') == ref:
                raise ValueError('Registration timestamp has a linked widget dependency')
            if f.get('layoutSpec'):
                layout = json.loads(f['layoutSpec']); revised = prune(layout) or {'children': []}
                if revised != layout: update['layoutSpec'] = json.dumps(revised, separators=(',', ':'))
            if f.get('sortColRefs'):
                sorting = json.loads(f['sortColRefs'])
                if not isinstance(sorting, list): raise ValueError('Invalid widget sort order')
                revised = [value for value in sorting if value not in (ref, -ref)]
                if revised != sorting: update['sortColRefs'] = json.dumps(revised)
            if update: actions.append(['UpdateRecord', '_grist_Views_section', row['id'], update])
        actions.append(['RemoveColumn', 'Experiments', 'RegisteredAt'])
    if apply and actions: client.apply(actions)
    return actions
