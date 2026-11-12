"""Keep an assessment editor in the Electrostatics tray, linked by row ID."""
import copy
import json
from .setup import _mcp
from .assessment_widget import assessment_html


def install_electrostatics_assessment(client, doc_id):
    call = lambda name, **kw: _mcp(client, name, {'doc_id':doc_id, **kw})
    views = client.records('_grist_Views')
    sections = client.records('_grist_Views_section')
    tables = client.records('_grist_Tables')
    pages = [r for r in views if r['fields'].get('name') == 'Electrostatics']
    experiments = [r for r in tables if r['fields'].get('tableId') == 'Experiments']
    if len(pages) != 1 or len(experiments) != 1:
        raise ValueError('Electrostatics page or Experiments table missing or ambiguous')
    page, table = pages[0]['id'], experiments[0]['id']
    on_page = [r for r in sections if r['fields'].get('parentId') == page]
    grids = [r for r in on_page if r['fields'].get('title') == 'Individual acquisitions']
    cards = [r for r in on_page if r['fields'].get('title') == 'Assess experiment']
    if len(grids) != 1 or len(cards) > 1:
        raise ValueError('Electrostatics acquisition selector or assessment editor ambiguous')
    grid = grids[0]
    if grid['fields'].get('tableRef') != table or grid['fields'].get('parentKey') != 'record':
        raise ValueError('Unexpected acquisition selector table or type')
    try:
        layout = json.loads(pages[0]['fields'].get('layoutSpec') or '{}')
    except (TypeError, ValueError):
        raise ValueError('Invalid Electrostatics layout') from None
    if (not isinstance(layout, dict) or not isinstance(layout.get('children'), list)
            or not isinstance(layout.get('collapsed', []), list)):
        raise ValueError('Unexpected Electrostatics layout; preserve for review')
    def leaves(node):
        if isinstance(node, list):
            return [wid for item in node for wid in leaves(item)]
        if isinstance(node, dict):
            return ([node['leaf']] if 'leaf' in node else []) + leaves(node.get('children', [])) + leaves(node.get('collapsed', []))
        return []
    ids = leaves(layout)
    page_ids = {r['id'] for r in on_page}
    if (any(not isinstance(wid, int) or isinstance(wid, bool) or wid <= 0 for wid in ids)
            or len(ids) != len(set(ids)) or not set(ids) <= page_ids or grid['id'] not in ids):
        raise ValueError('Unexpected Electrostatics layout widget references')
    card = cards[0] if cards else None
    card_id = card['id'] if card else None
    if card:
        cf = card['fields']
        if (cf.get('tableRef') != table or cf.get('parentKey') not in ('creator', 'custom')
                or cf.get('linkSrcSectionRef') not in (None, 0, grid['id'])
                or cf.get('linkSrcColRef') or cf.get('linkTargetColRef')):
            raise ValueError('Existing assessment editor has custom table, type or link')
        if any(r['fields'].get('linkSrcSectionRef') == card_id for r in sections):
            raise ValueError('Assessment editor controls another widget')
    # Follow the selector's upstream graph before adding the new edge.
    by_id = {r['id']:r for r in sections}
    seen = set()
    current = grid['id']
    while current:
        if current in seen or current == card_id or current not in by_id:
            raise ValueError('Assessment row-ID link would be cyclic or unresolved')
        seen.add(current)
        current = by_id[current]['fields'].get('linkSrcSectionRef') or 0
    if card:
        existing_html = call('grist_get_creator_html', widget_id=card_id).get('html') or ''
        if existing_html and 'data-nmr-assessment="1"' not in existing_html:
            raise ValueError('Existing editor has custom HTML; preserve for review')
        settings = call('grist_get_custom_widget_settings', widget_id=card_id)
        if settings.get('unmapped_required_columns'):
            raise ValueError('Existing editor has required column mappings')
    else:
        created = call('grist_add_page_widget', page_id=page, widget_options={
            'table_id':'Experiments', 'type':'creator', 'title':'Assess experiment'})
        card_id = created.get('id', created.get('widget_id'))
        if not isinstance(card_id, int) or isinstance(card_id, bool) or card_id <= 0:
            raise ValueError('Assessment editor creation failed')
        fresh = [r for r in client.records('_grist_Views_section') if r['id'] == card_id]
        if (len(fresh) != 1 or fresh[0]['fields'].get('parentId') != page
                or fresh[0]['fields'].get('tableRef') != table
                or fresh[0]['fields'].get('parentKey') not in ('creator', 'custom')):
            raise ValueError('Unexpected created assessment editor')
        card = fresh[0]
        existing_html = ''
        settings = call('grist_get_custom_widget_settings', widget_id=card_id)
    html = assessment_html(section_id=card_id)
    if existing_html != html:
        call('grist_set_creator_html', widget_id=card_id, html=html)
    if settings.get('access') != 'full':
        settings = call('grist_set_custom_widget_settings', widget_id=card_id, access='full')
    if settings.get('access') != 'full' or settings.get('unmapped_required_columns'):
        raise ValueError('Assessment widget settings incomplete')
    # Preserve every existing visible branch and unrelated tray entry. On initial
    # creation use the captured layout, before AddWidget inserted its default leaf.
    def remove(node):
        if isinstance(node, list):
            return [result for item in node if (result := remove(item)) is not None]
        if not isinstance(node, dict):
            return node
        if node.get('leaf') == card_id:
            return None
        result = {k:remove(v) if k in ('children', 'collapsed') else copy.deepcopy(v) for k,v in node.items()}
        if 'children' in result and not result['children'] and 'leaf' not in result:
            return None
        return result
    collapsed = layout.get('collapsed', [])
    already_collapsed = any(entry.get('leaf') == card_id for entry in collapsed if isinstance(entry, dict))
    desired = copy.deepcopy(layout) if already_collapsed else remove(layout)
    if not already_collapsed:
        desired.setdefault('collapsed', []).append({'leaf':card_id})
    actions = []
    cf = card['fields']
    link = {'linkSrcSectionRef':grid['id'], 'linkSrcColRef':0, 'linkTargetColRef':0}
    if any(cf.get(k) != v for k,v in link.items()):
        actions.append(['UpdateRecord','_grist_Views_section',card_id,link])
    # Creation may have changed layoutSpec, so compare against current metadata.
    current_views = client.records('_grist_Views') if not cards else views
    current_page = next(r for r in current_views if r['id'] == page)
    if json.loads(current_page['fields'].get('layoutSpec') or '{}') != desired:
        actions.append(['UpdateRecord','_grist_Views',page,{'layoutSpec':json.dumps(desired,separators=(',',':'))}])
    if actions:
        client.apply(actions)
    return {'page':page, 'grid':grid['id'], 'card':card_id}
