"""Per-client visibility control for the main experiment grid."""
import json

CORE_JS = r'''
const ExperimentVisibility=(()=>{
 function select(records,show){
  if(!Array.isArray(records))throw Error('Experiment list is unavailable.');
  const ids=[],seen=new Set();let hidden=0,conflicted=0;
  for(const r of records){
   if(!Number.isInteger(r.id)||r.id<=0||seen.has(r.id)||typeof r.ExperimentType!=='string'&&r.ExperimentType!=null||typeof r.ReconciliationReason!=='string'&&r.ReconciliationReason!=null||typeof r.ScanError!=='string'&&r.ScanError!=null)throw Error('Experiment list is invalid.');
   seen.add(r.id);const special=['Calibration','1D'].includes(r.ExperimentType),conflict=Boolean(r.ReconciliationReason||r.ScanError);
   if(special&&conflict)conflicted++;
   if(special&&!conflict&&!show)hidden++;else ids.push(r.id);
  }
  return {ids,hidden,conflicted};
 }
 return {select};
})();
'''
UI_JS = r'''
const $=id=>document.getElementById(id);let records=null,show=false,busy=false,sequence=0,queue=Promise.resolve();
function fail(){records=null;$('toggle').disabled=true;$('status').textContent='Could not load experiment list';}
async function update(){
 const token=++sequence;if(!records)return;busy=true;$('toggle').disabled=true;
 try{const result=ExperimentVisibility.select(records,show);queue=queue.catch(()=>{}).then(()=>grist.setSelectedRows(result.ids));await queue;if(token!==sequence)return;
 $('toggle').textContent=show?'Hide calibration / 1D':'Show calibration / 1D';
 $('status').textContent=(show?'Showing calibration and 1D':result.hidden+' calibration / 1D hidden')+(result.conflicted?' · '+result.conflicted+' with data warnings included':'')+' · Other filters still apply';
 }catch(e){if(token===sequence)fail();}finally{if(token===sequence){busy=false;if(records)$('toggle').disabled=false;}}
}
$('toggle').addEventListener('click',()=>{if(!busy&&records){show=!show;update();}});
grist.ready({requiredAccess:'read table',allowSelectBy:true});
const initialized=grist.setSelectedRows([]);
initialized.catch(fail);
grist.onRecords(async rows=>{try{await initialized;records=rows;await update();}catch(e){fail();}});
'''


def experiment_visibility_html(widget_id=None, grid_id=None):
    for value in (widget_id, grid_id):
        if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
            raise ValueError('widget and grid IDs must be positive integers')
    return '''<!doctype html><html lang="en" data-nmr-experiment-visibility="1" data-widget-id="''' + str(widget_id) + '''" data-grid-id="''' + str(grid_id) + '''"><head><meta name="viewport" content="width=device-width,initial-scale=1"><script src="https://docs.getgrist.com/grist-plugin-api.js"></script><style>body{margin:0;padding:5px 8px;font:13px system-ui,sans-serif;color:var(--grist-theme-text,#262633);background:var(--grist-theme-widget-bg,#fff)}button{font:inherit;font-weight:600;padding:7px 12px;border:1px solid #32865a;border-radius:4px;background:#e7f3eb;color:#244d33;cursor:pointer}button:disabled{opacity:.6}#status{margin-left:9px;color:#555}</style></head><body><button id="toggle" type="button" disabled>Show calibration / 1D</button><span id="status" role="status">Loading experiment list</span><script>''' + CORE_JS + '\n' + UI_JS + '</script></body></html>'
