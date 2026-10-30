"""Small Grist API adapter. Ambiguous mutations are reconciled by callers."""
import time
from urllib.parse import quote
import requests


class GristError(RuntimeError):
    pass


class GristClient:
    def __init__(self, config, token, session=None, sleep=time.sleep):
        self.base = config.api_base
        self.timeout = config.timeout
        self.session = session or requests.Session()
        self.session.headers.update({'Authorization':'Bearer '+token})
        self.sleep = sleep
        self.calls = 0

    def request(self, method, path, **kwargs):
        # GET is safe to retry. POST/PATCH retries belong to a reconciliation loop;
        # never blindly repeat a mutation after a timeout or server failure.
        attempts = 3 if method == 'GET' else 1
        for attempt in range(attempts):
            try:
                self.calls += 1
                response = self.session.request(method,self.base+path,timeout=self.timeout,
                                                allow_redirects=False,**kwargs)
            except requests.RequestException as error:
                if method == 'GET' and attempt+1 < attempts:
                    self.sleep(min(2**attempt,4))
                    continue
                raise GristError('Grist request failed ('+error.__class__.__name__+'); pending work retained') from None
            if response.status_code in (429,502,503,504) and attempt+1 < attempts:
                try: delay = min(max(float(response.headers.get('Retry-After',2**attempt)),0),10)
                except ValueError: delay = 2**attempt
                self.sleep(delay)
                continue
            if not 200 <= response.status_code < 300:
                # Bodies and URLs may contain private content. Never log them.
                raise GristError(f'Grist HTTP {response.status_code}; pending work retained')
            return response
        raise GristError('Grist request failed')

    def records(self,table):
        return self.request('GET','/tables/'+quote(table,safe='')+'/records').json()['records']

    def tables(self):
        return self.request('GET','/tables').json()['tables']

    def _write_records(self,method,table,records):
        # Grist PATCH requires the same field set on every record in a request.
        # Group POST as well: absent fields must stay absent, rather than filling
        # a union with nulls and accidentally clearing another row's human cells.
        groups = {}
        for index,record in enumerate(records):
            groups.setdefault(tuple(sorted(record['fields'])),[]).append((index,record))
        ordered = [None]*len(records)
        for group in groups.values():
            response = self.request(method,'/tables/'+quote(table,safe='')+'/records',
                json={'records':[record for _,record in group]})
            try:
                result = response.json()
            except ValueError:
                if method == 'POST':
                    raise GristError('Invalid Grist record creation response; reconcile pending work before retrying') from None
                result = None
            returned = result.get('records') if isinstance(result,dict) else None
            if method == 'POST' and (not isinstance(returned,list) or len(returned) != len(group)):
                raise GristError('Unexpected Grist record creation response; reconcile pending work before retrying')
            for position,(original,record) in enumerate(group):
                ordered[original] = returned[position] if isinstance(returned,list) and len(returned)==len(group) else {'id':record['id']}
        # All original work is acknowledged only after every subgroup succeeds.
        # A timeout keeps the complete upstream batch queued for stable-key GET
        # reconciliation, including any already committed earlier subgroups.
        return {'records':ordered}

    def add_records(self,table,fields):
        return self._write_records('POST',table,[{'fields':f} for f in fields])

    def add(self,table,fields):
        result = self.add_records(table,[fields])
        return result['records'][0]['id']

    def update_records(self,table,records):
        return self._write_records('PATCH',table,records)

    def apply(self,actions):
        return self.request('POST','/apply',json=actions).json()

    def attachments(self):
        result = self.request('GET','/attachments').json()
        return result.get('records',[]) if isinstance(result,dict) else result

    def upload(self,name,content):
        return self.request('POST','/attachments',files={'upload':(name,content,'text/csv')}).json()

    def download_attachment(self,attachment_id):
        return self.request('GET','/attachments/'+str(int(attachment_id))+'/download').content

    def download(self):
        return self.request('GET','/download').content

    def test_connection(self):
        return {'connected':True,'tables':len(self.tables())}

    def _mcp_rpc(self,method,params=None,notification=False):
        """Minimal Streamable HTTP MCP, using this same locally enrolled API key.

        No OAuth cache is read. SSE and JSON replies are supported; callers must
        review a tool's inputSchema before invoking a mutation.
        """
        import json
        from urllib.parse import urlsplit
        origin = urlsplit(self.base)
        url = origin.scheme+'://'+origin.netloc+'/api/mcp'
        headers = {'Accept':'application/json, text/event-stream','Content-Type':'application/json'}
        if getattr(self,'_mcp_session',None):
            headers['Mcp-Session-Id']=self._mcp_session
        if getattr(self,'_mcp_protocol',None):
            headers['MCP-Protocol-Version']=self._mcp_protocol
        payload = {'jsonrpc':'2.0','method':method}
        if params is not None: payload['params']=params
        if not notification:
            self._mcp_counter=getattr(self,'_mcp_counter',0)+1
            payload['id']=self._mcp_counter
        try:
            self.calls += 1
            response=self.session.request('POST',url,headers=headers,json=payload,
                timeout=self.timeout,allow_redirects=False)
        except requests.RequestException as error:
            raise GristError('Grist MCP request failed ('+error.__class__.__name__+'); reconcile tool state before retrying') from None
        if not 200 <= response.status_code < 300:
            raise GristError(f'Grist MCP HTTP {response.status_code}')
        session=response.headers.get('Mcp-Session-Id')
        if session: self._mcp_session=session
        if notification and not response.content:
            return None
        try:
            if 'text/event-stream' in response.headers.get('Content-Type',''):
                messages=[]
                for block in response.text.replace('\r\n','\n').split('\n\n'):
                    data='\n'.join(line[5:].lstrip() for line in block.splitlines() if line.startswith('data:'))
                    if data:
                        messages.append(json.loads(data))
                message=next(item for item in reversed(messages) if item.get('id')==payload.get('id'))
            else:
                message=response.json()
        except (ValueError,StopIteration,TypeError):
            raise GristError('Invalid Grist MCP response') from None
        if message.get('error'):
            raise GristError('Grist MCP returned a JSON-RPC error')
        return message.get('result')

    def _mcp_initialize(self):
        if getattr(self,'_mcp_ready',False):
            return
        result=self._mcp_rpc('initialize',{'protocolVersion':'2025-03-26',
            'capabilities':{},'clientInfo':{'name':'nmr-register','version':'0.1.0'}})
        self._mcp_protocol=result.get('protocolVersion','2025-03-26')
        self._mcp_rpc('notifications/initialized',notification=True)
        self._mcp_ready=True

    def mcp_tools(self):
        self._mcp_initialize()
        return self._mcp_rpc('tools/list',{})

    def mcp_call(self,name,arguments):
        self._mcp_initialize()
        return self._mcp_rpc('tools/call',{'name':name,'arguments':arguments})
