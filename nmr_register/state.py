"""Durable local observations and outbox with a single-process advisory lock."""
from contextlib import contextmanager
import fcntl
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import time


def serialized(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'))


class State:
    def __init__(self, home):
        self.home = Path(home)
        self.home.mkdir(parents=True,exist_ok=True,mode=0o700)
        self.db = sqlite3.connect(self.home/'state.sqlite3')
        os.chmod(self.home/'state.sqlite3',0o600)
        self.db.row_factory = sqlite3.Row
        self.db.executescript('''
          PRAGMA journal_mode=WAL;
          PRAGMA synchronous=FULL;
          CREATE TABLE IF NOT EXISTS observations(key TEXT PRIMARY KEY, fingerprint TEXT NOT NULL,
            signature TEXT NOT NULL, payload TEXT NOT NULL, first_seen REAL NOT NULL, last_seen REAL NOT NULL, context TEXT NOT NULL DEFAULT '{}');
          CREATE TABLE IF NOT EXISTS outbox(key TEXT PRIMARY KEY, payload TEXT NOT NULL, queued REAL NOT NULL);
          CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY,value TEXT NOT NULL);
          CREATE TABLE IF NOT EXISTS uploads(hash TEXT PRIMARY KEY, attachment INTEGER, state TEXT NOT NULL,
            filename TEXT NOT NULL, updated REAL NOT NULL);
        ''')

        columns = {r[1] for r in self.db.execute('PRAGMA table_info(observations)')}
        if 'context' not in columns:
            self.db.execute("ALTER TABLE observations ADD COLUMN context TEXT NOT NULL DEFAULT '{}'")
            self.db.commit()

    def close(self): self.db.close()

    @contextmanager
    def lock(self):
        fd = os.open(self.home/'collector.lock',os.O_RDWR|os.O_CREAT,0o600)
        try:
            try:
                fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB)
            except BlockingIOError:
                raise RuntimeError('Another collector is already running') from None
            yield self
        finally:
            os.close(fd)

    def get_meta(self,key,default=None):
        row = self.db.execute('SELECT value FROM meta WHERE key=?',(key,)).fetchone()
        return json.loads(row[0]) if row else default

    def set_meta(self,key,value):
        with self.db:
            self.db.execute('INSERT OR REPLACE INTO meta VALUES (?,?)',(key,serialized(value)))

    def observe(self, records, now=None):
        now = now or time.time()
        changed = 0
        with self.db:
            for incoming in records:
                context = {k:v for k,v in incoming.items() if k.startswith('_')}
                fields = {k:v for k,v in incoming.items() if not k.startswith('_')}
                key = fields['AcquisitionKey']
                old = self.db.execute('SELECT * FROM observations WHERE key=?',(key,)).fetchone()
                signature = hashlib.sha256(serialized({'fields':fields,'context':context}).encode()).hexdigest()
                stable = bool(old and old['signature'] == signature)
                reason = ''
                if old and old['fingerprint'] and fields['ObservedFingerprint'] and old['fingerprint'] != fields['ObservedFingerprint']:
                    reason = 'acquisition identity changed at existing path; reconcile before reusing assessment'
                elif old:
                    reason = json.loads(old['payload']).get('ReconciliationReason','')
                old_context = json.loads(old['context']) if old else {}
                old_files = {v[0]:tuple(v[3:]) for v in old_context.get('_RawFileStats',[])}
                current_files = {v[0]:tuple(v[3:]) for v in context.get('_RawFileStats',[])}
                if any(name in old_files and old_files[name] != identity for name,identity in current_files.items()):
                    reason = 'raw acquisition file replaced at existing path; reconcile before reusing assessment'
                matches = self.db.execute('SELECT key FROM observations WHERE fingerprint=? AND key<>?',
                                          (fields['ObservedFingerprint'],key)).fetchall() if fields['ObservedFingerprint'] else []
                if matches:
                    reason = 'matching acquisition metadata at another path; possible copy or rename: '+', '.join(r[0] for r in matches)
                fields['TransferStable'] = stable
                fields['ReconciliationReason'] = reason
                previous = json.loads(old['payload']) if old else None
                comparable = {k:v for k,v in (previous or {}).items() if k != 'LastChangedAt'}
                if fields != comparable:
                    fields['LastChangedAt'] = now
                    payload = serialized(fields)
                    self.db.execute('INSERT OR REPLACE INTO outbox VALUES (?,?,?)',(key,payload,now))
                    changed += 1
                else:
                    fields['LastChangedAt'] = previous['LastChangedAt']
                    payload = serialized(fields)
                self.db.execute('INSERT OR REPLACE INTO observations VALUES (?,?,?,?,?,?,?)',
                                (key,fields['ObservedFingerprint'],signature,payload,old['first_seen'] if old else now,now,serialized(context)))
            # Flag BOTH paths when a copy appears; never propagate assessments.
            for fingerprint, in self.db.execute("SELECT fingerprint FROM observations WHERE fingerprint != '' GROUP BY fingerprint HAVING count(*)>1").fetchall():
                rows = self.db.execute('SELECT key,payload FROM observations WHERE fingerprint=?',(fingerprint,)).fetchall()
                for row in rows:
                    fields = json.loads(row['payload'])
                    reason = 'matching acquisition metadata at another path; possible copy or rename: '+', '.join(r['key'] for r in rows if r['key'] != row['key'])
                    if fields.get('ReconciliationReason') != reason:
                        fields['ReconciliationReason'] = reason
                        fields['LastChangedAt'] = now
                        payload = serialized(fields)
                        self.db.execute('UPDATE observations SET payload=? WHERE key=?',(payload,row['key']))
                        self.db.execute('INSERT OR REPLACE INTO outbox VALUES (?,?,?)',(row['key'],payload,now))
                        changed += 1
        return changed

    def pending(self):
        return [(row['key'],json.loads(row['payload'])) for row in self.db.execute('SELECT * FROM outbox ORDER BY queued,key')]

    def acknowledge(self,keys):
        with self.db:
            self.db.executemany('DELETE FROM outbox WHERE key=?',[(key,) for key in keys])

    def bind_document(self,api_base):
        previous = self.get_meta('document_api_base')
        if previous and previous != api_base:
            raise ValueError('Local state belongs to a different Grist document; use a separate state folder')
        if not previous:
            self.set_meta('document_api_base',api_base)

    def status(self):
        return {'observations':self.db.execute('SELECT count(*) FROM observations').fetchone()[0],
                'pending':self.db.execute('SELECT count(*) FROM outbox').fetchone()[0],
                'last_scan':self.get_meta('last_scan'),'last_sync':self.get_meta('last_sync'),
                'last_error':self.get_meta('last_error','')}
