from datetime import datetime,timedelta,timezone
from pathlib import Path
import sqlite3
import tempfile
import unittest
from nmr_register.backup import validate_backup,create_backup,restore_copy,REQUIRED_TABLES


def make_doc(path,attachment=True,external=False):
    db=sqlite3.connect(path)
    for table in REQUIRED_TABLES:db.execute('CREATE TABLE '+table+'(id INTEGER PRIMARY KEY)')
    db.execute('CREATE TABLE _grist_Attachments(id INTEGER PRIMARY KEY,fileIdent TEXT,storageId TEXT)')
    db.execute('CREATE TABLE _gristsys_Files(ident TEXT PRIMARY KEY,data BLOB)')
    if attachment:
        db.execute('INSERT INTO _grist_Attachments VALUES (1,?,?)',('hash.csv','external' if external else ''))
        if not external:db.execute('INSERT INTO _gristsys_Files VALUES (?,?)',('hash.csv',b'original numbers csv'))
    db.commit();db.close()


class BackupTests(unittest.TestCase):
    def setUp(self):self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
    def tearDown(self):self.tmp.cleanup()

    def test_seven_daily_rotation_and_premigration_not_rotated(self):
        source=self.root/'source.grist';make_doc(source)
        class Client:
            def download(self):return source.read_bytes()
        client=Client();backups=self.root/'backups';start=datetime(2026,10,1,tzinfo=timezone.utc)
        pre=create_backup(client,backups,True,now=start)
        for day in range(9):create_backup(client,backups,now=start+timedelta(days=day))
        self.assertEqual(len(list(backups.glob('daily-*.grist'))),7)
        self.assertTrue(Path(pre['path']).exists())
        latest=backups/'daily-2026-10-09.grist'
        self.assertEqual(validate_backup(latest)['attachments'],1)
        restored=restore_copy(latest,self.root/'restored.grist')
        connection=sqlite3.connect(restored)
        self.assertEqual(connection.execute('SELECT data FROM _gristsys_Files').fetchone()[0],b'original numbers csv')
        connection.close()
        with self.assertRaises(ValueError):restore_copy(latest,restored)

    def test_invalid_download_never_overwrites_good_backup(self):
        source=self.root/'source.grist';make_doc(source)
        class Client:
            content=source.read_bytes()
            def download(self):return self.content
        client=Client();day=datetime(2026,10,1,tzinfo=timezone.utc)
        result=create_backup(client,self.root/'backups',now=day)
        before=Path(result['path']).read_bytes();client.content=b'html error page'
        with self.assertRaises(ValueError):create_backup(client,self.root/'backups',now=day)
        self.assertEqual(Path(result['path']).read_bytes(),before)
        self.assertFalse(list((self.root/'backups').glob('.nmr-backup-*')))

    def test_external_attachments_rejected(self):
        source=self.root/'external.grist';make_doc(source,external=True)
        with self.assertRaisesRegex(ValueError,'external attachments'):validate_backup(source)

    def test_header_and_missing_tables_rejected(self):
        invalid=self.root/'invalid.grist';invalid.write_bytes(b'not sqlite')
        with self.assertRaisesRegex(ValueError,'SQLite'):validate_backup(invalid)
        minimal=self.root/'minimal.grist';db=sqlite3.connect(minimal);db.execute('CREATE TABLE unrelated(id)');db.close()
        with self.assertRaisesRegex(ValueError,'missing required'):validate_backup(minimal)
