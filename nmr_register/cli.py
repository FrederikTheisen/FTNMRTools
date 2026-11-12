"""Desktop collector commands. Credentials are entered only through a local prompt."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import sys
from .backup import create_backup,validate_backup,restore_copy
from .config import Config,Source,load_config,save_config,default_home
from .credentials import enroll,read_token
from .grist import GristClient,GristError
from .state import State
from .sync import run


def parser():
    main = argparse.ArgumentParser(prog='nmr-register')
    main.add_argument('--config',type=Path,help='Non-secret configuration JSON')
    sub = main.add_subparsers(dest='command',required=True)
    sync = sub.add_parser('sync',help='Scan locally and flush changed observations')
    sync.add_argument('--dry-run',action='store_true',help='Scan only; no API or outbox changes')
    sub.add_parser('status',help='Show local collector state without contacting Grist')
    credential = sub.add_parser('credential',help='Enroll API key through a terminal or zenity')
    credential.add_argument('--test',action='store_true',help='Test existing key with a read-only request')
    config = sub.add_parser('config',help='Show or initialize non-secret configuration')
    config.add_argument('--init',action='store_true')
    config.add_argument('--doc-url',default=None)
    config.add_argument('--source',action='append',default=[],metavar='ID=ABSOLUTE_PATH')
    config.add_argument('--state-dir',type=Path)
    config.add_argument('--all-acquisitions',action='store_true',help='Explicitly broaden collection after pilot acceptance')
    config.add_argument('--results-manifest',type=Path)
    config.add_argument('--result-source-id',default='')
    backup = sub.add_parser('backup',help='Download and validate a self-contained .grist backup')
    backup.add_argument('--premigration',action='store_true')
    backup.add_argument('--directory',type=Path)
    backup.add_argument('--validate',type=Path)
    backup.add_argument('--restore-to',type=Path,help='With --validate, restore to a new local file only')
    setup = sub.add_parser('setup',help='Review or apply fresh register schema and layout')
    setup.add_argument('--apply',action='store_true')
    electro = sub.add_parser('electrostatics',help='Review or install the narrow Electrostatics page and reviewed pilot fields')
    electro.add_argument('--catalog',type=Path,required=True,help='Reviewed acquisition-catalog.csv')
    electro.add_argument('--apply',action='store_true',help='Back up, install the page and apply the reviewed pilot fields')
    migrate = sub.add_parser('import',help='Plan a selected pilot or full read-only migration')
    migrate.add_argument('--workbook',type=Path,required=True)
    migrate.add_argument('--reviews',type=Path,required=True)
    migrate.add_argument('--source-id',required=True)
    migrate.add_argument('--sample',action='append')
    migrate.add_argument('--target',action='append')
    migrate.add_argument('--catalog',type=Path)
    migrate.add_argument('--status-workbook',action='append',type=Path,default=[])
    migrate.add_argument('--overview',type=Path)
    migrate.add_argument('--full',action='store_true',help='Build full review plan only')
    migrate.add_argument('--apply',action='store_true',help='Apply selected pilot only, insert-only')
    results = sub.add_parser('results',help='Queue and attach explicitly selected CURRENT numbers CSVs')
    results.add_argument('--manifest',type=Path,required=True)
    results.add_argument('--source-id',required=True)
    results.add_argument('--dry-run',action='store_true')
    snapshot = sub.add_parser('review-snapshot',help='Export an explicit stable review snapshot without changing legacy inputs')
    snapshot.add_argument('--source-reviews',type=Path,required=True)
    snapshot.add_argument('--source-id',required=True)
    snapshot.add_argument('--output',type=Path,required=True)
    return main


def _client(config):
    return GristClient(config,read_token(config.home,config.doc_url,config.credential_file or None))


def install_electrostatics(client,config,catalog,apply=False):
    """Review or restart the narrow, pilot-only page and field installation."""
    from .migration import reviewed_pilot_patch,apply_reviewed_pilot_patch
    from .setup import (upgrade_acquisition_schema,upgrade_labelling_choices,upgrade_electrostatics_schema,
        build_layout,apply_electrostatics_ui,apply_experiments_ui,apply_experiment_registration_widget,apply_sample_registration_widget)
    from .samples_ui import apply_samples_ui
    before=reviewed_pilot_patch(client.records('Experiments'),catalog)
    if not apply:
        return {'dry_run':True,'pilot_patch':before,'scope':'Electrostatics page and reviewed existing pilot fields'}
    if before['issues']:
        raise ValueError('Pilot patch has unresolved identity/evidence issues; run without --apply for details')
    try:
        backup=create_backup(client,config.home/'backups',True)
        upgrade_acquisition_schema(client)
        upgrade_labelling_choices(client)
        upgrade_electrostatics_schema(client)
        build_layout(client,electrostatics_only=True)
        doc_id=config.api_base.rsplit('/',1)[1]
        layout=apply_electrostatics_ui(client,doc_id)
        form=apply_experiment_registration_widget(client,doc_id)
        experiments_layout=apply_experiments_ui(client,doc_id)
        sample_form=apply_sample_registration_widget(client,doc_id)
        samples_layout=apply_samples_ui(client,doc_id)
        # New Choice cells can read as empty strings instead of absent/None.
        # Recompute the guarded plan after schema creation and reject changes
        # to its reviewed field updates or evidence during installation.
        after=reviewed_pilot_patch(client.records('Experiments'),catalog)
        if (after['issues'] or after['records'] != before['records'] or
                after['catalog_sha256'] != before['catalog_sha256']):
            raise ValueError('Pilot fields or catalog changed during setup; review again before resuming')
        applied=apply_reviewed_pilot_patch(after,client,dry_run=False)
    except GristError as error:
        if '429' in str(error):
            raise GristError('Grist HTTP 429: site quota or rate limit reached; wait for access to recover, then rerun electrostatics --apply. Completed steps are restartable.') from None
        raise
    return {'dry_run':False,'backup':backup['path'],'layout':layout,
            'experiment_form':form,'experiments_layout':experiments_layout,'sample_form':sample_form,'samples_layout':samples_layout,'pilot_patch':applied}


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        if args.command == 'config' and args.init:
            from .migration import PILOT_TARGETS
            targets = None if args.all_acquisitions else list(PILOT_TARGETS)
            sources=[]
            for item in args.source:
                ident,sep,path = item.partition('=')
                if not sep: raise ValueError('--source expects ID=ABSOLUTE_PATH')
                sources.append(Source(ident,str(Path(path).resolve()),targets=targets))
            config = Config(sources=sources or [Source('nmrdata',str((Path.home()/'Documents/Projects/nmrdata').resolve()),targets=targets)],
                state_dir=str(args.state_dir) if args.state_dir else '',
                results_manifest=str(args.results_manifest) if args.results_manifest else '',
                result_source_id=args.result_source_id)
            if args.doc_url: config.doc_url=args.doc_url
            result = {'config':str(save_config(config,args.config))}
        elif args.command == 'backup' and args.validate:
            result = validate_backup(args.validate)
            if args.restore_to:
                result['restored_to']=str(restore_copy(args.validate,args.restore_to))
        elif args.command == 'setup' and not args.apply:
            from .schema import schema_actions
            result = {'dry_run':True,'actions':schema_actions()}
        else:
            config = load_config(args.config)
            if args.command == 'config':
                from dataclasses import asdict
                result = asdict(config)
            elif args.command == 'credential':
                if args.test:
                    result = _client(config).test_connection()
                else:
                    result = {'credential':str(enroll(config.home,lambda token:GristClient(config,token))),'connection':'read-only test passed'}
            else:
                state = State(config.home)
                try:
                    with state.lock():
                        if args.command == 'status':
                            result = state.status()
                        elif args.command == 'sync':
                            previous_sync = state.get_meta('last_sync')
                            client = None if args.dry_run else _client(config)
                            result = run(config,state,client,args.dry_run)
                            if not args.dry_run:
                                from .results import select_results,attach_results
                                selections = None
                                if config.results_manifest:
                                    if not config.result_source_id:
                                        raise ValueError('results_manifest requires result_source_id')
                                    source = next((s for s in config.sources if s.id == config.result_source_id),None)
                                    if source is None: raise ValueError('Unknown result_source_id')
                                    selections = select_results(config.results_manifest,config.result_source_id,config.max_result_bytes,source.targets)
                                result.update(attach_results(state,client,selections,config.max_result_bytes))
                                today=datetime.now(timezone.utc).date().isoformat()
                                if state.get_meta('last_backup_day')!=today:
                                    result['backup']=create_backup(client,config.home/'backups')
                                    state.set_meta('last_backup_day',today)
                        elif args.command == 'backup':
                            result = create_backup(_client(config),args.directory or config.home/'backups',args.premigration)
                        elif args.command == 'setup':
                            from .setup import create_schema,build_layout,build_registration_widgets
                            client = _client(config)
                            snapshot = create_backup(client,config.home/'backups',True,required_tables=('_grist_Tables','_grist_Tables_column'))
                            create_schema(client)
                            build_layout(client)
                            build_registration_widgets(client,config.api_base.rsplit('/',1)[1])
                            result = {'setup':'created','premigration':snapshot['path']}
                        elif args.command == 'electrostatics':
                            result = install_electrostatics(_client(config),config,args.catalog,args.apply)
                        elif args.command == 'import':
                            from .migration import build_pilot,build_full,apply_pilot,PILOT_SAMPLES,PILOT_TARGETS
                            if args.full and args.apply:
                                raise ValueError('Full migration application requires separate recorded pilot and Mac acceptance; CLI only builds its review plan')
                            plan = build_full(args.workbook,args.reviews,args.source_id) if args.full else build_pilot(args.workbook,args.reviews,args.sample or PILOT_SAMPLES,args.target or PILOT_TARGETS,args.source_id)
                            from .migration import reconcile_sources
                            plan = reconcile_sources(plan,catalog=args.catalog,status_workbooks=args.status_workbook,overview=args.overview)
                            if not args.apply:
                                result = plan
                            else:
                                client = _client(config)
                                snapshot=create_backup(client,config.home/'backups',True)
                                result = apply_pilot(plan,client,dry_run=False)
                                result['premigration']=snapshot['path']
                        elif args.command == 'review-snapshot':
                            from .legacy import review_snapshot
                            import os,tempfile
                            if args.output.resolve() == args.source_reviews.resolve():
                                raise ValueError('Snapshot output must differ from the source review TSV')
                            content=review_snapshot(args.source_reviews,_client(config).records('Experiments'),source_id=args.source_id)
                            args.output.parent.mkdir(parents=True,exist_ok=True,mode=0o700)
                            fd,tmp=tempfile.mkstemp(prefix='.nmr-reviews-',dir=args.output.parent)
                            try:
                                with os.fdopen(fd,'wb') as handle:
                                    handle.write(content);handle.flush();os.fsync(handle.fileno())
                                os.replace(tmp,args.output)
                            finally:
                                if os.path.exists(tmp):os.unlink(tmp)
                            result={'snapshot':str(args.output),'bytes':len(content),'source_unchanged':True}
                        elif args.command == 'results':
                            from .results import select_results,attach_results
                            source=next((s for s in config.sources if s.id==args.source_id),None)
                            if source is None: raise ValueError('Unknown result source ID')
                            selected=select_results(args.manifest,args.source_id,config.max_result_bytes,source.targets)
                            result = {'selected':len(selected),'dry_run':True} if args.dry_run else attach_results(state,_client(config),selected,config.max_result_bytes)
                        else:
                            raise ValueError('Unknown command')
                except (OSError,ValueError,RuntimeError,GristError) as error:
                    state.set_meta('last_error',str(error))
                    if args.command == 'sync' and not args.dry_run and 'previous_sync' in locals():
                        state.set_meta('last_sync',previous_sync)
                    if args.command == 'sync' and not args.dry_run and 'client' in locals() and client is not None:
                        from .sync import sync_scan_status
                        try: sync_scan_status(state,client,config,str(error))
                        except (GristError,OSError,ValueError,RuntimeError): pass
                    raise
                finally:
                    state.close()
        print(json.dumps(result,indent=2,default=str))
        return 0
    except (OSError,ValueError,RuntimeError,GristError) as error:
        # Transport errors have already been sanitized; local errors carry no key.
        print('nmr-register: '+str(error),file=sys.stderr)
        return 1
