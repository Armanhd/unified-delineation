#!/usr/bin/env python3
"""CAMELS-SPAT TDX preview batch: region prepare -> global ID audit -> basin DEM -> exports -> summary.

Two phases of basin processing are intentional: nested HRUs need ALL upstream
station Step05 watersheds from phase A before phase B starts. Every site is
independent: a failed Slurm task does not cancel other sites. No results are
approved for CWARHM/SUMMA production.

Copy to 1_delineation/common/09_tdx_batch.py alongside 09_global_tdx_id_audit.py.
Run `python 1_delineation/common/09_tdx_batch.py --help` for subcommands.
"""
from __future__ import annotations
import argparse
import csv
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
import traceback
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'batch_tdx'
NAME = re.compile(r'^[a-zA-Z0-9_-]+$')
SCRIPT = Path(__file__).resolve()
AUDIT = SCRIPT.with_name('09_global_tdx_id_audit.py')
DEFAULT_PYTHON = '/home/arman.haddadchi/.conda/envs/unified_delineation/bin/python'
TDX = ROOT / '1_delineation/tdx'
DEFAULT_CONTROLS = ROOT / '0_control_files/tdx_camels_spat'

def controls_dir(run):
    config = json_file(run / 'run_config.json', {})
    return Path(config.get('controls_dir', str(DEFAULT_CONTROLS)))

def manifest_path(run):
    return controls_dir(run) / 'batch_controls_manifest.csv'


def now():
    return datetime.now(timezone.utc).isoformat()


def read_csv(path):
    with Path(path).open(newline='', encoding='utf-8-sig') as fh:
        return list(csv.DictReader(fh))


def write_csv(path, rows, columns):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open('w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def json_file(path, default=None):
    try:
        return json.loads(Path(path).read_text())
    except (FileNotFoundError, ValueError, OSError):
        return default


def save_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(payload, indent=2) + '\n')
    os.replace(tmp, path)


def run_path(label):
    if not NAME.fullmatch(label):
        raise ValueError('Run label must use letters, digits, hyphen or underscore only')
    return BASE / label


def prep(args):
    basecontrols = Path(args.controls_dir).expanduser().resolve()
    manifest = basecontrols / 'batch_controls_manifest.csv'
    source = read_csv(manifest)
    only = set(args.sites.split(',')) if args.sites else set()
    rows = [r for r in source if str(r.get('source_input_ready','')).lower() == 'true'
            and (not only or r['site_id'] in only)]
    if args.limit:
        rows = rows[:args.limit]
    if not rows:
        raise RuntimeError('No source-input-ready sites selected')
    selected_ids = {r['site_id'] for r in rows}
    if only - selected_ids:
        print('WARNING: Not selected (not source ready or absent):',sorted(only-selected_ids))
    sites = []
    for r in rows:
        site, region = r['site_id'], str(r.get('tdx_id','')).strip()
        control = Path(r['control_path'])
        if not control.is_absolute():
            control = (ROOT / control).resolve()
        if not control.is_file():
            raise FileNotFoundError(f'Control missing: {control}')
        if not region:
            raise ValueError(f'Source-ready site without TDX region: {site}')
        sites.append({'site_id':site,'tdx_id':region,'control_path':str(control),
                      'snap_status':r.get('snap_status','')})
    region_sites = {}
    for row in sites:
        region_sites.setdefault(row['tdx_id'], row['site_id'])
    regions = [{'tdx_id':region,'example_site':site} for region,site in sorted(region_sites.items())]
    out = run_path(args.label)
    if out.exists() and (out/'sites.csv').exists() and not args.overwrite:
        raise FileExistsError(f'{out} exists; use --overwrite to regenerate the task lists')
    out.mkdir(parents=True, exist_ok=True)
    write_csv(out/'sites.csv', sites, ['site_id','tdx_id','control_path','snap_status'])
    write_csv(out/'regions.csv',regions,['tdx_id','example_site'])
    config = {'label':args.label,'created_at_utc':now(), 'project':str(ROOT),
              'selected_source_ready_sites':len(sites),'selected_regions':len(regions),
              'source_ready_manifest_count':sum(str(x.get('source_input_ready','')).lower()=='true' for x in source),
              'total_manifest_gauges':len(source),'source_inputs_unready':len(source)-sum(str(x.get('source_input_ready','')).lower()=='true' for x in source),
              'purpose':'PREVIEW_ONLY_NOT_PRODUCTION_APPROVED',
              'controls_dir':str(basecontrols)}
    save_json(out/'run_config.json',config)
    print(json.dumps(config,indent=2))
    print('REGIONS:',out/'regions.csv')
    print('SITES:',out/'sites.csv')


def status_path(run, phase, key):
    return run/'status'/f'{phase}_{key}.json'


def command_step(run, key, stage, cmd, stdout_json=None):
    logdir = run/'logs'/key
    logdir.mkdir(parents=True, exist_ok=True)
    stdout = logdir/f'{stage}.out'
    stderr = logdir/f'{stage}.err'
    print(f'  {key}: {stage} - ' + shlex.join(map(str, cmd)), flush=True)
    with stdout.open('w') as so, stderr.open('w') as se:
        proc = subprocess.run(list(map(str,cmd)),cwd=ROOT,stdout=so,stderr=se,check=False)
    if proc.returncode:
        excerpt = '\n'.join(stderr.read_text(errors='replace').splitlines()[-8:])
        raise RuntimeError(f'{stage} exited {proc.returncode}: {excerpt or stdout.read_text(errors="replace")[-500:]}')
    if stdout_json is not None:
        text = stdout.read_text(encoding='utf-8')
        # Step 02 is expected to print JSON only; never persist a malformed partial file.
        try:
            payload = json.loads(text)
        except ValueError as exc:
            raise RuntimeError(f'{stage} did not print pure JSON: {text[-500:]}') from exc
        save_json(stdout_json, payload)
    return {'stdout':str(stdout),'stderr':str(stderr)}


def site_config(row):
    return row['site_id'], Path(row['control_path'])


def region_stage(args, run, row):
    region, site = row['tdx_id'], row['example_site']
    control = controls_dir(run)/f'control_{site}.toml'
    diag = ROOT/'outputs'/site/'diagnostics'
    diag.mkdir(parents=True, exist_ok=True)
    command_step(run, 'region_'+region, '02_region',
                 [args.python,TDX/'02_select_tdx_region.py','--control',control],
                 diag/'tdx_region_selection.json')
    command_step(run, 'region_'+region, '04_index',
                 [args.python,TDX/'04_trace_tdx_upstream.py','--control',control,
                  '--build-topology-index','--index-only'])
    try:
        command_step(run, 'region_'+region, '04b_plan',
                     [args.python,TDX/'04b_plan_shared_nested_splits.py','--region',region])
    except RuntimeError as exc:
        # Keep this region's valid topology usable for basin DEM previews, while
        # the ID audit blocks final exports until the shared plan is repaired.
        return {'split_planner_ready':False,'split_planner_error':str(exc),
                'geometry_plan_ready':False}
    plan = ROOT/'indexes/gauge_catalogs/tdx_shared_network/nested_split_plans'/f'{region}_nested_summary.json'
    details = json_file(plan,{})
    count = int(details.get('candidate_parent_reaches_requiring_split',0))
    geometry_error = None
    if count:
        try:
            command_step(run, 'region_'+region, '04c_geometry',
                         [args.python,TDX/'04c_validate_shared_split_geometry.py','--region',region])
        except RuntimeError as exc:
            # Do not hold up ordinary basins because one planned nested geometry
            # fails. Nested exports in this region will be blocked separately.
            geometry_error = str(exc)
    return {'split_parents_planned':count,'split_planner_ready':True,
            'geometry_plan_ready':geometry_error is None,
            'geometry_plan_error':geometry_error,
            'geometry_skipped_no_splits':count==0,
            'planner_review_parents':details.get('parents_needing_review',0)}


def basin_stage(args, run, row):
    site, control = site_config(row)
    diag = ROOT/'outputs'/site/'diagnostics'
    root = ROOT/'outputs'/site
    diag.mkdir(parents=True, exist_ok=True)
    if (json_file(status_path(run,'region',row['tdx_id']),{}).get('status')!='COMPLETED'):
        raise RuntimeError('Regional topology index preparation failed or was not run')
    if not args.force and (diag/'tdx_refinement_summary.json').is_file() and \
       (root/'tdx_upstream_refined.gpkg').is_file() and \
       (diag/'tdx_upstream_summary.json').is_file():
        print(f'  {site}: reusing existing completed Steps 02-05 (use --force to rerun)',flush=True)
        return {'reused_step05':True}
    command_step(run,site,'02_region',[args.python,TDX/'02_select_tdx_region.py','--control',control],
                 diag/'tdx_region_selection.json')
    command_step(run,site,'03_outlet',[args.python,TDX/'03_select_tdx_outlet.py','--control',control,'--allow-review'])
    command_step(run,site,'04_upstream',[args.python,TDX/'04_trace_tdx_upstream.py','--control',control,
                                        '--require-topology-index','--require-shared-registry'])
    demcmd = [args.python,TDX/'05_refine_tdx_outlet.py','--control',control,
              '--dem-min-accumulation-cells',str(args.dem_min_accumulation_cells)]
    if args.max_connected_void_cells is not None:
        demcmd.extend(['--max-connected-void-cells',str(args.max_connected_void_cells)])
    command_step(run,site,'05_dem',demcmd)
    return {'reused_step05':False}


def choose_nested(site, region):
    import geopandas as gpd
    plan = ROOT/'indexes/gauge_catalogs/tdx_shared_network/nested_split_plans'/f'{region}_nested_split_parents.csv'
    if not plan.is_file():
        raise RuntimeError(f'Missing region 04b plan {plan}')
    parents = {r.get('parent_tdx_reach_id','') for r in read_csv(plan)}
    native = ROOT/'outputs'/site/'tdx_upstream_native.gpkg'
    frame = gpd.read_file(native,layer='upstream_reaches',columns=['reach_id'])
    actual = {str(x).removesuffix('.0') for x in frame['reach_id']}
    related = json_file(ROOT/'outputs'/site/'diagnostics/step04'/f'{site}_tdx_shared_gauge_topology.json',{})
    count = int(related.get('related_gauge_count',0))
    return bool(parents&actual or count), sorted(parents&actual), count


def export_stage(args,run,row):
    site, control = site_config(row)
    state = json_file(status_path(run,'basin',site),{})
    if state.get('status') != 'COMPLETED':
        raise RuntimeError('Steps 02-05 did not finish for this site in this run')
    audit = json_file(run/'id_audit/tdx_global_id_audit_summary.json',{})
    if row['tdx_id'] in audit.get('regions_blocked_for_candidate_export',[]):
        raise RuntimeError('Global ID audit blocked this TDX region; inspect id_audit')
    region_state = json_file(status_path(run,'region',row['tdx_id']),{})
    if not region_state.get('split_planner_ready',False):
        raise RuntimeError('Regional split planner was not ready')
    nested, split_parents, related_count = choose_nested(site,row['tdx_id'])
    if nested and not region_state.get('geometry_plan_ready',False):
        raise RuntimeError('Regional split geometry stage failed; nested export blocked')
    if nested:
        command_step(run,site,'05b_nested',[args.python,TDX/'05b_build_nested_hru_preview.py',
                                           '--control',control,'--allow-review','--overwrite'])
        command_step(run,site,'06b_nested_export',[args.python,TDX/'06b_export_nested_cwarhm_preview.py',
                                                  '--control',control,'--allow-review','--overwrite'])
        command_step(run,site,'07b_nested_validation',[args.python,TDX/'07b_validate_nested_cwarhm_preview.py',
                                                      '--control',control])
        result = json_file(ROOT/'outputs'/site/'validation/tdx_nested_preview/step07b_validation_summary.json',{})
    else:
        command_step(run,site,'06_export',[args.python,TDX/'06_export_cwarhm_shapefiles.py',
                                          '--control',control,'--allow-review','--overwrite'])
        command_step(run,site,'07_validate',[args.python,TDX/'07_validate_and_visualise.py',
                                            '--control',control,'--skip-maps'])
        result = json_file(ROOT/'outputs'/site/'validation/step07_validation_summary.json',{})
    if result.get('structural_validation_pass') is not True:
        raise RuntimeError('Structural validation did not pass (inspect Step07 summary)')
    return {'nested':nested,'split_parents_in_network':split_parents,
            'related_gauge_count':related_count,
            'validation_status':result.get('overall_status','UNKNOWN'),
            'structural_validation_pass':True,
            'production_export_authorized':False}


def task(args):
    run = run_path(args.label)
    rows = read_csv(run/('regions.csv' if args.phase=='region' else 'sites.csv'))
    index = args.offset + args.task_id
    if not (0 <= index < len(rows)):
        raise IndexError(f'Task index {index} out of range (length={len(rows)})')
    row = rows[index]
    key = row['tdx_id'] if args.phase=='region' else row['site_id']
    stat = status_path(run,args.phase,key)
    if not args.force and json_file(stat,{}).get('status') == 'COMPLETED':
        print(f'{args.phase} {key}: already completed in this run; skipping',flush=True)
        return 0
    payload = {'phase':args.phase,'key':key,'tdx_id':row['tdx_id'],
               'started_at_utc':now(),'status':'RUNNING','production_export_authorized':False}
    save_json(stat,payload)
    try:
        if args.phase=='region':
            details=region_stage(args,run,row)
        elif args.phase=='basin':
            details=basin_stage(args,run,row)
        else:
            details=export_stage(args,run,row)
        payload.update(details)
        payload['status']='COMPLETED'
    except Exception as exc:
        payload['status']='FAILED'
        payload['error']=str(exc)
        payload['traceback_last_lines']=traceback.format_exc().splitlines()[-8:]
    payload['ended_at_utc']=now()
    save_json(stat,payload)
    print(json.dumps(payload,indent=2),flush=True)
    return 0 if payload['status']=='COMPLETED' else 1


def submit_job(argv, dry_run):
    print(shlex.join(argv),flush=True)
    if dry_run:
        return 'DRYRUN'
    result=subprocess.run(argv,cwd=ROOT,capture_output=True,text=True,check=False)
    if result.returncode:
        raise RuntimeError('sbatch failed: '+result.stderr+' '+result.stdout)
    job=result.stdout.strip().split(';')[0]
    if not job.isdigit():
        raise RuntimeError('Unexpected sbatch --parsable response '+result.stdout)
    print('SUBMITTED:',job,flush=True)
    return job


def submit(args):
    run=run_path(args.label)
    sites=read_csv(run/'sites.csv')
    regions=read_csv(run/'regions.csv')
    if not sites or not regions:
        raise ValueError('Prepare run first')
    slogs=run/'slurm_logs'
    slogs.mkdir(parents=True,exist_ok=True)
    python=Path(args.python)
    if not args.dry_run and not python.is_file():
        raise FileNotFoundError(f'Python not found: {python}')
    def submit_array(phase, records, limit, dependency):
        jobs=[]
        for start in range(0,len(records),args.array_chunk):
            end=min(start+args.array_chunk,len(records))
            task_args=[str(python),str(SCRIPT),'task','--label',args.label,'--phase',phase,
                       '--offset',str(start),'--task-id','${SLURM_ARRAY_TASK_ID}',
                       '--python',str(python),'--dem-min-accumulation-cells',str(args.dem_min_accumulation_cells)]
            if args.max_connected_void_cells is not None:
                task_args.extend(['--max-connected-void-cells',str(args.max_connected_void_cells)])
            if args.force:
                task_args.append('--force')
            if args.preview_ordinary_slope_fallback:
                task_args.append('--preview-ordinary-slope-fallback')
            # Keep the Slurm array environment variable expandable by Bash.
            wrap=' '.join(x if x=='${SLURM_ARRAY_TASK_ID}' else shlex.quote(x) for x in task_args)
            cmd=['sbatch','--parsable',f'--job-name=tdx_{phase}',
                 f'--array=0-{end-start-1}%{limit}',
                 f'--time={args.time}',f'--mem={args.mem}', '--cpus-per-task=1',
                 f'--chdir={ROOT}',
                 f'--output={slogs}/{phase}_%A_%a.out',
                 f'--error={slogs}/{phase}_%A_%a.err']
            if dependency:
                cmd.append('--dependency=afterany:'+':'.join(dependency))
            cmd.extend(['--wrap',wrap])
            jobs.append(submit_job(cmd,args.dry_run))
        return jobs
    regionjobs=submit_array('region',regions,args.region_parallel,[])
    auditcmd=['sbatch','--parsable','--job-name=tdx_idaudit','--time=01:00:00',
              '--mem=8G','--cpus-per-task=1',f'--chdir={ROOT}',
              f'--output={slogs}/id_audit_%j.out',f'--error={slogs}/id_audit_%j.err',
              '--dependency=afterany:'+':'.join(regionjobs),
              '--wrap',shlex.join([str(python),str(AUDIT),'--project-root',str(ROOT),
                                   '--output-dir',str(run/'id_audit'),
                                   '--expected-region-count',str(args.expected_region_count)])]
    auditjob=submit_job(auditcmd,args.dry_run)
    # Phase A after audit: every watershed independently executes Steps 02-05.
    basinjobs=submit_array('basin',sites,args.basin_parallel,[auditjob])
    # Phase B needs EVERY upstream gauge DEM from Phase A; afterany avoids
    # losing independent basins if a few other array tasks fail.
    exportjobs=submit_array('export',sites,args.export_parallel,basinjobs)
    collectcmd=['sbatch','--parsable','--job-name=tdx_collect','--time=00:30:00',
                '--mem=2G','--cpus-per-task=1',f'--chdir={ROOT}',
                f'--output={slogs}/collect_%j.out',f'--error={slogs}/collect_%j.err',
                '--dependency=afterany:'+':'.join(exportjobs),
                '--wrap',shlex.join([str(python),str(SCRIPT),'collect','--label',args.label])]
    collectjob=submit_job(collectcmd,args.dry_run)
    result={'label':args.label,'dry_run':args.dry_run,'region_jobs':regionjobs,
            'id_audit_job':auditjob,'basin_jobs':basinjobs,
            'export_jobs':exportjobs,'collect_job':collectjob,
            'selected_sites':len(sites),'regions':len(regions),
            'preview_only':True,'production_export_authorized':False}
    if not args.dry_run:
        save_json(run/'submitted_jobs.json',result)
    print(json.dumps(result,indent=2))


def audit_cross_basin_ids(run, completed_rows):
    """Check exact nested IDs repeat consistently across independently exported basins.

    Distinguish a basin outlet (NextDownID=0) from a genuine contradictory edge.
    This is a candidate consistency check, not approval of observed hydrology.
    """
    by_reach={}
    errors=[]
    checked_files=0
    compared_reach_uses=0
    conflict_sites=set()
    for r in completed_rows:
        if r.get('export_status')!='COMPLETED' or str(r.get('nested_export')).lower()!='true':
            continue
        site,region=r['site_id'],r['tdx_id']
        mapping=ROOT/'outputs'/site/'cwarhm_exports/tdx_nested_preview/metadata/reach_id_mapping.csv'
        if not mapping.is_file():
            errors.append({'site_id':site,'tdx_id':region,'TDX_ID':'',
                           'issue':'NESTED_EXPORT_MAPPING_MISSING','previous_site':''})
            conflict_sites.add(site)
            continue
        checked_files+=1
        planfile=(ROOT/'indexes/gauge_catalogs/tdx_shared_network/nested_split_plans'/
                  f'{region}_nested_candidate_pieces.csv')
        known={v.get('candidate_comid_text') for v in read_csv(planfile)} if planfile.exists() else set()
        for entry in read_csv(mapping):
            token=entry.get('TDX_ID_exact','').strip()
            parent=entry.get('parent_native_TDX_ID','').strip()
            down=entry.get('NextDownID_exact','').strip()
            key=(region,token)
            compared_reach_uses+=1
            issue=None
            if '.' in token and token not in known:
                issue='FRACTIONAL_ID_ABSENT_FROM_REGIONAL_SPLIT_PLAN'
            old=by_reach.get(key)
            if old is not None:
                if old['parent']!=parent:
                    issue='SAME_ID_HAS_DIFFERENT_ORIGINAL_PARENT'
                elif down not in ('','0') and old['down'] not in ('','0') and down!=old['down']:
                    issue='SAME_ID_HAS_DIFFERENT_NONOUTLET_NEXTDOWN'
                if old['down'] in ('','0') and down not in ('','0'):
                    old['down']=down
                    old['site']=site
            else:
                by_reach[key]={'parent':parent,'down':down,'site':site}
            if issue:
                errors.append({'site_id':site,'tdx_id':region,'TDX_ID':token,
                               'issue':issue,'previous_site':old['site'] if old else ''})
                conflict_sites.add(site)
                if old:conflict_sites.add(old['site'])
    report={'nested_mapping_files_checked':checked_files,
            'reach_uses_checked':compared_reach_uses,
            'distinct_region_and_reach_ids':len(by_reach),
            'cross_basin_conflict_count':len(errors),
            'affected_sites':sorted(conflict_sites),
            'scope':'Only candidate nested exports that passed their per-basin validation',
            'global_production_approval':False}
    write_csv(run/'tdx_cross_basin_id_conflicts.csv', errors,
              ['site_id','tdx_id','TDX_ID','issue','previous_site'])
    save_json(run/'tdx_cross_basin_id_summary.json',report)
    return report


def collect(args):
    run=run_path(args.label)
    rows=[]
    audit=json_file(run/'id_audit/tdx_global_id_audit_summary.json',{})
    for entry in read_csv(run/'sites.csv'):
        site, region=entry['site_id'],entry['tdx_id']
        basin=json_file(status_path(run,'basin',site),{})
        export=json_file(status_path(run,'export',site),{})
        region_record=json_file(status_path(run,'region',region),{})
        if export.get('status')=='COMPLETED':
            state='STRUCTURAL_PASS_'+export.get('validation_status','REVIEW')
        elif export.get('status')=='FAILED':
            state='EXPORT_FAILED_OR_BLOCKED'
        elif basin.get('status')=='FAILED':
            state='BASIN_PREP_FAILED'
        elif region_record.get('status')=='FAILED':
            state='REGION_PREP_FAILED'
        else:
            state='INCOMPLETE'
        rows.append({'site_id':site,'tdx_id':region,
                     'snap_status':entry.get('snap_status',''),
                     'region_status':region_record.get('status','NOT_RUN'),
                     'step02_to05_status':basin.get('status','NOT_RUN'),
                     'export_status':export.get('status','NOT_RUN'),
                     'nested_export':export.get('nested',''),
                     'validation_status':export.get('validation_status',''),
                     'workflow_status':state,
                     'problem':export.get('error') or basin.get('error') or region_record.get('error',''),
                     'production_export_authorized':False})
    cross=audit_cross_basin_ids(run,rows)
    for row in rows:
        if row['site_id'] in cross['affected_sites']:
            row['workflow_status']='CROSS_BASIN_ID_CONFLICT'
            row['problem']=(row['problem']+'; ' if row['problem'] else '')+'See tdx_cross_basin_id_conflicts.csv'
    write_csv(run/'tdx_batch_results.csv',rows,list(rows[0]) if rows else ['site_id'])
    selected={r['site_id'] for r in rows}
    full_rows=list(rows)
    for source in read_csv(manifest_path(run)):
        if source['site_id'] in selected:
            continue
        ready=str(source.get('source_input_ready','')).lower()=='true'
        full_rows.append({'site_id':source['site_id'], 'tdx_id':source.get('tdx_id',''),
                          'snap_status':source.get('snap_status',''),
                          'region_status':'NOT_SELECTED','step02_to05_status':'NOT_RUN',
                          'export_status':'NOT_RUN','nested_export':'',
                          'validation_status':'',
                          'workflow_status':'SOURCE_READY_NOT_SELECTED' if ready else 'PENDING_SOURCE_OR_SNAP',
                          'problem':source.get('pending_flags',''),
                          'production_export_authorized':False})
    full_rows.sort(key=lambda row:row['site_id'])
    allfile=run/'tdx_batch_all_gauges_status.csv'
    write_csv(allfile,full_rows,list(rows[0]) if rows else list(full_rows[0]))
    summary={'label':args.label,'selected_sites':len(rows),
             'total_manifest_gauges':len(full_rows),
             'statuses':dict(sorted(Counter(r['workflow_status'] for r in rows).items())),
             'all_gauges_statuses':dict(sorted(Counter(r['workflow_status'] for r in full_rows).items())),
             'global_audit_coverage_complete':audit.get('full_global_id_collision_checked',False),
             'global_audit_blocked_regions':audit.get('regions_blocked_for_candidate_export',[]),
             'cross_basin_candidate_id_audit':cross,
             'results_csv':str(run/'tdx_batch_results.csv'),
             'all_gauges_csv':str(allfile),
             'production_export_authorized':False}
    save_json(run/'tdx_batch_results_summary.json',summary)
    print(json.dumps(summary,indent=2))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    sub=ap.add_subparsers(dest='action',required=True)
    p=sub.add_parser('prepare',help='Create ready-site and regional task manifests (no Slurm submission)')
    p.add_argument('--label',required=True)
    p.add_argument('--controls-dir',default=str(DEFAULT_CONTROLS),
                   help='Dataset-specific directory containing batch_controls_manifest.csv and per-site controls')
    p.add_argument('--sites',help='Comma-separated site IDs for a pilot')
    p.add_argument('--limit',type=int,help='Only first N source-ready sites')
    p.add_argument('--overwrite',action='store_true')
    p.set_defaults(func=prep)
    s=sub.add_parser('submit',help='Submit dependency-controlled Slurm arrays')
    s.add_argument('--label',required=True)
    s.add_argument('--dry-run',action='store_true')
    s.add_argument('--python',default=DEFAULT_PYTHON)
    s.add_argument('--region-parallel',type=int,default=2)
    s.add_argument('--basin-parallel',type=int,default=8)
    s.add_argument('--export-parallel',type=int,default=6)
    s.add_argument('--array-chunk',type=int,default=500)
    s.add_argument('--time',default='04:00:00')
    s.add_argument('--mem',default='32G')
    s.add_argument('--dem-min-accumulation-cells',type=int,default=100)
    s.add_argument('--max-connected-void-cells',type=int,default=None,
                   help='OPTIONAL exploratory DEM-void override; keep default safety threshold if omitted')
    s.add_argument('--expected-region-count',type=int,default=15)
    s.add_argument('--preview-ordinary-slope-fallback',action='store_true',
                   help='Deprecated compatibility flag: updated exporters automatically recover TDX slopes')
    s.add_argument('--force',action='store_true',help='Regenerate even if outputs exist')
    s.set_defaults(func=submit)
    t=sub.add_parser('task',help='One Slurm array task (internal)')
    t.add_argument('--label',required=True)
    t.add_argument('--phase',choices=('region','basin','export'),required=True)
    t.add_argument('--offset',type=int,default=0)
    t.add_argument('--task-id',type=int,required=True)
    t.add_argument('--python',default=DEFAULT_PYTHON)
    t.add_argument('--dem-min-accumulation-cells',type=int,default=100)
    t.add_argument('--max-connected-void-cells',type=int)
    t.add_argument('--preview-ordinary-slope-fallback',action='store_true')
    t.add_argument('--force',action='store_true')
    t.set_defaults(func=task)
    c=sub.add_parser('collect',help='Recalculate cross-site progress and failures')
    c.add_argument('--label',required=True)
    c.set_defaults(func=collect)
    args=ap.parse_args()
    if args.action=='submit' and (min(args.region_parallel,args.basin_parallel,args.export_parallel,args.array_chunk)<1):
        ap.error('Array sizes and parallel limits must be positive')
    sys.exit(args.func(args) or 0)


if __name__=='__main__':
    main()
