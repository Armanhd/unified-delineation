#!/usr/bin/env python3
"""Collect CAMELS-SPAT TDX stage readiness and QA statuses without running models.

Run after any stage, before or after Slurm arrays. Does not modify any basin data.
"""
from __future__ import annotations
import argparse
import csv
import json
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]


def load_json(path):
    try:
        return json.loads(path.read_text()) if path.exists() else None
    except (ValueError, OSError):
        return {'_parse_error':True}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project-root',type=Path,default=ROOT)
    ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    root=args.project_root.resolve()
    manifest=root/'0_control_files/tdx_camels_spat/batch_controls_manifest.csv'
    if not manifest.exists():
        raise FileNotFoundError(manifest)
    with manifest.open(newline='') as f:
        items=list(csv.DictReader(f))
    output=args.output or root/'indexes/gauge_catalogs/tdx_batch_qa_status.csv'
    output.parent.mkdir(parents=True,exist_ok=True)
    rows=[]
    for record in items:
        site=record['site_id']
        diag=root/'outputs'/site/'diagnostics'
        s04=load_json(diag/'tdx_upstream_summary.json')
        s05=load_json(diag/'tdx_refinement_summary.json')
        s06=load_json(diag/'step06_preflight.json')
        s07=load_json(root/'outputs'/site/'validation/step07_validation_summary.json')
        snap=record.get('snap_status','')
        source_ready=record.get('source_input_ready','').lower()=='true'
        nested=bool(s04 and (s04.get('shared_related_gauge_count') or 0)>0)
        export_status=(s06 or {}).get('status','NOT_RUN')
        geom_status=(s07 or {}).get('overall_status','NOT_RUN')
        if not source_ready:
            status='PENDING_SOURCE_OR_SNAP'
        elif nested or (s06 and 'PLANNED_SPLIT_PARENTS_PRESENT' in (s06.get('blocked_reasons') or [])):
            status='NESTED_PARTITION_PENDING'
        elif not s04:
            status='AWAITING_STEP04'
        elif not s05:
            status='AWAITING_STEP05'
        elif s05.get('qa_status')=='REVIEW':
            status='DEFERRED_DEM_REVIEW'
        elif not s06:
            status='AWAITING_STEP06'
        elif export_status=='BLOCKED':
            status='BLOCKED_STEP06'
        elif not s07:
            status='AWAITING_STEP07'
        elif geom_status=='BLOCKED':
            status='STRUCTURAL_QA_FAILED'
        else:
            status='PROVISIONAL_VALIDATED'
        rows.append({
            'site_id':site,'tdx_id':record.get('tdx_id',''),
            'snap_status':snap,'source_input_ready':source_ready,
            'step04_done':bool(s04),'step05_done':bool(s05),
            'step05_qa_status':(s05 or {}).get('qa_status',''),
            'step05_area_difference_pct':(s05 or {}).get('refined_area_difference_pct',''),
            'shared_related_gauge_count':(s04 or {}).get('shared_related_gauge_count',''),
            'step06_status':export_status,'step07_status':geom_status,
            'workflow_status':status,'final_export_authorized':False,
        })
    keys=list(rows[0]) if rows else []
    with output.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=keys)
        writer.writeheader();writer.writerows(rows)
    summary={'gauges':len(rows),'workflow_counts':dict(sorted(Counter(r['workflow_status'] for r in rows).items())),
             'source_input_ready':sum(r['source_input_ready'] for r in rows),
             'step04_done':sum(r['step04_done'] for r in rows),
             'step05_done':sum(r['step05_done'] for r in rows),
             'step07_done':sum(r['step07_status']!='NOT_RUN' for r in rows),
             'final_export_authorized':False,
             'status_csv':str(output)}
    summary_path=output.with_suffix('.summary.json')
    summary_path.write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    main()
