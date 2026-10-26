#!/usr/bin/env python3
"""Step 04b: plan deterministic, region-wide candidate TDX nested-gauge cuts.

The regional TDX topology index is read once. All mapped gauges in the same
region are considered together. Native LINKNO is retained for parent reaches
that do not need a shared nested-gauge split. A parent becomes a candidate for
splitting only if (a) it hosts two or more gauges, or (b) one of its gauges is
upstream of another gauge on a different TDX reach.

This planner DOES NOT cut stream geometries or HRU polygons and DOES NOT approve
fractional IDs for export. Step 03 orientation, cut placement, Copernicus
partition, full inventory collision checks, and the repaired TDX source regions
must be validated first.

Usage:
 python 1_delineation/tdx/04b_plan_shared_nested_splits.py --region 7020024600
 python 1_delineation/tdx/04b_plan_shared_nested_splits.py --region 7020024600 --strict
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REGISTRY = ROOT / 'indexes/gauge_catalogs/tdx_shared_network/tdx_gauge_reach_registry.csv'
TOPOLOGY = ROOT / 'indexes/tdx_regional_topology'
OUT = ROOT / 'indexes/gauge_catalogs/tdx_shared_network/nested_split_plans'

PAIR_COLS = ['tdx_id','upstream_site_id','downstream_site_id','upstream_parent_reach_id',
             'downstream_parent_reach_id','n_downstream_links','upstream_snap_status',
             'downstream_snap_status','qa_flags']
GROUP_COLS = ['tdx_id','parent_tdx_reach_id','cut_reason','n_gauges_on_parent',
              'n_nested_upstream_gauges','nested_upstream_site_ids','all_parent_gauge_ids',
              'n_candidate_pieces','candidate_piece_ids','group_status','qa_flags',
              'final_export_authorized']
PIECE_COLS = ['tdx_id','parent_tdx_reach_id','piece_index_from_line_start',
              'candidate_comid_text','candidate_integer_x10','from_fraction',
              'to_fraction','left_boundary_site','right_boundary_site',
              'approx_length_m','qa_flags','final_export_authorized']


def norm(value):
    value = str(value if value is not None else '').strip()
    if value.endswith('.0') and value[:-2].isdigit():
        value = value[:-2]
    return value or None


def number(v):
    try:
        out = float(v)
        return out if math.isfinite(out) else None
    except (ValueError, TypeError):
        return None


def load_csv(path):
    with path.open(newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        return list(reader), set(reader.fieldnames or [])


def write_csv(path, columns, rows):
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def collect_gauges(registry, region):
    rows, columns = load_csv(registry)
    required = {'site_id','tdx_id','parent_tdx_reach_id','status',
                'fraction_from_line_start','reach_length_m','group_status','qa_flags'}
    missing = required - columns
    if missing:
        raise ValueError(f'Registry missing columns: {sorted(missing)}')
    gauges = {}
    by_reach = defaultdict(list)
    for raw in rows:
        if norm(raw['tdx_id']) != region or raw['status'].upper() == 'FAILED':
            continue
        site = raw['site_id'].strip()
        reach = norm(raw['parent_tdx_reach_id'])
        if not site or not reach:
            continue
        if site in gauges:
            raise ValueError(f'Duplicate site in regional registry: {site}')
        row = dict(raw)
        row['site_id'] = site
        row['parent_tdx_reach_id'] = reach
        gauges[site] = row
        by_reach[reach].append(row)
    return gauges, by_reach


def load_topology(path):
    con = sqlite3.connect(path.resolve().as_uri() + '?mode=ro', uri=True)
    try:
        valid = con.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='links'").fetchone()
        if not valid:
            raise ValueError('Regional SQLite file lacks links table')
        nodes = {}
        for rid, ds in con.execute('SELECT rid, ds FROM links'):
            rid = norm(rid)
            if rid is not None:
                nodes[rid] = norm(ds)
        return nodes
    finally:
        con.close()


def find_pairs(gauges, by_reach, downstream, region, max_links):
    """Find all downstream nested gauges for each upstream reach (no SQL in loop)."""
    pairs = []
    nested_upstream = defaultdict(set)
    failures = []
    for start, source_sites in sorted(by_reach.items()):
        if start not in downstream:
            failures.append(f'Registry parent reach absent from index: {start}')
            continue
        current, visited = start, {start}
        for dist in range(1, max_links + 1):
            nxt = downstream.get(current)
            if nxt in (None, '0', '-1', current):
                break
            if nxt in visited:
                failures.append(f'Cycle on downstream path from {start}, at {nxt}')
                break
            visited.add(nxt)
            current = nxt
            if current in by_reach:
                for up in source_sites:
                    for down in by_reach[current]:
                        if up['site_id'] == down['site_id']:
                            continue
                        qa = []
                        for g in (up, down):
                            if g['status'].upper() != 'PROVISIONAL':
                                qa.append(f'{g["site_id"]}_SNAP_REVIEW')
                            if g['group_status'] != 'CANDIDATE_READY_FOR_TOPOLOGY':
                                qa.append(f'{g["site_id"]}_REGISTRY_REVIEW')
                        pairs.append({
                            'tdx_id':region,
                            'upstream_site_id':up['site_id'],
                            'downstream_site_id':down['site_id'],
                            'upstream_parent_reach_id':start,
                            'downstream_parent_reach_id':current,
                            'n_downstream_links':dist,
                            'upstream_snap_status':up['status'],
                            'downstream_snap_status':down['status'],
                            'qa_flags':';'.join(sorted(set(qa))),
                        })
                        nested_upstream[start].add(up['site_id'])
        else:
            failures.append(f'Max links {max_links} reached from {start}; topology not fully searched')
    return pairs, nested_upstream, failures


def build_plan(region, gauges, by_reach, downstream, pairs, nested_upstream,
               min_endpoint_m, min_separation_m):
    plan, pieces = [], []
    used_codes = defaultdict(list)
    for parent in sorted(by_reach, key=lambda s:(not s.isdigit(), int(s) if s.isdigit() else s)):
        members = by_reach[parent]
        is_multi = len(members) > 1
        is_nested = bool(nested_upstream.get(parent))
        if not (is_multi or is_nested):
            continue  # Native ID must remain untouched for this parent.
        reasons = []
        if is_multi:
            reasons.append('MULTIPLE_GAUGES_SAME_PARENT')
        if is_nested:
            reasons.append('UPSTREAM_OF_ANOTHER_GAUGE')
        errors = set()
        if parent not in downstream:
            errors.add('PARENT_NOT_IN_TOPOLOGY')
        if not parent.isdigit():
            errors.add('NONINTEGER_PARENT_LINKNO')
        for member in members:
            if member['status'].upper() != 'PROVISIONAL':
                errors.add('SNAP_REVIEW')
            if member['group_status'] != 'CANDIDATE_READY_FOR_TOPOLOGY':
                errors.add('GROUP_REVIEW')
            if member.get('qa_flags'):
                errors.add('REGISTRY_QA_FLAGS')
        fracs = [(number(g.get('fraction_from_line_start')),g) for g in members]
        if any(f is None or not (0 < f < 1) for f, _ in fracs):
            errors.add('MISSING_OR_ENDPOINT_FRACTION')
        fracs = sorted(((f,g) for f,g in fracs if f is not None),
                       key=lambda item:(item[0],item[1]['site_id']))
        lengths = [number(g.get('reach_length_m')) for g in members]
        if any(x is None or x <= 0 for x in lengths):
            errors.add('MISSING_REACH_LENGTH')
            ref_length = None
        else:
            ref_length = sum(lengths)/len(lengths)
            if (max(lengths)-min(lengths))/ref_length > 0.05:
                errors.add('REACH_LENGTH_INCONSISTENT')
        if ref_length is not None:
            if any(min(f,1-f)*ref_length < min_endpoint_m for f,_ in fracs):
                errors.add('CUT_NEAR_ENDPOINT')
            if any((b[0]-a[0])*ref_length < min_separation_m for a,b in zip(fracs,fracs[1:])):
                errors.add('CUTS_TOO_CLOSE')
        count = len(members) + 1
        if count > 9:
            errors.add('MORE_THAN_NINE_SEGMENTS')
        ids = [f'{parent}.{i}' for i in range(1,count+1)] if count <= 9 and parent.isdigit() else []
        # WARNING: compatibility check is only within THIS region; full global
        # region inventory must be examined before any final ID authorization.
        for idtext in ids:
            code = int(parent)*10 + int(idtext.rsplit('.',1)[1])
            used_codes[code].append(parent)
            if str(code) in downstream:
                errors.add('ALTERNATIVE_MIXED_X10_COLLISION_REVIEW')
        if len(fracs) == len(members) and len(ids) == count:
            boundaries = [0.] + [f for f,_ in fracs] + [1.]
            for j, (low, high) in enumerate(zip(boundaries, boundaries[1:])):
                item = {
                    'tdx_id':region,
                    'parent_tdx_reach_id':parent,
                    'piece_index_from_line_start':j+1,
                    'candidate_comid_text':ids[j],
                    'candidate_integer_x10':int(parent)*10 + j+1,
                    'from_fraction':round(low,12), 'to_fraction':round(high,12),
                    'left_boundary_site':fracs[j-1][1]['site_id'] if j else '',
                    'right_boundary_site':fracs[j][1]['site_id'] if j<len(fracs) else '',
                    'approx_length_m':round(ref_length*(high-low),3) if ref_length else '',
                    'qa_flags':'', 'final_export_authorized':'false',
                }
                pieces.append(item)
        plan.append({
            'tdx_id':region,'parent_tdx_reach_id':parent,
            'cut_reason':';'.join(reasons),
            'n_gauges_on_parent':len(members),
            'n_nested_upstream_gauges':len(nested_upstream.get(parent,set())),
            'nested_upstream_site_ids':';'.join(sorted(nested_upstream.get(parent,set()))),
            'all_parent_gauge_ids':';'.join(g['site_id'] for _,g in fracs),
            'n_candidate_pieces':count,
            'candidate_piece_ids':';'.join(ids),
            'group_status':'PENDING_VALIDATION',
            'qa_flags':'',
            'final_export_authorized':'false',
            '_issues':errors,
        })
    for code, parents in used_codes.items():
        if len(set(parents))>1:
            for row in plan:
                if row['parent_tdx_reach_id'] in parents:
                    row['_issues'].add('PROPOSED_ID_COLLISION')
    by_parent = {}
    for row in plan:
        row['qa_flags'] = ';'.join(sorted(row.pop('_issues')))
        row['group_status'] = 'REVIEW_REQUIRED' if row['qa_flags'] else 'READY_FOR_ORIENTATION_AND_GEOMETRY_CHECK'
        by_parent[row['parent_tdx_reach_id']] = row
    for item in pieces:
        item['qa_flags'] = by_parent[item['parent_tdx_reach_id']]['qa_flags']
    return plan, pieces


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--region', required=True, help='10-digit TDX region identifier')
    p.add_argument('--registry', type=Path, default=REGISTRY)
    p.add_argument('--topology-index', type=Path, default=None)
    p.add_argument('--output-dir', type=Path, default=OUT)
    p.add_argument('--max-links', type=int, default=100000)
    p.add_argument('--min-endpoint-m', type=float, default=30.)
    p.add_argument('--min-separation-m', type=float, default=30.)
    p.add_argument('--strict', action='store_true', help='Fail if the plan contains any QA issues')
    args = p.parse_args()
    region = norm(args.region)
    if not region or not region.isdigit():
        p.error('--region must be a numeric TDX identifier')
    if args.max_links < 1 or args.min_endpoint_m < 0 or args.min_separation_m < 0:
        p.error('Invalid distance limits')
    topology = args.topology_index or TOPOLOGY / f'{region}_topology.sqlite'
    if not topology.is_file():
        raise FileNotFoundError(f'Build regional topology first: {topology}')
    gauges, by_reach = collect_gauges(args.registry, region)
    if not gauges:
        raise ValueError(f'No mapped gauges in registry for region {region}')
    downstream = load_topology(topology)
    pairs, nested_upstream, failures = find_pairs(gauges, by_reach, downstream, region, args.max_links)
    plan, pieces = build_plan(region, gauges, by_reach, downstream, pairs, nested_upstream,
                              args.min_endpoint_m, args.min_separation_m)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    prefix = f'{region}_nested'
    write_csv(args.output_dir / f'{prefix}_gauge_pairs.csv', PAIR_COLS, sorted(
        pairs, key=lambda r:(r['upstream_site_id'], int(r['n_downstream_links']), r['downstream_site_id'])))
    write_csv(args.output_dir / f'{prefix}_split_parents.csv', GROUP_COLS, plan)
    write_csv(args.output_dir / f'{prefix}_candidate_pieces.csv', PIECE_COLS, pieces)
    problems = {r['parent_tdx_reach_id']:r['qa_flags'] for r in plan if r['qa_flags']}
    summary = {
        'tdx_id':region,'regional_topology_reaches':len(downstream),
        'assigned_gauges':len(gauges),'gauge_parent_reaches':len(by_reach),
        'different_parent_nested_pairs':len(pairs),
        'candidate_parent_reaches_requiring_split':len(plan),
        'candidate_pieces':len(pieces),
        'parents_needing_review':len(problems),
        'flagged_parent_reaches':problems,
        'topology_warnings':failures[:50],
        'regional_native_id_collision_checked':True,
        'candidate_export_id_encoding':'EXACT_DECIMAL_TDX_ID_1DP',
        'mixed_x10_collisions_are_advisory':True,
        'full_global_native_id_collision_checked':False,
        'stream_geometry_split_performed':False,
        'hru_polygon_split_performed':False,
        'fractional_ids_finalized':False,
        'final_export_authorized':False,
        'note':'Global candidate planning in one region; source geometries are unchanged. Full cross-region consistency and physical cut QA remain pending.',
    }
    (args.output_dir / f'{prefix}_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    print('Files:', args.output_dir / f'{prefix}_split_parents.csv')
    if args.strict and (problems or failures):
        raise RuntimeError('Nested-split plan requires QA review; see summary JSON')

if __name__ == '__main__':
    main()
