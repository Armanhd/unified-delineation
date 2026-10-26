#!/usr/bin/env python3
"""Audit proposed TDX nested reach IDs across ALL available regional SQLite indexes.

This is a read-only CANDIDATE audit. It verifies exact decimal IDs and integer-x10
codes, source native LINKNO collisions (regional and cross-regional), parent IDs,
regional split-plan consistency and inventory coverage. It never authorizes output
for CWARHM/SUMMA. Regional LINKNO values can repeat in independent regions, so
cross-region clashes are reported separately from within-region fatal clashes.

Usage:
 python 1_delineation/common/09_global_tdx_id_audit.py
 python 1_delineation/common/09_global_tdx_id_audit.py --expected-region-count 15
"""
from __future__ import annotations
import argparse
import csv
import json
import re
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REGEX_ID = re.compile(r'^(0|[1-9][0-9]*)\.([1-9])$')


def read_csv(path):
    if not path.is_file():
        return []
    with path.open(newline='', encoding='utf-8-sig') as fh:
        return list(csv.DictReader(fh))


def write_csv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=columns, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def index_signature(con, region):
    row = con.execute("SELECT value FROM metadata WHERE key='signature'").fetchone()
    if row is None:
        return 'MISSING_SIGNATURE'
    try:
        signature = json.loads(row[0])
        if str(signature.get('tdx_id', region)) != region:
            return 'REGION_MISMATCH'
        for kind in ('stream', 'basin'):
            name = signature.get(kind + '_path')
            if not name:
                continue
            path = Path(name)
            if not path.is_file():
                return 'SOURCE_MISSING_' + kind.upper()
            st = path.stat()
            if (kind + '_size') in signature and st.st_size != signature[kind + '_size']:
                return 'SOURCE_REPLACED_' + kind.upper()
            if (kind + '_mtime_ns') in signature and st.st_mtime_ns != signature[kind + '_mtime_ns']:
                return 'SOURCE_CHANGED_' + kind.upper()
        return 'PASS'
    except (ValueError, OSError, sqlite3.Error):
        return 'INVALID_SIGNATURE'


def norm(value):
    s = str(value or '').strip()
    return s[:-2] if s.endswith('.0') and s[:-2].isdigit() else s


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project-root', type=Path, default=ROOT)
    ap.add_argument('--output-dir', type=Path)
    ap.add_argument('--expected-region-count', type=int, default=15,
                    help='Expected TDX regions in North America (coverage diagnostic; no fictitious IDs)')
    ap.add_argument('--expected-regions', default='', help='Comma-separated mandatory TDX region IDs')
    ap.add_argument('--strict', action='store_true', help='Exit nonzero for fatal collisions or missing inventory')
    args = ap.parse_args()
    root = args.project_root.resolve()
    output = args.output_dir or (root / 'indexes/gauge_catalogs/tdx_shared_network/global_id_audit')
    output.mkdir(parents=True, exist_ok=True)
    plans = root / 'indexes/gauge_catalogs/tdx_shared_network/nested_split_plans'
    indexroot = root / 'indexes/tdx_regional_topology'
    manifest = root / '0_control_files/tdx_camels_spat/batch_controls_manifest.csv'
    manifest_rows = read_csv(manifest)
    ready_regions = {norm(r.get('tdx_id')) for r in manifest_rows
                     if str(r.get('source_input_ready', '')).lower() == 'true' and norm(r.get('tdx_id'))}
    mandatory = {x.strip() for x in args.expected_regions.split(',') if x.strip()}
    expected_from_manifest = {norm(r.get('tdx_id')) for r in manifest_rows if norm(r.get('tdx_id'))}
    inventory = sorted(indexroot.glob('*_topology.sqlite'))
    indexfiles = {p.name.removesuffix('_topology.sqlite'): p for p in inventory}
    planfiles = {p.name.removesuffix('_nested_candidate_pieces.csv'): p
                 for p in plans.glob('*_nested_candidate_pieces.csv')}
    regions = sorted(set(indexfiles) | set(planfiles) | ready_regions | mandatory)
    rows = []
    by_region = defaultdict(list)
    regions_blocked = set()
    region_warnings = defaultdict(list)
    native_scan = {}
    planned_parents = {}
    candidate_by_region = defaultdict(list)
    code_by_region = defaultdict(dict)
    for region in regions:
        parents_file = plans / f'{region}_nested_split_parents.csv'
        pieces_file = plans / f'{region}_nested_candidate_pieces.csv'
        if region in ready_regions and not parents_file.is_file():
            region_warnings[region].append('REGIONAL_PLAN_MISSING')
            regions_blocked.add(region)
        parents_rows = read_csv(parents_file)
        parents_map = {norm(r.get('parent_tdx_reach_id')): r for r in parents_rows}
        planned_parents[region] = parents_map
        if not indexfiles.get(region):
            region_warnings[region].append('REGIONAL_TOPOLOGY_INDEX_MISSING')
            if region in ready_regions:
                regions_blocked.add(region)
        for r in read_csv(pieces_file):
            text = norm(r.get('candidate_comid_text'))
            par = norm(r.get('parent_tdx_reach_id'))
            rec = {'tdx_id':region, 'parent_tdx_reach_id':par, 'candidate_comid_text':text,
                   'candidate_integer_x10':norm(r.get('candidate_integer_x10')),
                   'qa_flags': norm(r.get('qa_flags')), 'checks':[], 'regional_status':'PASS'}
            match = REGEX_ID.fullmatch(text)
            if not match or match.group(1) != par:
                rec['checks'].append('INVALID_DECIMAL_ID_OR_PARENT')
            else:
                code = int(match.group(1)) * 10 + int(match.group(2))
                rec['integer_code'] = code
                if rec['candidate_integer_x10'] != str(code):
                    rec['checks'].append('X10_MAPPING_MISMATCH')
                if code in code_by_region[region]:
                    rec['checks'].append('DUPLICATE_CANDIDATE_CODE_IN_REGION')
                    code_by_region[region][code]['checks'].append('DUPLICATE_CANDIDATE_CODE_IN_REGION')
                else:
                    code_by_region[region][code] = rec
            if par not in parents_map:
                rec['checks'].append('CANDIDATE_PARENT_ABSENT_FROM_SPLIT_PLAN')
            elif text not in parents_map[par].get('candidate_piece_ids','').split(';'):
                rec['checks'].append('CANDIDATE_NOT_LISTED_IN_PARENT_PLAN')
            candidate_by_region[region].append(rec)
            rows.append(rec)
    # Candidate x10 integers are only a lossless *diagnostic* representation.
    # CWARHM nested previews actually export decimal native/fractional IDs.
    # Thus x10 == some native LINKNO is not an actual ID collision in that
    # export representation. Keep such findings as encoding warnings.
    # Only query actual native IDs referenced by candidate integer-x10 codes.
    # The SQLite native reach indexes are PK-indexed; no huge multi-region load.
    all_codes = sorted({str(r['integer_code']) for r in rows if 'integer_code' in r})
    all_parent_ids = sorted({r['parent_tdx_reach_id'] for r in rows if r['parent_tdx_reach_id']})
    native_hits = defaultdict(set)
    parent_hits = defaultdict(set)
    for region, filename in sorted(indexfiles.items()):
        try:
            con = sqlite3.connect(filename.resolve().as_uri() + '?mode=ro', uri=True)
            try:
                signature_status = index_signature(con, region)
                native_scan[region] = {'status':signature_status,
                                       'reach_count':con.execute('SELECT COUNT(*) FROM links').fetchone()[0]}
                if signature_status != 'PASS':
                    region_warnings[region].append('INDEX_' + signature_status)
                    if region in ready_regions:
                        regions_blocked.add(region)
                for values, dest in ((all_codes, native_hits), (all_parent_ids, parent_hits)):
                    for i in range(0, len(values), 450):
                        chunk = values[i:i + 450]
                        if not chunk:
                            continue
                        found = con.execute('SELECT rid FROM links WHERE rid IN (' +
                                            ','.join('?' for _ in chunk) + ')', chunk)
                        for (rid,) in found:
                            dest[str(rid)].add(region)
            finally:
                con.close()
        except (sqlite3.Error, OSError, ValueError) as exc:
            native_scan[region] = {'status':'UNREADABLE_INDEX', 'error':str(exc)}
            region_warnings[region].append('INDEX_UNREADABLE')
            if region in ready_regions:
                regions_blocked.add(region)
    for rec in rows:
        region = rec['tdx_id']
        par = rec['parent_tdx_reach_id']
        code = str(rec.get('integer_code',''))
        if region not in parent_hits.get(par,set()):
            rec['checks'].append('ORIGINAL_PARENT_MISSING_FROM_REGIONAL_NATIVE_INDEX')
        same_region_collision = region in native_hits.get(code,set())
        if same_region_collision:
            rec['checks'].append('ALTERNATIVE_MIXED_X10_COLLIDES_WITH_NATIVE_LINKNO_IN_REGION_REVIEW')
        cross = sorted(native_hits.get(code,set()) - {region})
        rec['cross_region_native_id_hits'] = ';'.join(cross)
        if cross:
            rec['checks'].append('CROSS_REGION_NATIVE_LINKNO_EQUAL_TO_X10_CODE_REVIEW')
        rec['checks'] = sorted(set(rec['checks']))
        advisory_checks = {
            'CROSS_REGION_NATIVE_LINKNO_EQUAL_TO_X10_CODE_REVIEW',
            'ALTERNATIVE_MIXED_X10_COLLIDES_WITH_NATIVE_LINKNO_IN_REGION_REVIEW',
        }
        fatal = [x for x in rec['checks'] if x not in advisory_checks]
        rec['regional_status'] = 'BLOCKED' if fatal else 'REVIEW' if rec['checks'] else 'PASS'
        rec['qa_flags'] = ';'.join(filter(None,(rec['qa_flags'],';'.join(rec['checks']))))
        rec['checks'] = ';'.join(rec['checks'])
        if fatal:
            regions_blocked.add(region)
        by_region[region].append(rec)
    # Cross-region fractional duplicates matter for a globally unified ID space,
    # even if independent per-basin CWARHM domains can legitimately reuse LINKNO.
    candidate_code_regions = defaultdict(set)
    for rec in rows:
        if rec.get('integer_code') is not None:
            candidate_code_regions[rec['integer_code']].add(rec['tdx_id'])
    cross_region_candidate_codes = {str(code):sorted(rs) for code,rs in candidate_code_regions.items() if len(rs)>1}
    cross_region_code_count = len(cross_region_candidate_codes)
    # Check plan matches exact candidate inventory (includes 0-split regions).
    for region, parents in planned_parents.items():
        actual = {r['candidate_comid_text'] for r in candidate_by_region[region]}
        for parent, plan in parents.items():
            expected = set(filter(None,plan.get('candidate_piece_ids','').split(';')))
            if not expected.issubset(actual):
                region_warnings[region].append('PARENT_PLAN_PIECES_MISSING:' + parent)
                regions_blocked.add(region)
    missing_ready = sorted(ready_regions - set(indexfiles))
    missing_mandatory = sorted(mandatory - set(indexfiles))
    manifest_no_region = sum(not norm(r.get('tdx_id')) for r in manifest_rows)
    # The *count* of indexes is not proof of geographic completeness: some
    # may cover continents outside North America. Require the explicit, known
    # set of North American region IDs before claiming FULL continental audit.
    inventory_coverage_complete = (bool(mandatory) and
                                   len(mandatory) >= args.expected_region_count and
                                   not missing_mandatory and not missing_ready and
                                   all(x['status']=='PASS' for x in native_scan.values()) and
                                   manifest_no_region == 0)
    # Do not conflate an audit of available index files with a complete audit of
    # all 1426 gauges or all officially expected North American regions.
    summary = {
        'scope':'ALL_AVAILABLE_TDX_REGIONAL_SQLITE_TOPOLOGY_INDEXES',
        'candidate_reach_id_count':len(rows),
        'candidate_parent_count':len({(r['tdx_id'],r['parent_tdx_reach_id']) for r in rows}),
        'candidate_status_counts':dict(Counter(r['regional_status'] for r in rows)),
        'regional_native_index_count':len(indexfiles),
        'regional_native_reach_total':sum(v.get('reach_count',0) for v in native_scan.values()),
        'expected_north_american_region_count':args.expected_region_count,
        'explicit_expected_region_ids_supplied':sorted(mandatory),
        'source_ready_regions':sorted(ready_regions),
        'missing_source_ready_topology_indexes':missing_ready,
        'missing_mandatory_indexes':missing_mandatory,
        'manifest_gauges_without_assigned_region':manifest_no_region,
        'missing_source_ready_plan_files':sorted(r for r in ready_regions if not (plans/f'{r}_nested_split_parents.csv').is_file()),
        'regions_blocked_for_candidate_export':sorted(regions_blocked),
        'regions_with_warnings':{r:sorted(set(v)) for r,v in region_warnings.items() if v},
        'cross_region_candidate_code_reuse_count':cross_region_code_count,
        'cross_region_candidate_code_reuse_sample':dict(list(sorted(cross_region_candidate_codes.items()))[:30]),
        'full_north_american_inventory_coverage_confirmed':inventory_coverage_complete,
        'full_global_id_collision_checked':inventory_coverage_complete,
        'global_unique_reach_id_namespace_approved':False,
        'candidate_export_id_encoding':'EXACT_DECIMAL_TDX_ID_1DP',
        'alternative_mixed_x10_collisions_are_advisory':True,
        'final_production_export_authorized':False,
        'interpretation':'Audit uses exact decimal TDX_ID for preview export. Mixed x10-vs-native LINKNO overlaps are advisory only. Real decimal candidate duplicates, missing parents, missing plans/indexes and invalid candidate mappings remain blockers. Continental production approval is not granted.',
        'regional_index_status':native_scan,
    }
    cols = ['tdx_id','parent_tdx_reach_id','candidate_comid_text','candidate_integer_x10',
            'regional_status','qa_flags','cross_region_native_id_hits']
    write_csv(output/'tdx_candidate_id_audit.csv',rows,cols)
    (output/'tdx_global_id_audit_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k != 'regional_index_status'},indent=2))
    print('SUMMARY:',output/'tdx_global_id_audit_summary.json')
    if args.strict and (regions_blocked or not inventory_coverage_complete):
        raise SystemExit(2)


if __name__ == '__main__':
    main()
