#!/usr/bin/env python3
"""TDX Step 02f: construct a deterministic, GLOBAL candidate cut registry.

Input is the completed Step 02e consolidated gauge-snap CSV. The script is
provider-neutral: it requires only standardized gauge IDs and snapped TDX
attributes; it NEVER opens or modifies TDX source GeoPackages.

The same parent TDX reach is processed ONCE, using the positions of ALL gauges
snapped to it, rather than being split separately for each basin.

IMPORTANT: Candidate numeric suffixes are *proposals*, not approved IDs.
Full upstream/downstream topology, line orientation, headwater decisions,
whole-TDX ID collision checking and Copernicus polygon splitting remain required
before a final CWARHM shapefile export. For uncut reaches retain native LINKNO.

Usage:
 python 1_delineation/tdx/02f_build_shared_reach_registry.py
 python 1_delineation/tdx/02f_build_shared_reach_registry.py --strict
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SNAPS = ROOT / 'indexes/gauge_catalogs/tdx_snapping/global_500m/tdx_gauge_snaps_all.csv'
DEFAULT_CATALOG = ROOT / 'indexes/gauge_catalogs/camels-spat_gauge_catalog_validated.csv'
DEFAULT_OUT = ROOT / 'indexes/gauge_catalogs/tdx_shared_network'

GAUGE_COLUMNS = [
    'site_id', 'status', 'warnings', 'tdx_id', 'parent_tdx_reach_id',
    'downstream_reach_id', 'snap_distance_m', 'expected_area_km2',
    'along_reach_m', 'reach_length_m', 'fraction_from_line_start',
    'n_gauges_on_parent', 'candidate_piece_before_gauge',
    'candidate_piece_after_gauge', 'group_status', 'qa_flags',
    'needs_headwater_or_snap_review', 'final_export_authorized',
]
GROUP_COLUMNS = [
    'tdx_id', 'parent_tdx_reach_id', 'n_gauges', 'site_ids',
    'n_candidate_pieces', 'candidate_piece_ids', 'group_status',
    'qa_flags', 'mean_reach_length_m', 'min_gauge_separation_m',
    'final_export_authorized',
]
PIECE_COLUMNS = [
    'tdx_id', 'parent_tdx_reach_id', 'piece_index_from_line_start',
    'candidate_comid_text', 'candidate_integer_x10',
    'from_fraction', 'to_fraction', 'approx_length_m',
    'left_boundary_site', 'right_boundary_site',
    'orientation', 'group_status', 'final_export_authorized',
]


def rows_from_csv(path: Path):
    with path.open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f'Empty CSV header: {path}')
        return list(reader)


def write_csv(path, rows, fields):
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def as_number(value, field, site):
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(f'{site}: invalid {field} = {value!r}') from None
    if not math.isfinite(number):
        raise ValueError(f'{site}: non-finite {field} = {value!r}')
    return number


def native_id(value):
    value = str(value if value is not None else '').strip()
    if value.endswith('.0'):
        value = value[:-2]
    if not value.isdecimal():
        return None
    return str(int(value))


def yes_no(value):
    return str(value or 'true').strip().lower() in ('true', 'yes', '1')


def grouped_registry(snaps, args):
    groups = defaultdict(list)
    gauge_records = {}
    known_native_ids = set()
    status_counts = Counter()
    for row in snaps:
        site = row['site_id'].strip()
        if site in gauge_records:
            raise ValueError(f'Duplicate gauge in snap results: {site}')
        status = row.get('status', '').strip()
        status_counts[status] += 1
        region, parent = native_id(row.get('tdx_id')), native_id(row.get('reach_id'))
        base = {
            'site_id': site, 'status': status, 'warnings': row.get('warnings', ''),
            'tdx_id': region or '', 'parent_tdx_reach_id': parent or '',
            'downstream_reach_id': native_id(row.get('downstream_reach_id')) or '',
            'snap_distance_m': row.get('snap_distance_m', ''),
            'expected_area_km2': row.get('expected_area_km2', ''),
            'along_reach_m': row.get('along_reach_m', ''),
            'reach_length_m': row.get('reach_length_m', ''),
            'fraction_from_line_start': '',
            'n_gauges_on_parent': 0,
            'candidate_piece_before_gauge': '',
            'candidate_piece_after_gauge': '',
            'group_status': 'NO_PARENT_REACH',
            'qa_flags': 'UNRESOLVED_SNAP' if status == 'FAILED' else '',
            'needs_headwater_or_snap_review': 'true' if status != 'PROVISIONAL' else 'false',
            'final_export_authorized': 'false',
        }
        gauge_records[site] = base
        if region and parent:
            known_native_ids.add(parent)
            issues = []
            try:
                along = as_number(row.get('along_reach_m'), 'along_reach_m', site)
                length = as_number(row.get('reach_length_m'), 'reach_length_m', site)
                if length <= 0 or along < 0 or along > length + 1e-5:
                    issues.append('INVALID_REACH_POSITION')
                    frac = None
                else:
                    frac = max(0., min(1., along / length))
            except ValueError:
                issues.append('INVALID_REACH_POSITION')
                length, frac = None, None
            groups[(region, parent)].append({
                'site': site, 'status': status, 'length': length,
                'fraction': frac, 'issues': issues,
            })
        elif status != 'FAILED':
            base['qa_flags'] = 'MISSING_PARENT_REACH'
            base['needs_headwater_or_snap_review'] = 'true'
    return groups, gauge_records, known_native_ids, status_counts


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--snaps', type=Path, default=DEFAULT_SNAPS)
    ap.add_argument('--catalog', type=Path, default=DEFAULT_CATALOG)
    ap.add_argument('--output-dir', type=Path, default=DEFAULT_OUT)
    ap.add_argument('--min-endpoint-m', type=float, default=30.)
    ap.add_argument('--min-cut-separation-m', type=float, default=30.)
    ap.add_argument('--max-length-difference-pct', type=float, default=5.)
    ap.add_argument('--strict', action='store_true',
                    help='Fail unless every snap is PROVISIONAL and every group passes QA. '
                         'Still does not approve export or finalize reach IDs.')
    args = ap.parse_args()
    if min(args.min_endpoint_m, args.min_cut_separation_m, args.max_length_difference_pct) < 0:
        ap.error('Minimum distances and length tolerance must be >= 0')

    snaps = rows_from_csv(args.snaps)
    catalog = rows_from_csv(args.catalog)
    if not snaps or not catalog:
        raise ValueError('Catalogue and snaps must not be empty')
    for key in ('site_id', 'tdx_id', 'reach_id', 'status', 'along_reach_m', 'reach_length_m'):
        if key not in snaps[0]:
            raise ValueError(f'Missing snap column: {key}')
    sites = [r['site_id'].strip() for r in catalog]
    if len(sites) != len(set(sites)):
        raise ValueError('Duplicate site ID in catalogue')
    if len(snaps) != len(sites) or {r['site_id'].strip() for r in snaps} != set(sites):
        raise ValueError('Consolidated snaps do not match the complete catalogue')

    process_target = {r['site_id'].strip(): yes_no(r.get('process_basin')) for r in catalog}
    groups, gauge_records, known_native_ids, counts = grouped_registry(snaps, args)
    pieces = []
    group_rows = []
    proposed_codes = defaultdict(list)
    for (region, parent), members in sorted(groups.items(), key=lambda kv: (int(kv[0][0]), int(kv[0][1]))):
        members.sort(key=lambda x: (float('inf') if x['fraction'] is None else x['fraction'], x['site']))
        issues = set()
        for member in members:
            issues.update(member['issues'])
            if member['status'] != 'PROVISIONAL':
                issues.add('UNRESOLVED_MEMBER_SNAP')
        lengths = [m['length'] for m in members if m['length'] is not None and m['length'] > 0]
        reference_length = float(median(lengths)) if lengths else None
        if not lengths or len(lengths) != len(members):
            issues.add('MISSING_REACH_LENGTH')
        elif (max(lengths) - min(lengths)) / reference_length * 100. > args.max_length_difference_pct:
            issues.add('REACH_LENGTH_INCONSISTENT')
        valid = [m for m in members if m['fraction'] is not None]
        if len(valid) != len(members):
            issues.add('INVALID_REACH_POSITION')
        gap = None
        if reference_length:
            for m in valid:
                if min(m['fraction'], 1. - m['fraction']) * reference_length < args.min_endpoint_m:
                    issues.add('GAUGE_TOO_CLOSE_TO_REACH_ENDPOINT')
            gaps = [
                (b['fraction'] - a['fraction']) * reference_length
                for a, b in zip(valid, valid[1:])
            ]
            gap = min(gaps) if gaps else None
            if gap is not None and gap < args.min_cut_separation_m:
                issues.add('COINCIDENT_OR_NEARBY_GAUGE_CUTS')
        n_parts = len(members) + 1
        if n_parts > 9:
            issues.add('MORE_THAN_NINE_ONE_DECIMAL_SEGMENTS')
        if len(members) == 1:
            # Still a candidate: a different downstream basin may need the
            # split reach although only one gauge was assigned this parent.
            pass

        # Numeric suffix convention is deterministic from LINE START to END,
        # not upstream to downstream until river orientation is confirmed.
        ids = [f'{parent}.{i}' for i in range(1, n_parts + 1)] if n_parts <= 9 else []
        for s in ids:
            code = int(parent) * 10 + int(s.rsplit('.', 1)[-1])
            proposed_codes[code].append((region, parent, s))
            if str(code) in known_native_ids:
                issues.add('POTENTIAL_NATIVE_ID_X10_COLLISION')
            if code > 2**53 - 1:
                issues.add('ID_EXCEEDS_SAFE_FLOAT_INTEGER')
        # Postpone emitting pieces for malformed input: only emit valid ranges.
        if not issues.intersection({'INVALID_REACH_POSITION', 'MISSING_REACH_LENGTH', 'MORE_THAN_NINE_ONE_DECIMAL_SEGMENTS'}):
            cuts = [0.] + [m['fraction'] for m in members] + [1.]
            for i in range(n_parts):
                a, b = cuts[i], cuts[i + 1]
                item = {
                    'tdx_id': region,
                    'parent_tdx_reach_id': parent,
                    'piece_index_from_line_start': i + 1,
                    'candidate_comid_text': ids[i],
                    'candidate_integer_x10': int(parent) * 10 + i + 1,
                    'from_fraction': round(a, 12),
                    'to_fraction': round(b, 12),
                    'approx_length_m': round(reference_length * (b - a), 3),
                    'left_boundary_site': members[i-1]['site'] if i > 0 else '',
                    'right_boundary_site': members[i]['site'] if i < len(members) else '',
                    'orientation': 'UNRESOLVED',
                    'group_status': '',
                    'final_export_authorized': 'false',
                }
                pieces.append(item)
        for i, member in enumerate(members):
            gr = gauge_records[member['site']]
            gr['n_gauges_on_parent'] = len(members)
            gr['fraction_from_line_start'] = (
                round(member['fraction'], 12) if member['fraction'] is not None else ''
            )
            gr['candidate_piece_before_gauge'] = ids[i] if ids else ''
            gr['candidate_piece_after_gauge'] = ids[i+1] if ids else ''
        group_rows.append({
            'tdx_id': region, 'parent_tdx_reach_id': parent,
            'n_gauges': len(members), 'site_ids': ';'.join(m['site'] for m in members),
            'n_candidate_pieces': n_parts, 'candidate_piece_ids': ';'.join(ids),
            'group_status': '', 'qa_flags': '',
            'mean_reach_length_m': round(reference_length, 3) if reference_length else '',
            'min_gauge_separation_m': round(gap, 3) if gap is not None else '',
            'final_export_authorized': 'false',
            '_issues': issues,
        })

    # Check candidate suffixes against other candidate suffixes across regions.
    for code, owners in proposed_codes.items():
        if len(owners) > 1:
            for row in group_rows:
                if any(row['tdx_id'] == r and row['parent_tdx_reach_id'] == p for r, p, _ in owners):
                    row['_issues'].add('CROSS_REGION_X10_ID_COLLISION')

    group_status_by_key = {}
    for row in group_rows:
        issues = row.pop('_issues')
        status = 'CANDIDATE_NEEDS_REVIEW' if issues else 'CANDIDATE_READY_FOR_TOPOLOGY'
        row['group_status'] = status
        row['qa_flags'] = ';'.join(sorted(issues))
        group_status_by_key[(row['tdx_id'], row['parent_tdx_reach_id'])] = (status, row['qa_flags'])
    for item in pieces:
        item['group_status'] = group_status_by_key[(item['tdx_id'], item['parent_tdx_reach_id'])][0]
    for row in gauge_records.values():
        key = (row['tdx_id'], row['parent_tdx_reach_id'])
        if key in group_status_by_key:
            row['group_status'], row['qa_flags'] = group_status_by_key[key]
            if row['group_status'] != 'CANDIDATE_READY_FOR_TOPOLOGY':
                row['needs_headwater_or_snap_review'] = 'true'
        elif row['status'] == 'FAILED':
            row['group_status'] = 'NO_PARENT_REACH'
            row['qa_flags'] = 'UNRESOLVED_SNAP'
            row['needs_headwater_or_snap_review'] = 'true'

    group_status_counts = Counter(g['group_status'] for g in group_rows)
    blocked = [r for r in gauge_records.values() if r['needs_headwater_or_snap_review'] == 'true']
    signature = hashlib.sha256(args.snaps.read_bytes() + b'\0' + args.catalog.read_bytes()).hexdigest()
    summary = {
        'status': 'CANDIDATE_ONLY_NOT_APPROVED_FOR_EXPORT',
        'input_gauges': len(snaps),
        'target_basins': sum(process_target.values()),
        'snap_status_counts': dict(sorted(counts.items())),
        'unique_gauge_parent_reaches': len(group_rows),
        'parents_with_multiple_gauges': sum(int(r['n_gauges']) > 1 for r in group_rows),
        'parents_with_single_gauge': sum(int(r['n_gauges']) == 1 for r in group_rows),
        'candidate_piece_count': len(pieces),
        'parent_group_status_counts': dict(sorted(group_status_counts.items())),
        'gauges_needing_snap_or_parent_review': len(blocked),
        'fingerprint_sha256_of_input_snaps_and_catalog': signature,
        'whole_network_upstream_downstream_nesting_tested': False,
        'native_tdx_global_ID_collision_audit_complete': False,
        'fractional_ids_finalized': False,
        'source_gpkg_files_modified': False,
        'notes': [
            'Fractions .1 to .9 are draft names measured from the original LINE START, not flow direction.',
            'A parent with one gauge is included because it may be nested inside another downstream basin.',
            'Do not use candidate fractional IDs in CWARHM until orientation, full topology and polygon partitions pass validation.',
            'Only basins requiring no split keep their original integer TDX reach IDs.',
            'Rerun this registry from the complete, corrected gauge snaps after TDX source restoration.',
        ],
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / 'tdx_gauge_reach_registry.csv',
              [gauge_records[k] for k in sorted(gauge_records)], GAUGE_COLUMNS)
    write_csv(args.output_dir / 'tdx_shared_parent_groups.csv', group_rows, GROUP_COLUMNS)
    write_csv(args.output_dir / 'tdx_candidate_split_pieces.csv', pieces, PIECE_COLUMNS)
    (args.output_dir / 'tdx_shared_registry_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))
    print('OUTPUT:', args.output_dir)

    if args.strict and (blocked or summary['parent_group_status_counts'].get('CANDIDATE_NEEDS_REVIEW')):
        raise SystemExit('STRICT QA FAILED: unresolved gauge or parent reach groups; candidate outputs are retained')


if __name__ == '__main__':
    main()
