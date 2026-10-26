#!/usr/bin/env python3
"""Consolidate Step 02b Slurm array output; audit gauge coverage and candidate IDs.

Run from project root (script can live in 1_delineation/tdx/).
Only output directories are written; original batches and TDX data stay unchanged.

Example:
  python 1_delineation/tdx/02e_consolidate_gauge_snaps.py
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from collections import Counter
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
SNAP_COLUMNS = [
    'site_id', 'latitude', 'longitude', 'expected_area_km2', 'tdx_id',
    'reach_id', 'downstream_reach_id', 'snapped_lat', 'snapped_lon',
    'snap_distance_m', 'along_reach_m', 'reach_length_m', 'end_distance_m',
    'tdx_upstream_area_km2', 'tdx_downstream_area_km2', 'stream_order',
    'n_candidates', 'n_regions_with_candidates', 'status', 'warnings',
]
CANDIDATE_COLUMNS = [
    'site_id', 'rank', 'tdx_id', 'reach_id', 'downstream_reach_id',
    'distance_m', 'snapped_lat', 'snapped_lon', 'along_reach_m',
    'reach_length_m', 'end_distance_m', 'upstream_area_km2',
    'downstream_area_km2', 'stream_order', 'area_outside_range_pct',
]
STATES = {'PROVISIONAL', 'REVIEW', 'FAILED'}


def read_csv(path: Path, required: list[str]) -> list[dict[str, str]]:
    with path.open(newline='', encoding='utf-8-sig') as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames or any(c not in reader.fieldnames for c in required):
            raise ValueError(f'Incorrect/missing columns in {path}')
        return list(reader)


def write_csv(path: Path, rows: list[dict[str, str]], fields: list[str]):
    with path.open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--catalog', type=Path, default=PROJECT / 'indexes/gauge_catalogs/camels-spat_gauge_catalog_validated.csv')
    p.add_argument('--batch-dir', type=Path, default=PROJECT / 'indexes/gauge_catalogs/tdx_snapping/batch_500m')
    p.add_argument('--output-dir', type=Path, default=PROJECT / 'indexes/gauge_catalogs/tdx_snapping/global_500m')
    p.add_argument('--batch-size', type=int, default=20)
    p.add_argument('--expected-radius-m', type=float, default=500)
    args = p.parse_args()
    if args.batch_size < 1:
        p.error('batch-size must be positive')

    catalog = read_csv(args.catalog, ['site_id'])
    expected_ids = sorted(row['site_id'].strip() for row in catalog)
    if len(set(expected_ids)) != len(expected_ids) or any(not x for x in expected_ids):
        p.error('Catalogue has missing or duplicate site IDs')
    if not expected_ids:
        p.error('Catalogue is empty')

    n = len(expected_ids)
    expected_count = math.ceil(n / args.batch_size)
    problems: list[str] = []
    merged_snaps: list[dict[str, str]] = []
    merged_candidates: list[dict[str, str]] = []
    batch_report = []
    expected_files = set()

    for index in range(expected_count):
        start = index * args.batch_size
        stop = min(start + args.batch_size, n)
        suffix = f'{start:04d}_{stop:04d}'
        paths = {
            'snaps': args.batch_dir / f'tdx_gauge_snaps_{suffix}.csv',
            'candidates': args.batch_dir / f'tdx_snap_candidates_{suffix}.csv',
            'summary': args.batch_dir / f'tdx_snap_summary_{suffix}.json',
        }
        expected_files.update(paths.values())
        absent = [str(path) for path in paths.values() if not path.is_file()]
        if absent:
            problems.append(f'Batch {suffix}: missing files: {absent}')
            continue
        try:
            snaps = read_csv(paths['snaps'], SNAP_COLUMNS)
            candidates = read_csv(paths['candidates'], CANDIDATE_COLUMNS)
            summary = json.loads(paths['summary'].read_text(encoding='utf-8'))
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            problems.append(f'Batch {suffix}: unreadable content: {exc}')
            continue

        expected_slice = expected_ids[start:stop]
        ids = [row['site_id'].strip() for row in snaps]
        if ids != expected_slice:
            problems.append(
                f'Batch {suffix}: gauge IDs or order differ from catalogue; '
                f'expected {len(expected_slice)} and found {len(ids)} rows'
            )
        if summary.get('n_gauges') != len(snaps):
            problems.append(f'Batch {suffix}: summary gauge count mismatch')
        if Counter(row['status'] for row in snaps) != Counter(summary.get('status_counts', {})):
            problems.append(f'Batch {suffix}: status counts inconsistent with summary')
        if summary.get('search_radius_m') != args.expected_radius_m:
            problems.append(f'Batch {suffix}: search radius differs from expected {args.expected_radius_m} m')
        if Path(summary.get('catalog', '')).name != args.catalog.name:
            problems.append(f'Batch {suffix}: source catalogue filename mismatch')

        candidate_by_id: dict[str, list[dict[str, str]]] = {}
        for c in candidates:
            candidate_by_id.setdefault(c['site_id'], []).append(c)
            if c['site_id'] not in set(expected_slice):
                problems.append(f'Batch {suffix}: unknown candidate site {c["site_id"]}')

        for row in snaps:
            site = row['site_id']
            if row['status'] not in STATES:
                problems.append(f'Batch {suffix}: {site} has unexpected status {row["status"]}')
            ranked = candidate_by_id.get(site, [])
            if row['status'] == 'FAILED':
                if row['reach_id']:
                    problems.append(f'Batch {suffix}: {site} is FAILED but has a selected reach ID')
                if ranked:
                    problems.append(f'Batch {suffix}: {site} is FAILED but has candidate reaches')
            else:
                if not row['reach_id'] or not row['tdx_id']:
                    problems.append(f'Batch {suffix}: {site} has status {row["status"]} but no TDX reach')
                first = next((c for c in ranked if c['rank'] == '1'), None)
                if first is None or (first['tdx_id'], first['reach_id']) != (row['tdx_id'], row['reach_id']):
                    problems.append(f'Batch {suffix}: {site} selected reach disagrees with rank-1 candidate')

        merged_snaps.extend(snaps)
        merged_candidates.extend(candidates)
        batch_report.append({'batch': suffix, 'n_gauges': len(snaps), 'status_counts': dict(Counter(s['status'] for s in snaps))})

    unexpected = sorted(path.name for pattern in ('tdx_gauge_snaps_*.csv', 'tdx_snap_candidates_*.csv', 'tdx_snap_summary_*.json')
                        for path in args.batch_dir.glob(pattern) if path not in expected_files)
    if unexpected:
        problems.append(f'Unexpected batch files in input directory: {unexpected[:15]}')

    actual_ids = [r['site_id'] for r in merged_snaps]
    counts = Counter(actual_ids)
    duplicates = sorted(k for k, v in counts.items() if v != 1)
    missing = sorted(set(expected_ids) - set(actual_ids))
    extra = sorted(set(actual_ids) - set(expected_ids))
    if duplicates:
        problems.append(f'Duplicate gauge IDs: {duplicates[:20]}')
    if missing:
        problems.append(f'Missing gauge IDs: {missing[:20]}')
    if extra:
        problems.append(f'Unexpected gauge IDs: {extra[:20]}')
    if len(merged_snaps) != n:
        problems.append(f'Expected {n} merged rows, got {len(merged_snaps)}')

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = sorted(merged_snaps, key=lambda r: r['site_id'])
    write_csv(args.output_dir / 'tdx_gauge_snaps_all.csv', rows, SNAP_COLUMNS)
    write_csv(args.output_dir / 'tdx_snap_candidates_all.csv', merged_candidates, CANDIDATE_COLUMNS)
    for status, label in (('PROVISIONAL', 'provisional'), ('REVIEW', 'review'), ('FAILED', 'failed')):
        write_csv(args.output_dir / f'tdx_gauge_snaps_{label}.csv', [r for r in rows if r['status'] == status], SNAP_COLUMNS)

    status_counts = Counter(row['status'] for row in rows)
    summary = {
        'validation': 'PASS' if not problems else 'FAIL',
        'catalog': str(args.catalog),
        'batch_directory': str(args.batch_dir),
        'expected_gauges': n,
        'merged_gauge_rows': len(rows),
        'expected_batches': expected_count,
        'processed_batches': len(batch_report),
        'candidate_rows': len(merged_candidates),
        'status_counts': {s: status_counts.get(s, 0) for s in ('PROVISIONAL', 'REVIEW', 'FAILED')},
        'missing_gauge_ids': missing,
        'duplicate_gauge_ids': duplicates,
        'unexpected_gauge_ids': extra,
        'problems': problems,
        'notes': [
            'A PASS validates batch integrity only, NOT hydrological correctness of provisional gauge snapping.',
            'Do not use REVIEW or FAILED outlets for automatic reach splitting until resolved.',
            'The output is provisional and the original TDX source data remain unchanged.',
        ],
        'batches': batch_report,
    }
    (args.output_dir / 'tdx_global_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({key: val for key, val in summary.items() if key != 'batches'}, indent=2))
    print(f'Outputs: {args.output_dir}')
    return 0 if not problems else 1


if __name__ == '__main__':
    sys.exit(main())
