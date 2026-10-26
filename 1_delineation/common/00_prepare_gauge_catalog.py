
#!/usr/bin/env python3
"""Validate a portable gauge CSV and prepare a deterministic gauge registry.

Required: site_id, latitude, longitude
Optional: expected_area_km2, process_basin (default true)
Coordinates use EPSG:4326. Input files are never modified.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter
from pathlib import Path


def numeric(value, field, site):
    try:
        result = float(value)
    except (TypeError, ValueError):
        raise ValueError(f'{site}: invalid {field}: {value!r}') from None
    if not math.isfinite(result):
        raise ValueError(f'{site}: {field} is not finite')
    return result


def yes_no(value, site):
    value = str(value or 'true').strip().lower()
    if value in ('true', '1', 'yes'):
        return True
    if value in ('false', '0', 'no'):
        return False
    raise ValueError(f'{site}: invalid process_basin: {value!r}')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--catalog', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, default=Path('indexes/gauge_catalogs'))
    args = p.parse_args()

    with args.catalog.open(newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        required = {'site_id', 'latitude', 'longitude'}
        missing = required - set(reader.fieldnames or [])
        if missing:
            p.error(f'Missing CSV columns: {sorted(missing)}')
        records = list(reader)

    if not records:
        raise ValueError('Gauge catalogue is empty')

    gauges = {}
    for line, row in enumerate(records, start=2):
        site = (row['site_id'] or '').strip()
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*', site):
            raise ValueError(f'Line {line}: unsafe or empty site_id: {site!r}')
        if site in gauges:
            raise ValueError(f'Line {line}: duplicate site_id: {site}')

        lat = numeric(row['latitude'], 'latitude', site)
        lon = numeric(row['longitude'], 'longitude', site)
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            raise ValueError(f'{site}: coordinates outside EPSG:4326 bounds')

        area_text = (row.get('expected_area_km2') or '').strip()
        area = numeric(area_text, 'expected_area_km2', site) if area_text else None
        if area is not None and area <= 0:
            raise ValueError(f'{site}: expected_area_km2 must be positive')

        gauges[site] = {
            'site_id': site,
            'latitude': lat,
            'longitude': lon,
            'expected_area_km2': area if area is not None else '',
            'process_basin': yes_no(row.get('process_basin'), site),
        }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = args.catalog.stem
    catalog_out = args.output_dir / f'{stem}_validated.csv'
    summary_out = args.output_dir / f'{stem}_summary.json'

    with catalog_out.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(next(iter(gauges.values()))))
        writer.writeheader()
        writer.writerows(gauges[k] for k in sorted(gauges))

    summary = {
        'source': str(args.catalog.resolve()),
        'n_gauges': len(gauges),
        'n_target_basins': sum(r['process_basin'] for r in gauges.values()),
        'n_reference_only_gauges': sum(not r['process_basin'] for r in gauges.values()),
        'n_with_expected_area': sum(r['expected_area_km2'] != '' for r in gauges.values()),
        'prefix_counts': dict(sorted(Counter(k.split('_')[0] for k in gauges).items())),
        'nesting_status': 'NOT_YET_DETERMINED',
    }
    summary_out.write_text(json.dumps(summary, indent=2) + '\n')

    print(json.dumps(summary, indent=2))
    print(f'Validated catalogue: {catalog_out}')
    print(f'Summary: {summary_out}')


if __name__ == '__main__':
    main()
