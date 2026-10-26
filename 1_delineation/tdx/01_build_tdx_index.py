#!/usr/bin/env python3
"""Build a portable, metadata-based index of TDX hydrofabric regions.

Only GeoPackage metadata is read. No full feature tables are loaded.
The CSV is portable; the local health report is intentionally separate.
"""
import argparse
import csv
import json
import sqlite3
from pathlib import Path


def inspect(path):
    if not path.exists():
        return {'status': 'MISSING', 'error': 'File not found'}
    try:
        # Immutable avoids creating journal sidecars on shared read-only data.
        uri = path.resolve().as_uri() + '?mode=ro&immutable=1'
        with sqlite3.connect(uri, uri=True, timeout=5) as db:
            rows = db.execute('''SELECT c.table_name, c.srs_id,
                c.min_x, c.min_y, c.max_x, c.max_y,
                g.geometry_type_name
                FROM gpkg_contents c LEFT JOIN gpkg_geometry_columns g
                ON c.table_name = g.table_name
                WHERE c.data_type='features' ORDER BY c.table_name''').fetchall()
            if not rows:
                return {'status': 'NO_FEATURE_LAYER', 'error': 'No feature layers'}
            layers = [{'name':r[0], 'srs_id':r[1], 'bounds':list(r[2:6]),
                       'geometry_type':r[6]} for r in rows]
            return {'status':'METADATA_OK', 'layers':layers, 'error':''}
    except (sqlite3.Error, OSError, ValueError) as exc:
        return {'status':'METADATA_ERROR', 'error':str(exc)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tdx-dir', type=Path, required=True)
    parser.add_argument('--index-dir', type=Path, required=True)
    args = parser.parse_args()
    root = args.tdx_dir.expanduser().resolve()
    if not root.is_dir():
        parser.error(f'TDX directory not found: {root}')
    args.index_dir.mkdir(parents=True, exist_ok=True)
    paths = sorted(root.glob('*.gpkg'))
    regions = {}
    for p in paths:
        for suffix, kind in [('-basins.gpkg','basins'), ('-streamnet.gpkg','streamnet')]:
            if p.name.endswith(suffix):
                region = p.name[:-len(suffix)]
                regions.setdefault(region, {})[kind] = p
                break
    portable, health = [], []
    for region in sorted(regions):
        entry = {'tdx_id':region}
        checks = {}
        for kind in ('basins', 'streamnet'):
            p = regions[region].get(kind)
            result = inspect(p) if p else {'status':'MISSING', 'error':'Missing counterpart'}
            checks[kind] = result
            entry[f'{kind}_filename'] = p.name if p else f'{region}-{kind}.gpkg'
            # Keep all discovered layers to avoid assuming a fixed name.
            entry[f'{kind}_layers_json'] = json.dumps(result.get('layers', []), separators=(',', ':'))
            health.append({'tdx_id':region, 'type':kind,
                           'filename':entry[f'{kind}_filename'],
                           'status':result['status'], 'error':result.get('error','')})
        # A region's spatial envelope can be derived from either feature layer.
        # Bounds are retained in native CRS, never mislabeled as longitude/latitude.
        candidates = [layer for kind in ('basins','streamnet')
                      for layer in checks[kind].get('layers', [])
                      if all(v is not None for v in layer['bounds'])]
        if candidates and len({str(x['srs_id']) for x in candidates}) == 1:
            entry['bbox_srs_id'] = candidates[0]['srs_id']
            entry['min_x'] = min(x['bounds'][0] for x in candidates)
            entry['min_y'] = min(x['bounds'][1] for x in candidates)
            entry['max_x'] = max(x['bounds'][2] for x in candidates)
            entry['max_y'] = max(x['bounds'][3] for x in candidates)
        else:
            entry.update({k:'' for k in ('bbox_srs_id','min_x','min_y','max_x','max_y')})
        portable.append(entry)
    fields = ['tdx_id','basins_filename','streamnet_filename',
              'basins_layers_json','streamnet_layers_json','bbox_srs_id',
              'min_x','min_y','max_x','max_y']
    with (args.index_dir/'tdx_regions.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(portable)
    with (args.index_dir/'tdx_local_health.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['tdx_id','type','filename','status','error'])
        writer.writeheader()
        writer.writerows(health)
    print(f'Regions discovered: {len(portable)}')
    print(f'GeoPackages discovered: {len(paths)}')
    print(f'Metadata OK: {sum(r["status"] == "METADATA_OK" for r in health)}')
    print(f'Metadata issues: {sum(r["status"] != "METADATA_OK" for r in health)}')
    print(f'Portable index: {args.index_dir / "tdx_regions.csv"}')
    print(f'Local health report: {args.index_dir / "tdx_local_health.csv"}')
    print('NOTE: Metadata bounds may not equal exact regional footprints.')


if __name__ == '__main__':
    main()
