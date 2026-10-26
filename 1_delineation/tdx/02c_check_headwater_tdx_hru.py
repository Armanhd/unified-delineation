#!/usr/bin/env python3
"""Step 02c: identify the original TDX subbasin covering a headwater gauge.

Read-only diagnostic. Does not invent IDs, alter TDX features, or generate a
DEM-derived channel. Identifies the correct native parent before DEM work.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import fiona
import geopandas as gpd
from pyproj import CRS, Transformer
from shapely.geometry import Point, shape
from shapely.ops import transform

PROJECT = Path(__file__).resolve().parents[2]
TDX_ROOT = Path('/work/comphyd_lab/data/geospatial-data/TDX_hydrofabrics/TDX_hydrofabrics')
CATALOG = PROJECT / 'indexes/gauge_catalogs/camels-spat_gauge_catalog_validated.csv'
REGION_INDEX = PROJECT / 'indexes/tdx_regions.csv'
SNAP_DIR = PROJECT / 'indexes/gauge_catalogs/tdx_snapping'


def row_for(path: Path, field: str, value: str):
    with path.open(newline='', encoding='utf-8-sig') as f:
        items = [r for r in csv.DictReader(f) if r.get(field) == value]
    if len(items) != 1:
        raise ValueError(f'Expected one {field}={value} in {path}; found {len(items)}')
    return items[0]


def normalize_id(value):
    text = str(value or '').strip()
    return text[:-2] if text.endswith('.0') else text


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--site', required=True)
    p.add_argument('--catalog', type=Path, default=CATALOG)
    p.add_argument('--region-index', type=Path, default=REGION_INDEX)
    p.add_argument('--snap-dir', type=Path, default=SNAP_DIR)
    p.add_argument('--tdx-root', type=Path, default=TDX_ROOT)
    p.add_argument('--radius-m', type=float, default=2000.0)
    p.add_argument('--output-dir', type=Path, default=SNAP_DIR / 'headwater_diagnostics')
    args = p.parse_args()
    if args.radius_m <= 0:
        p.error('--radius-m must be > 0')

    catalog = row_for(args.catalog, 'site_id', args.site)
    snap_file = args.snap_dir / f'tdx_gauge_snaps_{args.site}.csv'
    snap = row_for(snap_file, 'site_id', args.site)
    if not snap['tdx_id'] or not snap['reach_id']:
        raise ValueError(f'{args.site}: no selected TDX region/reach in {snap_file}')
    region = row_for(args.region_index, 'tdx_id', snap['tdx_id'])
    tdx_file = args.tdx_root / region['basins_filename']
    if not tdx_file.is_file():
        raise FileNotFoundError(tdx_file)

    lat, lon = float(catalog['latitude']), float(catalog['longitude'])
    gauge = Point(lon, lat)
    metric = gpd.GeoSeries([gauge], crs='EPSG:4326').estimate_utm_crs()
    if metric is None:
        raise RuntimeError('Cannot estimate local metric CRS')
    to_metric = Transformer.from_crs(4326, metric, always_xy=True).transform
    gauge_metric = transform(to_metric, gauge)
    search_metric = gauge_metric.buffer(args.radius_m).envelope
    layers = [x['name'] for x in json.loads(region['basins_layers_json']) if x.get('name')]
    if not layers:
        raise ValueError(f'No TDX basin layer in region {snap["tdx_id"]}')

    candidates = []
    for layer in layers:
        with fiona.open(tdx_file, layer=layer) as src:
            if not src.crs and not src.crs_wkt:
                raise ValueError(f'Missing CRS in {tdx_file}:{layer}')
            native = CRS.from_user_input(src.crs_wkt or src.crs)
            bbox_transform = Transformer.from_crs(metric, native, always_xy=True).transform
            native_to_metric = Transformer.from_crs(native, metric, always_xy=True).transform
            native_to_geo = Transformer.from_crs(native, 4326, always_xy=True).transform
            bbox = transform(bbox_transform, search_metric).bounds
            if 'streamID' not in src.schema['properties']:
                raise KeyError(f'Missing streamID in {tdx_file}:{layer}')
            for feature in src.filter(bbox=bbox):
                if not feature['geometry']:
                    continue
                river_id = normalize_id(feature['properties'].get('streamID'))
                if not river_id:
                    continue
                geom = shape(feature['geometry'])
                if geom.is_empty:
                    continue
                poly_metric = transform(native_to_metric, geom)
                dist = float(poly_metric.distance(gauge_metric))
                if dist > args.radius_m:
                    continue
                contains = bool(poly_metric.covers(gauge_metric))
                poly_geo = transform(native_to_geo, geom)
                candidates.append({
                    'site_id': args.site,
                    'tdx_id': snap['tdx_id'],
                    'layer': layer,
                    'hru_reach_id': river_id,
                    'gauge_in_polygon': contains,
                    'distance_to_polygon_m': round(dist, 3),
                    'native_hru_area_km2': round(poly_metric.area / 1e6, 6),
                    'polygon_bounds_wgs84': ','.join(f'{v:.7f}' for v in poly_geo.bounds),
                })

    candidates.sort(key=lambda r: (not r['gauge_in_polygon'], r['distance_to_polygon_m'], r['hru_reach_id']))
    containing = sorted({r['hru_reach_id'] for r in candidates if r['gauge_in_polygon']})
    selected_river = normalize_id(snap['reach_id'])
    parent_id = containing[0] if len(containing) == 1 else None
    warnings = []
    if not containing:
        warnings.append('GAUGE_NOT_INSIDE_ANY_NEARBY_TDX_HRU')
    if len(containing) > 1:
        warnings.append('MULTIPLE_TDX_HRUS_COVER_GAUGE')
    if parent_id and parent_id != selected_river:
        warnings.append('GAUGE_HRU_ID_DIFFERS_FROM_NEAREST_TDX_REACH_ID')

    result = {
        'site_id': args.site,
        'gauge_lat': lat,
        'gauge_lon': lon,
        'expected_area_km2': float(catalog['expected_area_km2']),
        'tdx_region': snap['tdx_id'],
        'nearest_tdx_reach_id': selected_river,
        'nearest_reach_distance_m': float(snap['snap_distance_m']),
        'containing_hru_ids': containing,
        'provisional_parent_tdx_id': parent_id,
        'status': 'PARENT_HRU_IDENTIFIED' if parent_id else 'REVIEW',
        'warnings': warnings,
        'next_step': 'Verify Copernicus D8 catchment overlaps this TDX HRU before extending the reach',
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    json_file = args.output_dir / f'{args.site}_tdx_hru_diagnostic.json'
    csv_file = args.output_dir / f'{args.site}_tdx_hru_candidates.csv'
    json_file.write_text(json.dumps(result, indent=2) + '\n')
    with csv_file.open('w', newline='') as f:
        fields = ['site_id', 'tdx_id', 'layer', 'hru_reach_id', 'gauge_in_polygon',
                  'distance_to_polygon_m', 'native_hru_area_km2', 'polygon_bounds_wgs84']
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(candidates)
    print(json.dumps(result, indent=2))
    print(f'HRU candidates: {csv_file}')
    print(f'Diagnostic JSON: {json_file}')


if __name__ == '__main__':
    main()
