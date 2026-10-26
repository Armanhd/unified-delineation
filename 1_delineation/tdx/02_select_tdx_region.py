#!/usr/bin/env python3
"""Select a TDX-Hydro region for a gauge, using a portable regional index.

Run standalone with --control or --lat/--lon. Import select_tdx_region() from
future per-basin workflows. Never writes to the source TDX dataset.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import tomllib
from pathlib import Path
from typing import Any

import fiona
from pyproj import CRS, Transformer
from shapely.geometry import Point, shape
from shapely.ops import transform


def _layers(raw: str | None) -> list[dict[str, Any]]:
    if not raw or str(raw).lower() == 'nan':
        return []
    value = json.loads(raw)
    return value if isinstance(value, list) else []


def _preferred_layer(layers: list[dict[str, Any]], kind: str) -> dict[str, Any] | None:
    preferred = 'polygon' if kind == 'basins' else 'line'
    for layer in layers:
        if preferred in str(layer.get('geometry_type', '')).lower():
            return layer
    return layers[0] if layers else None


def _contains_bbox(row: dict[str, str], lon: float, lat: float) -> bool:
    try:
        xmin, ymin, xmax, ymax = (float(row[k]) for k in ('min_x','min_y','max_x','max_y'))
        srs = str(row.get('bbox_srs_id', '')).strip()
        if srs and srs.lower() != 'nan' and int(float(srs)) != 4326:
            p = Transformer.from_crs('EPSG:4326', CRS.from_epsg(int(float(srs))), always_xy=True).transform(Point(lon, lat))
            x, y = p.x, p.y
        else:
            if not srs or srs.lower() == 'nan':
                return True  # unknown CRS: do not reject prematurely
            x, y = lon, lat
        return xmin <= x <= xmax and ymin <= y <= ymax
    except (ValueError, KeyError, TypeError):
        return True  # unknown bounds: inspect actual geometry


def _gauge_in_basin(path: Path, layer: str, lon: float, lat: float, search_deg: float = 0.001) -> tuple[bool, int, str]:
    """Read only features in a small spatial filter; confirm polygon covers gauge."""
    with fiona.open(path, layer=layer) as src:
        if not src.crs and not src.crs_wkt:
            raise ValueError('Basin layer has no coordinate reference system')
        crs = CRS.from_user_input(src.crs_wkt or src.crs)
        convert = Transformer.from_crs('EPSG:4326', crs, always_xy=True).transform
        pt = transform(convert, Point(lon, lat))
        # Query a small bounding box in the native CRS; never read all features.
        delta = search_deg
        corners = [transform(convert, Point(lon + dx, lat + dy)) for dx in (-delta, delta) for dy in (-delta, delta)]
        bbox = (min(p.x for p in corners), min(p.y for p in corners),
                max(p.x for p in corners), max(p.y for p in corners))
        tested = 0
        for feature in src.filter(bbox=bbox):
            if not feature['geometry']:
                continue
            tested += 1
            geom = shape(feature['geometry'])
            if geom.covers(pt):
                return True, tested, crs.to_string()
        return False, tested, crs.to_string()


def select_tdx_region(*, lat: float, lon: float, tdx_dir: str | Path,
                      index_csv: str | Path, diagnostics_dir: str | Path | None = None,
                      site_id: str = 'site') -> dict[str, Any]:
    """Return a confirmed region or raise an actionable error.

    Ambiguous multiple confirmed regions are never silently resolved.
    """
    if not (-90 <= lat <= 90 and -180 <= lon <= 180):
        raise ValueError('Invalid gauge latitude/longitude')
    tdx_dir, index_csv = Path(tdx_dir).expanduser(), Path(index_csv).expanduser()
    if not index_csv.is_file():
        raise FileNotFoundError(f'Regional index not found: {index_csv}')
    if not tdx_dir.is_dir():
        raise NotADirectoryError(f'TDX source directory not found: {tdx_dir}')
    with index_csv.open(newline='') as fh:
        index_rows = list(csv.DictReader(fh))
    candidates = [r for r in index_rows if _contains_bbox(r, lon, lat)]
    records: list[dict[str, Any]] = []
    confirmed: list[dict[str, Any]] = []
    for row in candidates:
        region = row.get('tdx_id', '')
        basin_name = row.get('basins_filename', '')
        stream_name = row.get('streamnet_filename', '')
        basin_path, stream_path = tdx_dir / basin_name, tdx_dir / stream_name
        record: dict[str, Any] = {'tdx_id': region, 'status': '', 'details': '',
                                  'basins_filename': basin_name, 'streamnet_filename': stream_name}
        try:
            basin_layers = _layers(row.get('basins_layers_json'))
            stream_layers = _layers(row.get('streamnet_layers_json'))
            basin_layer = _preferred_layer(basin_layers, 'basins')
            stream_layer = _preferred_layer(stream_layers, 'streamnet')
            if not basin_path.is_file() or not stream_path.is_file():
                raise FileNotFoundError('Missing basin or stream GeoPackage')
            if basin_layer is None or stream_layer is None:
                raise ValueError('Missing layer metadata in regional index')
            covers, tested, crs = _gauge_in_basin(basin_path, basin_layer['name'], lon, lat)
            record['features_tested'] = tested
            record['basins_crs'] = crs
            if covers:
                # Check that indexed stream layer really exists without scanning geometries.
                with fiona.open(stream_path, layer=stream_layer['name']):
                    pass
                record['status'] = 'CONFIRMED'
                confirmed.append({'tdx_id': region,
                                  'basins_file': str(basin_path),
                                  'streamnet_file': str(stream_path),
                                  'basins_layer': basin_layer['name'],
                                  'streamnet_layer': stream_layer['name'],
                                  'lat': lat, 'lon': lon})
            else:
                record['status'] = 'NOT_CONTAINED'
        except Exception as exc:
            record['status'] = 'ERROR'
            record['details'] = f'{type(exc).__name__}: {exc}'
        records.append(record)
    if diagnostics_dir is not None:
        directory = Path(diagnostics_dir)
        directory.mkdir(parents=True, exist_ok=True)
        fields = ['tdx_id','status','details','basins_filename','streamnet_filename','features_tested','basins_crs']
        with (directory / f'{site_id}_tdx_region_candidates.csv').open('w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(records)
    if len(confirmed) == 1:
        result = confirmed[0]
        if diagnostics_dir is not None:
            (Path(diagnostics_dir) / f'{site_id}_tdx_region_selection.json').write_text(json.dumps(result, indent=2) + '\n')
        return result
    errors = [r for r in records if r['status'] == 'ERROR']
    if len(confirmed) > 1:
        raise RuntimeError(f'Multiple TDX regions contain gauge {site_id}: '
                           f'{[r["tdx_id"] for r in confirmed]}. Review candidate diagnostics.')
    if errors:
        raise RuntimeError(f'No confirmed TDX region for {site_id}; {len(errors)} candidate file(s) could not be read. '
                           'See tdx_region_candidates.csv; do not assume the gauge is outside TDX coverage.')
    raise RuntimeError(f'No TDX basin polygon covers gauge {site_id} ({lat}, {lon}); '
                       f'{len(candidates)} regional bounding-box candidates tested.')


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--control', type=Path, help='Basin TOML control file')
    parser.add_argument('--lat', type=float)
    parser.add_argument('--lon', type=float)
    parser.add_argument('--site-id')
    parser.add_argument('--tdx-dir', type=Path)
    parser.add_argument('--index', type=Path, default=Path('indexes/tdx_regions.csv'))
    parser.add_argument('--diagnostics-dir', type=Path)
    args = parser.parse_args()
    config: dict[str, Any] = {}
    if args.control:
        with args.control.open('rb') as fh:
            config = tomllib.load(fh)
    domain = config.get('domain', {})
    inp = config.get('input', {})
    out = config.get('output', {})
    lat = args.lat if args.lat is not None else domain.get('lat')
    lon = args.lon if args.lon is not None else domain.get('lon')
    site_id = args.site_id or domain.get('site_id', 'site')
    tdx_dir = args.tdx_dir or inp.get('tdx_directory')
    if lat is None or lon is None or not tdx_dir:
        parser.error('Provide --lat, --lon and --tdx-dir, or a TOML with [domain] lat/lon and [input] tdx_directory')
    if config.get('domain', {}).get('delineation_method', 'TDX').upper() != 'TDX':
        parser.error('This selector currently supports only delineation_method = "TDX"')
    diagnostics = args.diagnostics_dir
    if diagnostics is None and out.get('output_directory'):
        diagnostics = Path(out['output_directory']) / site_id / out.get('diagnostics_subdirectory', 'diagnostics') / 'step02'
    try:
        result = select_tdx_region(lat=float(lat), lon=float(lon), site_id=site_id,
                                   tdx_dir=tdx_dir, index_csv=args.index,
                                   diagnostics_dir=diagnostics)
    except Exception as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 1
    print(json.dumps(result, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
