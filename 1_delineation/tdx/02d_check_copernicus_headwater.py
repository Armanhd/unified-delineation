#!/usr/bin/env python3
"""Diagnose Copernicus GLO-30 D8 headwaters within an original TDX HRU.

Read-only with respect to TDX / original DEM and Steps 01-07.
Given a gauge with completed 02b and 02c diagnostics, find plausible D8
outlet pixels near the OBSERVED gauge coordinates and compare their upstream
catchments with the identified original TDX HRU and reference drainage area.

This does NOT yet change/export a river reach. Retain original TDX IDs.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import tempfile
from contextlib import ExitStack
from pathlib import Path

import fiona
import geopandas as gpd
import numpy as np
import rasterio
from rasterio.features import shapes
from rasterio.merge import merge
from rasterio.fill import fillnodata
from scipy.ndimage import distance_transform_edt
from rasterio.warp import transform as transform_coords
from rasterio.warp import transform_bounds
from shapely.geometry import Point, shape
from shapely.ops import unary_union

PROJECT = Path(__file__).resolve().parents[2]
REGIONS = PROJECT / 'indexes/tdx_regions.csv'
CATALOG = PROJECT / 'indexes/gauge_catalogs/camels-spat_gauge_catalog_validated.csv'
DIAG_DIR = PROJECT / 'indexes/gauge_catalogs/tdx_snapping/headwater_diagnostics'
TDX_ROOT = Path('/work/comphyd_lab/data/geospatial-data/TDX_hydrofabrics/TDX_hydrofabrics')
DEM_ROOT = Path('/work/comphyd_lab/data/geospatial-data/copernicus/glo-30')
DIRMAP = (64, 128, 1, 2, 4, 8, 16, 32)  # pysheds convention


def csv_row(file: Path, field: str, key: str):
    with file.open(newline='', encoding='utf-8-sig') as f:
        records = [r for r in csv.DictReader(f) if str(r.get(field)) == str(key)]
    if len(records) != 1:
        raise ValueError(f'Expected one {field}={key} in {file}; found {len(records)}')
    return records[0]


def normalize_id(v):
    s = str(v or '').strip()
    return s[:-2] if s.endswith('.0') else s


def copernicus_tiles(root: Path, wgs_bounds):
    west, south, east, north = map(float, wgs_bounds)
    if not (west >= -180 and east <= 180 and south >= -90 and north <= 90 and west < east):
        raise ValueError('Bounds not supported (antimeridian/polar handling required)')
    tiles = []
    for lat in range(math.floor(south), math.ceil(north)):
        for lon in range(math.floor(west), math.ceil(east)):
            ns = f'N{lat:02d}' if lat >= 0 else f'S{-lat:02d}'
            ew = f'E{lon:03d}' if lon >= 0 else f'W{-lon:03d}'
            name = f'Copernicus_DSM_COG_10_{ns}_00_{ew}_00_DEM'
            tile = root / name / f'{name}.tif'
            if not tile.is_file():
                raise FileNotFoundError(f'Copernicus GLO-30 tile missing: {tile}')
            tiles.append(tile)
    return tiles


def original_hru(region, parent_id, gauge, metric, radius_m, tdx_root):
    path = tdx_root / region['basins_filename']
    if not path.is_file():
        raise FileNotFoundError(path)
    layers = [x['name'] for x in json.loads(region['basins_layers_json']) if x.get('name')]
    search = gpd.GeoSeries([gauge], crs='EPSG:4326').to_crs(metric).buffer(radius_m)
    matches = []
    for layer in layers:
        with fiona.open(path, layer=layer) as src:
            native = src.crs_wkt or src.crs
            if not native:
                raise RuntimeError(f'No CRS: {path}:{layer}')
            bounds = tuple(gpd.GeoSeries(search, crs=metric).to_crs(native).total_bounds)
            for feat in src.filter(bbox=bounds):
                if not feat['geometry']:
                    continue
                if normalize_id(feat['properties'].get('streamID')) != parent_id:
                    continue
                geom = shape(feat['geometry'])
                if not geom.is_valid:
                    geom = geom.buffer(0)
                if geom.is_empty:
                    continue
                candidate = gpd.GeoSeries([geom], crs=native).to_crs('EPSG:4326').iloc[0]
                if candidate.covers(gauge):
                    matches.append(candidate)
    if len(matches) != 1:
        raise RuntimeError(f'Expected exactly one matching TDX parent HRU; found {len(matches)}')
    return matches[0]


def cell_area_m2(affine, crs, lat):
    if crs.is_geographic:
        return (abs(affine.a) * 111320 * math.cos(math.radians(lat))
                * abs(affine.e) * 110574)
    return abs(affine.a * affine.e)


def write_json(path: Path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--site', required=True)
    ap.add_argument('--catalog', type=Path, default=CATALOG)
    ap.add_argument('--region-index', type=Path, default=REGIONS)
    ap.add_argument('--diagnostics-dir', type=Path, default=DIAG_DIR)
    ap.add_argument('--tdx-root', type=Path, default=TDX_ROOT)
    ap.add_argument('--dem-directory', type=Path, default=DEM_ROOT)
    ap.add_argument('--buffer-m', type=float, default=2500,
                    help='Metric buffer around complete parent HRU for DEM conditioning')
    ap.add_argument('--snap-radius-m', type=float, default=120,
                    help='Max distance of candidate D8 outlet cell from observed gauge')
    ap.add_argument('--max-candidate-watersheds', type=int, default=5)
    ap.add_argument('--max-dem-cells', type=int, default=15000000)
    ap.add_argument('--max-void-percent', type=float, default=1.0,
                    help='Fail if missing DEM cells exceed this fraction of the window, in percent')
    ap.add_argument('--max-void-depth-cells', type=float, default=5.0,
                    help='Fail if a missing pixel is farther than this many pixels from valid DEM')
    args = ap.parse_args()
    if args.buffer_m <= 0 or args.snap_radius_m <= 0 or args.max_candidate_watersheds < 1:
        ap.error('Buffer and radius must be >0; candidate count must be >=1')
    if not (0 < args.max_void_percent <= 100) or args.max_void_depth_cells < 1:
        ap.error('Void percent must be in (0, 100]; void depth must be >= 1 pixel')
    try:
        from pysheds.grid import Grid
    except ImportError as exc:
        raise RuntimeError('pysheds is required (use unified_delineation conda environment)') from exc
    if not hasattr(np, 'in1d'):
        np.in1d = np.isin  # compatibility with older pysheds

    site = args.site
    catalog = csv_row(args.catalog, 'site_id', site)
    diagnostic = args.diagnostics_dir / f'{site}_tdx_hru_diagnostic.json'
    parent_report = json.loads(diagnostic.read_text())
    if parent_report['status'] != 'PARENT_HRU_IDENTIFIED':
        raise ValueError(f'{site}: 02c parent HRU is not uniquely identified')
    parent_id = str(parent_report['provisional_parent_tdx_id'])
    region_id = str(parent_report['tdx_region'])
    region = csv_row(args.region_index, 'tdx_id', region_id)
    gauge = Point(float(catalog['longitude']), float(catalog['latitude']))
    expected = float(catalog['expected_area_km2'])
    metric = gpd.GeoSeries([gauge], crs='EPSG:4326').estimate_utm_crs()
    if metric is None:
        raise RuntimeError('Failed to estimate local UTM CRS')
    parent = original_hru(region, parent_id, gauge, metric, args.buffer_m, args.tdx_root)
    parent_metric = gpd.GeoSeries([parent], crs=4326).to_crs(metric).iloc[0]
    parent_area = float(parent_metric.area / 1e6)

    # Local working DEM covers parent polygon and a user-controlled buffer.
    bounds_metric = parent_metric.buffer(args.buffer_m)
    bounds_wgs = gpd.GeoSeries([bounds_metric], crs=metric).to_crs(4326).total_bounds
    if args.dem_directory.is_file():
        tiles = [args.dem_directory]
    else:
        tiles = copernicus_tiles(args.dem_directory, bounds_wgs)
    print(f'{site}: original TDX HRU={parent_id}, parent area={parent_area:.4f} km2')
    print(f'Copernicus tiles: {len(tiles)}', flush=True)

    with ExitStack() as stack:
        srcs = [stack.enter_context(rasterio.open(t)) for t in tiles]
        dem_crs = srcs[0].crs
        if dem_crs is None or any(ds.crs != dem_crs for ds in srcs):
            raise ValueError('Selected DEM tiles lack a common CRS')
        dem_bounds = transform_bounds('EPSG:4326', dem_crs, *bounds_wgs, densify_pts=21)
        mosaic, affine = merge(srcs, bounds=dem_bounds, indexes=1, masked=True)
        raster = np.ma.asarray(mosaic[0])
        if raster.size > args.max_dem_cells:
            raise RuntimeError(f'DEM window {raster.size:,} cells exceeds --max-dem-cells')
        invalid = np.ma.getmaskarray(raster) | ~np.isfinite(np.ma.getdata(raster))
        raw = np.asarray(np.ma.getdata(raster), dtype='float32').copy()

        # Record the original coverage gaps BEFORE limited interpolation.
        # Original DEM tiles are never modified. Large/deep gaps cause a hard failure.
        void_count = int(invalid.sum())
        void_percent = 100.0 * void_count / invalid.size
        max_void_depth = 0.0
        void_raster_path = None
        if void_count:
            args.diagnostics_dir.mkdir(parents=True, exist_ok=True)
            void_raster_path = args.diagnostics_dir / f'{site}_copernicus_dem_voids.tif'
            with rasterio.open(
                void_raster_path, 'w', driver='GTiff', height=raw.shape[0],
                width=raw.shape[1], count=1, dtype='uint8', crs=dem_crs,
                transform=affine, nodata=255, compress='LZW',
            ) as dst:
                dst.write(invalid.astype('uint8'), 1)

            if void_count == invalid.size:
                raise RuntimeError('The selected DEM window has no valid elevation cells')
            distances = distance_transform_edt(invalid)
            max_void_depth = float(np.max(distances[invalid]))
            print(f'DEM missing pixels: {void_count:,}/{invalid.size:,} '
                  f'({void_percent:.4f}%), maximum gap depth={max_void_depth:.2f} pixels',
                  flush=True)
            if (void_percent > args.max_void_percent or
                    max_void_depth > args.max_void_depth_cells):
                raise RuntimeError(
                    f'DEM gaps exceed safe local-interpolation limits; '
                    f'void={void_percent:.4f}% (max {args.max_void_percent}%), '
                    f'depth={max_void_depth:.2f} pixels '
                    f'(max {args.max_void_depth_cells}). '
                    f'Inspect void mask: {void_raster_path}'
                )

            # GDAL fillnodata interpolates locally from valid neighbours.
            # Do not pretend interpolated pixels are original elevations:
            # original invalid mask is retained for watershed QA below.
            raw = fillnodata(
                np.where(invalid, 0.0, raw).astype('float32'),
                mask=(~invalid).astype('uint8'),
                max_search_distance=max(1.0, math.ceil(args.max_void_depth_cells)),
                smoothing_iterations=0,
            )
            if not np.isfinite(raw).all():
                raise RuntimeError(f'Interpolated DEM still contains invalid values: {void_raster_path}')
            print('Only small local DEM gaps interpolated in working copy; '
                  'original gaps retained for diagnostics.', flush=True)
        else:
            print('DEM missing pixels: 0', flush=True)

    gx, gy = transform_coords('EPSG:4326', dem_crs, [gauge.x], [gauge.y])
    gx, gy = float(gx[0]), float(gy[0])
    if dem_crs.is_geographic:
        meters_per_x = 111320.0 * math.cos(math.radians(gauge.y))
        meters_per_y = 110574.0
    else:
        # For projected non-metric DEM CRSs this would need explicit unit conversion.
        meters_per_x = meters_per_y = 1.0
    gr, gc = rasterio.transform.rowcol(affine, gx, gy)
    if not (0 <= gr < raw.shape[0] and 0 <= gc < raw.shape[1]):
        raise ValueError('Gauge outside DEM working window')
    if invalid[gr, gc]:
        raise RuntimeError(
            'Observed gauge falls on an originally missing Copernicus DEM cell. '
            'Do not rely on interpolated elevations at the gauge; review the DEM void '
            f'mask: {void_raster_path}'
        )

    print(f'Conditioning DEM window {raw.shape[0]} x {raw.shape[1]} ...', flush=True)
    with tempfile.TemporaryDirectory(prefix=f'tdx_headwater_{site}_') as temp:
        tif = Path(temp) / 'local_dem.tif'
        profile = dict(driver='GTiff', count=1, height=raw.shape[0], width=raw.shape[1],
                       dtype='float32', crs=dem_crs, transform=affine)
        with rasterio.open(tif, 'w', **profile) as ds:
            ds.write(raw, 1)
        grid = Grid.from_raster(str(tif))
        dem = grid.read_raster(str(tif))
        dem_conditioned = grid.resolve_flats(grid.fill_depressions(grid.fill_pits(dem)))
        direction = grid.flowdir(dem_conditioned, dirmap=DIRMAP)
        accum = grid.accumulation(direction, dirmap=DIRMAP)
        print('D8 direction and accumulation computed.', flush=True)

        max_dr = math.ceil(args.snap_radius_m / abs(affine.e * meters_per_y)) + 1
        max_dc = math.ceil(args.snap_radius_m / abs(affine.a * meters_per_x)) + 1
        area_per_cell = cell_area_m2(affine, dem_crs, gauge.y)
        candidates = []
        for row in range(max(0, gr-max_dr), min(raw.shape[0], gr+max_dr+1)):
            for col in range(max(0, gc-max_dc), min(raw.shape[1], gc+max_dc+1)):
                x, y = rasterio.transform.xy(affine, row, col, offset='center')
                dist = math.hypot((x-gx)*meters_per_x, (y-gy)*meters_per_y)
                value = float(accum[row, col])
                if dist > args.snap_radius_m or not math.isfinite(value) or value < 1:
                    continue
                approx_area = value * area_per_cell / 1e6
                # Candidate ranking combines proximity and reference-area agreement.
                area_mismatch = abs(math.log(max(approx_area, 1e-8) / expected))
                score = area_mismatch + 0.35 * (dist / args.snap_radius_m)
                candidates.append((score, dist, row, col, x, y, approx_area, value))
        if not candidates:
            raise RuntimeError('No valid Copernicus D8 outlet cells close to observed gauge')
        candidates.sort(key=lambda t: (t[0], t[1], t[2], t[3]))
        shortlisted = candidates[:args.max_candidate_watersheds]
        output_rows = []
        geometries = []
        for rank, (score, dist, row, col, x, y, accum_area, accum_cells) in enumerate(shortlisted, 1):
            catch = grid.catchment(x=int(col), y=int(row), fdir=direction,
                                   dirmap=DIRMAP, xytype='index')
            mask = np.asarray(catch, dtype=bool)
            mask[row, col] = True
            touches_edge = bool(mask[0].any() or mask[-1].any() or mask[:, 0].any() or mask[:, -1].any())
            polygons = [shape(geom) for geom, val in shapes(mask.astype('uint8'),
                        mask=mask, transform=affine) if val == 1]
            poly = unary_union(polygons)
            if not poly.is_valid:
                poly = poly.buffer(0)
            poly_geo = gpd.GeoSeries([poly], crs=dem_crs).to_crs(4326).iloc[0]
            poly_metric = gpd.GeoSeries([poly_geo], crs=4326).to_crs(metric).iloc[0]
            area = float(poly_metric.area / 1e6)
            overlap = float(poly_metric.intersection(parent_metric).area / 1e6)
            outside = max(0., area - overlap)
            output_rows.append({
                'site_id': site, 'parent_tdx_id': parent_id, 'rank': rank,
                'snap_distance_m': round(dist, 3),
                'snapped_lon': round(float(gpd.GeoSeries([Point(x, y)], crs=dem_crs).to_crs(4326).iloc[0].x), 9),
                'snapped_lat': round(float(gpd.GeoSeries([Point(x, y)], crs=dem_crs).to_crs(4326).iloc[0].y), 9),
                'd8_accumulation_cells': int(accum_cells),
                'approx_accumulation_area_km2': round(accum_area, 6),
                'd8_watershed_area_km2': round(area, 6),
                'reference_area_km2': round(expected, 6),
                'area_error_pct': round(100 * (area-expected)/expected, 3),
                'parent_hru_area_km2': round(parent_area, 6),
                'inside_parent_hru_km2': round(overlap, 6),
                'outside_parent_hru_km2': round(outside, 6),
                'percent_d8_inside_parent': round(100 * overlap/area, 3) if area > 0 else 0.,
                'touches_dem_window_edge': touches_edge,
                'original_dem_void_cells_in_watershed': int(np.count_nonzero(mask & invalid)),
                'ranking_score': round(score, 6),
            })
            geometries.append(poly_geo)
            print(f'  D8 candidate {rank}: distance={dist:.1f} m, area={area:.3f} km2, '
                  f'within original TDX HRU={100*overlap/area:.1f}%, edge={touches_edge}, '
                  f'filled DEM cells in catchment={int(np.count_nonzero(mask & invalid))}', flush=True)

    args.diagnostics_dir.mkdir(parents=True, exist_ok=True)
    out_csv = args.diagnostics_dir / f'{site}_copernicus_d8_candidates.csv'
    out_json = args.diagnostics_dir / f'{site}_copernicus_d8_summary.json'
    out_gpkg = args.diagnostics_dir / f'{site}_copernicus_d8_diagnostic.gpkg'
    import pandas as pd
    pd.DataFrame(output_rows).to_csv(out_csv, index=False)
    if out_gpkg.exists():
        out_gpkg.unlink()
    gpd.GeoDataFrame({'rank': list(range(1, len(geometries)+1)), 'site_id': site,
                      'parent_id': parent_id}, geometry=geometries,
                     crs=4326).to_file(out_gpkg, layer='candidate_d8_watersheds', driver='GPKG')
    gpd.GeoDataFrame({'site_id': [site], 'parent_id': [parent_id]},
                     geometry=[parent], crs=4326).to_file(out_gpkg,
                     layer='original_tdx_parent_hru', driver='GPKG', mode='a')
    gpd.GeoDataFrame({'site_id': [site]}, geometry=[gauge], crs=4326).to_file(
        out_gpkg, layer='observed_gauge', driver='GPKG', mode='a')
    best = output_rows[0]
    result = {
        'site_id': site, 'tdx_region': region_id, 'original_tdx_id': parent_id,
        'reference_area_km2': expected, 'original_tdx_hru_area_km2': parent_area,
        'dem_tiles': [str(t) for t in tiles], 'dem_cell_count': int(raw.size),
        'original_dem_void_cells': void_count,
        'original_dem_void_percent': round(void_percent, 6),
        'max_original_dem_void_depth_cells': round(max_void_depth, 3),
        'dem_void_cells_interpolated_in_working_copy': void_count,
        'original_dem_void_mask_tif': str(void_raster_path) if void_raster_path else None,
        'snap_radius_m': args.snap_radius_m, 'dem_buffer_m': args.buffer_m,
        'candidate_count_within_radius': len(candidates),
        'd8_candidates_evaluated': len(output_rows), 'best_candidate': best,
        'status': 'DIAGNOSTIC_ONLY',
        'notes': [
            'No original TDX reach or subbasin geometry or identifiers were modified.',
            'D8 catchment and gauge-to-cell choice must be reviewed before reach export.',
            'If a candidate catchment contains interpolated DEM cells, review their positions '
            'and their influence on D8 topology before accepting the delineation.',
            'If gauge lies upstream of TDX channel start, DEM-derived gauge-upstream reach '
            'will replace downstream original TDX geometry in this basin export, keeping original TDX ID.',
            'Fractional IDs are reserved for actual nested-gauge splits, not headwater extension.'
        ],
        'files': {'candidates_csv': str(out_csv), 'diagnostic_gpkg': str(out_gpkg)}
    }
    write_json(out_json, result)
    print(json.dumps({'site': site, 'tdx_id': parent_id, 'best_candidate': best,
                      'csv': str(out_csv), 'gpkg': str(out_gpkg), 'json': str(out_json)}, indent=2))


if __name__ == '__main__':
    main()
