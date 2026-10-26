#!/usr/bin/env python3
"""Step 07: read-only validation and visualization of Step 06 TDX exports.

Usage:
  python 1_delineation/tdx/07_validate_and_visualize.py --control 0_control_files/control_USGS_10343500.toml
"""
import argparse
import json
import math
import tomllib
from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from shapely.geometry import Point

VARIANTS = ('tdx_native', 'tdx_clipped')


def args_parse():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--control', required=True, type=Path)
    p.add_argument('--dpi', type=int, default=300)
    p.add_argument('--skip-maps', action='store_true', help='Batch mode: validate geometry without rendering figures')
    return p.parse_args()


def union(gdf):
    return gdf.geometry.union_all()


def make_valid(gdf, name):
    if gdf.empty:
        raise ValueError(f'{name} is empty')
    return bool(gdf.geometry.notna().all() and (~gdf.geometry.is_empty).all() and gdf.geometry.is_valid.all())


def load_layer(base, site, variant, layer):
    path = base / variant / 'shapefiles' / layer / f'{site}_{layer}.shp'
    if not path.is_file():
        raise FileNotFoundError(path)
    return gpd.read_file(path)


def projected_area(gdf, crs):
    return union(gdf.to_crs(crs)).area / 1e6


def validate(site, variant, frames, proj_crs, obs_area, outlet_id, gauge_point):
    hru = frames['catchment']; river = frames['river_network']; basin = frames['river_basins']
    req_hru = {'COMID', 'GRU_ID', 'HRU_ID', 'HRU_area', 'hru_to_seg', 'center_lat', 'center_lon'}
    req_river = {'COMID', 'NextDownID', 'slope', 'lengthkm'}
    for layer, required, df in [('catchment', req_hru, hru), ('river_network', req_river, river)]:
        missing = sorted(required - set(df.columns))
        if missing:
            raise ValueError(f'{variant} {layer} missing columns: {missing}')
    basin_ids = [int(v) for v in basin['COMID']] if 'COMID' in basin else []
    hru_ids = [int(v) for v in hru['HRU_ID']]
    river_ids = [int(v) for v in river['COMID']]
    id_set = set(river_ids)
    next_down = {int(r.COMID): int(r.NextDownID) for r in river.itertuples()}
    outlets = sorted(k for k, v in next_down.items() if v in (0, -1))
    invalid_down = sorted((k, v) for k, v in next_down.items() if v not in (0, -1) and v not in id_set)
    cycles = []
    for start in id_set:
        seen = set(); cur = start
        while cur in next_down and next_down[cur] not in (0, -1):
            if cur in seen:
                cycles.append(start); break
            seen.add(cur); cur = next_down[cur]
    hp = hru.to_crs(proj_crs)
    rp = river.to_crs(proj_crs)
    bp = basin.to_crs(proj_crs)
    gauge = gpd.GeoSeries([gauge_point], crs='EPSG:4326').to_crs(proj_crs).iloc[0]
    sum_hru_area = float(hp.geometry.area.sum() / 1e6)
    union_hru_area = float(union(hp).area / 1e6)
    basin_area = float(union(bp).area / 1e6)
    overlap_area = max(0., sum_hru_area - union_hru_area)
    missing_basin_area = float(union(bp).difference(union(hp)).area / 1e6)
    outside_basin_area = float(union(hp).difference(union(bp)).area / 1e6)
    hru_by_id = {int(row.HRU_ID): row.geometry for row in hp.itertuples()}
    river_in_hru = {}
    for row in rp.itertuples():
        rid = int(row.COMID)
        geom = hru_by_id.get(rid)
        river_in_hru[str(rid)] = round(float(row.geometry.intersection(geom).length / row.geometry.length), 6) if geom is not None and row.geometry.length else None
    slopes = pd.to_numeric(river['slope'], errors='coerce')
    lengths = pd.to_numeric(river['lengthkm'], errors='coerce')
    area_values = pd.to_numeric(hru['HRU_area'], errors='coerce')
    area_error = np.abs(area_values.to_numpy(dtype=float) - hp.geometry.area.to_numpy(dtype=float))
    gauge_to_outlet = None
    outlet_rows = rp[rp['COMID'].astype(int) == outlet_id]
    if len(outlet_rows):
        gauge_to_outlet = float(outlet_rows.geometry.iloc[0].distance(gauge))
    checks = {
        'valid_HRU_geometries': make_valid(hru, 'HRU'),
        'valid_river_geometries': make_valid(river, 'river'),
        'valid_basin_geometries': make_valid(basin, 'basin'),
        'unique_HRU_IDs': len(set(hru_ids)) == len(hru_ids),
        'river_basins_have_per_reach_COMIDs': len(basin_ids)==len(river_ids) and set(basin_ids)==set(river_ids),
        'unique_river_basin_COMIDs':len(basin_ids)==len(set(basin_ids)),
        'unique_river_COMIDs': len(id_set) == len(river_ids),
        'HRU_to_reach_mapping_valid': set(hru['hru_to_seg'].astype(int)).issubset(id_set),
        'HRU_IDs_equal_river_COMIDs': set(hru_ids) == id_set,
        'exactly_one_outlet': len(outlets) == 1,
        'expected_outlet_reach': outlets == [outlet_id],
        'valid_downstream_references': not invalid_down,
        'acyclic_network': not cycles,
        'positive_finite_slopes': bool(np.isfinite(slopes).all() and (slopes > 0).all()),
        'positive_finite_lengths': bool(np.isfinite(lengths).all() and (lengths > 0).all()),
        'HRU_area_matches_geometry_within_1pct': bool(np.isfinite(area_error).all() and np.all(area_error <= np.maximum(1, hp.geometry.area.to_numpy() * .01))),
        'HRU_overlaps_below_0_01_km2': overlap_area < .01,
        'basin_coverage_difference_below_0_01_km2': missing_basin_area < .01 and outside_basin_area < .01,
        'outlet_reach_passes_gauge_within_50m': gauge_to_outlet is not None and gauge_to_outlet <= 50.,
    }
    return {
        'site_id': site, 'variant': variant, 'HRU_count': len(hru), 'reach_count': len(river),
        'HRU_IDs': sorted(hru_ids), 'river_COMIDs': sorted(river_ids), 'outlet_reach_ids': outlets,
        'HRU_area_sum_km2': round(sum_hru_area, 6), 'HRU_union_area_km2': round(union_hru_area, 6),
        'basin_area_km2': round(basin_area, 6), 'observed_area_km2': obs_area,
        'basin_area_difference_from_observed_pct': round((basin_area / obs_area - 1) * 100, 3) if obs_area > 0 else None,
        'HRU_overlap_area_km2': round(overlap_area, 7), 'basin_not_covered_by_HRUs_km2': round(missing_basin_area, 7),
        'HRUs_outside_basin_km2': round(outside_basin_area, 7), 'gauge_distance_to_outlet_reach_m': round(gauge_to_outlet, 3) if gauge_to_outlet is not None else None,
        'min_river_length_fraction_inside_own_HRU': min((v for v in river_in_hru.values() if v is not None), default=None),
        'river_length_fraction_inside_own_HRU': river_in_hru,
        'river_slopes_by_COMID': {str(int(k)): float(v) if pd.notna(v) else None for k, v in zip(river['COMID'], slopes)},
        'invalid_downstream_references': invalid_down, 'cycle_start_reaches': cycles,
        'checks': checks, 'all_checks_pass': all(checks.values())
    }


def annotate(ax, hru, river, proj_crs):
    hp = hru.to_crs(proj_crs); rp = river.to_crs(proj_crs)
    for row in hp.itertuples():
        p = row.geometry.representative_point()
        ax.annotate(f'HRU {int(row.HRU_ID)}', (p.x, p.y), xytext=(0, 12), textcoords='offset points',
                    fontsize=9, fontweight='bold', color='#244b36', ha='center', va='bottom',
                    bbox=dict(boxstyle='round,pad=.23', fc='white', ec='#68a37d', alpha=.93), zorder=8)
    for row in rp.itertuples():
        p = row.geometry.interpolate(.5, normalized=True)
        ax.annotate(f'R {int(row.COMID)}', (p.x, p.y), xytext=(0, -13), textcoords='offset points',
                    fontsize=8, fontweight='bold', color='#842d26', ha='center', va='top',
                    bbox=dict(boxstyle='round,pad=.2', fc='white', ec='#cf8a80', alpha=.95), zorder=9)


def style_axis(ax, layers, gauge, proj_crs, title, native_boundary=None, d8_boundary=None, label=True):
    hp = layers['catchment'].to_crs(proj_crs)
    rp = layers['river_network'].to_crs(proj_crs)
    hp.plot(ax=ax, color='#cce5d0', edgecolor='#568b63', linewidth=1.1, alpha=.65, zorder=2)
    rp.plot(ax=ax, color='#b73329', linewidth=2.1, zorder=4)
    if native_boundary is not None:
        native_boundary.to_crs(proj_crs).boundary.plot(ax=ax, color='#777777', linestyle='--', linewidth=1.5, zorder=3)
    if d8_boundary is not None:
        d8_boundary.to_crs(proj_crs).boundary.plot(ax=ax, color='#2675bc', linestyle=':', linewidth=1.7, zorder=5)
    gp = gpd.GeoSeries([gauge], crs='EPSG:4326').to_crs(proj_crs).iloc[0]
    ax.scatter([gp.x], [gp.y], marker='x', s=100, color='black', linewidths=2.2, zorder=10)
    if label:
        annotate(ax, layers['catchment'], layers['river_network'], proj_crs)
    ax.set_title(title, fontsize=13, weight='bold')
    ax.set_aspect('equal', adjustable='datalim')
    ax.set_xlabel('Easting (m)'); ax.set_ylabel('Northing (m)')
    ax.ticklabel_format(axis='both', style='plain', useOffset=False)
    ax.tick_params(labelsize=8)
    ax.grid(alpha=.15)
    ax.margins(.12)


def legend(ax, include_comparison=False):
    items = [Patch(facecolor='#cce5d0', edgecolor='#568b63', label='HRUs'),
             Line2D([0], [0], color='#b73329', lw=2, label='TDX river reaches'),
             Line2D([0], [0], color='black', marker='x', linestyle='None', markersize=8, label='Outlet gauge')]
    if include_comparison:
        items.extend([Line2D([0], [0], color='#777777', ls='--', label='Native boundary'),
                      Line2D([0], [0], color='#2675bc', ls=':', label='Copernicus D8 boundary')])
    ax.legend(handles=items, loc='best', fontsize=8, framealpha=.95)


def save_variant_map(path, site, variant, layers, gauge, proj_crs, native_boundary, d8_boundary, dpi):
    fig, ax = plt.subplots(figsize=(11, 9))
    style_axis(ax, layers, gauge, proj_crs, f'{site} | {variant}',
               native_boundary=native_boundary if variant == 'tdx_clipped' else None,
               d8_boundary=d8_boundary, label=True)
    legend(ax, include_comparison=(variant == 'tdx_clipped'))
    fig.tight_layout(); fig.savefig(path, dpi=dpi, bbox_inches='tight'); plt.close(fig)


def save_comparison(path, site, frames, gauge, proj_crs, native_boundary, d8_boundary, dpi):
    fig, axes = plt.subplots(1, 2, figsize=(19, 9))
    for ax, variant in zip(axes, VARIANTS):
        style_axis(ax, frames[variant], gauge, proj_crs, variant,
                   native_boundary=native_boundary if variant == 'tdx_clipped' else None,
                   d8_boundary=d8_boundary, label=True)
        legend(ax, include_comparison=variant == 'tdx_clipped')
    fig.suptitle(f'{site} | TDX native vs clipped (HRU and river IDs)', fontsize=15, weight='bold')
    fig.tight_layout(); fig.savefig(path, dpi=dpi, bbox_inches='tight'); plt.close(fig)


def main():
    args = args_parse()
    with args.control.open('rb') as f:
        config = tomllib.load(f)
    site = config['domain']['site_id']
    proj_crs = config.get('delineation', {}).get('projected_crs', 'EPSG:26910')
    obs_area = float(config['domain'].get('observed_area_km2', 0))
    output_root = Path(config['output']['output_directory']) / site
    export_root = output_root / 'cwarhm_exports'
    validation_root = output_root / 'validation'
    summary_path = export_root / 'step06_summary.json'
    if not summary_path.exists():
        raise FileNotFoundError(f'Step 06 summary not found: {summary_path}')
    summary06 = json.loads(summary_path.read_text())
    outlet_id = int(summary06['outlet_reach_id'])
    gauge = Point(float(summary06['outlet_gauge_lon']), float(summary06['outlet_gauge_lat']))
    frames = {variant: {layer: load_layer(export_root, site, variant, layer)
                        for layer in ('catchment', 'river_basins', 'river_network')} for variant in VARIANTS}
    native_boundary = frames['tdx_native']['river_basins']
    refined_path = output_root / 'tdx_upstream_refined.gpkg'
    d8_boundary = None
    if refined_path.is_file():
        import fiona
        if 'refinement_boundary' in fiona.listlayers(refined_path):
            d8_boundary = gpd.read_file(refined_path, layer='refinement_boundary')
    results = {}
    for variant in VARIANTS:
        target = validation_root / variant
        target.mkdir(parents=True, exist_ok=True)
        results[variant] = validate(site, variant, frames[variant], proj_crs, obs_area, outlet_id, gauge)
        (target / 'validation_summary.json').write_text(json.dumps(results[variant], indent=2) + '\n')
        if not args.skip_maps:
            save_variant_map(target / 'validation_map.jpg', site, variant, frames[variant], gauge,
                             proj_crs, native_boundary, d8_boundary, args.dpi)
    structural_pass = all(results[v]['all_checks_pass'] for v in VARIANTS)
    preflight = summary06.get('preflight', {})
    hydro_qa = preflight.get('status', 'UNKNOWN')
    comparison = {
                  'site_id':site,
                  'structural_validation_pass':structural_pass,
                  'hydrological_QA_status':hydro_qa,
                  'production_export_authorized':False,
                  'overall_status':('BLOCKED' if not structural_pass else
                                    'REVIEW' if hydro_qa in ('REVIEW','BLOCKED','UNKNOWN') else
                                    'PROVISIONAL'),
                  'review_flags':preflight.get('review_flags', []),
                  'variant_results': results,
                  'same_HRU_IDs': results['tdx_native']['HRU_IDs'] == results['tdx_clipped']['HRU_IDs'],
                  'same_river_COMIDs': results['tdx_native']['river_COMIDs'] == results['tdx_clipped']['river_COMIDs'],
                  'native_area_km2': results['tdx_native']['basin_area_km2'],
                  'clipped_area_km2': results['tdx_clipped']['basin_area_km2']}
    if d8_boundary is not None:
        d8p = d8_boundary.to_crs(proj_crs)
        cp = frames['tdx_clipped']['catchment'].to_crs(proj_crs)
        d8_geom = union(d8p); clipped_geom = union(cp)
        comparison['Copernicus_D8_area_km2'] = round(d8_geom.area / 1e6, 6)
        comparison['D8_coverage_by_clipped_HRUs_pct'] = round(100 * clipped_geom.intersection(d8_geom).area / d8_geom.area, 4)
        comparison['D8_uncovered_area_km2'] = round(d8_geom.difference(clipped_geom).area / 1e6, 6)
    if not args.skip_maps:
        save_comparison(validation_root / 'native_vs_clipped_comparison.jpg', site, frames, gauge,
                        proj_crs, native_boundary, d8_boundary, args.dpi)
    (validation_root / 'step07_validation_summary.json').write_text(json.dumps(comparison, indent=2) + '\n')
    print(json.dumps({'site_id': site, 'outputs': str(validation_root),
                      'checks_pass': {v: results[v]['all_checks_pass'] for v in VARIANTS},
                      'overall_status':comparison['overall_status'],
                      'production_export_authorized':False,
                      'same_HRU_IDs': comparison['same_HRU_IDs'],
                      'same_river_COMIDs': comparison['same_river_COMIDs'],
                      'D8_coverage_by_clipped_HRUs_pct': comparison.get('D8_coverage_by_clipped_HRUs_pct'),
                      'figures': ([] if args.skip_maps else
                                  [str(validation_root / v / 'validation_map.jpg') for v in VARIANTS] +
                                  [str(validation_root / 'native_vs_clipped_comparison.jpg')])}, indent=2))


if __name__ == '__main__':
    main()
