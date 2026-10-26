#!/usr/bin/env python3
"""Step 06: export tdx_native and tdx_clipped CWARHM-style shapefiles.

Original TDX reach IDs are never renumbered. Both exports terminate the outlet
river at the gauge; native polygon boundaries are retained in the native export.
Clipping provenance is written to metadata, not encoded in identifiers.
"""
from __future__ import annotations

import argparse
import json
import sys
import csv
import tomllib
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point, LineString, MultiLineString
from shapely.ops import substring, linemerge

# Share one TDX slope policy across ordinary and nested exporters.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'common'))
from tdx_slope_resolver import TDXSlopeResolver


def normalize_id(value):
    if value is None or pd.isna(value):
        raise ValueError('Missing TDX reach ID')
    value = str(value).strip()
    if value.endswith('.0') and value[:-2].lstrip('-').isdigit():
        value = value[:-2]
    if not value.lstrip('-').isdigit():
        raise ValueError(f'TDX reach ID is not integer-like: {value!r}')
    return int(value)


def pick_column(df, candidates, required=False):
    columns = {c.lower(): c for c in df.columns}
    for candidate in candidates:
        if candidate.lower() in columns:
            return columns[candidate.lower()]
    if required:
        raise KeyError(f'Expected one of {candidates}; actual columns: {list(df.columns)}')
    return None


def linestring(geom):
    if isinstance(geom, LineString):
        return geom
    if isinstance(geom, MultiLineString):
        merged = linemerge(geom)
        if isinstance(merged, LineString):
            return merged
    raise ValueError('Outlet reach is not a single connected LineString; inspect geometry')


def outlet_upstream_piece(line, gauge, orientation):
    """Return portion from gauge toward upstream end, with native orientation."""
    line = linestring(line)
    distance = line.project(gauge)
    if not (0 < distance < line.length):
        raise ValueError('Gauge projection lies at endpoint of outlet reach; cannot clip safely')
    if orientation == 'downstream_to_upstream':
        piece = substring(line, distance, line.length)
    elif orientation == 'upstream_to_downstream':
        piece = substring(line, 0, distance)
    else:
        raise ValueError('Invalid outlet reach orientation')
    if not isinstance(piece, LineString) or piece.length <= 0:
        raise RuntimeError('Clipped outlet reach is empty')
    return piece, distance


def projected_area_km2(df, metric):
    return float(df.to_crs(metric).geometry.union_all().area / 1e6)


def write_shape(gdf, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f'Refusing to overwrite existing shapefile: {path}. Use --overwrite.')
    gdf.to_file(path, driver='ESRI Shapefile', index=False, encoding='UTF-8')


def remove_shape(path):
    for ext in ('.shp', '.shx', '.dbf', '.prj', '.cpg', '.qix', '.fix'):
        candidate = path.with_suffix(ext)
        if candidate.exists():
            candidate.unlink()


def retrieve_original_slopes(cfg, reach_ids):
    """Query just selected LINKNO rows from original TDX stream GeoPackage.

    Never load the entire regional TDX stream network into memory.
    """
    inputs = cfg.get('input', {})
    tdx_dir = Path(inputs.get('tdx_directory', ''))
    region = str(inputs.get('tdx_id', '')).strip()
    filename = inputs.get('tdx_streams_file') or f'{region}-streamnet.gpkg'
    source = tdx_dir / filename
    if not source.exists():
        raise FileNotFoundError(f'Original TDX stream GeoPackage not found: {source}')
    try:
        import fiona
        layers = fiona.listlayers(source)
    except Exception as e:
        raise RuntimeError(f'Cannot list TDX stream layers in {source}: {e}') from e
    configured = inputs.get('tdx_streams_layer', '')
    layer = configured if configured else next((x for x in layers if 'streamnet' in x.lower()), None)
    if not layer:
        raise ValueError(f'Cannot identify stream layer; available: {layers}')
    field = cfg.get('tdx_fields', {}).get('reach_id', 'LINKNO')
    slope_field = cfg.get('tdx_fields', {}).get('slope', 'Slope')
    ids = sorted(int(x) for x in reach_ids)
    where = f'"{field}" IN ({",".join(map(str, ids))})'
    try:
        subset = gpd.read_file(source, layer=layer, where=where, ignore_geometry=True)
    except Exception as e:
        raise RuntimeError(f'Filtered TDX lookup failed (will not read full stream layer): {e}') from e
    rid_col = pick_column(subset, [field, 'LINKNO'], required=True)
    sl_col = pick_column(subset, [slope_field, 'Slope', 'slope'], required=True)
    subset['_rid'] = subset[rid_col].map(normalize_id)
    if subset['_rid'].duplicated().any():
        raise ValueError('Duplicate LINKNO in original TDX slope query')
    missing = set(ids) - set(subset['_rid'])
    if missing:
        raise ValueError(f'Original TDX reaches missing for slope recovery: {sorted(missing)}')
    result = {}
    for _, r in subset.iterrows():
        val = pd.to_numeric(r[sl_col], errors='coerce')
        result[int(r['_rid'])] = float(val) if pd.notna(val) and np.isfinite(val) and val > 0 else None
    return result, {'source_file':str(source), 'source_layer':layer,
                    'original_id_field':rid_col, 'original_slope_field':sl_col,
                    'missing_or_nonpositive_ids':[k for k,v in result.items() if v is None]}


def spatial_diagnostics(hrus, reaches, outlet_id, metric, variant, tolerance_m=50.0):
    """Quantify river/HRU correspondence and physical junction gaps; do not silently snap."""
    by_hru = {int(row.COMID): row.geometry for _, row in hrus.iterrows()}
    by_river = {int(row.COMID): row for _, row in reaches.iterrows()}
    records = []
    for rid, reach in sorted(by_river.items()):
        line = reach.geometry
        poly = by_hru[rid]
        covered = float(line.intersection(poly.buffer(1.0)).length)
        ratio = min(1.0, covered / line.length)
        next_id = int(reach.NextDownID)
        if next_id:
            # TDX line orientation is downstream -> upstream; check both endpoints
            # to avoid assuming upstream source geometry has identical orientation.
            down_line = by_river[next_id].geometry
            junction_gap = min(Point(line.coords[0]).distance(down_line),
                               Point(line.coords[-1]).distance(down_line))
        else:
            junction_gap = None
        records.append({'variant':variant, 'COMID':rid, 'NextDownID':next_id,
                        'river_length_m':round(float(line.length),3),
                        'length_inside_own_HRU_m':round(covered,3),
                        'own_HRU_length_fraction':round(ratio,6),
                        'junction_gap_m':round(float(junction_gap),3) if junction_gap is not None else None,
                        'junction_within_tolerance':bool(junction_gap <= tolerance_m) if junction_gap is not None else None})
    return records



def read_json_optional(path):
    if path.is_file():
        return json.loads(path.read_text())
    return {}


def planned_split_parents(cfg, region):
    """Read only existing 04b plan; missing plan means not independently verified."""
    root = Path(__file__).resolve().parents[2]
    path = (root / 'indexes/gauge_catalogs/tdx_shared_network/nested_split_plans'
            / f'{region}_nested_split_parents.csv')
    if not path.is_file():
        return None, path
    with path.open(newline='') as handle:
        rows = list(csv.DictReader(handle))
    return {str(r['parent_tdx_reach_id']) for r in rows}, path


def qa_preflight(cfg, site, root, diag, outlet, native, refined, rivers):
    step05 = read_json_optional(diag / 'tdx_refinement_summary.json')
    step04 = read_json_optional(diag / 'step04' / f'{site}_tdx_shared_gauge_topology.json')
    region = str(outlet.get('tdx_id') or cfg.get('input',{}).get('tdx_id', ''))
    parents, plan_path = planned_split_parents(cfg, region)
    river_ids = {str(int(x)) for x in rivers['_rid']}
    overlapping = sorted(river_ids & parents) if parents is not None else []
    all_native = set(native['_rid']); all_refined = set(refined['_rid']); all_reaches = set(rivers['_rid'])
    blocked = []
    review = []
    if not step05:
        blocked.append('STEP05_SUMMARY_MISSING')
    if not step04 or not step04.get('shared_relationships_analyzed'):
        blocked.append('STEP04_SHARED_ANALYSIS_MISSING')
    elif step04.get('related_gauge_count', 0):
        blocked.append('NESTED_GAUGES_REQUIRE_SHARED_HRU_PARTITION')
    if parents is None:
        blocked.append('REGIONAL_04B_SPLIT_PLAN_MISSING')
    if overlapping:
        blocked.append('PLANNED_SPLIT_PARENTS_PRESENT')
    if all_native != all_reaches:
        blocked.append('NATIVE_HRU_REACH_IDS_DIFFER')
    if all_refined != all_reaches:
        blocked.append('REFINED_HRU_REACH_IDS_DIFFER')
    if not step05.get('qa_status') == 'PROVISIONAL':
        review.extend(step05.get('qa_flags') or ['STEP05_NOT_APPROVED'])
    if outlet.get('review_required'):
        review.append('STEP03_REVIEW_REQUIRED')
    if cfg.get('nested_gauges', {}).get('enabled') or cfg.get('headwater', {}).get('enabled'):
        blocked.append('FEATURE_FLAGS_ENABLED_BUT_NOT_INTEGRATED')
    status = 'BLOCKED' if blocked else 'REVIEW' if review else 'PROVISIONAL'
    result = {
        'site_id': site, 'status': status, 'preview_only': True,
        'production_export_authorized': False,
        'blocked_reasons': sorted(set(blocked)), 'review_flags': sorted(set(review)),
        'step05_qa_status': step05.get('qa_status'),
        'step05_area_difference_pct': step05.get('refined_area_difference_pct'),
        'shared_related_gauge_count': step04.get('related_gauge_count'),
        'split_plan_file': str(plan_path),
        'split_parents_in_upstream_network': overlapping,
        'native_HRU_count': len(native), 'refined_HRU_count': len(refined),
        'native_reach_count': len(rivers),
    }
    (diag / 'step06_preflight.json').write_text(json.dumps(result, indent=2)+'\n')
    return result

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--control', type=Path, required=True)
    ap.add_argument('--junction-tolerance-m', type=float, default=50.0)
    ap.add_argument('--overwrite', action='store_true', help='Replace Step 06 preview files only')
    ap.add_argument('--allow-review', action='store_true', help='Write REVIEW previews (never production-ready); does not override blocking QA')
    ap.add_argument('--preflight-only', action='store_true', help='Check dependencies and report readiness without writing shapefiles')
    args = ap.parse_args()
    with args.control.open('rb') as f:
        cfg = tomllib.load(f)
    site = cfg['domain']['site_id']
    root = Path(cfg['output']['output_directory']).expanduser() / site
    native_path = root / 'tdx_upstream_native.gpkg'
    refined_path = root / 'tdx_upstream_refined.gpkg'
    diag = root / cfg['output'].get('diagnostics_subdirectory', 'diagnostics')
    outlet_selection = json.loads((diag / 'tdx_outlet_selection.json').read_text())
    outlet_id = normalize_id(outlet_selection['reach_id'])
    gauge = Point(float(outlet_selection['snapped_lon']), float(outlet_selection['snapped_lat']))
    settings = cfg.get('delineation', {})
    metric = settings.get('projected_crs', 'EPSG:26910')
    export_crs = settings.get('output_crs', 'EPSG:4326')
    orientation = outlet_selection.get('line_orientation')
    if orientation not in ('downstream_to_upstream', 'upstream_to_downstream'):
        raise ValueError('Step 03 must verify and record outlet line_orientation before Step 06')
    clip_river = bool(settings.get('clip_outlet_reach', True))
    if not clip_river:
        raise ValueError('Step 06 requires clip_outlet_reach=true to terminate river at gauge')

    # Read inputs before creating any output. Step 04/05 GeoPackages remain immutable.
    native = gpd.read_file(native_path, layer='native_subbasins')
    refined = gpd.read_file(refined_path, layer='refined_subbasins')
    rivers = gpd.read_file(native_path, layer='upstream_reaches')
    if any(x.empty for x in (native, refined, rivers)):
        raise ValueError('Input polygon or reach layer is empty')
    if any(x.crs is None for x in (native, refined, rivers)):
        raise ValueError('All source layers require a CRS')

    poly_id = pick_column(native, ['reach_id'], required=True)
    refined_id = pick_column(refined, ['reach_id'], required=True)
    river_id = pick_column(rivers, ['reach_id'], required=True)
    native['_rid'] = native[poly_id].map(normalize_id)
    refined['_rid'] = refined[refined_id].map(normalize_id)
    rivers['_rid'] = rivers[river_id].map(normalize_id)
    for label, df in [('native HRUs',native),('refined HRUs',refined),('river reaches',rivers)]:
        if df['_rid'].duplicated().any():
            raise ValueError(f'Duplicate original reach IDs in {label}')
    native_ids = set(native['_rid']); refined_ids = set(refined['_rid']); river_ids = set(rivers['_rid'])
    if outlet_id not in native_ids or outlet_id not in refined_ids:
        raise ValueError('Outlet ID absent from one of the HRU datasets')
    # Preflight is written even for unresolved or REVIEW basins.
    preflight = qa_preflight(cfg, site, root, diag, outlet_selection, native, refined, rivers)
    print('Step 06 preflight:', json.dumps(preflight, indent=2), flush=True)
    if args.preflight_only:
        return
    if preflight['status'] == 'BLOCKED':
        raise RuntimeError('Basin blocked pending shared splits/HRU partition or missing dependencies; see step06_preflight.json')
    if preflight['status'] == 'REVIEW' and not args.allow_review:
        raise RuntimeError('Basin has hydrological QA flags. Use --allow-review for *preview-only* exports, or leave for later.')
    if refined_ids != native_ids:
        raise ValueError('Refined HRU IDs differ from native IDs. Resolve before CWARHM export.')

    # Use metric CRS for lengths/areas and gauge projection.
    rivers_m = rivers.to_crs(metric).copy()
    native_m = native.to_crs(metric).copy()
    refined_m = refined.to_crs(metric).copy()
    gauge_m = gpd.GeoSeries([gauge], crs='EPSG:4326').to_crs(metric).iloc[0]
    outlet_idx = rivers_m.index[rivers_m['_rid'] == outlet_id][0]
    full_line = linestring(rivers_m.at[outlet_idx, 'geometry'])
    clipped_line, along = outlet_upstream_piece(full_line, gauge_m, orientation)
    gauge_offset_m = float(full_line.distance(gauge_m))
    if gauge_offset_m > float(cfg.get('outlet', {}).get('search_radius_m', 500)):
        raise ValueError(f'Gauge is {gauge_offset_m:.1f} m from selected reach')
    rivers_m.at[outlet_idx, 'geometry'] = clipped_line

    # Build reach attributes, retaining original numeric TDX IDs.
    downstream_field = pick_column(rivers, ['downstream_id', cfg.get('tdx_fields', {}).get('downstream_id', ''), 'DSLINKNO', 'NextDownID', 'nextdownid'])
    if downstream_field is None:
        raise KeyError(f'Cannot determine downstream connectivity; available fields: {list(rivers.columns)}')
    slope_field = pick_column(rivers, [cfg.get('tdx_fields', {}).get('slope', ''), 'Slope', 'slope'])
    slope_resolver = TDXSlopeResolver(cfg)
    try:
        recovered_slopes, slope_substitutions, slope_provenance = slope_resolver.resolve(river_ids)
    finally:
        slope_resolver.close()
    # Any source-slope substitution is a documented preview-only QA condition.
    slope_substitutions_path = diag / 'step06_slope_substitutions.csv'
    if slope_substitutions:
        pd.DataFrame(slope_substitutions).to_csv(slope_substitutions_path, index=False)
    reach_records = []
    for _, r in rivers_m.iterrows():
        rid = int(r['_rid'])
        if rid == outlet_id:
            next_id = 0
        elif downstream_field:
            raw = r[downstream_field]
            next_id = normalize_id(raw) if pd.notna(raw) else 0
        else:
            raise KeyError('Missing downstream connectivity field')
        if rid != outlet_id and next_id not in river_ids:
            raise ValueError(f'Reach {rid} points to {next_id}, not in selected network')
        geom = linestring(r.geometry)
        length_m = float(geom.length)
        if length_m <= 0:
            raise ValueError(f'Reach {rid} has zero length')
        # The original parent TDX slope is authoritative; split/shortened geometry
        # never changes the slope. Missing original slopes inherit from a valid
        # connected original TDX upstream reach, then downstream as last resort.
        slope = recovered_slopes[rid]
        reach_records.append({'COMID':rid,'NextDownID':next_id,'slope':slope,
                              'lengthkm':length_m/1000,'length_m':length_m,
                              'geometry':geom})
    reaches = gpd.GeoDataFrame(reach_records, crs=metric)
    if sum(reaches['NextDownID'].eq(0)) != 1:
        raise ValueError('Expected exactly one terminal outlet reach')
    reach_lookup = dict(zip(reaches.COMID, reaches.NextDownID))
    for rid in reach_lookup:
        seen = set(); current = rid
        while current != 0:
            if current in seen:
                raise ValueError(f'River network contains cycle starting at {rid}')
            seen.add(current)
            if current not in reach_lookup:
                raise ValueError(f'Reach {rid} routes outside exported network')
            current = reach_lookup[current]
    # NO FRACTIONAL IDS: this preview path forbids all globally planned split reaches.
    # No reach is renumbered; HRU_ID = native TDX reach ID (one HRU per native subbasin).
    # GRU_ID also uses the native reach ID. If a future CWARHM scheme needs separate GRUs,
    # change the mapping in an explicit, versioned export, not by silently renumbering.
    def make_hrus(polys):
        rows = []
        for _, r in polys.iterrows():
            rid = int(r['_rid'])
            geom = r.geometry
            if geom is None or geom.is_empty or geom.area <= 0:
                raise ValueError(f'Invalid or empty HRU polygon {rid}')
            center = gpd.GeoSeries([geom.representative_point()], crs=metric).to_crs('EPSG:4326').iloc[0]
            rows.append({'COMID':rid, 'GRU_ID':rid, 'HRU_ID':rid,
                         'HRU_area':float(geom.area), 'area':float(geom.area)/1e6,
                         'center_lat':float(center.y), 'center_lon':float(center.x),
                         'hru_to_seg':rid, 'geometry':geom})
        result = gpd.GeoDataFrame(rows, crs=metric)
        if result.geometry.union_all().area <= 0:
            raise ValueError('Invalid HRU union')
        sum_area = float(result.geometry.area.sum())
        union_area = float(result.geometry.union_all().area)
        if sum_area - union_area > max(1.0, union_area * 1e-7):
            raise ValueError('HRU polygons overlap materially')
        return result

    export_root = root / 'cwarhm_exports'
    native_area_by_id = {int(r['_rid']): float(r.geometry.area) for _, r in native_m.iterrows()}
    native_river_lengths = {int(r['_rid']): float(linestring(r.geometry).length) for _, r in rivers.to_crs(metric).iterrows()}
    summary = {'site_id':site, 'outlet_reach_id':outlet_id,
               'preflight':preflight,'preview_only':True,
               'production_export_authorized':False,
               'outlet_gauge_lon':gauge.x,'outlet_gauge_lat':gauge.y,
               'outlet_river_clipped':True,'original_outlet_length_m':round(full_line.length,3),
               'exported_outlet_length_m':round(clipped_line.length,3),
               'gauge_to_reach_distance_m':round(gauge_offset_m,3),
               'outlet_geometry_orientation':orientation,
               'id_policy':'All HRU_ID, GRU_ID, COMID, hru_to_seg retain original TDX reach IDs; no fractional/new IDs',
               'slope_provenance':slope_provenance,
               'slope_substitutions':slope_substitutions,
               'slope_substitutions_csv':str(slope_substitutions_path) if slope_substitutions else None,
               'exports':{}}
    for variant, polygons in [('tdx_native',native_m),('tdx_clipped',refined_m)]:
        hrus = make_hrus(polygons)
        # CWARHM river_basins MUST contain one polygon per reach, NOT a
        # dissolved single watershed. We retain the exact HRU/river mapping.
        river_basins = hrus[['COMID', 'geometry']].copy()
        river_basins['area_km2'] = hrus['area'].values
        watershed_area = float(hrus.geometry.union_all().area/1e6)
        folder = export_root / variant / 'shapefiles'
        files = {
            'catchment':folder/'catchment'/f'{site}_catchment.shp',
            'river_basins':folder/'river_basins'/f'{site}_river_basins.shp',
            'river_network':folder/'river_network'/f'{site}_river_network.shp',
        }
        if not args.overwrite:
            existing = [str(p) for p in files.values() if p.exists()]
            if existing:
                raise FileExistsError('Exports already exist. Use --overwrite: '+', '.join(existing))
        else:
            for p in files.values():
                remove_shape(p)
        # catchment = HRUs; river_basins = one polygon per mapped river reach.
        # All shapefiles use configured output CRS.
        write_shape(hrus.to_crs(export_crs), files['catchment'])
        write_shape(river_basins.to_crs(export_crs), files['river_basins'])
        write_shape(reaches.to_crs(export_crs), files['river_network'])
        metadata_dir = export_root / variant / 'metadata'
        metadata_dir.mkdir(parents=True, exist_ok=True)
        spatial_rows = spatial_diagnostics(hrus, reaches, outlet_id, metric, variant, args.junction_tolerance_m)
        pd.DataFrame(spatial_rows).to_csv(metadata_dir/'river_hru_spatial_diagnostics.csv',index=False)
        rows = []
        for _, h in hrus.iterrows():
            rid = int(h.COMID)
            area0 = native_area_by_id[rid]
            area1 = float(h.HRU_area)
            len0 = native_river_lengths[rid]
            len1 = float(reaches.loc[reaches.COMID==rid,'length_m'].iloc[0])
            rows.append({'original_TDX_ID':rid,'is_outlet':rid==outlet_id,
                         'HRU_polygon_clipped':bool(area1 < area0 - 1.0),
                         'reach_clipped_at_gauge':rid==outlet_id,
                         'native_HRU_area_m2':round(area0,3),
                         'exported_HRU_area_m2':round(area1,3),
                         'native_reach_length_m':round(len0,3),
                         'exported_reach_length_m':round(len1,3),
                         'area_removed_m2':round(max(0,area0-area1),3),
                         'length_removed_m':round(max(0,len0-len1),3),
                         'polygon_source':'native TDX' if variant=='tdx_native' else 'TDX clipped to Copernicus D8 watershed'})
        pd.DataFrame(rows).sort_values('original_TDX_ID').to_csv(metadata_dir/'clipping_provenance.csv',index=False)
        meta = {'variant':variant,'site_id':site,'original_TDX_ID_preserved':True,
                'outlet_reach_id':outlet_id,'outlet_reach_clipped':True,
                'native_polygon_preserved':variant=='tdx_native',
                'refined_polygon_source':'Copernicus GLO-30 conditioned DEM D8' if variant=='tdx_clipped' else None,
                'HRU_count':len(hrus),'reach_count':len(reaches),
                'watershed_area_km2':round(watershed_area,6),
                'river_basins_are_per_reach':True,
                'preview_only':True,'production_export_authorized':False,
                'qa_status':preflight['status'],
                'review_flags':sorted(set(preflight['review_flags'] +
                     (['SUBSTITUTED_SOURCE_TDX_SLOPES_REVIEW'] if slope_substitutions else []))),
                'files':{k:str(v) for k,v in files.items()},
                'spatial_diagnostics_file':str(metadata_dir/'river_hru_spatial_diagnostics.csv'),
                'min_river_length_fraction_in_own_HRU':min(x['own_HRU_length_fraction'] for x in spatial_rows),
                'max_upstream_junction_gap_m':max((x['junction_gap_m'] or 0) for x in spatial_rows),
                'notes':[('Native HRU boundaries are preserved.' if variant=='tdx_native' else 'HRU boundaries are clipped to Copernicus D8.'),
                 'The outlet river is truncated at the gauge in both variants.',
                 'PREVIEW ONLY: no basin has final production export authorization at this stage.',
                 'Original source GeoPackages are not modified.',
                 'Clipping flags and original sizes are in clipping_provenance.csv.']}
        (metadata_dir/'export_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
        summary['exports'][variant] = meta
    (export_root/'step06_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    try:
        main()
    except Exception as exc:
        print(f'ERROR: {type(exc).__name__}: {exc}',file=sys.stderr)
        sys.exit(1)
