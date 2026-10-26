#!/usr/bin/env python3
"""Export Step 05b candidate nested HRUs/reaches to *preview-only* CWARHM shapefiles.

This is deliberately separate from the existing Step 06 native/clipped exporter.
It neither modifies original TDX GeoPackages nor authorizes model production.

Fractional IDs are stored as (a) exact text TDX_IDs and (b) numerical DBF
COMID/HRU_ID/NextDownID fields compatible with the one-decimal TDX convention.
A companion manifest retains exact string IDs and the lossless integer x10 code.

Usage:
 python 1_delineation/tdx/06b_export_nested_cwarhm_preview.py \\
   --control 0_control_files/tdx_camels_spat/control_CAN_02HC025.toml --preflight-only
 python 1_delineation/tdx/06b_export_nested_cwarhm_preview.py \\
   --control 0_control_files/tdx_camels_spat/control_CAN_02HC025.toml --allow-review
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import re
import statistics
import sqlite3
import sys
import tomllib
from decimal import Decimal
from pathlib import Path

import geopandas as gpd
from shapely.geometry import Point, LineString, MultiLineString, GeometryCollection
from shapely.ops import linemerge, substring, unary_union
from shapely.validation import explain_validity
from shapely import make_valid, set_precision

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'common'))
from tdx_slope_resolver import TDXSlopeResolver

ID_RE = re.compile(r'^(0|[1-9][0-9]*)(?:\.([1-9]))?$')


def id_token(value):
    """Accept text decimal IDs only if exactly native integer or .1..9."""
    text = str(value).strip()
    m = ID_RE.fullmatch(text)
    if m is None:
        raise ValueError(f'Unsupported TDX ID {text!r}; only native integer and .1...9 IDs are allowed')
    return text


def scale10(value):
    value = id_token(value)
    m = ID_RE.fullmatch(value)
    return int(m.group(1)) * 10 + (int(m.group(2)) if m.group(2) else 0)


def mixed_scaled_id(value):
    """Compatibility check if downstream CWARHM scales only fractional IDs."""
    value = id_token(value)
    return scale10(value) if '.' in value else int(value)


def numeric_id(value):
    return float(Decimal(id_token(value)))


def read_json(path):
    if not path.is_file():
        raise FileNotFoundError(f'Required summary missing: {path}')
    return json.loads(path.read_text())


def column(df, name):
    found = next((c for c in df.columns if c.lower() == name.lower()), None)
    if not found:
        raise KeyError(f'Missing {name} field; found {list(df.columns)}')
    return found


def only_line(geom):
    if isinstance(geom, LineString):
        return geom
    if isinstance(geom, MultiLineString):
        merged = linemerge(geom)
        if isinstance(merged, LineString):
            return merged
    raise ValueError('A river reach must be a single continuous LineString')


def polygonal_part(geom):
    """Keep only polygonal components of make_valid output, never buffer(0)."""
    if geom is None or geom.is_empty:
        return GeometryCollection()
    if geom.geom_type in ('Polygon', 'MultiPolygon'):
        return geom
    if geom.geom_type == 'GeometryCollection':
        parts = []
        for piece in geom.geoms:
            if piece.geom_type in ('Polygon', 'MultiPolygon') and not piece.is_empty:
                parts.append(piece)
        return unary_union(parts) if parts else GeometryCollection()
    return GeometryCollection()


def validate_and_repair_hrus(hru, id_column, diag, max_area_change_m2=100.0,
                             max_relative_area_change=0.001):
    """Report each invalid HRU and repair only very small, bounded errors.

    A zero-area or non-polygonal geometry cannot be repaired here. Repairs are
    explicitly flagged as REVIEW; no production authorization is conferred.
    """
    records = []
    repaired = []
    failures = []
    for idx, row in hru.iterrows():
        rid = str(row[id_column])
        geom = row.geometry
        kind = 'None' if geom is None else geom.geom_type
        original_area = float(geom.area) if geom is not None else 0.0
        reason = ('NULL_GEOMETRY' if geom is None else
                  'EMPTY_GEOMETRY' if geom.is_empty else
                  explain_validity(geom) if not geom.is_valid else
                  'VALID')
        status = 'PASS'
        area_change = 0.0
        result = geom
        if geom is None or geom.is_empty or original_area <= 0:
            status = 'BLOCKED_EMPTY_OR_ZERO_AREA'
        elif geom.geom_type not in ('Polygon', 'MultiPolygon'):
            status = 'BLOCKED_NON_POLYGON'
        elif not geom.is_valid:
            candidate = polygonal_part(make_valid(geom))
            if candidate.is_empty or not candidate.is_valid or candidate.area <= 0:
                status = 'BLOCKED_REPAIR_UNUSABLE'
            else:
                area_change = abs(float(candidate.area) - original_area)
                if (area_change > max_area_change_m2 or
                        area_change > max_relative_area_change * original_area):
                    status = 'BLOCKED_REPAIR_AREA_CHANGED'
                else:
                    result = candidate
                    status = 'REPAIRED_REVIEW'
        records.append({
            'reach_id': rid, 'original_type': kind,
            'original_validity': reason, 'original_area_m2': round(original_area, 7),
            'area_change_m2': round(area_change, 7), 'status': status,
        })
        if status.startswith('BLOCKED'):
            failures.append((rid, reason, status))
        if status == 'REPAIRED_REVIEW':
            repaired.append(rid)
            hru.at[idx, 'geometry'] = result
    report = diag / 'step06b_hru_geometry_qa.csv'
    report.parent.mkdir(parents=True, exist_ok=True)
    with report.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    if failures:
        detail = '; '.join(f'{rid}: {reason} [{status}]' for rid, reason, status in failures[:10])
        raise ValueError(f'Step05b HRU geometry check failed for {len(failures)} reach(es): '
                         f'{detail}. Inspect {report}. Do not use --allow-review to bypass this.')
    return repaired, report


def source_slopes(cfg, parents):
    """Read only relevant slopes via original GPKG attributes (not full geometries)."""
    inp = cfg['input']
    gpkg = Path(inp['tdx_directory']) / inp['tdx_streams_file']
    if not gpkg.is_file():
        raise FileNotFoundError(f'TDX slope source missing: {gpkg}')
    layer = inp.get('tdx_streams_layer')
    if not layer:
        import fiona
        layer = next((s for s in fiona.listlayers(gpkg) if 'streamnet' in s.lower()), None)
    if not layer:
        raise ValueError('TDX source stream layer not identified')
    field = cfg.get('tdx_fields', {}).get('reach_id', 'LINKNO')
    slope_field = cfg.get('tdx_fields', {}).get('slope', 'Slope')
    def quote(x):
        return '"' + str(x).replace('"','""') + '"'
    result = {}
    con = sqlite3.connect(gpkg.resolve().as_uri() + '?mode=ro', uri=True)
    try:
        all_columns = {r[1] for r in con.execute(f'PRAGMA table_info({quote(layer)})')}
        if field not in all_columns or slope_field not in all_columns:
            raise ValueError(f'TDX source missing {field} or {slope_field}')
        # SQLite parameter limit could be 999; bounded batches avoid it.
        ids = sorted({int(p) for p in parents})
        for i in range(0,len(ids),500):
            batch = ids[i:i+500]
            query = (f'SELECT {quote(field)}, {quote(slope_field)} FROM {quote(layer)} '
                     f'WHERE {quote(field)} IN ({",".join("?" for _ in batch)})')
            for orig, slope in con.execute(query,batch):
                rid = str(int(orig))
                val = float(slope) if slope is not None else math.nan
                result[rid] = val if math.isfinite(val) and val > 0 else None
    finally:
        con.close()
    missing = sorted(set(map(str, parents)) - set(result))
    if missing:
        raise ValueError(f'Cannot recover original TDX slopes for parents: {missing[:12]}')
    return result, {'source_file':str(gpkg),'source_layer':layer,
                    'missing_or_nonpositive_parent_slopes':[p for p,v in result.items() if v is None]}


def remove_shapefile(path):
    for ext in ['.shp','.shx','.dbf','.prj','.cpg','.qix','.fix']:
        p = path.with_suffix(ext)
        if p.exists():
            p.unlink()


def write_shapefile(frame,path):
    path.parent.mkdir(parents=True,exist_ok=True)
    frame.to_file(path,driver='ESRI Shapefile',index=False,encoding='UTF-8')


def fix_collapsed_polygon_components(geometry, crs):
    """Discard collapsed non-polygon pieces without discarding valid area.

    make_valid alone can leave invalid MultiPolygons after GPKG -> SHP
    serialization or projection. A tiny precision grid in the export CRS
    prevents zero-width polygon components from reappearing on read-back.
    """
    if geometry is None or geometry.is_empty:
        raise ValueError('Null or empty exported HRU')
    if geometry.geom_type not in ('Polygon','MultiPolygon'):
        raise ValueError(f'Exported HRU is not a polygon: {geometry.geom_type}')
    candidate = polygonal_part(make_valid(geometry))
    if candidate.is_empty:
        raise ValueError('No polygonal area remains after make_valid')
    # 1e-10 degrees is <= ~0.012 mm; for projected coordinates use 0.0001 m.
    is_geo = bool(crs.is_geographic)
    grid = 1e-10 if is_geo else 0.0001
    candidate = polygonal_part(set_precision(candidate, grid, mode='valid_output'))
    if candidate.is_empty or not candidate.is_valid or candidate.area <= 0:
        raise ValueError('Exported polygon cannot be repaired without data loss')
    return candidate


def verify_and_repair_written_polygons(files, site, metric_crs, diag,
                                       max_area_change_m2=100.0,
                                       max_relative_area_change=0.001):
    """Reopen actual SHP files, repair collapsed rings and verify again.

    Both the catchment and river_basins datasets are rewritten together from
    one geometry per exact TDX_ID to avoid cross-layer inconsistencies.
    Fixes are provisional and reported. No significant area changes permitted.
    """
    catch = gpd.read_file(files['catchment'])
    basins = gpd.read_file(files['river_basins'])
    if len(catch) != len(basins) or len(catch) == 0:
        raise ValueError('Shapefile re-read: HRU/basin counts differ')
    if catch['TDX_ID'].duplicated().any() or basins['TDX_ID'].duplicated().any():
        raise ValueError('Shapefile re-read: duplicate TDX_ID')
    if set(catch['TDX_ID']) != set(basins['TDX_ID']):
        raise ValueError('Shapefile re-read: HRU/basin TDX_ID mismatch')
    if catch.crs != basins.crs:
        raise ValueError('Shapefile re-read: catchment and basin CRS mismatch')
    report = []
    repairs = []
    original_union = None
    try:
        original_union = catch.to_crs(metric_crs).geometry.union_all()
    except Exception:
        # Invalid rings can prevent union: individual area checks still apply.
        pass
    # On the first pass, validate each layer. For each broken ID, repair a
    # common geometry and propagate it to both layers.
    for rid in sorted(set(catch['TDX_ID'])):
        i = catch.index[catch['TDX_ID'] == rid][0]
        j = basins.index[basins['TDX_ID'] == rid][0]
        gc = catch.at[i, 'geometry']
        gb = basins.at[j, 'geometry']
        if (gc is not None and not gc.is_empty and gc.is_valid and
                gb is not None and not gb.is_empty and gb.is_valid):
            # Basins and catchment originate from the same HRU.
            continue
        reasons = ';'.join(explain_validity(g) if g is not None else 'NULL'
                           for g in (gc,gb))
        source = gc if gc is not None and not gc.is_empty else gb
        fixed = fix_collapsed_polygon_components(source, catch.crs)
        source_metric = gpd.GeoSeries([source], crs=catch.crs).to_crs(metric_crs).iloc[0]
        fixed_metric = gpd.GeoSeries([fixed], crs=catch.crs).to_crs(metric_crs).iloc[0]
        before = float(source_metric.area)
        after = float(fixed_metric.area)
        delta = abs(after-before)
        if (delta > max_area_change_m2 or
                delta > max_relative_area_change * max(before,1.0)):
            raise ValueError(f'Post-SHP fix for {rid} changed area by {delta:.3f} m2')
        catch.at[i, 'geometry'] = fixed
        basins.at[j, 'geometry'] = fixed
        repairs.append(str(rid))
        report.append({'reach_id':str(rid),'original_error':reasons,
                       'before_area_m2':round(before,6),
                       'after_area_m2':round(after,6),
                       'abs_area_change_m2':round(delta,6),
                       'status':'REPAIRED_REVIEW'})
    report_path = diag / 'step06b_post_shapefile_geometry_qa.csv'
    with report_path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=['reach_id','original_error',
            'before_area_m2','after_area_m2','abs_area_change_m2','status'])
        writer.writeheader()
        writer.writerows(report)
    if repairs:
        # Rewriting a generated preview only, never the Step05b or TDX input.
        for key,frame in [('catchment',catch),('river_basins',basins)]:
            remove_shapefile(files[key])
            write_shapefile(frame,files[key])
    # Mandatory second read: do not report success if serialization recreates
    # the same defects.
    check_c = gpd.read_file(files['catchment'])
    check_b = gpd.read_file(files['river_basins'])
    for key,gdf in [('catchment',check_c),('river_basins',check_b)]:
        bad = gdf.geometry.isna() | gdf.geometry.is_empty | ~gdf.geometry.is_valid
        if bad.any():
            details = [(str(r['TDX_ID']),explain_validity(r.geometry))
                       for _,r in gdf.loc[bad].head(8).iterrows()]
            raise ValueError(f'{key} contains invalid geometry after shapefile '
                             f'repair: {details}')
    if original_union is not None:
        new_union = check_c.to_crs(metric_crs).geometry.union_all()
        changed = float(original_union.symmetric_difference(new_union).area)
        if changed > max(200.0, 1e-5 * max(float(original_union.area), 1.0)):
            raise ValueError(f'Post-export union changed excessively: {changed:.3f} m2')
    return repairs, report_path


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--control',type=Path,required=True)
    ap.add_argument('--preflight-only',action='store_true')
    ap.add_argument('--allow-review',action='store_true')
    ap.add_argument('--overwrite',action='store_true',help='Overwrite only nested preview shapefiles')
    ap.add_argument('--gauge-tolerance-m',type=float,default=50.)
    ap.add_argument('--preview-slope-fallback', choices=('none', 'neighbor-median'), default='none',
                    help='PREVIEW ONLY: for missing/nonpositive TDX slopes, use median of nearby valid reach slopes; flag replacements. Default blocks.')
    args=ap.parse_args()
    with args.control.open('rb') as f:
        cfg=tomllib.load(f)
    site=cfg['domain']['site_id']
    root=Path(cfg['output']['output_directory'])/site
    diag=root/cfg['output'].get('diagnostics_subdirectory','diagnostics')
    input_gpkg=root/'tdx_nested_partition_preview.gpkg'
    summary05b=read_json(diag/'tdx_nested_partition_summary.json')
    outlet=read_json(diag/'tdx_outlet_selection.json')
    if summary05b.get('site_id')!=site or outlet.get('site_id',site)!=site:
        raise ValueError('Control, Step 03 and Step 05b site IDs differ')
    if not summary05b.get('output_written') or not input_gpkg.is_file():
        raise ValueError('Step 05b must write the candidate GeoPackage first')
    if summary05b.get('blockers') or summary05b.get('status')=='BLOCKED':
        raise ValueError(f'Step 05b has blockers: {summary05b.get("blockers")}')
    if summary05b.get('production_export_authorized'):
        raise ValueError('Unexpected Step 05b production authorization; this path always emits previews')
    metric=cfg.get('delineation',{}).get('projected_crs','EPSG:32619')
    export_crs=cfg.get('delineation',{}).get('output_crs','EPSG:4326')
    hru=gpd.read_file(input_gpkg,layer='nested_refined_subbasins').to_crs(metric)
    reach=gpd.read_file(input_gpkg,layer='nested_reaches').to_crs(metric)
    boundary=gpd.read_file(input_gpkg,layer='refinement_boundary').to_crs(metric)
    if hru.empty or reach.empty or boundary.empty:
        raise ValueError('Step 05b contains an empty layer')
    id_hru=column(hru,'reach_id')
    id_reach=column(reach,'reach_id')
    ds_col=column(reach,'downstream_id')
    parent_hru=column(hru,'parent_reach_id')
    parent_reach=column(reach,'parent_reach_id')
    hru['_id']=hru[id_hru].map(id_token)
    reach['_id']=reach[id_reach].map(id_token)
    hru['_parent']=hru[parent_hru].map(id_token)
    reach['_parent']=reach[parent_reach].map(id_token)
    reach['_down']=reach[ds_col].map(id_token)
    ids=list(reach['_id'])
    hids=list(hru['_id'])
    if len(set(ids))!=len(ids) or len(set(hids))!=len(hids):
        raise ValueError('Duplicate reach or HRU IDs')
    if set(ids)!=set(hids):
        raise ValueError('Step 05b reach IDs and HRU IDs differ')
    scaled={scale10(s) for s in ids}
    if len(scaled)!=len(ids):
        raise ValueError('IDs collide after CWARHM all-ID x10 conversion')
    # Mixed encoding is NOT used in exported COMID/NextDownID fields. A native
    # 26821 and fractional 2682.1 remain distinct decimal numeric values.
    # Record mixed-encoding risks for provenance, never block a decimal export.
    mixed_codes=[mixed_scaled_id(s) for s in ids]
    mixed_collisions=len(mixed_codes)-len(set(mixed_codes))
    mapping=dict(zip(reach['_id'],reach['_down']))
    if sorted(x for x,d in mapping.items() if d=='0')!=[summary05b['proposed_outlet_piece_id']]:
        raise ValueError('Outlet reach differs from Step 05b summary, or multiple outlets exist')
    for rid,ds in mapping.items():
        if ds!='0' and ds not in mapping:
            raise ValueError(f'NextDown ID {ds} for {rid} missing from exported network')
    for rid in mapping:
        seen=set();current=rid
        while current!='0':
            if current in seen:
                raise ValueError(f'Cycle detected starting from {rid}')
            seen.add(current)
            current=mapping[current]
    if any(int(p) < 0 for p in hru['_parent']) or any(int(p) < 0 for p in reach['_parent']):
        raise ValueError('Invalid native parent reach ID')
    repaired_hrus, geometry_qa_csv = validate_and_repair_hrus(hru, '_id', diag)
    if not all(row.geometry is not None and not row.geometry.is_empty and row.geometry.is_valid and row.geometry.length>0 for row in reach.itertuples()):
        raise ValueError('Invalid/empty Step05b reach geometry')
    union=hru.geometry.union_all()
    bgeom=boundary.geometry.union_all()
    area_sum=float(hru.geometry.area.sum())
    overlap=max(0.,(area_sum-union.area)/1e6)
    coverage=float(union.symmetric_difference(bgeom).area/1e6)
    if overlap>max(0.002,union.area/1e6*0.001) or coverage>max(0.002,bgeom.area/1e6*0.001):
        raise ValueError(f'Candidate HRU coverage invalid: overlap={overlap:.6f} km2, symmetric diff={coverage:.6f} km2')
    source_kind_col=column(reach,'source_kind')
    source_sl_col=column(reach,'slope')
    slope_resolver = TDXSlopeResolver(cfg)
    try:
        slopes, slope_substitutions, slopes_provenance = slope_resolver.resolve(set(reach['_parent']))
    finally:
        slope_resolver.close()
    gauge=gpd.GeoSeries([Point(float(outlet['snapped_lon']),float(outlet['snapped_lat']))],crs='EPSG:4326').to_crs(metric).iloc[0]
    outlet_id=id_token(summary05b['proposed_outlet_piece_id'])
    parent_outlet=id_token(summary05b['original_outlet_parent_tdx_id'])
    outlet_row=reach.loc[reach['_id']==outlet_id].iloc[0]
    outlet_line=only_line(outlet_row.geometry)
    gauge_distance=float(outlet_line.distance(gauge))
    if gauge_distance>args.gauge_tolerance_m:
        raise ValueError(f'Gauge {gauge_distance:.1f}m from proposed outlet reach (> {args.gauge_tolerance_m}m)')
    outlet_original_length=float(outlet_line.length)
    if outlet_id==parent_outlet:
        # Unsplit outlet still has its full TDX reach; truncate at the gauge.
        orientation=outlet.get('line_orientation')
        if orientation not in ('upstream_to_downstream','downstream_to_upstream'):
            raise ValueError('Step 03 outlet orientation not verified')
        at=float(outlet_line.project(gauge))
        if not (0<at<outlet_line.length):
            raise ValueError('Gauge at outlet reach endpoint; cannot truncate safely')
        cut=(substring(outlet_line,at,outlet_line.length) if orientation=='downstream_to_upstream'
             else substring(outlet_line,0,at))
        if not isinstance(cut,LineString) or cut.length<=0:
            raise ValueError('Invalid clipped outlet LineString')
        reach.loc[reach['_id']==outlet_id,'geometry']=cut
    # Split outlet piece has already been terminated at the global gauge cut.
    # Verify an endpoint is reasonably near the snapped gauge.
    else:
        endpoint_gap=min(Point(outlet_line.coords[0]).distance(gauge),
                         Point(outlet_line.coords[-1]).distance(gauge))
        if endpoint_gap>args.gauge_tolerance_m:
            raise ValueError(f'Split outlet piece does not end at gauge; gap={endpoint_gap:.1f}m')
    gauge_distance_after=float(reach.loc[reach['_id']==outlet_id].geometry.iloc[0].distance(gauge))
    hru_records=[]
    for _,row in hru.iterrows():
        sid=row['_id']; poly=row.geometry
        center=gpd.GeoSeries([poly.representative_point()],crs=metric).to_crs('EPSG:4326').iloc[0]
        hru_records.append({'COMID':numeric_id(sid),'TDX_ID':sid,'GRU_ID':numeric_id(sid),
                            'HRU_ID':numeric_id(sid),'HRU_area':float(poly.area),
                            'area':float(poly.area)/1e6,'center_lat':float(center.y),
                            'center_lon':float(center.x),'hru_to_seg':numeric_id(sid),
                            'parent_id':int(row['_parent']),'geometry':poly})
    hru_export=gpd.GeoDataFrame(hru_records,crs=metric)
    # Every native or split reach inherits its original parent's TDX slope.
    # If invalid, recover using positive ORIGINAL source slopes, prioritizing
    # upstream main reaches and then downstream main reaches.
    candidate_slopes = {row['_id']: slopes[int(row['_parent'])]
                        for _, row in reach.iterrows()}
    slope_fallback_csv = diag / 'step06b_preview_slope_substitutions.csv'
    if slope_substitutions:
        # Include native pieces affected by a substituted parent source slope.
        parent_to_pieces = {}
        for _, row in reach.iterrows():
            parent_to_pieces.setdefault(int(row['_parent']), []).append(row['_id'])
        for record in slope_substitutions:
            record['affected_reach_ids'] = ';'.join(parent_to_pieces.get(
                record['parent_tdx_reach_id'], []))
        with slope_fallback_csv.open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(slope_substitutions[0]))
            writer.writeheader()
            writer.writerows(slope_substitutions)

    reach_records = []
    for _, row in reach.iterrows():
        sid = row['_id']
        parent = row['_parent']
        down = row['_down']
        line = only_line(row.geometry)
        val = candidate_slopes[sid]
        reach_records.append({'COMID': numeric_id(sid), 'TDX_ID': sid,
                              'NextDownID': numeric_id(down),
                              'slope': float(val), 'lengthkm': line.length / 1000,
                              'length_m': float(line.length), 'parent_id': int(parent),
                              'geometry': line})
    reach_export = gpd.GeoDataFrame(reach_records, crs=metric)
    river_basins = hru_export[['COMID', 'TDX_ID', 'parent_id', 'geometry']].copy()
    river_basins['area_km2'] = hru_export['area'].values
    review = list(summary05b.get('review_flags', []))
    if slope_substitutions:
        review.append('PREVIEW_ONLY_SUBSTITUTED_MISSING_TDX_SLOPES')
    if repaired_hrus:
        review.append('STEP06B_MINOR_HRU_GEOMETRY_REPAIR')
    # Existing Step05b findings are accepted for provisional preview, not production.
    status='REVIEW' if review or summary05b.get('status')=='REVIEW' else 'PROVISIONAL'
    manifest={
        'site_id':site,'stage':'STEP06B_NESTED_EXPORT_PREVIEW',
        'preview_only':True,'production_export_authorized':False,
        'status':status,'review_flags':sorted(set(review)),
        'source_05b_geopackage':str(input_gpkg),
        'hru_geometry_qa_csv':str(geometry_qa_csv),
        'hru_geometry_repairs':repaired_hrus,
        'hru_count':len(hru_export),'reach_count':len(reach_export),
        'split_reach_count':sum('.' in x for x in ids),
        'outlet_tdx_id':outlet_id,'outlet_numeric_id':numeric_id(outlet_id),
        'original_outlet_parent_tdx_id':parent_outlet,
        'outlet_original_length_m':round(outlet_original_length,3),
        'outlet_export_length_m':round(float(reach_export.loc[reach_export.TDX_ID==outlet_id,'length_m'].iloc[0]),3),
        'gauge_distance_to_export_outlet_m':round(gauge_distance_after,3),
        'overlap_km2':round(overlap,9),'boundary_symmetric_difference_km2':round(coverage,9),
        'slope_provenance':slopes_provenance,
        'preview_slope_fallback_policy': 'original_parent_then_upstream_then_downstream',
        'preview_slope_substitutions': slope_substitutions,
        'preview_slope_substitutions_csv': str(slope_fallback_csv) if slope_substitutions else None,
        'full_global_id_collision_checked':False,
        'headwater_extension_geometry_generated':False,
        'cross_basin_consistency_approved':False,
        'post_shapefile_geometries_validated': False,
        'numeric_id_policy':'EXACT_DECIMAL_TDX_ID_1DP: native LINKNO and one-decimal split IDs in numeric DBF fields; x10 stored in mapping only',
        'alternative_mixed_x10_collision_count':mixed_collisions,
        'files':{},
    }
    outroot=root/'cwarhm_exports'/'tdx_nested_preview'
    files={
        'catchment':outroot/'shapefiles/catchment'/f'{site}_catchment.shp',
        'river_basins':outroot/'shapefiles/river_basins'/f'{site}_river_basins.shp',
        'river_network':outroot/'shapefiles/river_network'/f'{site}_river_network.shp',
    }
    manifest['files']={k:str(v) for k,v in files.items()}
    meta=outroot/'metadata'
    (diag/'step06b_preflight.json').write_text(json.dumps(manifest,indent=2)+'\n')
    if args.preflight_only:
        print(json.dumps(manifest,indent=2));return
    if status=='REVIEW' and not args.allow_review:
        raise RuntimeError('Step 05b is REVIEW. Use --allow-review for preview-only export.')
    if any(p.exists() for p in files.values()) and not args.overwrite:
        raise FileExistsError('Nested preview already exists; rerun with --overwrite')
    if args.overwrite:
        for p in files.values():remove_shapefile(p)
    for frame, key in [(hru_export,'catchment'),(river_basins,'river_basins'),(reach_export,'river_network')]:
        write_shapefile(frame.to_crs(export_crs),files[key])
    post_repairs, post_qa_csv = verify_and_repair_written_polygons(files,site,metric,diag)
    manifest['post_shapefile_geometry_repairs'] = post_repairs
    manifest['post_shapefile_geometry_qa_csv'] = str(post_qa_csv)
    manifest['post_shapefile_geometries_validated'] = True
    if post_repairs:
        manifest['review_flags'] = sorted(set(manifest['review_flags'] +
            ['STEP06B_POST_SHAPEFILE_COLLAPSED_COMPONENT_REPAIR']))
        manifest['status'] = 'REVIEW'
    meta.mkdir(parents=True,exist_ok=True)
    with (meta/'reach_id_mapping.csv').open('w',newline='') as f:
        writer=csv.writer(f)
        writer.writerow(['TDX_ID_exact','TDX_ID_x10_integer','CWARHM_mixed_scaled_ID','COMID_numeric','parent_native_TDX_ID','NextDownID_exact','NextDownID_x10_integer','NextDownID_mixed_scaled'])
        for _,r in reach.iterrows():
            writer.writerow([r['_id'],scale10(r['_id']),mixed_scaled_id(r['_id']),numeric_id(r['_id']),r['_parent'],r['_down'],scale10(r['_down']),mixed_scaled_id(r['_down']) if r['_down']!='0' else 0])
    (meta/'export_metadata.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2))

if __name__=='__main__':
    try:main()
    except Exception as exc:
        print(f'ERROR: {type(exc).__name__}: {exc}',file=sys.stderr)
        sys.exit(1)
