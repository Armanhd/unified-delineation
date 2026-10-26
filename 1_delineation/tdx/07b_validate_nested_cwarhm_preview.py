#!/usr/bin/env python3
"""Validate Step06b nested preview CWARHM shapefiles, preserving .1..9 reach IDs.

Independent of legacy Step07 which casts IDs to int and therefore must NOT be
used as-is on fractionally split reaches.

Only structural compatibility is tested. Hydrological REVIEW / missing global
native ID collision checks remain unapproved.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import sys
import tomllib
from pathlib import Path
import re

import geopandas as gpd
import numpy as np
from shapely.geometry import Point


def readjson(path):
    if not path.is_file():
        raise FileNotFoundError(f'Cannot find {path}')
    return json.loads(path.read_text())


def f_id(value):
    """Parse DBF numeric ID using integer x10 with tolerance (not int(COMID))."""
    if value is None or not math.isfinite(float(value)) or float(value)<0:
        raise ValueError(f'Invalid DBF reach ID {value}')
    scaled=round(float(value)*10)
    if abs(float(value)*10-scaled)>1.e-5:
        raise ValueError(f'ID cannot be represented with 1 decimal: {value}')
    return int(scaled)


EXACT_ID_RE = re.compile(r'^(0|[1-9][0-9]*)(?:\.([1-9]))?$')


def expected_x10_from_text(value):
    token = str(value).strip()
    match = EXACT_ID_RE.fullmatch(token)
    if match is None:
        raise ValueError(f'Invalid exact TDX_ID: {token!r}')
    return int(match.group(1)) * 10 + int(match.group(2) or 0)


def string_ids(frame, attr='TDX_ID'):
    if attr not in frame:
        raise ValueError(f'Missing exact TDX_ID text field in {list(frame.columns)}')
    return list(frame[attr].astype(str).str.strip())


def check_shape(frame, kind):
    if frame.crs is None or frame.empty:
        raise ValueError(f'{kind} missing CRS or empty')
    if frame.geometry.is_empty.any() or frame.geometry.isna().any():
        raise ValueError(f'{kind} contains null/empty geometry')
    if not frame.geometry.is_valid.all():
        raise ValueError(f'{kind} contains invalid geometries')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--control',type=Path,required=True)
    ap.add_argument('--min-split-reach-in-hru',type=float,default=0.5)
    ap.add_argument('--junction-tolerance-m',type=float,default=50.)
    ap.add_argument('--outlet-tolerance-m',type=float,default=50.)
    args=ap.parse_args()
    with args.control.open('rb') as f:cfg=tomllib.load(f)
    site=cfg['domain']['site_id']
    root=Path(cfg['output']['output_directory'])/site
    exp=root/'cwarhm_exports/tdx_nested_preview'
    meta=readjson(exp/'metadata/export_metadata.json')
    if meta.get('site_id')!=site or not meta.get('preview_only') or meta.get('production_export_authorized'):
        raise ValueError('Nested metadata inconsistent or not preview-only')
    frames={k:gpd.read_file(v) for k,v in meta['files'].items()}
    for k,frame in frames.items():check_shape(frame,k)
    hru=frames['catchment'];bas=frames['river_basins'];river=frames['river_network']
    H=string_ids(hru);B=string_ids(bas);R=string_ids(river)
    ids_ok=(len(H)==len(set(H)) and len(B)==len(set(B)) and len(R)==len(set(R)) and set(H)==set(B)==set(R))
    if not ids_ok:raise ValueError('Duplicate or inconsistent string TDX IDs')
    expected={s:expected_x10_from_text(s) for s in R}
    numeric_ok=True
    for frame in (hru,bas,river):
        for row in frame.itertuples(index=False):
            if f_id(getattr(row,'COMID'))!=expected[str(row.TDX_ID).strip()]:
                numeric_ok=False
    extra_numeric=True
    for _,row in hru.iterrows():
        rid=expected[str(row['TDX_ID']).strip()]
        if any(f_id(row[k])!=rid for k in ('GRU_ID','HRU_ID','hru_to_seg')):
            extra_numeric=False
    downstream={str(r['TDX_ID']).strip():f_id(r['NextDownID']) for _,r in river.iterrows()}
    reverse={v:k for k,v in expected.items()}
    outlets=[key for key,ds in downstream.items() if ds==0]
    invalid=[(key,ds) for key,ds in downstream.items() if ds!=0 and ds not in reverse]
    cycles=[]
    for start in downstream:
        seen=set();at=start
        while at is not None:
            if at in seen:
                cycles.append(start);break
            seen.add(at)
            ds=downstream[at]
            at=reverse.get(ds) if ds else None
            if ds and at is None:break
    if len({f_id(x) for x in river.COMID})!=len(R):
        numeric_ok=False
    hp=hru.to_crs(cfg['delineation']['projected_crs'])
    rp=river.to_crs(cfg['delineation']['projected_crs'])
    bp=bas.to_crs(cfg['delineation']['projected_crs'])
    hru_by={str(r['TDX_ID']).strip():r.geometry for _,r in hp.iterrows()}
    river_by={str(r['TDX_ID']).strip():r.geometry for _,r in rp.iterrows()}
    union=hp.geometry.union_all()
    overlap=max(0.,(float(hp.geometry.area.sum())-union.area)/1e6)
    b_union=bp.geometry.union_all()
    basin_matching=float(b_union.symmetric_difference(union).area/1e6)
    source=root/'tdx_nested_partition_preview.gpkg'
    refined_boundary=gpd.read_file(source,layer='refinement_boundary').to_crs(hp.crs).geometry.union_all()
    refined_mismatch=float(refined_boundary.symmetric_difference(union).area/1e6)
    areas=np.asarray(hp.geometry.area,dtype=float)
    stored=np.asarray(hru.HRU_area,dtype=float)
    stored_area_error=np.max(np.abs(areas-stored)/np.maximum(areas,1.))
    slopes=np.asarray(river.slope,dtype=float)
    lengths=np.asarray(river.lengthkm,dtype=float)
    length_geom=np.asarray(rp.geometry.length,dtype=float)/1000
    length_error=float(np.max(np.abs(lengths-length_geom)/np.maximum(length_geom,1.e-9)))
    coverage=[];gaps=[]
    for sid,geom in river_by.items():
        own=hru_by[sid]
        ratio=float(geom.intersection(own.buffer(30.)).length)/max(geom.length,1.e-10)
        coverage.append((sid,ratio))
        ds=downstream[sid]
        if ds:
            other=river_by[reverse[ds]]
            gap=min(Point(geom.coords[0]).distance(other),Point(geom.coords[-1]).distance(other))
            gaps.append((sid,float(gap)))
    split_coverage=[(sid,frac) for sid,frac in coverage if '.' in sid]
    too_far=[(sid,round(d,2)) for sid,d in gaps if d>args.junction_tolerance_m]
    too_far_split=[(sid,round(frac,4)) for sid,frac in split_coverage if frac<args.min_split_reach_in_hru]
    outlet_id=meta['outlet_tdx_id']
    outlet_sel=readjson(root/'diagnostics/tdx_outlet_selection.json')
    gauge=gpd.GeoSeries([Point(float(outlet_sel['snapped_lon']),float(outlet_sel['snapped_lat']))],crs='EPSG:4326').to_crs(hp.crs).iloc[0]
    gauge_dist=float(river_by[outlet_id].distance(gauge)) if outlet_id in river_by else math.inf
    checks={
      'matching_unique_exact_TDX_IDs':ids_ok,
      'numerical_COMIDs_preserve_fractional_IDs':numeric_ok,
      'HRU_GRU_and_mapping_IDs_match':extra_numeric,
      'exactly_one_expected_outlet':outlets==[outlet_id],
      'valid_downstream_references':not invalid,
      'acyclic_network':not cycles,
      'positive_finite_slopes':bool(np.isfinite(slopes).all() and (slopes>0).all()),
      'positive_finite_lengths':bool(np.isfinite(lengths).all() and (lengths>0).all()),
      'stored_length_matches_geometry':length_error<0.01,
      'stored_HRU_area_matches_geometry':stored_area_error<0.01,
      'HRU_overlap_below_0_002_km2':overlap<0.002,
      'river_basins_match_HRUs':basin_matching<0.002,
      'HRUs_cover_step05b_boundary':refined_mismatch<max(0.002,refined_boundary.area/1e6*0.001),
      'split_reach_inside_own_HRU_at_least_fraction':not too_far_split,
      'outlet_reach_at_gauge':gauge_dist<=args.outlet_tolerance_m,
    }
    checks={name:bool(flag) for name,flag in checks.items()}
    structural_pass=all(checks.values())
    warnings=[]
    if too_far:warnings.append(f'{len(too_far)} junction gaps exceed {args.junction_tolerance_m:g}m')
    if meta.get('review_flags'):warnings.append('Original Step05b hydrological REVIEW flags retained')
    if not meta.get('full_global_id_collision_checked'):warnings.append('Full North American TDX ID collision check pending')
    if not meta.get('headwater_extension_geometry_generated'):warnings.append('DEM-derived headwater reach extension not implemented')
    result={'site_id':site,'stage':'STEP07B_NESTED_PREVIEW_VALIDATION',
            'structural_validation_pass':structural_pass,
            'hydrological_QA_status':meta['status'],
            'overall_status':'BLOCKED' if not structural_pass else 'REVIEW' if meta['status']=='REVIEW' else 'PROVISIONAL',
            'production_export_authorized':False,
            'checks':checks,
            'reach_count':len(river),'hru_count':len(hru),'fractional_ID_count':sum('.' in s for s in R),
            'outlet_TDX_ID':outlet_id,'gauge_to_exported_outlet_m':round(gauge_dist,4),
            'HRU_overlap_km2':round(overlap,8),'HRU_boundary_mismatch_km2':round(refined_mismatch,8),
            'river_basins_difference_km2':round(basin_matching,8),
            'max_junction_gap_m':max((d for _,d in gaps),default=0.),
            'junctions_over_tolerance':too_far[:30],
            'split_reaches_below_coverage_threshold':too_far_split,
            'hydrological_review_flags':meta.get('review_flags',[]),
            'warnings':warnings,
            'full_global_id_collision_checked':False,
            'final_export_authorized':False,
           }
    dest=root/'validation/tdx_nested_preview/step07b_validation_summary.json'
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if not structural_pass:sys.exit(2)

if __name__=='__main__':
    try:main()
    except Exception as exc:
        print(f'ERROR: {type(exc).__name__}: {exc}',file=sys.stderr)
        sys.exit(1)
