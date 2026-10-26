#!/usr/bin/env python3
"""Step 04c: preview geometry cuts and connectivity rewiring for a TDX region.

Read-only TDX stream GeoPackages, immutable 04b split plan, and immutable Step04
regional topology index. Only a handful of parent reach geometries are read.

Outputs are CANDIDATES, NOT an approved hydrofabric. Native polygons are NOT
partitioned; no SUMMA/mizuRoute IDs are finalized. Headwaters, Copernicus D8,
and final whole-continent ID collision checks remain for later stages.

Examples:
  python 1_delineation/tdx/04c_validate_shared_split_geometry.py --region 7020024600
  python 1_delineation/tdx/04c_validate_shared_split_geometry.py --region 7020024600 --max-gap-m 150
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sqlite3
from collections import defaultdict
from pathlib import Path

import geopandas as gpd
from pyproj import CRS, Transformer
from shapely.ops import substring, transform
from shapely import wkb

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INDEX_ROOT = ROOT / 'indexes/tdx_regional_topology'
DEFAULT_PLAN_ROOT = ROOT / 'indexes/gauge_catalogs/tdx_shared_network/nested_split_plans'
DEFAULT_OUTPUT_ROOT = ROOT / 'indexes/gauge_catalogs/tdx_shared_network/validated_geometry'

PARENT_COLUMNS = ['tdx_id','parent_tdx_reach_id','planning_status','planning_qa_flags',
                  'geometry_status','geometry_qa_flags','orientation','metric_crs',
                  'native_length_m','piece_length_sum_m','downstream_parent','downstream_geom_gap_m',
                  'geometry_gap_margin_m','n_pieces','geometry_written','final_export_authorized']
EDGE_COLUMNS = ['source_reach_id','source_original_parent','source_kind',
                'new_downstream_id','target_original_parent','edge_role','qa_flags',
                'final_export_authorized']
GAUGE_COLUMNS = ['site_id','parent_tdx_reach_id','fraction_from_line_start',
                 'upstream_piece_id','downstream_piece_id','orientation',
                 'geometry_qa_flags','final_export_authorized']

def norm(value):
    if value is None:
        return None
    s = str(value).strip()
    if not s or s.lower() in ('nan', 'none', 'null'):
        return None
    if s.endswith('.0') and s[:-2].lstrip('-').isdigit():
        return s[:-2]
    return s


def qname(name):
    return '"' + str(name).replace('"', '""') + '"'


def open_ro(path):
    p = Path(path).resolve()
    if not p.is_file():
        raise FileNotFoundError(p)
    return sqlite3.connect(p.as_uri() + '?mode=ro', uri=True)


def read_csv(path, needed):
    with Path(path).open(newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        missing = set(needed) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f'{path}: missing columns {sorted(missing)}')
        return list(reader)


def save_csv(path, columns, rows):
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def gpkg_geometry(blob):
    if blob is None:
        return None
    data = bytes(blob)
    if len(data) < 8 or data[:2] != b'GP':
        raise ValueError('Invalid GeoPackage geometry header')
    env = (data[3] >> 1) & 7
    sizes = {0:0, 1:32, 2:48, 3:48, 4:64}
    if env not in sizes:
        raise ValueError(f'Unsupported GeoPackage envelope type {env}')
    return wkb.loads(data[8+sizes[env]:])


def source_metadata(index, region):
    row = index.execute("SELECT value FROM metadata WHERE key='signature'").fetchone()
    if not row:
        raise ValueError('Regional topology index missing signature; rebuild with Step 04')
    sig = json.loads(row[0])
    if str(sig.get('tdx_id')) != region:
        raise ValueError(f'Region index has {sig.get("tdx_id")}, not {region}')
    for key, pathkey in (('stream','stream_path'), ('basin','basin_path')):
        path = Path(sig[pathkey])
        if not path.is_file():
            raise FileNotFoundError(path)
        stat = path.stat()
        if stat.st_size != sig[f'{key}_size'] or stat.st_mtime_ns != sig[f'{key}_mtime_ns']:
            raise RuntimeError(f'Stale topology index: source {path} changed. Rebuild Step 04 regional index.')
    return sig


def line_info(con, layer):
    result = con.execute(
        'SELECT column_name,srs_id FROM gpkg_geometry_columns WHERE table_name=?', (layer,)).fetchone()
    if result is None:
        raise ValueError(f'No geometry layer metadata for {layer}')
    col, epsg = result
    fields = {r[1] for r in con.execute(f'PRAGMA table_info({qname(layer)})')}
    if 'fid' not in fields:
        raise ValueError('Source layer missing fid key')
    if not epsg or epsg <= 0:
        raise ValueError(f'Unknown source CRS (srs_id={epsg})')
    return col, CRS.from_epsg(epsg), fields


def choose_metric_crs(line, native_crs):
    centroid = transform(Transformer.from_crs(native_crs, 'EPSG:4326', always_xy=True).transform, line.centroid)
    lon, lat = centroid.x, centroid.y
    if not (-80 < lat < 84):
        raise ValueError(f'Invalid UTM latitude {lat}')
    zone = min(60,max(1,int(math.floor((lon + 180.) / 6.)) + 1))
    return CRS.from_epsg((32600 if lat >= 0 else 32700) + zone)


def num(value, label):
    try:
        v = float(value)
        if not math.isfinite(v):
            raise ValueError()
        return v
    except (ValueError, TypeError):
        raise ValueError(f'Invalid {label}: {value!r}') from None


def identify_orientation(line_m, downstream_m, max_gap_m, min_margin_m):
    if downstream_m is None:
        return None, None, None, ['DOWNSTREAM_GEOMETRY_NOT_FOUND']
    if downstream_m.is_empty:
        return None, None, None, ['DOWNSTREAM_GEOMETRY_EMPTY']
    from shapely.geometry import Point
    start_gap = float(Point(line_m.coords[0]).distance(downstream_m))
    end_gap = float(Point(line_m.coords[-1]).distance(downstream_m))
    best = min(start_gap, end_gap)
    margin = abs(start_gap - end_gap)
    issues = []
    if best > max_gap_m:
        issues.append('DISCONNECTED_DOWNSTREAM_GEOMETRY')
    if margin < min_margin_m:
        issues.append('AMBIGUOUS_FLOW_ORIENTATION')
    orientation = 'downstream_to_upstream' if start_gap < end_gap else 'upstream_to_downstream'
    if issues:
        orientation = None
    return orientation, best, margin, issues


def split_metric_line(line_m, pieces, eps=1e-8):
    if line_m.geom_type != 'LineString' or line_m.is_empty or line_m.length <= 0:
        raise ValueError('Parent stream geometry is not a valid LineString')
    parts = []
    end_prev = 0.
    for piece in pieces:
        start = num(piece['from_fraction'], 'from_fraction')
        end = num(piece['to_fraction'], 'to_fraction')
        if abs(start - end_prev) > eps or end - start <= eps or end > 1.+eps:
            raise ValueError('Candidate fractions are noncontiguous or invalid')
        part = substring(line_m, max(0., min(1., start)), max(0., min(1., end)), normalized=True)
        if part.is_empty or part.geom_type != 'LineString' or part.length < 0.01:
            raise ValueError('Candidate cut produced empty or sub-centimetre line')
        parts.append(part)
        end_prev = end
    if abs(end_prev - 1.) > eps:
        raise ValueError('Candidate fractions do not cover full parent reach')
    return parts


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--region', required=True)
    ap.add_argument('--topology-index', type=Path)
    ap.add_argument('--plan-dir', type=Path, default=DEFAULT_PLAN_ROOT)
    ap.add_argument('--output-dir', type=Path, default=DEFAULT_OUTPUT_ROOT)
    ap.add_argument('--max-gap-m', type=float, default=200.)
    ap.add_argument('--min-orientation-margin-m', type=float, default=5.)
    ap.add_argument('--min-cut-length-m', type=float, default=1.)
    ap.add_argument('--strict', action='store_true', help='Nonzero exit for unresolved geometry or planner QA')
    args = ap.parse_args()
    region = args.region.strip()
    if not region.isdigit():
        ap.error('--region must be numeric')
    if min(args.max_gap_m, args.min_orientation_margin_m, args.min_cut_length_m) < 0:
        ap.error('Distance controls must be nonnegative')
    index_path = args.topology_index or (DEFAULT_INDEX_ROOT / f'{region}_topology.sqlite')
    plan_file = args.plan_dir / f'{region}_nested_split_parents.csv'
    piece_file = args.plan_dir / f'{region}_nested_candidate_pieces.csv'
    parents = read_csv(plan_file, ['parent_tdx_reach_id','group_status','qa_flags','n_candidate_pieces'])
    pieces = read_csv(piece_file, ['parent_tdx_reach_id','candidate_comid_text','from_fraction','to_fraction','left_boundary_site','right_boundary_site','candidate_integer_x10'])
    if not parents:
        raise ValueError('Empty split parent plan')
    parent_rows = {norm(r['parent_tdx_reach_id']): r for r in parents}
    if len(parent_rows) != len(parents):
        raise ValueError('Duplicate split parent row')
    by_parent = defaultdict(list)
    for piece in pieces:
        parent = norm(piece['parent_tdx_reach_id'])
        if parent not in parent_rows:
            raise ValueError(f'Piece for unplanned parent {parent}')
        by_parent[parent].append(piece)
    for parent, record in parent_rows.items():
        if len(by_parent[parent]) != int(record['n_candidate_pieces']):
            raise ValueError(f'{parent}: candidate piece count mismatch')
        by_parent[parent].sort(key=lambda x:num(x['from_fraction'], 'from_fraction'))
    print(f'Loaded {len(parents)} split parents and {len(pieces)} candidate pieces', flush=True)
    topo = open_ro(index_path)
    source = None
    try:
        sig = source_metadata(topo, region)
        source = open_ro(sig['stream_path'])
        geom_col, native_crs, src_cols = line_info(source, sig['stream_layer'])
        slope_col = 'Slope' if 'Slope' in src_cols else None
        fid_by_rid = {}
        ds_by_rid = {}
        for rid, fid, ds in topo.execute('SELECT rid,fid,ds FROM links'):
            rid = norm(rid)
            fid_by_rid[rid] = fid
            ds_by_rid[rid] = norm(ds)
        all_native = set(fid_by_rid)
        def fetch(reach):
            fid = fid_by_rid.get(reach)
            if fid is None:
                return None, None
            cols = qname(geom_col) + (','+qname(slope_col) if slope_col else '')
            r = source.execute(f'SELECT {cols} FROM {qname(sig["stream_layer"])} WHERE fid=?', (fid,)).fetchone()
            if not r:
                return None, None
            return gpkg_geometry(r[0]), (r[1] if slope_col else None)
        result_rows = []
        segment_rows = []
        geom_records = []
        orientations = {}
        cached = {}
        global_id_codes = set()
        emitted_piece_ids = set()
        for parent in sorted(parent_rows, key=lambda s:int(s) if s.isdigit() else s):
            pr = parent_rows[parent]
            qa = set(filter(None, pr.get('qa_flags','').split(';')))
            geometry_issues = set()
            metric_crs = ''
            dsgap = None; margin = None; original_len = None; split_len = None
            orientation = None
            native, slope = fetch(parent)
            orig_ds = ds_by_rid.get(parent)
            if native is None or native.is_empty or native.geom_type != 'LineString':
                geometry_issues.add('INVALID_OR_MISSING_PARENT_GEOMETRY')
            else:
                metric = choose_metric_crs(native, native_crs)
                metric_crs = metric.to_string()
                trans = Transformer.from_crs(native_crs, metric, always_xy=True).transform
                line_m = transform(trans, native)
                original_len = float(line_m.length)
                cached[parent] = (native, slope, line_m, metric)
                ds_geom, _ = fetch(orig_ds)
                downstream_m = transform(trans, ds_geom) if ds_geom is not None and not ds_geom.is_empty else None
                orientation, dsgap, margin, issues = identify_orientation(
                    line_m, downstream_m, args.max_gap_m, args.min_orientation_margin_m)
                geometry_issues.update(issues)
                if not geometry_issues:
                    try:
                        cuts = split_metric_line(line_m, by_parent[parent])
                        if min(seg.length for seg in cuts) < args.min_cut_length_m:
                            raise ValueError('Candidate segment shorter than minimum')
                        split_len = sum(seg.length for seg in cuts)
                        if abs(split_len-original_len)>max(0.05, original_len*1.e-5):
                            raise ValueError('Split piece lengths do not sum to original')
                        back = Transformer.from_crs(metric,native_crs,always_xy=True).transform
                        for piece, geom in zip(by_parent[parent], cuts):
                            piece_id = piece['candidate_comid_text']
                            code = int(piece['candidate_integer_x10'])
                            # Compare actual exported IDs, not a hypothetical mixed
                            # integer-x10 encoding.  For example, native LINKNO 26821
                            # and split reach 2682.1 are different exact IDs.
                            # Keep integer_x10 for diagnostics only.
                            actual_id = norm(piece_id)
                            if actual_id in all_native or actual_id in emitted_piece_ids:
                                raise ValueError('Duplicate exact candidate reach ID: ' + str(actual_id))
                            emitted_piece_ids.add(actual_id)
                            global_id_codes.add(code)
                            segment_rows.append({'tdx_id':region,'parent_reach_id':parent,
                                'reach_id':piece_id,'integer_x10':code,
                                'piece_index':int(piece['piece_index_from_line_start']) if piece.get('piece_index_from_line_start') else len(segment_rows),
                                'from_frac':num(piece['from_fraction'],'from_fraction'),
                                'to_frac':num(piece['to_fraction'],'to_fraction'),
                                'length_m':float(geom.length), 'slope':float(slope) if slope is not None else None,
                                'orientation':orientation,'planning_qa_flags':pr.get('qa_flags',''),
                                'geometry':transform(back,geom)})
                    except Exception as exc:
                        geometry_issues.add('SPLIT_VALIDATION_FAILED:'+str(exc))
            orientations[parent] = orientation
            result_rows.append({
                'tdx_id':region,'parent_tdx_reach_id':parent,
                'planning_status':pr['group_status'], 'planning_qa_flags':pr.get('qa_flags',''),
                'geometry_status':'PASS' if not geometry_issues else 'REVIEW',
                'geometry_qa_flags':';'.join(sorted(geometry_issues)),
                'orientation':orientation or '', 'metric_crs':metric_crs,
                'native_length_m':round(original_len,3) if original_len is not None else '',
                'piece_length_sum_m':round(split_len,3) if split_len is not None else '',
                'downstream_parent':orig_ds or '',
                'downstream_geom_gap_m':round(dsgap,3) if dsgap is not None else '',
                'geometry_gap_margin_m':round(margin,3) if margin is not None else '',
                'n_pieces':len(by_parent[parent]), 'geometry_written':not geometry_issues,
                'final_export_authorized':False})
        # Do not publish invalid split pieces; avoid partial geometry for a parent.
        geometry_bad = {row['parent_tdx_reach_id'] for row in result_rows if row['geometry_status'] != 'PASS'}
        segment_rows = [r for r in segment_rows if r['parent_reach_id'] not in geometry_bad]
        valid = {r['parent_reach_id'] for r in segment_rows}
        pieces_by_parent = {p:by_parent[p] for p in valid}
        ids_by_parent = {p:[q['candidate_comid_text'] for q in pieces_by_parent[p]] for p in valid}
        upstreammost = {}
        downstreammost = {}
        for parent, ids in ids_by_parent.items():
            ori = orientations[parent]
            upstreammost[parent] = ids[-1] if ori=='downstream_to_upstream' else ids[0]
            downstreammost[parent] = ids[0] if ori=='downstream_to_upstream' else ids[-1]
        edges = []
        for parent in sorted(valid):
            ids = ids_by_parent[parent]
            ori = orientations[parent]
            ds_parent = ds_by_rid.get(parent)
            ds_target = upstreammost.get(ds_parent, ds_parent)
            for i, piece_id in enumerate(ids):
                ds_id = (ids[i-1] if i>0 else ds_target) if ori=='downstream_to_upstream' else (ids[i+1] if i+1<len(ids) else ds_target)
                edge_qa = '' if ds_parent not in parent_rows or ds_parent in valid else 'DOWNSTREAM_PARENT_SPLIT_UNRESOLVED'
                edges.append({'source_reach_id':piece_id,'source_original_parent':parent,
                    'source_kind':'SPLIT_PIECE','new_downstream_id':ds_id or '',
                    'target_original_parent':ds_parent or '',
                    'edge_role':'PIECE_LINK','qa_flags':edge_qa,'final_export_authorized':False})
        for orig, ds in ds_by_rid.items():
            if orig in parent_rows or ds not in parent_rows:
                continue
            if ds in valid:
                edges.append({'source_reach_id':orig,'source_original_parent':orig,
                    'source_kind':'UNCHANGED_UPSTREAM_REACH','new_downstream_id':upstreammost[ds],
                    'target_original_parent':ds,'edge_role':'REWIRED_TO_SPLIT_PARENT',
                    'qa_flags':'','final_export_authorized':False})
            else:
                edges.append({'source_reach_id':orig,'source_original_parent':orig,
                    'source_kind':'UNCHANGED_UPSTREAM_REACH','new_downstream_id':'',
                    'target_original_parent':ds,'edge_role':'PENDING_SPLIT_PARENT_GEOMETRY',
                    'qa_flags':'TARGET_SPLIT_PARENT_UNRESOLVED','final_export_authorized':False})
        # Capture every proposed gauge cut as upstream/downstream piece IDs.
        gauge_rows = []
        for parent, plist in pieces_by_parent.items():
            ori = orientations[parent]
            for i, piece in enumerate(plist[:-1]):
                site = piece['right_boundary_site']
                if not site:
                    continue
                left = piece['candidate_comid_text']
                right = plist[i+1]['candidate_comid_text']
                gauge_rows.append({'site_id':site,'parent_tdx_reach_id':parent,
                    'fraction_from_line_start':piece['to_fraction'],
                    'upstream_piece_id':right if ori=='downstream_to_upstream' else left,
                    'downstream_piece_id':left if ori=='downstream_to_upstream' else right,
                    'orientation':ori,
                    'geometry_qa_flags':next((r['geometry_qa_flags'] for r in result_rows if r['parent_tdx_reach_id']==parent),''),
                    'final_export_authorized':False})
        # Ensure each emitted rewired source has exactly one NextDown candidate.
        if len({e['source_reach_id'] for e in edges}) != len(edges):
            raise RuntimeError('Duplicate source reach in rewiring table')
        if len({r['reach_id'] for r in segment_rows}) != len(segment_rows):
            raise RuntimeError('Duplicate split piece ID')
        output = args.output_dir
        output.mkdir(parents=True, exist_ok=True)
        prefix = f'{region}_nested'
        save_csv(output / f'{prefix}_geometry_qa.csv', PARENT_COLUMNS, result_rows)
        save_csv(output / f'{prefix}_candidate_nextdown.csv', EDGE_COLUMNS, edges)
        save_csv(output / f'{prefix}_gauge_cut_nodes.csv', GAUGE_COLUMNS, gauge_rows)
        out_gpkg = output / f'{prefix}_split_reaches_preview.gpkg'
        if segment_rows:
            gg = gpd.GeoDataFrame(segment_rows, geometry='geometry', crs=native_crs).to_crs('EPSG:4326')
            gg.to_file(out_gpkg, layer='candidate_split_reaches', driver='GPKG')
        elif out_gpkg.is_file():
            out_gpkg.unlink()  # Remove stale preview after a no-valid-geometry rerun.
        planning_bad = [r['parent_tdx_reach_id'] for r in result_rows if r['planning_qa_flags']]
        geom_bad = [r['parent_tdx_reach_id'] for r in result_rows if r['geometry_status'] != 'PASS']
        summary = {
            'tdx_id':region,
            'candidate_parents':len(parents),'candidate_pieces':len(pieces),
            'planning_review_parents':len(planning_bad),
            'geometry_pass_parents':len(parents)-len(geom_bad),
            'geometry_review_parents':len(geom_bad),
            'split_geometry_preview_pieces':len(segment_rows),
            'candidate_nextdown_rows':len(edges),
            'gauge_cut_node_rows':len(gauge_rows),
            'regional_native_id_collision_checked':True,
            'actual_reach_id_encoding':'EXACT_DECIMAL_TDX_ID_1DP',
            'full_global_native_id_collision_checked':False,
            'hru_polygon_partition_performed':False,
            'copernicus_d8_refinement_performed':False,
            'source_gpkg_modified':False,
            'fractional_ids_finalized':False,
            'final_export_authorized':False,
            'output_geopackage':str(out_gpkg) if segment_rows else None,
            'planning_review_parent_ids':planning_bad,
            'geometry_review_parent_ids':geom_bad,
            'note':'Read-only source, diagnostic split geometry and NextDown candidates only. All exports require further QA and HRU partition.',
        }
        (output / f'{prefix}_geometry_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        print(json.dumps(summary,indent=2))
        if args.strict and (planning_bad or geom_bad or any(e['qa_flags'] for e in edges)):
            raise RuntimeError('Plan contains unresolved QA; see output files')
    finally:
        if source is not None:
            source.close()
        topo.close()

if __name__=='__main__':
    main()
