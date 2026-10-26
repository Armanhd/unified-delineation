#!/usr/bin/env python3
"""Step 05b: PREVIEW ONLY of DEM-guided nested TDX HRU partitions.

Inputs: Step04 native GeoPackage, Step05 refined GeoPackage for target and for
all gauges cutting TDX reaches in its traced upstream network; Step04c regional
split preview GeoPackage, candidate NextDown edges and gauge cut-node CSV.

The upstream-gauge DEM watershed is a *mask* that partitions its original TDX
parent HRU at the gauge. This is NOT an independent terrain-derived subreach
watershed, not an approved CWARHM export, and not a source-data modification.
Unresolved DEM alignment, nested-mask inconsistency and area QA are retained.

Example:
 python 1_delineation/tdx/05b_build_nested_hru_preview.py \
   --control 0_control_files/tdx_camels_spat/control_CAN_02HC025.toml \
   --preflight-only
 python 1_delineation/tdx/05b_build_nested_hru_preview.py \
   --control 0_control_files/tdx_camels_spat/control_CAN_02HC025.toml \
   --allow-review
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import tomllib
from collections import defaultdict, deque
from pathlib import Path

import geopandas as gpd
import numpy as np
from shapely.geometry import GeometryCollection
from shapely.ops import unary_union
from shapely import make_valid


def normalize_id(value):
    """Keep fractions as text; never cast fractional IDs through float."""
    if value is None:
        return ''
    s = str(value).strip()
    if s.endswith('.0') and s[:-2].lstrip('-').isdigit():
        return s[:-2]
    return s


def poly(geometry):
    """Return valid polygonal part (can be an empty GeometryCollection)."""
    if geometry is None or geometry.is_empty:
        return GeometryCollection()
    if not geometry.is_valid:
        geometry = make_valid(geometry)
    if geometry.geom_type in ('Polygon', 'MultiPolygon'):
        return geometry
    if geometry.geom_type == 'GeometryCollection':
        parts = [g for g in geometry.geoms if g.geom_type in ('Polygon', 'MultiPolygon') and not g.is_empty]
        return unary_union(parts) if parts else GeometryCollection()
    raise ValueError(f'Expected polygonal geometry, got {geometry.geom_type}')


def read_csv(path, required=()):
    with Path(path).open(newline='', encoding='utf-8-sig') as f:
        r = csv.DictReader(f)
        missing = set(required).difference(r.fieldnames or ())
        if missing:
            raise ValueError(f'{path} missing columns {sorted(missing)}')
        return list(r)


def read_geo(path, layer, metric):
    frame = gpd.read_file(path, layer=layer)
    if frame.crs is None:
        raise ValueError(f'{path}::{layer} is missing CRS')
    return frame.to_crs(metric)


def load_site_boundary(output_root, site, metric):
    path = output_root/site/'tdx_upstream_refined.gpkg'
    summary = output_root/site/'diagnostics'/'tdx_refinement_summary.json'
    if not path.is_file() or not summary.is_file():
        raise FileNotFoundError(f'{site}: Step05 refined GPKG and summary required ({path})')
    frame = read_geo(path,'refinement_boundary',metric)
    if frame.empty:
        raise ValueError(f'{site}: empty refinement boundary')
    return poly(frame.geometry.union_all()), json.loads(summary.read_text())


def valid_network(reaches, outlet_id):
    """Validate candidate reach -> NextDown mapping; outlet explicitly terminates."""
    if outlet_id not in reaches:
        raise ValueError('Selected outlet ID absent from candidate reach network')
    for rid, record in reaches.items():
        ds = record['downstream_id']
        if rid == outlet_id:
            if ds != '0':
                raise ValueError('Outlet must terminate with downstream ID 0')
        elif ds not in reaches:
            raise ValueError(f'Candidate reach {rid} points outside network: {ds}')
    for start in reaches:
        visited = set()
        current = start
        while current != '0':
            if current in visited:
                raise ValueError(f'Cycle in proposed network starting at {start}')
            visited.add(current)
            current = reaches[current]['downstream_id']
    return True


def build_reach_graph(native_reaches, split_reaches, edges, parent_to_pieces, outlet_parent, site, cuts):
    """Preserve Step04 native drainage connectivity when inserting regional split pieces.

    Step04c candidate edges can reference split pieces outside this basin, or an
    incomplete regional candidate topology. Blindly replacing each native edge
    with those candidates can disconnect nearly all upstream reaches. Rebuild
    connections from the basin's original Step04 downstream IDs and use the
    regional geometry only to define the order of pieces within a split parent.

    Never invent hydrological connections: every inter-parent edge must already
    exist in the original Step04 network. Return diagnostics of edge differences.
    """
    originals = {}
    for row in native_reaches.itertuples():
        rid = normalize_id(row.reach_id)
        if rid in originals:
            raise ValueError(f'Duplicate native reach {rid}')
        originals[rid] = {
            'reach_id': rid, 'parent_id': rid, 'geometry': row.geometry,
            'downstream_id': normalize_id(row.downstream_id),
            'source_kind': 'NATIVE',
        }
    if outlet_parent not in originals:
        raise ValueError(f'Original outlet parent {outlet_parent} absent from Step04')

    generated = {}
    split_flow_order = {}
    for rid, original in originals.items():
        if rid not in parent_to_pieces:
            generated[rid] = dict(original)
            continue
        piece_ids = list(parent_to_pieces[rid])
        if not piece_ids or len(set(piece_ids)) != len(piece_ids):
            raise ValueError(f'Invalid or duplicate split pieces for parent {rid}')
        for pid in piece_ids:
            if pid not in split_reaches:
                raise ValueError(f'Planned piece {pid} lacks geometry')
            if normalize_id(split_reaches[pid].get('parent_id')) != rid:
                raise ValueError(f'Split piece {pid} has incorrect parent ID')
        orientation = split_reaches[piece_ids[0]].get('orientation')
        if orientation not in ('downstream_to_upstream', 'upstream_to_downstream'):
            raise ValueError(f'Unknown split line orientation for {rid}: {orientation}')
        if any(split_reaches[pid].get('orientation') != orientation for pid in piece_ids):
            raise ValueError(f'Inconsistent split orientation in parent {rid}')
        # Regional piece_index is from line start. First entry in flow order
        # is upstream-most; last is downstream-most.
        flow_order = (list(reversed(piece_ids)) if orientation == 'downstream_to_upstream'
                      else piece_ids)
        split_flow_order[rid] = flow_order
        for pid in piece_ids:
            if pid in generated or pid in originals:
                raise ValueError(f'Candidate split ID collides with native reach: {pid}')
            record = dict(split_reaches[pid])
            record.update(reach_id=pid, parent_id=rid, source_kind='SPLIT_PIECE')
            generated[pid] = record

    def upstream_entry(parent):
        """Reach/piece entered by upstream flow at the upstream end of a parent."""
        return split_flow_order[parent][0] if parent in split_flow_order else parent

    # Reconstruct inter-parent connectivity solely from original Step04 IDs.
    # Original upstream reach X -> original downstream reach Y becomes
    # X (or X's downstream-most piece) -> Y's upstream-most piece.
    repaired_edge_rows = []
    for rid, original in originals.items():
        original_down = original['downstream_id']
        if original_down in ('', '0', '-1', 'None'):
            target = '0'
        elif original_down in originals:
            target = upstream_entry(original_down)
        else:
            # The outlet may normally route outside the selected watershed.
            # Other out-of-basin links must not be silently invented.
            if rid == outlet_parent:
                target = '0'
            else:
                raise ValueError(f'Native reach {rid} routes to absent parent {original_down}')
        if rid in split_flow_order:
            flow = split_flow_order[rid]
            for upstream_piece, downstream_piece in zip(flow, flow[1:]):
                generated[upstream_piece]['downstream_id'] = downstream_piece
            generated[flow[-1]]['downstream_id'] = target
        else:
            generated[rid]['downstream_id'] = target

    # Compare candidate Step04c edges with reconstructed Step04 connectivity.
    # This is advisory, since the originals are authoritative for this basin.
    for rid, record in generated.items():
        if rid in edges:
            candidate = normalize_id(edges[rid])
            if candidate != record['downstream_id']:
                repaired_edge_rows.append({
                    'reach_id': rid,
                    'step04c_candidate_downstream': candidate,
                    'restored_step04_downstream': record['downstream_id'],
                })

    own_cuts = [r for r in cuts if r['site_id'] == site and
                normalize_id(r['parent_tdx_reach_id']) == outlet_parent]
    if len(own_cuts) > 1:
        raise ValueError(f'Duplicate cut-node assignment for outlet {site}')
    outlet_id = normalize_id(own_cuts[0]['upstream_piece_id']) if own_cuts else outlet_parent
    if outlet_id not in generated:
        raise ValueError(f'Cannot identify outlet reach/piece {outlet_id}')

    # Keep only candidates routing upstream of this specific gauge cut.
    incoming = defaultdict(list)
    for rid, record in generated.items():
        incoming[record['downstream_id']].append(rid)
    reachable = set()
    todo = deque([outlet_id])
    while todo:
        rid = todo.popleft()
        if rid in reachable:
            continue
        reachable.add(rid)
        todo.extend(incoming.get(rid, []))
    kept = {rid: dict(generated[rid]) for rid in reachable}
    kept[outlet_id]['downstream_id'] = '0'
    valid_network(kept, outlet_id)

    # Check against the original Step04 *parent* network. Any parent draining
    # into this gauge in Step04 must remain represented in the nested preview.
    expected_parents = set()
    for rid in originals:
        seen = set()
        at = rid
        while at in originals and at not in seen:
            if at == outlet_parent:
                expected_parents.add(rid)
                break
            seen.add(at)
            at = originals[at]['downstream_id']
    retained_parents = {record['parent_id'] for record in kept.values()}
    missing_parents = expected_parents - retained_parents
    if missing_parents:
        raise ValueError('Step05b connectivity lost original upstream parents: '
                         + ', '.join(sorted(missing_parents)[:15]))
    details = {
        'original_native_reach_count': len(originals),
        'expected_upstream_parent_count': len(expected_parents),
        'retained_parent_count': len(retained_parents),
        'retained_candidate_reach_count': len(kept),
        'step04c_edges_different_from_native': repaired_edge_rows,
        'step04c_edge_difference_count': len(repaired_edge_rows),
        'method': 'RECONSTRUCT_FROM_STEP04_NATIVE_PARENT_TOPOLOGY',
    }
    return kept, outlet_id, set(originals), details


def piece_partition(parent_geom, planned_pieces, cuts, cut_watersheds, orientation):
    """Partition parent polygon using hydrologically upstream masks at gauge cuts.

    Piece ordering follows original line start; line orientation controls how
    watershed differences are translated into consecutive reach polygons.
    """
    if len(planned_pieces) != len(cuts)+1:
        raise ValueError('One more reach piece than gauge cuts required')
    if orientation not in ('downstream_to_upstream','upstream_to_downstream'):
        raise ValueError(f'Unknown original reach orientation: {orientation}')
    parent_geom = poly(parent_geom)
    if parent_geom.is_empty:
        raise ValueError('Empty native parent polygon')
    masks = [poly(cut_watersheds[r['site_id']].intersection(parent_geom)) for r in cuts]
    # Increasing along-line fractions -> upstream masks shrink for D2U,
    # and increase for U2D. Discrepancies indicate inconsistent DEM flow networks.
    non_nested_m2 = 0.0
    for a,b in zip(masks[:-1], masks[1:]):
        violation = b.difference(a) if orientation=='downstream_to_upstream' else a.difference(b)
        non_nested_m2 += float(violation.area)
    parts=[]
    if orientation=='downstream_to_upstream':
        for i in range(len(planned_pieces)):
            if i==0:
                geom = parent_geom.difference(masks[0])
            elif i==len(planned_pieces)-1:
                geom = masks[-1]
            else:
                geom = masks[i-1].difference(masks[i])
            parts.append(poly(geom))
    else:
        for i in range(len(planned_pieces)):
            if i==0:
                geom = masks[0]
            elif i==len(planned_pieces)-1:
                geom = parent_geom.difference(masks[-1])
            else:
                geom = masks[i].difference(masks[i-1])
            parts.append(poly(geom))
    # Resolve D8 overlap inconsistencies deterministically while recording them.
    used = GeometryCollection()
    disjoint = []
    for geom in parts:
        geom = poly(geom.difference(used))
        disjoint.append(geom)
        used = poly(used.union(geom))
    missing_m2 = float(parent_geom.difference(used).area)
    overlap_m2 = max(0., sum(g.area for g in parts)-float(used.area))
    return dict(zip(planned_pieces,disjoint)), {
        'non_nested_mask_area_km2':non_nested_m2/1e6,
        'coverage_gap_km2':missing_m2/1e6,
        'overlap_resolved_km2':overlap_m2/1e6,
        'parent_area_km2':float(parent_geom.area)/1e6}


def prepare(control, project, allow_review, preflight_only, overwrite, max_mask_mismatch_pct):
    with control.open('rb') as f:
        cfg = tomllib.load(f)
    site = cfg['domain']['site_id']
    region = str(cfg['input']['tdx_id'])
    metric = cfg['delineation']['projected_crs']
    output_root = Path(cfg['output']['output_directory']).expanduser().resolve()
    root = output_root/site
    diag = root/cfg['output'].get('diagnostics_subdirectory','diagnostics')
    dest_summary = diag/'step05b'/f'{site}_nested_partition_summary.json'
    dest_summary.parent.mkdir(parents=True,exist_ok=True)
    result = {'site_id':site,'tdx_id':region,'stage':'05b_NESTED_HRU_CANDIDATE',
              'preview_only':True,'production_export_authorized':False,
              'full_global_id_collision_checked':False,
              'headwater_extension_geometry_generated':False,
              'status':'BLOCKED','blockers':[],'review_flags':[],
              'note':'Candidate nested-D8 HRUs and proposed reach links. Do not feed to CWARHM/SUMMA until full global ID and hydrological QA.'}
    blockers = result['blockers']; reviews = result['review_flags']
    def block(message):
        blockers.append(message)
    geometry_root = project/'indexes/gauge_catalogs/tdx_shared_network/validated_geometry'
    stem = f'{region}_nested'
    split_gpkg = geometry_root/f'{stem}_split_reaches_preview.gpkg'
    edge_file = geometry_root/f'{stem}_candidate_nextdown.csv'
    gauge_file = geometry_root/f'{stem}_gauge_cut_nodes.csv'
    parent_qa_file = geometry_root/f'{stem}_geometry_qa.csv'
    required_files = [root/'tdx_upstream_native.gpkg',root/'tdx_upstream_refined.gpkg',
                      diag/'tdx_refinement_summary.json', split_gpkg,edge_file,gauge_file,parent_qa_file]
    for p in required_files:
        if not p.is_file():
            block(f'MISSING_INPUT:{p}')
    if blockers:
        dest_summary.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result,indent=2))
        return 2

    native = read_geo(root/'tdx_upstream_native.gpkg','upstream_reaches',metric)
    parent_hrus = read_geo(root/'tdx_upstream_refined.gpkg','refined_subbasins',metric)
    if native.empty or parent_hrus.empty:
        block('EMPTY_STEP04_OR_STEP05_LAYER')
    original_ids = {normalize_id(x) for x in native['reach_id']}
    parent_hrus['reach_id'] = parent_hrus.reach_id.map(normalize_id)
    if parent_hrus.reach_id.duplicated().any():
        block('DUPLICATED_NATIVE_HRU_REACH_ID')
    if set(parent_hrus.reach_id)!=original_ids:
        block('STEP04_STEP05_PARENT_IDS_DIFFER')
    parent_polys = {r.reach_id:poly(r.geometry) for r in parent_hrus.itertuples()}
    geometry_rows = read_csv(parent_qa_file,['parent_tdx_reach_id','planning_qa_flags','geometry_status'])
    geometry_qa = {normalize_id(x['parent_tdx_reach_id']):x for x in geometry_rows}
    cuts = read_csv(gauge_file,['site_id','parent_tdx_reach_id','fraction_from_line_start','upstream_piece_id'])
    split = read_geo(split_gpkg,'candidate_split_reaches',metric)
    split['reach_id'] = split['reach_id'].map(normalize_id)
    split_by_id = {r.reach_id:{'geometry':r.geometry,'orientation':r.orientation,
                              'slope':float(r.slope) if r.slope is not None and math.isfinite(float(r.slope)) else None,
                              'parent_id':normalize_id(r.parent_reach_id),
                              'length_m':float(r.length_m)} for r in split.itertuples()}
    if len(split_by_id)!=len(split):
        block('DUPLICATE_SPLIT_REACH_ID')
    pieces_by_parent = defaultdict(list)
    for r in split.itertuples():
        pieces_by_parent[normalize_id(r.parent_reach_id)].append((int(r.piece_index),r.reach_id))
    pieces_by_parent = {p:[x for _,x in sorted(v)] for p,v in pieces_by_parent.items()}
    edges = read_csv(edge_file,['source_reach_id','new_downstream_id','qa_flags'])
    edge_map = {normalize_id(e['source_reach_id']):normalize_id(e['new_downstream_id']) for e in edges}
    if len(edge_map)!=len(edges):
        block('DUPLICATE_NEXTDOWN_SOURCE')
    for e in edges:
        if e['qa_flags']:
            if normalize_id(e.get('source_original_parent','')) in original_ids:
                reviews.append('STEP04C_EDGE_REVIEW:'+normalize_id(e['source_reach_id']))
    site_cut_file = diag/'step04'/f'{site}_tdx_shared_gauge_topology.json'
    if not site_cut_file.is_file():
        block('MISSING_SITE_SHARED_GAUGE_RELATIONSHIP_JSON')
    else:
        meta = json.loads(site_cut_file.read_text())
        if not meta.get('shared_relationships_analyzed'):
            block('SITE_SHARED_RELATIONSHIPS_NOT_ANALYZED')
    if blockers:
        dest_summary.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result,indent=2))
        return 2
    original_outlet = normalize_id(json.loads((diag/'tdx_outlet_selection.json').read_text())['reach_id'])
    present_cut_parents = set(pieces_by_parent)&original_ids
    for p in sorted(present_cut_parents):
        qa = geometry_qa.get(p)
        if qa is None or qa['geometry_status']!='PASS':
            block('INVALID_PARENT_SPLIT_GEOMETRY:'+p)
        elif qa['planning_qa_flags']:
            reviews.append('STEP04B_PLANNING_REVIEW:'+p)
    if blockers:
        dest_summary.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result,indent=2))
        return 2
    network, outlet_id, _, network_qa = build_reach_graph(
        native, split_by_id, edge_map, pieces_by_parent, original_outlet, site, cuts)
    result['network_connectivity_qa'] = network_qa
    # Validate original reach paths: every planned parent piece is represented by the
    # candidate graph and exact native reach IDs remain for unchanged reaches.
    used_parents = {r['parent_id'] for r in network.values() if r['source_kind']=='SPLIT_PIECE'}
    result['original_outlet_parent_tdx_id'] = original_outlet
    result['proposed_outlet_piece_id'] = outlet_id
    result['planned_parents_in_basin'] = sorted(present_cut_parents)
    result['split_parents_in_selected_network'] = sorted(used_parents)

    cut_by_parent = defaultdict(list)
    for cut in cuts:
        p = normalize_id(cut['parent_tdx_reach_id'])
        if p in used_parents:
            cut_by_parent[p].append(cut)
    for group in cut_by_parent.values():
        group.sort(key=lambda r:float(r['fraction_from_line_start']))

    required_sites = sorted({r['site_id'] for group in cut_by_parent.values() for r in group})
    boundaries = {}
    watershed_metadata = {}
    for other in required_sites:
        try:
            boundary, meta = load_site_boundary(output_root, other, metric)
            boundaries[other] = boundary
            watershed_metadata[other] = meta
            if meta.get('qa_status') != 'PROVISIONAL':
                reviews.append('UPSTREAM_DEM_REVIEW:'+other)
        except (FileNotFoundError, ValueError) as exc:
            block('UPSTREAM_STEP05_MISSING_OR_BAD:'+other+':'+str(exc))
    current_boundary, own_meta = load_site_boundary(output_root,site,metric)
    if own_meta.get('qa_status') != 'PROVISIONAL':
        reviews.append('TARGET_DEM_REVIEW:'+site)
    result['required_upstream_gauge_watersheds'] = required_sites
    result['upstream_watershed_metadata'] = {s:{'qa_status':watershed_metadata[s].get('qa_status'),
                                               'area_difference_pct':watershed_metadata[s].get('refined_area_difference_pct')}
                                            for s in sorted(watershed_metadata)}
    if blockers:
        dest_summary.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result,indent=2))
        return 2
    mask_percentages={}
    for station,shape in boundaries.items():
        outside = float(shape.difference(current_boundary).area)
        percent = 100*outside/max(shape.area,1.e-10)
        mask_percentages[station] = round(percent,4)
        if percent>max_mask_mismatch_pct:
            reviews.append(f'UPSTREAM_DEM_OUTSIDE_TARGET_{station}:{percent:.1f}%')
    result['upstream_watershed_outside_target_pct'] = mask_percentages

    hru_geoms = {}
    partition_qa = {}
    split_parents = {rid:row['parent_id'] for rid,row in network.items() if row['source_kind']=='SPLIT_PIECE'}
    for parent,orig_polygon in parent_polys.items():
        if parent not in pieces_by_parent:
            if parent in network:
                hru_geoms[parent] = orig_polygon
            continue
        if parent not in used_parents:
            continue
        piece_ids = pieces_by_parent[parent]
        parent_cuts = cut_by_parent[parent]
        orientation = split_by_id[piece_ids[0]]['orientation']
        candidate_polys,qa = piece_partition(orig_polygon,piece_ids,parent_cuts,boundaries,orientation)
        partition_qa[parent] = {k:round(v,8) for k,v in qa.items()}
        for pid,geometry in candidate_polys.items():
            if pid in network:
                if geometry.is_empty or geometry.area<=0:
                    block('EMPTY_REQUIRED_SPLIT_HRU:'+pid)
                else:
                    hru_geoms[pid] = geometry
        if qa['non_nested_mask_area_km2']>max(0.001,qa['parent_area_km2']*0.001):
            reviews.append('INCONSISTENT_UPSTREAM_MASKS:'+parent)
        if qa['coverage_gap_km2']>max(0.001,qa['parent_area_km2']*0.0001):
            reviews.append('PARENT_PARTITION_COVERAGE_GAP:'+parent)
    result['per_split_parent_partition_qa'] = partition_qa
    if set(hru_geoms)!=set(network):
        block('REACH_HRU_IDS_DIFFER:' + ','.join(sorted(set(network)^set(hru_geoms))[:10]))
    if blockers:
        dest_summary.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result,indent=2))
        return 2

    # Partition must cover the already-refined current basin, without overlaps.
    union = poly(unary_union(list(hru_geoms.values())))
    gap_km2 = float(current_boundary.symmetric_difference(union).area)/1e6
    overlap_km2 = max(0.,(sum(g.area for g in hru_geoms.values())-union.area)/1e6)
    result['coverage_gap_km2'] = round(gap_km2,8)
    result['overlap_km2'] = round(overlap_km2,8)
    result['hru_count'] = len(hru_geoms)
    result['reach_count'] = len(network)
    result['split_reach_count'] = len(split_parents)
    result['original_native_reach_count'] = len(original_ids)
    if gap_km2 > max(0.002,current_boundary.area/1e6*0.001):
        reviews.append(f'NESTED_HRU_BOUNDARY_MISMATCH:{gap_km2:.4f}_km2')
    if overlap_km2>max(0.002,current_boundary.area/1e6*0.001):
        reviews.append(f'NESTED_HRU_OVERLAP:{overlap_km2:.4f}_km2')

    # Check split reach-line segment is spatially close to its own HRU (source
    # TDX line can deviate from DEM drainage path; review, never silently move).
    spatial_rows=[]
    for rid,row in network.items():
        geom=row['geometry']
        if geom.is_empty or geom.length<=0:
            block('INVALID_REACH_GEOMETRY:'+rid)
            continue
        buffer = hru_geoms[rid].buffer(30.)
        length_inside = float(geom.intersection(buffer).length)
        ratio = length_inside/max(float(geom.length),1.e-8)
        spatial_rows.append({'reach_id':rid,'parent_reach_id':row['parent_id'],
                             'source_kind':row['source_kind'],'nextdown_id':row['downstream_id'],
                             'own_hru_line_fraction_within_30m':round(ratio,6),
                             'hru_area_km2':round(float(hru_geoms[rid].area)/1e6,8)})
        if row['source_kind']=='SPLIT_PIECE' and ratio<0.5:
            reviews.append('SPLIT_LINE_HRU_GEOMETRY_MISMATCH:'+rid)
    result['minimum_split_reach_own_HRU_overlap_fraction'] = min(
        (r['own_hru_line_fraction_within_30m'] for r in spatial_rows if r['source_kind']=='SPLIT_PIECE'),
        default=None)
    if blockers:
        dest_summary.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result,indent=2))
        return 2
    result['review_flags']=sorted(set(reviews))
    result['status']='REVIEW' if reviews else 'PROVISIONAL'
    preview = root/'tdx_nested_partition_preview.gpkg'
    result['candidate_output_gpkg']=str(preview)
    if preflight_only:
        result['output_written']=False
    elif reviews and not allow_review:
        result['output_written']=False
        result['note'] += ' Use --allow-review to write REVIEW preview.'
    else:
        if preview.exists():
            if not overwrite:
                raise FileExistsError(f'Existing preview {preview}; use --overwrite')
            preview.unlink()
        hru_rows=[]
        for rid in sorted(hru_geoms):
            hru_rows.append({'reach_id':rid,'parent_reach_id':split_parents.get(rid,rid),
                             'partition_source':'DEM_GAUGE_MASK' if rid in split_parents else 'STEP05_PARENT',
                             'area_km2':float(hru_geoms[rid].area)/1e6, 'geometry':hru_geoms[rid]})
        reach_rows=[]
        for rid,row in sorted(network.items()):
            reach_rows.append({'reach_id':rid,'parent_reach_id':row['parent_id'],
                               'downstream_id':row['downstream_id'],
                               'source_kind':row['source_kind'],
                               'slope':row.get('slope'),
                               'length_m':float(row['geometry'].length),'geometry':row['geometry']})
        gpd.GeoDataFrame(hru_rows,geometry='geometry',crs=metric).to_crs('EPSG:4326').to_file(
            preview,layer='nested_refined_subbasins',driver='GPKG')
        gpd.GeoDataFrame(reach_rows,geometry='geometry',crs=metric).to_crs('EPSG:4326').to_file(
            preview,layer='nested_reaches',driver='GPKG',mode='a')
        gpd.GeoDataFrame([{'site_id':site,'geometry':current_boundary}],crs=metric).to_crs('EPSG:4326').to_file(
            preview,layer='refinement_boundary',driver='GPKG',mode='a')
        rows_file=diag/'step05b'/f'{site}_nested_spatial_qa.csv'
        with rows_file.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(spatial_rows[0]))
            writer.writeheader()
            writer.writerows(spatial_rows)
        result['output_written']=True
        result['spatial_qa_csv']=str(rows_file)
    dest_summary.write_text(json.dumps(result,indent=2)+'\n')
    (diag/'tdx_nested_partition_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    return 0


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--control',type=Path,required=True)
    parser.add_argument('--preflight-only',action='store_true')
    parser.add_argument('--allow-review',action='store_true',help='Write REVIEW diagnostic outputs; never production exports')
    parser.add_argument('--overwrite',action='store_true')
    parser.add_argument('--max-mask-mismatch-pct',type=float,default=5.)
    args=parser.parse_args()
    if args.max_mask_mismatch_pct<0 or args.max_mask_mismatch_pct>100:
        parser.error('Mask mismatch tolerance must lie between 0 and 100 percent')
    project=Path(__file__).resolve().parents[2]
    return prepare(args.control.resolve(), project, args.allow_review,
                   args.preflight_only,args.overwrite,args.max_mask_mismatch_pct)

if __name__=='__main__':
    try:
        sys.exit(main())
    except Exception as e:
        print(f'ERROR: {type(e).__name__}: {e}',file=sys.stderr)
        sys.exit(1)
