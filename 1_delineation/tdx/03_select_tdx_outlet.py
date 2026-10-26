#!/usr/bin/env python3
"""Step 03: Select and snap a gauge to a TDX-Hydro stream reach.

Read-only TDX access. Standalone TOMLs retain the original nearest-reach
selection. Batch TOMLs with [gauge_match] select the preassigned original TDX
reach (without splitting it), and verify it occurs among nearby candidates.
Failed global snaps are never used; REVIEW is diagnostic-only with an explicit
--allow-review flag. Neither shared nested-reach topology nor DEM headwater
extension is implemented by this script.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import tomllib
from pathlib import Path

import fiona
from pyproj import CRS, Transformer
from shapely.geometry import Point, box, shape
from shapely.ops import transform as transform_geom
from shapely.geometry import mapping
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def _number(value):
    try:
        if value is None or str(value).strip() == '':
            return None
        n = float(value)
        return n if math.isfinite(n) else None
    except (TypeError, ValueError):
        return None


def _identity(value):
    if value is None:
        return None
    s = str(value).strip()
    return s[:-2] if s.endswith('.0') and s[:-2].isdigit() else s


def _write_csv(path, records, fields):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(records)


def _metric_crs(lat, lon):
    if not (-80 <= lat <= 84):
        raise ValueError('Automatic local UTM projection is unavailable outside 80 S to 84 N; supply [delineation] projected_crs.')
    zone = min(60, max(1, int((lon + 180) // 6) + 1))
    return CRS.from_epsg((32600 if lat >= 0 else 32700) + zone)


def _area_km2(value, units):
    if value is None:
        return None
    factors = {'km2': 1.0, 'm2': 1e-6, 'ha': 0.01}
    return value * factors[units]


def _orientation(stream_file, stream_layer, selected, metric, reach_field, downstream_field):
    """Infer line digitization using the NEXT downstream reach's endpoint.

    Returns unknown unless a unique, spatially connected downstream line exists.
    """
    ds_id = selected.get('downstream_reach_id')
    if not ds_id or ds_id in ('0', '-1'):
        return 'unknown', 'No downstream reach identifier'
    geom = selected['_geom']
    if geom.geom_type != 'LineString':
        return 'unknown', 'Multipart selected reach'
    endpoints = [Point(geom.coords[0]), Point(geom.coords[-1])]
    native_to_metric = None
    with fiona.open(stream_file, layer=stream_layer) as src:
        native = CRS.from_user_input(src.crs_wkt or src.crs)
        native_to_metric = Transformer.from_crs(native, metric, always_xy=True).transform
        to_native = Transformer.from_crs(metric, native, always_xy=True).transform
        # Search a small neighborhood of the selected line's endpoints only.
        bounds = transform_geom(to_native, geom.envelope.buffer(100)).bounds
        matches = []
        for feat in src.filter(bbox=bounds):
            if _identity(feat['properties'].get(reach_field)) == ds_id and feat['geometry']:
                matches.append(transform_geom(native_to_metric, shape(feat['geometry'])))
    if len(matches) != 1:
        return 'unknown', f'Found {len(matches)} downstream reach geometries nearby'
    downstream = matches[0]
    d0, d1 = [downstream.distance(pt) for pt in endpoints]
    # Only accept a clear connection at one end, within 30 m.
    if min(d0, d1) > 30 or abs(d0 - d1) < 30:
        return 'unknown', f'Downstream connection ambiguous: {d0:.1f}, {d1:.1f} m'
    return ('upstream_to_downstream' if d1 < d0 else 'downstream_to_upstream'), 'Verified against downstream reach geometry'


def _make_preview(path, *, stream_file, stream_layer, basins_file, basins_layer,
                  lat, lon, metric, radius_m, selected, rows, site_id):
    """Local hydrofabric preview only; NOT a delineated upstream watershed."""
    to_metric = Transformer.from_crs('EPSG:4326', metric, always_xy=True).transform
    center = transform_geom(to_metric, Point(lon, lat))
    # Wider than snapping radius to show hydrological context.
    half = max(1500., radius_m * 4)
    window = box(center.x-half, center.y-half, center.x+half, center.y+half)
    fig, ax = plt.subplots(figsize=(10, 9))
    counts = {'basins': 0, 'streams': 0}
    if basins_file and Path(basins_file).is_file() and basins_layer:
        try:
            with fiona.open(basins_file, layer=basins_layer) as src:
                native = CRS.from_user_input(src.crs_wkt or src.crs)
                to_native = Transformer.from_crs(metric, native, always_xy=True).transform
                to_m = Transformer.from_crs(native, metric, always_xy=True).transform
                for feat in src.filter(bbox=transform_geom(to_native, window).bounds):
                    if not feat['geometry']:
                        continue
                    geom = transform_geom(to_m, shape(feat['geometry']))
                    if not geom.intersects(window):
                        continue
                    geom = geom.intersection(window)
                    polys = list(geom.geoms) if geom.geom_type == 'MultiPolygon' else ([geom] if geom.geom_type == 'Polygon' else [])
                    for poly in polys:
                        x,y = poly.exterior.xy
                        ax.fill(x,y,facecolor='#c7dfbd',edgecolor='#8aaa83',linewidth=.45,alpha=.5,zorder=1)
                        counts['basins'] += 1
        except Exception as exc:
            print(f'WARNING: Could not draw basin polygons: {exc}', file=sys.stderr)
    with fiona.open(stream_file, layer=stream_layer) as src:
        native = CRS.from_user_input(src.crs_wkt or src.crs)
        to_native = Transformer.from_crs(metric, native, always_xy=True).transform
        to_m = Transformer.from_crs(native, metric, always_xy=True).transform
        for feat in src.filter(bbox=transform_geom(to_native, window).bounds):
            if not feat['geometry']:
                continue
            geom = transform_geom(to_m, shape(feat['geometry']))
            if not geom.intersects(window):
                continue
            lines = list(geom.geoms) if geom.geom_type == 'MultiLineString' else ([geom] if geom.geom_type == 'LineString' else [])
            for line in lines:
                clip = line.intersection(window)
                pieces = list(clip.geoms) if hasattr(clip,'geoms') else [clip]
                for piece in pieces:
                    if piece.geom_type == 'LineString' and not piece.is_empty:
                        x,y=piece.xy
                        is_selected = _identity(feat['properties'].get(selected['_reach_field'])) == selected['reach_id']
                        ax.plot(x,y,color='#d14b35' if is_selected else '#3177a3',linewidth=2.8 if is_selected else 1.05,zorder=3 if is_selected else 2)
                        counts['streams'] += 1
    snapped = Point(selected['_snap_x'], selected['_snap_y'])
    ax.scatter([center.x],[center.y],marker='x',s=130,color='black',linewidth=2.3,zorder=6,label='Original gauge (TOML)')
    ax.scatter([snapped.x],[snapped.y],marker='o',s=95,facecolor='#ffbf45',edgecolor='black',zorder=7,label='Snapped gauge')
    ax.plot([center.x,snapped.x],[center.y,snapped.y],linestyle=':',color='black',linewidth=1.5,zorder=5)
    ax.annotate(f"LINKNO {selected['reach_id']}\nSnap {selected['distance_m']:.1f} m",(snapped.x,snapped.y),xytext=(14,16),textcoords='offset points',fontsize=9,bbox={'facecolor':'white','alpha':.8,'edgecolor':'none'})
    ax.set_xlim(center.x-half,center.x+half)
    ax.set_ylim(center.y-half,center.y+half)
    ax.set_aspect('equal')
    ax.set_xlabel(f'Easting (m; {metric.to_string()})')
    ax.set_ylabel('Northing (m)')
    ax.set_title(f'{site_id} | TDX outlet snapping preview\nLocal TDX subbasins and rivers — NOT the delineated catchment/HRUs')
    handles,labels=ax.get_legend_handles_labels()
    handles=[Patch(facecolor='#c7dfbd',edgecolor='#8aaa83',label='TDX local subbasin polygons'),
             plt.Line2D([0],[0],color='#3177a3',label='TDX river network'),
             plt.Line2D([0],[0],color='#d14b35',linewidth=3,label='Selected outlet reach')]+handles
    ax.legend(handles=handles,loc='lower left',fontsize=8,framealpha=.9)
    ax.grid(alpha=.15)
    fig.text(.5,.012,'Preview extent is local only. Full upstream catchment and HRUs will be generated in Step 04.',ha='center',fontsize=9)
    fig.savefig(path,dpi=220,bbox_inches='tight')
    plt.close(fig)
    return counts


def select_tdx_outlet(*, region, lat, lon, site_id, diagnostics_dir,
                      radius_m=500.0, max_candidates=10, observed_area_km2=None,
                      reach_id_field='LINKNO', downstream_id_field='DSLINKNO',
                      upstream_area_field='USContArea', downstream_area_field='DSContArea',
                      order_field='strmOrder', area_units='km2', manual_reach_id=None,
                      projected_crs=None, area_warning_pct=20.0,
                      preferred_reach_id=None):
    """Select a TDX reach, estimate gauge area, and save diagnostic preview."""
    if not (-90 <= lat <= 90 and -180 <= lon <= 180):
        raise ValueError('Invalid gauge latitude/longitude')
    if radius_m <= 0 or max_candidates < 1:
        raise ValueError('search_radius_m must be positive and max_candidates >= 1')
    if area_units not in ('km2','m2','ha'):
        raise ValueError('upstream_area_units must be km2, m2 or ha')
    diagnostics = Path(diagnostics_dir)
    diagnostics.mkdir(parents=True,exist_ok=True)
    stream_file = Path(region['streamnet_file'])
    stream_layer = region['streamnet_layer']
    if not stream_file.is_file():
        raise FileNotFoundError(stream_file)
    metric = CRS.from_user_input(projected_crs) if projected_crs else _metric_crs(lat,lon)
    if not metric.is_projected or metric.axis_info[0].unit_name.lower() not in ('metre','meter','metres','meters'):
        raise ValueError('projected_crs must use metres')
    to_metric = Transformer.from_crs('EPSG:4326',metric,always_xy=True).transform
    to_geo = Transformer.from_crs(metric,'EPSG:4326',always_xy=True).transform
    gauge = transform_geom(to_metric,Point(lon,lat))
    search_box = gauge.buffer(radius_m).envelope
    rows=[]
    with fiona.open(stream_file,layer=stream_layer) as src:
        native=CRS.from_user_input(src.crs_wkt or src.crs)
        to_native=Transformer.from_crs(metric,native,always_xy=True).transform
        to_m=Transformer.from_crs(native,metric,always_xy=True).transform
        bbox=transform_geom(to_native,search_box).bounds
        for feature in src.filter(bbox=bbox):
            if not feature['geometry']:
                continue
            props=dict(feature['properties'])
            rid=_identity(props.get(reach_id_field))
            if not rid:
                continue
            geom=transform_geom(to_m,shape(feature['geometry']))
            if geom.is_empty or geom.geom_type not in ('LineString','MultiLineString'):
                continue
            distance=gauge.distance(geom)
            if distance>radius_m:
                continue
            parts=list(geom.geoms) if geom.geom_type=='MultiLineString' else [geom]
            part=min(parts,key=lambda g:g.distance(gauge))
            along=part.project(gauge)
            snapped=part.interpolate(along)
            snapgeo=transform_geom(to_geo,snapped)
            up=_number(props.get(upstream_area_field))
            down=_number(props.get(downstream_area_field))
            row={'reach_id':rid,'distance_m':round(distance,3),'snapped_lon':round(snapgeo.x,9),
                 'snapped_lat':round(snapgeo.y,9),'stream_order':_number(props.get(order_field)),
                 'upstream_area_raw':up,'downstream_area_raw':down,
                 'upstream_area_km2':_area_km2(up,area_units),
                 'downstream_area_km2':_area_km2(down,area_units),
                 'along_reach_m':round(along,3),'reach_length_m':round(part.length,3),
                 'downstream_reach_id':_identity(props.get(downstream_id_field)),
                 '_geom':part,'_snap_x':snapped.x,'_snap_y':snapped.y,'_reach_field':reach_id_field}
            if any(r['reach_id']==rid for r in rows):
                raise ValueError(f'Duplicate reach ID near gauge: {rid}')
            rows.append(row)
    rows.sort(key=lambda r:(r['distance_m'],r['reach_id']))
    fields=['rank','reach_id','distance_m','snapped_lon','snapped_lat','stream_order',
            'upstream_area_raw','downstream_area_raw','upstream_area_km2','downstream_area_km2',
            'orientation','orientation_evidence','area_estimate_km2','area_difference_pct',
            'along_reach_m','reach_length_m','selected']
    if not rows:
        _write_csv(diagnostics/'outlet_candidates.csv',[],fields)
        raise RuntimeError(f'No TDX stream within {radius_m:g} m')
    shortlist=rows[:max_candidates]
    if preferred_reach_id is not None:
        if manual_reach_id is not None and str(manual_reach_id) not in ('', '0', 'None'):
            if _identity(manual_reach_id) != _identity(preferred_reach_id):
                raise ValueError('Manual reach and global gauge-match reach disagree')
        matches = [r for r in rows if r['reach_id'] == _identity(preferred_reach_id)]
        if len(matches) != 1:
            _write_csv(diagnostics/'outlet_candidates.csv', [], fields)
            raise ValueError(
                f'Global assigned TDX reach {preferred_reach_id} not found within '
                f'{radius_m:g} m for {site_id}; do not silently choose another reach'
            )
        chosen = matches[0]
        mode = 'global_gauge_match_provisional'
        if chosen not in shortlist:
            shortlist.append(chosen)
    elif manual_reach_id is not None and str(manual_reach_id) not in ('','0','None'):
        matches=[r for r in rows if r['reach_id']==_identity(manual_reach_id)]
        if not matches:
            raise ValueError(f'Manual reach {manual_reach_id} not found within {radius_m:g} m')
        chosen=matches[0]
        mode='manual'
        if chosen not in shortlist:
            shortlist.append(chosen)
    else:
        chosen=rows[0]
        mode='automatic_nearest_provisional'
    # Orientation verification is done only for shortlist to keep HPC runs efficient.
    for r in shortlist:
        orientation,evidence=_orientation(stream_file,stream_layer,r,metric,reach_id_field,downstream_id_field)
        r['orientation']=orientation
        r['orientation_evidence']=evidence
        up,down=r['upstream_area_km2'],r['downstream_area_km2']
        if up is not None and down is not None and r['reach_length_m']>0 and orientation!='unknown':
            frac=r['along_reach_m']/r['reach_length_m']
            frac=frac if orientation=='upstream_to_downstream' else 1-frac
            r['area_estimate_km2']=up+frac*(down-up)
        else:
            r['area_estimate_km2']=None
        r['area_difference_pct']=(100*abs(r['area_estimate_km2']-observed_area_km2)/observed_area_km2
            if r['area_estimate_km2'] is not None and observed_area_km2 and observed_area_km2>0 else None)
    for i,r in enumerate(shortlist,1):
        r['rank']=i
        r['selected']=r is chosen
    _write_csv(diagnostics/'outlet_candidates.csv',shortlist,fields)
    warnings=[]
    if chosen['orientation']=='unknown':
        warnings.append('River digitization orientation not verified; gauge-area interpolation withheld')
    if chosen['area_difference_pct'] is not None and chosen['area_difference_pct']>area_warning_pct:
        warnings.append(f'Drainage-area discrepancy {chosen["area_difference_pct"]:.1f}% exceeds {area_warning_pct:g}%')
    if len(rows)>1 and rows[1]['distance_m']-rows[0]['distance_m']<25 and mode!='manual':
        warnings.append('Multiple nearby reaches: review topology')
    if chosen['area_estimate_km2'] is None:
        warnings.append('No verified interpolated drainage-area estimate available')
    preview=diagnostics/'tdx_outlet_preview.jpg'
    try:
        counts=_make_preview(preview,stream_file=stream_file,stream_layer=stream_layer,
                             basins_file=region.get('basins_file'),basins_layer=region.get('basins_layer'),
                             lat=lat,lon=lon,metric=metric,radius_m=radius_m,
                             selected=chosen,rows=shortlist,site_id=site_id)
    except Exception as exc:
        warnings.append(f'Preview generation failed: {type(exc).__name__}: {exc}')
        counts={}
    result={'site_id':site_id,'tdx_id':region['tdx_id'],'selection_method':mode,
            'review_required':bool(warnings) or mode!='manual',
            'reach_id':chosen['reach_id'],'gauge_lat':lat,'gauge_lon':lon,
            'snapped_lat':chosen['snapped_lat'],'snapped_lon':chosen['snapped_lon'],
            'snap_distance_m':chosen['distance_m'],'metric_crs':metric.to_string(),
            'observed_area_km2':observed_area_km2,'tdx_area_units_assumed':area_units,
            'tdx_upstream_area_km2':chosen['upstream_area_km2'],
            'tdx_downstream_area_km2':chosen['downstream_area_km2'],
            'line_orientation':chosen['orientation'],
            'orientation_evidence':chosen['orientation_evidence'],
            'tdx_area_estimate_km2':chosen['area_estimate_km2'],
            'area_estimation_method':'linear_interpolation_diagnostic_only' if chosen['area_estimate_km2'] is not None else 'unavailable',
            'area_difference_pct':chosen['area_difference_pct'],
            'stream_order':chosen['stream_order'],'candidate_count_within_radius':len(rows),
            'warnings':warnings,'streamnet_file':str(stream_file),'streamnet_layer':stream_layer,
            'preview_jpg':str(preview) if preview.is_file() else None,
            'preview_features_drawn':counts}
    (diagnostics/'tdx_outlet_selection.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--control', type=Path, required=True)
    parser.add_argument('--region-json', type=Path, help='Region-selection JSON (defaults to per-site diagnostics)')
    parser.add_argument('--diagnostics-dir', type=Path)
    parser.add_argument('--allow-review', action='store_true',
                        help='Permit REVIEW global matches for diagnostic-only outlet selection; not export approval')
    args = parser.parse_args()
    with args.control.open('rb') as f:
        config = tomllib.load(f)
    domain, outlet, inp = (config.get(k, {}) for k in ('domain', 'outlet', 'input'))
    output, delineation, fields = (config.get(k, {}) for k in ('output', 'delineation', 'tdx_fields'))
    site_id = domain['site_id']
    gauge_match = config.get('gauge_match')
    for section in ('nested_gauges', 'headwater'):
        if config.get(section, {}).get('enabled', False):
            parser.error(f'[{section}] enabled=true is not implemented in Step 03; disable until the shared workflow is integrated')
    preferred_reach_id = None
    match_status = None
    if gauge_match is not None:
        match_status = str(gauge_match.get('snap_status', '')).strip().upper()
        if match_status == 'FAILED':
            parser.error(f'{site_id}: global gauge snapping FAILED; Step 03 cannot select an outlet')
        if match_status not in ('PROVISIONAL', 'REVIEW'):
            parser.error(f'{site_id}: invalid global match status {match_status!r}')
        if match_status == 'REVIEW' and not args.allow_review:
            parser.error(f'{site_id}: global match requires review; use --allow-review for DIAGNOSTIC testing only')
        preferred_reach_id = _identity(gauge_match.get('parent_tdx_reach_id'))
        region_id = _identity(gauge_match.get('tdx_region'))
        input_region = _identity(inp.get('tdx_id'))
        if not preferred_reach_id or preferred_reach_id in ('0', '-1'):
            parser.error(f'{site_id}: missing or invalid parent TDX reach ID')
        if not region_id or region_id != input_region:
            parser.error(f'{site_id}: gauge_match region and [input] tdx_id disagree')
    diagnostics = args.diagnostics_dir or (Path(output['output_directory']) / site_id /
                                             output.get('diagnostics_subdirectory', 'diagnostics'))
    region_json = args.region_json or diagnostics / 'tdx_region_selection.json'
    if not region_json.is_file():
        parser.error(f'Missing {region_json}. Run region-selection step first.')
    region = json.loads(region_json.read_text())
    if abs(float(region.get('lat')) - float(domain.get('lat'))) > 1e-9 or abs(float(region.get('lon')) - float(domain.get('lon'))) > 1e-9:
        parser.error('Gauge coordinates differ from saved region selection; rerun region-selection step')
    if gauge_match is not None:
        if _identity(region.get('tdx_id')) != _identity(gauge_match.get('tdx_region')):
            parser.error(f'{site_id}: Step 02 region differs from the global gauge-match region')
        for cfg_field, region_field in (
            ('tdx_streams_file', 'streamnet_file'), ('tdx_basins_file', 'basins_file'),
        ):
            wanted = str(inp.get(cfg_field, '')).strip()
            actual = Path(region.get(region_field, '')).name
            if not wanted or wanted != actual:
                parser.error(f'{site_id}: Step 02 {region_field} does not match TOML {cfg_field}')
        for cfg_field, region_field in (
            ('tdx_streams_layer', 'streamnet_layer'), ('tdx_basins_layer', 'basins_layer'),
        ):
            wanted = str(inp.get(cfg_field, '')).strip()
            if wanted and wanted != str(region.get(region_field, '')).strip():
                parser.error(f'{site_id}: Step 02 {region_field} differs from TOML {cfg_field}')
    try:
        result = select_tdx_outlet(
            region=region, lat=float(domain['lat']), lon=float(domain['lon']), site_id=site_id,
            diagnostics_dir=diagnostics, radius_m=float(outlet.get('search_radius_m', 500)),
            max_candidates=int(outlet.get('max_candidates', 10)),
            observed_area_km2=_number(domain.get('observed_area_km2')),
            reach_id_field=fields.get('reach_id', 'LINKNO'),
            upstream_area_field=fields.get('upstream_area', 'USContArea'),
            downstream_area_field=fields.get('downstream_area', 'DSContArea'),
            order_field=fields.get('stream_order', 'strmOrder'),
            downstream_id_field=fields.get('downstream_id', 'DSLINKNO'),
            area_units=outlet.get('upstream_area_units', 'km2'),
            manual_reach_id=outlet.get('manual_reach_id') if outlet.get('selection_method', 'automatic') == 'manual' else None,
            projected_crs=delineation.get('projected_crs'),
            area_warning_pct=float(outlet.get('area_difference_warning_pct', 20)),
            preferred_reach_id=preferred_reach_id)
    except Exception as exc:
        print(f'ERROR: {type(exc).__name__}: {exc}', file=sys.stderr)
        return 1
    if gauge_match is not None:
        result['global_snap_status'] = match_status
        result['global_parent_tdx_reach_id'] = preferred_reach_id
        result['shared_group_status'] = gauge_match.get('shared_group_status')
        result['shared_topology_applied'] = False
        result['headwater_extension_applied'] = False
        result['final_export_authorized'] = False
        if match_status == 'REVIEW':
            result['warnings'].append('Global gauge-match REVIEW: diagnostic only, not export-authorized')
        if gauge_match.get('requires_review', False):
            result['warnings'].append('Control manifest flags unresolved shared-network and/or source-data checks')
        result['review_required'] = True
        # Update the JSON with provenance used by downstream steps, but do not
        # manufacture any split IDs or alter original TDX source data.
        (Path(diagnostics) / 'tdx_outlet_selection.json').write_text(
            json.dumps(result, indent=2) + '\n'
        )
    print(json.dumps(result, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
