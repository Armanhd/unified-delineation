#!/usr/bin/env python3
"""Step 02b: snap a gauge catalogue to original TDX-Hydro reaches.

Read-only TDX access. Spatially filter each regional GeoPackage; do not load
whole river layers. Snaps are PROVISIONAL: this script does NOT split reaches,
assign fractional IDs, trace upstream networks, or generate basin TOMLs.

Examples:
  python 1_delineation/tdx/02b_snap_gauge_catalog.py --site USA_10343500
  python 1_delineation/tdx/02b_snap_gauge_catalog.py --start 0 --stop 20
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path

import fiona
from pyproj import CRS, Transformer
from shapely.geometry import Point, box, shape
from shapely.ops import transform

PROJECT = Path(__file__).resolve().parents[2]
DEFAULT_CATALOG = PROJECT / "indexes/gauge_catalogs/camels-spat_gauge_catalog_validated.csv"
DEFAULT_INDEX = PROJECT / "indexes/tdx_regions.csv"
DEFAULT_TDX = Path("/work/comphyd_lab/data/geospatial-data/TDX_hydrofabrics/TDX_hydrofabrics")
DEFAULT_OUTPUT = PROJECT / "indexes/gauge_catalogs/tdx_snapping"

SNAP_COLUMNS = [
    "site_id", "latitude", "longitude", "expected_area_km2", "tdx_id",
    "reach_id", "downstream_reach_id", "snapped_lat", "snapped_lon",
    "snap_distance_m", "along_reach_m", "reach_length_m", "end_distance_m",
    "tdx_upstream_area_km2", "tdx_downstream_area_km2", "stream_order",
    "n_candidates", "n_regions_with_candidates", "status", "warnings",
]
CANDIDATE_COLUMNS = [
    "site_id", "rank", "tdx_id", "reach_id", "downstream_reach_id",
    "distance_m", "snapped_lat", "snapped_lon", "along_reach_m",
    "reach_length_m", "end_distance_m", "upstream_area_km2",
    "downstream_area_km2", "stream_order", "area_outside_range_pct",
]


def write_csv(path, rows, columns):
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def ident(value):
    if value is None:
        return ""
    val = str(value).strip()
    if val.endswith(".0") and val[:-2].isdigit():
        val = val[:-2]
    return val


def local_crs(lat, lon):
    if -80 <= lat <= 84:
        zone = min(60, max(1, int((lon + 180) // 6) + 1))
        return CRS.from_epsg((32600 if lat >= 0 else 32700) + zone)
    # Polar stereographic at latitudes beyond UTM coverage.
    return CRS.from_epsg(3413 if lat > 0 else 3031)


def index_bounds(row):
    return tuple(float(row[k]) for k in ("min_x", "min_y", "max_x", "max_y"))


def index_contains(row, lat, lon, radius_m):
    """Conservative bbox prefilter; actual reach snap is spatially exact."""
    try:
        crs = CRS.from_epsg(int(row.get("bbox_srs_id") or 4326))
        if crs.is_geographic:
            # Larger-than-needed padding to allow near-boundary snap candidates.
            lat_pad = radius_m / 100000.0
            lon_pad = radius_m / max(1000.0, 100000.0 * abs(math.cos(math.radians(lat))))
            x, y = lon, lat
            xpad, ypad = lon_pad, lat_pad
        else:
            x, y = Transformer.from_crs(4326, crs, always_xy=True).transform(lon, lat)
            xpad = ypad = radius_m * 1.5
        xmin, ymin, xmax, ymax = index_bounds(row)
        return xmin - xpad <= x <= xmax + xpad and ymin - ypad <= y <= ymax + ypad
    except Exception:
        return False


def stream_layers(row):
    layers = json.loads(row.get("streamnet_layers_json") or "[]")
    return [entry["name"] for entry in layers if entry.get("name")]


def region_candidates(row, gauge, metric, tdx_root, radius_m, area_units):
    """Read only stream geometries near the gauge, from an indexed TDX GPKG."""
    file = tdx_root / row["streamnet_filename"]
    if not file.is_file():
        raise FileNotFoundError(file)
    layer_names = stream_layers(row)
    if not layer_names:
        raise ValueError(f"No indexed stream layer for region {row['tdx_id']}")

    to_geo = Transformer.from_crs(metric, 4326, always_xy=True).transform
    search = gauge.buffer(radius_m).envelope
    matches = []
    factors = {"km2": 1.0, "m2": 1e-6, "ha": 0.01}

    for layer in layer_names:
        with fiona.open(file, layer=layer) as src:
            if src.crs_wkt or src.crs:
                native_crs = CRS.from_user_input(src.crs_wkt or src.crs)
            else:
                raise ValueError(f"Missing CRS in {file}:{layer}")
            to_native = Transformer.from_crs(metric, native_crs, always_xy=True).transform
            to_metric = Transformer.from_crs(native_crs, metric, always_xy=True).transform
            bbox = transform(to_native, search).bounds
            for feature in src.filter(bbox=bbox):
                if not feature["geometry"]:
                    continue
                prop = feature["properties"]
                rid = ident(prop.get("LINKNO"))
                if not rid:
                    continue
                geom = transform(to_metric, shape(feature["geometry"]))
                if geom.is_empty or geom.geom_type not in ("LineString", "MultiLineString"):
                    continue
                lines = list(geom.geoms) if geom.geom_type == "MultiLineString" else [geom]
                line = min(lines, key=lambda part: part.distance(gauge))
                distance = float(line.distance(gauge))
                if distance > radius_m or line.length <= 0:
                    continue
                along = float(line.project(gauge))
                snapped = line.interpolate(along)
                snap_lon, snap_lat = transform(to_geo, snapped).coords[0]
                up = number(prop.get("USContArea"))
                down = number(prop.get("DSContArea"))
                matches.append({
                    "tdx_id": str(row["tdx_id"]), "reach_id": rid,
                    "downstream_reach_id": ident(prop.get("DSLINKNO")),
                    "distance_m": round(distance, 3),
                    "snapped_lat": round(snap_lat, 9),
                    "snapped_lon": round(snap_lon, 9),
                    "along_reach_m": round(along, 3),
                    "reach_length_m": round(line.length, 3),
                    "end_distance_m": round(min(along, line.length - along), 3),
                    "upstream_area_km2": up * factors[area_units] if up is not None else None,
                    "downstream_area_km2": down * factors[area_units] if down is not None else None,
                    "stream_order": number(prop.get("strmOrder")),
                })
    # A reach should occur once per region. Repeated IDs would make splitting unsafe.
    seen = set()
    for item in matches:
        key = (item["tdx_id"], item["reach_id"])
        if key in seen:
            raise ValueError(f"Duplicate TDX reach ID in candidate region: {key}")
        seen.add(key)
    return matches


def outside_area_range_pct(area, candidate):
    """Area plausibility only; does not assume line digitization direction."""
    a, b = candidate["upstream_area_km2"], candidate["downstream_area_km2"]
    if area is None or area <= 0 or a is None or b is None:
        return None
    delta = max(min(a, b) - area, area - max(a, b), 0)
    return round(100.0 * delta / area, 3)


def snap_one(site, regions, args):
    lat, lon = float(site["latitude"]), float(site["longitude"])
    area = number(site.get("expected_area_km2"))
    gauge_crs = local_crs(lat, lon)
    to_metric = Transformer.from_crs(4326, gauge_crs, always_xy=True).transform
    gauge = transform(to_metric, Point(lon, lat))
    regions = [r for r in regions if index_contains(r, lat, lon, args.radius_m)]
    if args.tdx_id:
        regions = [r for r in regions if str(r["tdx_id"]) == args.tdx_id]
    rows, region_errors = [], []
    for region in regions:
        try:
            rows.extend(region_candidates(region, gauge, gauge_crs,
                                          args.tdx_root, args.radius_m, args.area_units))
        except Exception as exc:
            region_errors.append(f"{region['tdx_id']}: {type(exc).__name__}: {exc}")

    for r in rows:
        r["area_outside_range_pct"] = outside_area_range_pct(area, r)
    rows.sort(key=lambda r: (r["distance_m"], r["tdx_id"], int(r["reach_id"]) if r["reach_id"].isdigit() else r["reach_id"]))
    warnings = []
    if region_errors:
        warnings.append("Regional data errors: " + " | ".join(region_errors[:3]))
    if not regions:
        warnings.append("NO_REGION_COVERAGE")
    elif not rows:
        warnings.append(f"NO_REACH_WITHIN_{args.radius_m:g}_M")
    if rows:
        best = rows[0]
        if best["distance_m"] > args.distance_warning_m:
            warnings.append("LARGE_SNAP_DISTANCE")
        if len(rows) > 1 and rows[1]["distance_m"] - rows[0]["distance_m"] < args.ambiguity_m:
            warnings.append("AMBIGUOUS_NEARBY_REACHES")
        if best["end_distance_m"] < args.endpoint_warning_m:
            warnings.append("GAUGE_NEAR_REACH_ENDPOINT")
        mismatch = best["area_outside_range_pct"]
        if mismatch is not None and mismatch > args.area_warning_pct:
            warnings.append("REFERENCE_AREA_OUTSIDE_TDX_REACH_RANGE")
        n_regions = len({r["tdx_id"] for r in rows})
        if n_regions > 1:
            warnings.append("MULTIPLE_TDX_REGIONS_NEAR_GAUGE")
        data = {
            "site_id": site["site_id"], "latitude": lat, "longitude": lon,
            "expected_area_km2": area if area is not None else "",
            "tdx_id": best["tdx_id"], "reach_id": best["reach_id"],
            "downstream_reach_id": best["downstream_reach_id"],
            "snapped_lat": best["snapped_lat"], "snapped_lon": best["snapped_lon"],
            "snap_distance_m": best["distance_m"],
            "along_reach_m": best["along_reach_m"],
            "reach_length_m": best["reach_length_m"],
            "end_distance_m": best["end_distance_m"],
            "tdx_upstream_area_km2": best["upstream_area_km2"],
            "tdx_downstream_area_km2": best["downstream_area_km2"],
            "stream_order": best["stream_order"],
            "n_candidates": len(rows), "n_regions_with_candidates": n_regions,
            "status": "REVIEW" if warnings else "PROVISIONAL",
            "warnings": "; ".join(warnings),
        }
    else:
        data = dict.fromkeys(SNAP_COLUMNS, "")
        data.update({
            "site_id": site["site_id"], "latitude": lat, "longitude": lon,
            "expected_area_km2": area if area is not None else "",
            "n_candidates": 0, "n_regions_with_candidates": 0,
            "status": "FAILED", "warnings": "; ".join(warnings)
        })
    candidates = [{"site_id": site["site_id"], "rank": i + 1, **r}
                  for i, r in enumerate(rows[:args.max_candidates])]
    return data, candidates


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--catalog", type=Path, default=DEFAULT_CATALOG)
    p.add_argument("--region-index", type=Path, default=DEFAULT_INDEX)
    p.add_argument("--tdx-root", type=Path, default=DEFAULT_TDX)
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--site", help="Process only this gauge (recommended first test)")
    p.add_argument("--start", type=int, default=0, help="First sorted gauge index, inclusive")
    p.add_argument("--stop", type=int, help="Last sorted gauge index, exclusive")
    p.add_argument("--tdx-id", help="Optional restriction to one TDX region")
    p.add_argument("--radius-m", type=float, default=500.0)
    p.add_argument("--max-candidates", type=int, default=10)
    p.add_argument("--area-units", choices=("km2", "m2", "ha"), default="m2")
    p.add_argument("--area-warning-pct", type=float, default=20.0)
    p.add_argument("--distance-warning-m", type=float, default=250.0)
    p.add_argument("--ambiguity-m", type=float, default=25.0)
    p.add_argument("--endpoint-warning-m", type=float, default=30.0)
    args = p.parse_args()
    if args.radius_m <= 0 or args.max_candidates < 1 or args.start < 0:
        p.error("radius-m and max-candidates must be positive; start must be >= 0")
    if args.site is None and args.stop is None:
        p.error("Provide --site for a test or --start/--stop for a batch; do not scan all gauges on a login node")

    with args.catalog.open(newline="", encoding="utf-8-sig") as fh:
        catalog = sorted(csv.DictReader(fh), key=lambda r: r["site_id"])
    with args.region_index.open(newline="", encoding="utf-8-sig") as fh:
        regions = list(csv.DictReader(fh))
    if not regions:
        p.error("The TDX region index is empty")

    if args.site:
        chosen = [x for x in catalog if x["site_id"] == args.site]
        if not chosen:
            p.error(f"Gauge {args.site} was not found in catalogue")
        suffix = args.site
    else:
        end = args.stop if args.stop is not None else len(catalog)
        if end <= args.start:
            p.error("stop must exceed start")
        chosen = catalog[args.start:end]
        suffix = "all" if args.start == 0 and end >= len(catalog) else f"{args.start:04d}_{end:04d}"
    if not chosen:
        p.error("No gauges selected")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    snaps, candidates = [], []
    for i, site in enumerate(chosen, 1):
        result, rows = snap_one(site, regions, args)
        snaps.append(result)
        candidates.extend(rows)
        print(f"[{i}/{len(chosen)}] {site['site_id']}: "
              f"{result['status']} reach={result['reach_id']} "
              f"region={result['tdx_id']} dist={result['snap_distance_m']}", flush=True)

    snap_path = args.output_dir / f"tdx_gauge_snaps_{suffix}.csv"
    candidate_path = args.output_dir / f"tdx_snap_candidates_{suffix}.csv"
    summary_path = args.output_dir / f"tdx_snap_summary_{suffix}.json"
    write_csv(snap_path, snaps, SNAP_COLUMNS)
    write_csv(candidate_path, candidates, CANDIDATE_COLUMNS)
    summary = {
        "catalog": str(args.catalog.resolve()),
        "region_index": str(args.region_index.resolve()),
        "tdx_root": str(args.tdx_root.resolve()),
        "n_gauges": len(snaps), "status_counts": dict(Counter(s["status"] for s in snaps)),
        "n_distinct_selected_regions": len({s["tdx_id"] for s in snaps if s["tdx_id"]}),
        "search_radius_m": args.radius_m,
        "note": "Provisional gauge snapping only. No TDX river IDs or geometries modified.",
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print("\n" + json.dumps(summary, indent=2))
    print(f"SNAPS:      {snap_path}")
    print(f"CANDIDATES: {candidate_path}")
    print(f"SUMMARY:    {summary_path}")


if __name__ == "__main__":
    main()
