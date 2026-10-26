#!/usr/bin/env python3
"""Create ArcGIS-ready QA GIS from completed, reviewed, and failed TDX basins.

Inputs: batch_tdx/LABEL/tdx_batch_all_gauges_status.csv (run 'collect' first),
validated CAMELS-SPAT catalog, optional per-site diagnostics and GIS outputs.

One review GeoPackage includes basin_boundaries, river_network, original_gauges,
and selected_outlets. Failed/unmatched gauges remain visible as original points.
The script never changes delineation outputs and never invents a missing outlet,
river, or basin polygon. It writes an attribute-only CSV for *every* status row.

Designed for Slurm: bounded per-basin GIS reads, progress prints, and write-to-disk
by chunks instead of holding the complete continental network in memory.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import traceback
from collections import Counter
from pathlib import Path

import geopandas as gpd
import pandas as pd
from shapely.geometry import Point

VARIANTS = ("tdx_nested_preview", "tdx_clipped", "tdx_native")
GIS_LAYERS = ("basin_boundaries", "river_network", "original_gauges", "selected_outlets")
FIELD_NAMES = [
    "site_id", "tdx_id", "wf_status", "snap_stat", "val_status", "str_pass",
    "has_basin", "has_river", "has_outlet", "geom_src", "river_src",
    "orig_lon", "orig_lat", "sel_lon", "sel_lat", "snap_m", "reach_id",
    "problem", "gis_error",
]


def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def finite_coord(value, limit):
    try:
        x = float(value)
    except (ValueError, TypeError):
        return None
    return x if math.isfinite(x) and abs(x) <= limit and x not in (-999, -9999) else None


def first_value(*values):
    for value in values:
        if value is not None and str(value).strip().lower() not in ("", "nan", "none", "null"):
            return value
    return ""


def json_optional(path):
    if not path.is_file():
        return {}, ""
    try:
        d = json.loads(path.read_text())
        if not isinstance(d, dict):
            return {}, "Not a JSON object: " + str(path)
        return d, ""
    except Exception as exc:
        return {}, f"Cannot read {path.name}: {type(exc).__name__}: {exc}"


def one_frame(path, layer=None, *, wanted=None):
    """Read a candidate GIS layer; reject missing CRS/invalid geometry types."""
    gdf = gpd.read_file(path, layer=layer)
    if gdf.crs is None or gdf.empty:
        raise ValueError("Empty GIS layer or missing CRS")
    gdf = gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty].copy()
    if wanted == "polygon":
        gdf = gdf[gdf.geometry.geom_type.isin(("Polygon", "MultiPolygon"))]
    elif wanted == "line":
        gdf = gdf[gdf.geometry.geom_type.isin(("LineString", "MultiLineString"))]
    if gdf.empty:
        raise ValueError(f"No {wanted} geometries")
    return gdf.to_crs("EPSG:4326")


def collect_boundary(project, site):
    base = project / "outputs" / site
    choices = []
    for variant in VARIANTS:
        folder = base / "cwarhm_exports" / variant / "shapefiles" / "catchment"
        choices.extend((p, None, "EXPORT_" + variant) for p in sorted(folder.glob("*.shp")))
    choices += [
        (base / "tdx_upstream_refined.gpkg", "refinement_boundary", "STEP05_BOUNDARY"),
        (base / "tdx_upstream_refined.gpkg", "refined_subbasins", "PARTIAL_STEP05"),
        (base / "tdx_upstream_native.gpkg", "native_subbasins", "PARTIAL_STEP04"),
    ]
    problems = []
    for path, layer, source in choices:
        if not path.is_file():
            continue
        try:
            frame = one_frame(path, layer=layer, wanted="polygon")
            combined = frame.geometry.union_all()
            if combined.is_empty:
                raise ValueError("Empty polygon union")
            return combined, source, problems
        except Exception as exc:
            problems.append(f"{path.name}:{type(exc).__name__}:{exc}")
    return None, "NO_POLYGON", problems


def collect_rivers(project, site):
    base = project / "outputs" / site
    choices = []
    for variant in VARIANTS:
        folder = base / "cwarhm_exports" / variant / "shapefiles" / "river_network"
        choices.extend((p, None, "EXPORT_" + variant) for p in sorted(folder.glob("*.shp")))
    choices += [
        (base / "tdx_nested_partition_preview.gpkg", "nested_reaches", "NESTED_STEP05B"),
        (base / "tdx_upstream_refined.gpkg", "upstream_reaches", "STEP05_REACHES"),
        (base / "tdx_upstream_native.gpkg", "upstream_reaches", "STEP04_REACHES"),
    ]
    problems = []
    for path, layer, source in choices:
        if not path.is_file():
            continue
        try:
            frame = one_frame(path, layer=layer, wanted="line")
            # Keep actual original/split identifiers for ArcGIS reach inspection.
            if "reach_id" not in frame.columns:
                for alt in ("COMID", "comid", "LINKNO"):
                    if alt in frame.columns:
                        frame["reach_id"] = frame[alt].astype(str)
                        break
            frame["reach_id"] = frame.get("reach_id", pd.Series("", index=frame.index)).astype(str)
            if "downstream_id" not in frame.columns:
                for alt in ("NextDownID", "DSLINKNO", "nextdownid"):
                    if alt in frame.columns:
                        frame["downstream_id"] = frame[alt].astype(str)
                        break
            frame["downstream_id"] = frame.get("downstream_id", pd.Series("", index=frame.index)).astype(str)
            return frame[["reach_id", "downstream_id", "geometry"]].copy(), source, problems
        except Exception as exc:
            problems.append(f"{path.name}:{type(exc).__name__}:{exc}")
    return None, "NO_RIVER", problems


def safe_row(text, max_len=2000):
    return str(text if text is not None else "")[:max_len]


def batch_append(output, layer, features, *, exists):
    if not features:
        return exists
    gdf = gpd.GeoDataFrame(features, geometry="geometry", crs="EPSG:4326")
    gdf.to_file(output, layer=layer, driver="GPKG", mode="a" if exists else "w", index=False)
    return True


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path.cwd())
    p.add_argument("--label", required=True)
    p.add_argument("--catalog", type=Path, help="Validated gauge catalog, with site_id, latitude, longitude")
    p.add_argument("--output-dir", type=Path)
    p.add_argument("--chunk-size", type=int, default=20, help="Number of basins per GeoPackage flush")
    p.add_argument("--progress-every", type=int, default=25)
    args = p.parse_args()
    if args.chunk_size < 1 or args.progress_every < 1:
        p.error("--chunk-size and --progress-every must be positive")
    root = args.root.resolve()
    base = root / "batch_tdx" / args.label
    status_path = base / "tdx_batch_all_gauges_status.csv"
    if not status_path.is_file():
        p.error(f"Missing {status_path}. Run the batch collect step first.")
    candidates = [
        root / "indexes/gauge_catalogs/camels-spat_gauge_catalog_validated.csv",
        root / "0_control_files/tdx_camels_spat/camels-spat_gauge_catalog.csv",
        root / "indexes/gauge_catalogs/camels-spat_gauge_catalog.csv",
    ]
    catalog = args.catalog.resolve() if args.catalog else next((x for x in candidates if x.is_file()), None)
    if catalog is None or not catalog.is_file():
        p.error("Missing catalog; use --catalog /absolute/path/to/catalog.csv")
    stations = {row["site_id"]: row for row in read_csv(catalog)}
    statuses = read_csv(status_path)
    out = (args.output_dir or base / "arcgis_review").resolve()
    out.mkdir(parents=True, exist_ok=True)
    gpkg = out / "camels_spat_tdx_review.gpkg"
    if gpkg.exists():
        gpkg.unlink()  # Only regenerated review artifact, never original basin data.

    counts = Counter()
    sources = Counter()
    cache = {key: [] for key in GIS_LAYERS}
    wrote = {key: False for key in GIS_LAYERS}
    status_csv = out / "all_basins_review_status.csv"
    mismatch_csv = out / "gis_read_errors.csv"
    with status_csv.open("w", newline="", encoding="utf-8") as status_stream, \
         mismatch_csv.open("w", newline="", encoding="utf-8") as mismatch_stream:
        csv_writer = csv.DictWriter(status_stream, fieldnames=FIELD_NAMES)
        csv_writer.writeheader()
        error_writer = csv.DictWriter(mismatch_stream, fieldnames=["site_id", "gis_error"])
        error_writer.writeheader()

        for ix, status in enumerate(statuses, start=1):
            site = str(status.get("site_id", "")).strip()
            if not site or "/" in site or ".." in site or "\\" in site:
                raise ValueError(f"Unsafe/missing site ID in status CSV: {site!r}")
            site_path = root / "outputs" / site
            gauge = stations.get(site, {})
            orig_lon = finite_coord(gauge.get("longitude"), 180)
            orig_lat = finite_coord(gauge.get("latitude"), 90)
            # Selection JSON exists only after Step 03; use no fabricated outlet for P01.
            selection_paths = [
                site_path / "diagnostics" / "tdx_outlet_selection.json",
                site_path / "diagnostics" / "step03" / f"{site}_tdx_outlet_selection.json",
            ]
            selected_path = next((x for x in selection_paths if x.is_file()), selection_paths[0])
            selected, json_error = json_optional(selected_path)
            sel_lon = finite_coord(selected.get("snapped_lon"), 180)
            sel_lat = finite_coord(selected.get("snapped_lat"), 90)
            errors = [json_error] if json_error else []
            try:
                boundary, boundary_src, issues = collect_boundary(root, site)
                errors += issues
            except Exception as exc:
                boundary, boundary_src = None, "NO_POLYGON"
                errors.append("Basin read failure: " + safe_row(exc))
            try:
                rivers, river_src, issues = collect_rivers(root, site)
                errors += issues
            except Exception as exc:
                rivers, river_src = None, "NO_RIVER"
                errors.append("River read failure: " + safe_row(exc))
            workflow = first_value(status.get("workflow_status"), status.get("status"), "UNKNOWN")
            problem = first_value(status.get("problem"), status.get("error"), status.get("error_message"))
            meta = {
                "site_id": site,
                "tdx_id": safe_row(first_value(status.get("tdx_id"), selected.get("tdx_id")), 30),
                "wf_status": safe_row(workflow, 80),
                "snap_stat": safe_row(first_value(status.get("snap_status"), selected.get("global_snap_status")), 80),
                "val_status": safe_row(first_value(status.get("validation_status"), status.get("validation")), 80),
                "str_pass": safe_row(first_value(status.get("structural_validation_pass")), 20),
                "has_basin": int(boundary is not None),
                "has_river": int(rivers is not None),
                "has_outlet": int(sel_lon is not None and sel_lat is not None),
                "geom_src": boundary_src,
                "river_src": river_src,
                "orig_lon": orig_lon, "orig_lat": orig_lat,
                "sel_lon": sel_lon, "sel_lat": sel_lat,
                "snap_m": finite_coord(selected.get("snap_distance_m"), 1e12),
                "reach_id": safe_row(first_value(selected.get("reach_id"), gauge.get("reach_id")), 80),
                "problem": safe_row(problem, 2000),
                "gis_error": safe_row(" | ".join(errors), 2000),
            }
            csv_writer.writerow(meta)
            if errors:
                error_writer.writerow({"site_id": site, "gis_error": meta["gis_error"]})
            if boundary is not None:
                cache["basin_boundaries"].append({**meta, "geometry": boundary})
                counts["boundaries"] += 1
            if rivers is not None:
                # Per-reach rows maintain reach IDs and downstream IDs for tracing.
                for r in rivers.itertuples():
                    cache["river_network"].append({**meta, "reach_id": safe_row(r.reach_id, 80),
                                                   "nextdown": safe_row(r.downstream_id, 80),
                                                   "geometry": r.geometry})
                counts["river_reaches"] += len(rivers)
                counts["sites_with_rivers"] += 1
            if orig_lon is not None and orig_lat is not None:
                cache["original_gauges"].append({**meta, "geometry": Point(orig_lon, orig_lat)})
                counts["original_gauges"] += 1
            else:
                counts["sites_without_original_coords"] += 1
            if sel_lon is not None and sel_lat is not None:
                cache["selected_outlets"].append({**meta, "geometry": Point(sel_lon, sel_lat)})
                counts["selected_outlets"] += 1
            counts["stations"] += 1
            sources[boundary_src] += 1
            if ix % args.chunk_size == 0 or ix == len(statuses):
                for layer in GIS_LAYERS:
                    if cache[layer]:
                        wrote[layer] = batch_append(gpkg, layer, cache[layer], exists=gpkg.exists())
                        cache[layer].clear()
            if ix % args.progress_every == 0 or ix == len(statuses):
                status_stream.flush()
                print(f"Processed {ix}/{len(statuses)} sites; {counts['boundaries']} boundaries, "
                      f"{counts['river_reaches']} reaches, {counts['original_gauges']} gauges", flush=True)

    summary = {
        "batch_label": args.label,
        "stations": counts["stations"],
        "original_gauges": counts["original_gauges"],
        "selected_outlets": counts["selected_outlets"],
        "basins_with_polygon": counts["boundaries"],
        "sites_with_rivers": counts["sites_with_rivers"],
        "river_reaches": counts["river_reaches"],
        "basins_without_polygon": counts["stations"] - counts["boundaries"],
        "polygon_source_counts": dict(sources),
        "gpkg_layers_written": [key for key in GIS_LAYERS if wrote[key]],
        "source_catalog": str(catalog),
        "geo_package": str(gpkg) if gpkg.exists() else None,
        "status_csv": str(status_csv),
        "gis_error_csv": str(mismatch_csv),
        "note": "Original gauges are included for failed/unmatched sites. No missing outlet, polygon, or river is fabricated. Partial GIS outputs are explicitly source-labelled; REVIEW does not mean approved.",
    }
    (out / "review_layers_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
