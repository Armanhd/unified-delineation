#!/usr/bin/env python3
"""03a - Generate portable, basin-specific TDX control TOMLs.

The original single-basin template (e.g. Sagehen) is read but NEVER modified.
Generate controls for one or many gauges with source-specific TDX region files,
new shared-gauge/headwater configuration, and a readiness manifest.

IMPORTANT: Controls are CONFIGURATION DRAFTS, not approved delineations.
Current Steps 03-07 do NOT implement shared nesting or headwater extension.
Both new feature flags are therefore written disabled until integrated/tested.

Default is dry-run.  --write creates controls.  Existing controls require
--overwrite. Missing/corrupt region metadata produces a flagged control with
empty source fields, NOT a dangerous fallback to the template's TDX region.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import tomllib
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CATALOG = ROOT / 'indexes/gauge_catalogs/camels-spat_gauge_catalog_validated.csv'
DEFAULT_SNAPS = ROOT / 'indexes/gauge_catalogs/tdx_snapping/global_500m/tdx_gauge_snaps_all.csv'
DEFAULT_REGISTRY = ROOT / 'indexes/gauge_catalogs/tdx_shared_network/tdx_gauge_reach_registry.csv'
DEFAULT_INDEX = ROOT / 'indexes/tdx_regions.csv'
DEFAULT_CONTROLS = ROOT / '0_control_files/tdx_camels_spat'


def read_csv(path: Path):
    with path.open(newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames:
            raise ValueError(f'CSV has no header: {path}')
        return list(reader)


def unique_by_id(rows, label):
    output = {}
    for row in rows:
        site = str(row.get('site_id') or '').strip()
        if not site or site in output:
            raise ValueError(f'{label}: missing or duplicate site_id {site!r}')
        output[site] = row
    return output


def is_enabled(value):
    value = str(value if value is not None else 'true').strip().lower()
    if value in ('true', 'yes', '1'):
        return True
    if value in ('false', 'no', '0'):
        return False
    raise ValueError(f'Invalid process_basin value: {value!r}')


def number(value, label):
    try:
        out = float(value)
    except (TypeError, ValueError):
        raise ValueError(f'{label}: invalid numeric value {value!r}') from None
    if not math.isfinite(out):
        raise ValueError(f'{label}: non-finite numeric value {value!r}')
    return out


def utm_crs(lat, lon):
    if not -80 <= lat <= 84:
        raise ValueError(f'Cannot auto-select UTM at latitude {lat}')
    zone = max(1, min(60, math.floor((lon + 180) / 6) + 1))
    return f'EPSG:{(32600 if lat >= 0 else 32700) + zone}'


def put(text, table, key, value):
    """Upsert a scalar into a standard TOML table, preserving template text."""
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f'Invalid non-finite value in [{table}] {key}')
    encoded = json.dumps(value, ensure_ascii=True)
    replacement = f'{key} = {encoded}\n'
    lines = text.splitlines(keepends=True)
    table_pattern = re.compile(r'^\s*\[([^\[\]]+)\]\s*(?:#.*)?$')
    key_pattern = re.compile(r'^\s*' + re.escape(key) + r'\s*=')
    start = None
    end = len(lines)
    for i, line in enumerate(lines):
        match = table_pattern.match(line.rstrip('\r\n'))
        if match:
            if start is not None:
                end = i
                break
            if match.group(1).strip() == table:
                start = i
    if start is None:
        return text.rstrip('\n') + f'\n\n[{table}]\n' + replacement
    matches = [i for i in range(start + 1, end) if key_pattern.match(lines[i])]
    if len(matches) > 1:
        raise ValueError(f'Multiple occurrences of [{table}] {key}')
    if matches:
        lines[matches[0]] = replacement
    else:
        lines.insert(end, replacement)
    return ''.join(lines)


def find_preferred_layer(raw, kind):
    if not raw or raw.strip() in ('[]', 'nan'):
        return ''
    try:
        layers = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise ValueError(f'Malformed indexed layer metadata: {raw[:90]}') from exc
    if not isinstance(layers, list):
        raise ValueError('Indexed layers must be a JSON list')
    preferred = 'polygon' if kind == 'basins' else 'line'
    for layer in layers:
        if preferred in str(layer.get('geometry_type', '')).lower():
            return str(layer.get('name', '')).strip()
    return str(layers[0].get('name', '')).strip() if layers else ''


def normalized_integer_id(value):
    raw = str(value or '').strip()
    if raw.endswith('.0') and raw[:-2].isdigit():
        raw = raw[:-2]
    if not raw.isdecimal():
        return ''
    return str(int(raw))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--template', type=Path, required=True)
    parser.add_argument('--catalog', type=Path, default=DEFAULT_CATALOG)
    parser.add_argument('--snaps', type=Path, default=DEFAULT_SNAPS)
    parser.add_argument('--registry', type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument('--region-index', type=Path, default=DEFAULT_INDEX)
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_CONTROLS)
    parser.add_argument('--start', type=int, default=0)
    parser.add_argument('--stop', type=int)
    parser.add_argument('--method', choices=('full', 'dem_refined'), default='dem_refined')
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()

    if args.start < 0 or (args.stop is not None and args.stop <= args.start):
        parser.error('Invalid --start / --stop')
    for path in (args.template, args.catalog, args.snaps, args.region_index):
        if not path.is_file():
            parser.error(f'Required input missing: {path}')
    template = args.template.read_text(encoding='utf-8')
    source = tomllib.loads(template)
    for section in ('domain', 'input', 'output', 'outlet', 'delineation'):
        if section not in source:
            parser.error(f'Template missing [{section}]')
    if not source['input'].get('tdx_directory'):
        parser.error('Template missing [input] tdx_directory')
    if not source['output'].get('output_directory'):
        parser.error('Template missing [output] output_directory')

    catalog = unique_by_id(read_csv(args.catalog), 'catalogue')
    snap_rows = unique_by_id(read_csv(args.snaps), 'snapping results')
    if set(catalog) != set(snap_rows):
        raise ValueError('Gauge catalogue and consolidated snapping results have different site IDs')

    index_rows = read_csv(args.region_index)
    index = {}
    for entry in index_rows:
        rid = normalized_integer_id(entry.get('tdx_id'))
        if not rid or rid in index:
            raise ValueError(f'Duplicate or invalid TDX region index key: {entry.get("tdx_id")!r}')
        index[rid] = entry

    registry = {}
    if args.registry.is_file():
        registry = unique_by_id(read_csv(args.registry), 'shared registry')
        if set(catalog) != set(registry):
            raise ValueError('Shared registry and catalogue IDs differ; regenerate 02f for the current snaps')

    selected = [catalog[k] for k in sorted(catalog)][args.start:args.stop]
    if not selected:
        parser.error('Selected gauge range is empty')

    records = []
    generated = []
    for row in selected:
        site = str(row['site_id']).strip()
        if not is_enabled(row.get('process_basin')):
            continue
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*', site):
            raise ValueError(f'Unsafe site ID: {site!r}')
        lat = number(row['latitude'], f'{site} latitude')
        lon = number(row['longitude'], f'{site} longitude')
        area = number(row['expected_area_km2'], f'{site} drainage area')
        if not (-90 <= lat <= 90 and -180 <= lon <= 180 and area > 0):
            raise ValueError(f'{site}: invalid coordinate or area')

        snap = snap_rows[site]
        reg = registry.get(site, {})
        snap_status = str(snap.get('status', '')).strip()
        snap_region = normalized_integer_id(snap.get('tdx_id'))
        reg_region = normalized_integer_id(reg.get('tdx_id'))
        if reg_region and snap_region and reg_region != snap_region:
            raise ValueError(f'{site}: mismatch between 02f and 02e regional assignments')
        region = reg_region or snap_region
        entry = index.get(region)
        reach = normalized_integer_id(reg.get('parent_tdx_reach_id')) or normalized_integer_id(snap.get('reach_id'))
        flags = []
        if snap_status != 'PROVISIONAL':
            flags.append('SNAP_' + (snap_status or 'UNKNOWN'))
        if not region:
            flags.append('NO_ASSIGNED_TDX_REGION')
        if region and entry is None:
            flags.append('REGION_NOT_IN_INDEX')
        if entry is not None:
            basin_name = str(entry.get('basins_filename') or '').strip()
            stream_name = str(entry.get('streamnet_filename') or '').strip()
            basin_layer = find_preferred_layer(entry.get('basins_layers_json', ''), 'basins')
            stream_layer = find_preferred_layer(entry.get('streamnet_layers_json', ''), 'streamnet')
            if not basin_layer or not stream_layer:
                flags.append('REGION_LAYER_METADATA_INCOMPLETE')
            if not basin_name or not stream_name:
                flags.append('REGION_FILENAME_MISSING')
        else:
            basin_name = stream_name = basin_layer = stream_layer = ''
        if not reach:
            flags.append('NO_MATCHED_PARENT_REACH')
        if not registry:
            flags.append('SHARED_REGISTRY_NOT_GENERATED')
        elif reg.get('group_status') not in ('CANDIDATE_READY_FOR_TOPOLOGY',):
            flags.append('SHARED_GROUP_NEEDS_REVIEW')
        if reg.get('final_export_authorized', '').strip().lower() == 'true':
            raise ValueError(f'{site}: candidate registry unexpectedly claims final export authorization')

        changes = [
            ('domain', 'site_id', site),
            ('domain', 'site_name', site),
            ('domain', 'lat', lat),
            ('domain', 'lon', lon),
            ('domain', 'observed_area_km2', area),
            ('domain', 'delineation_method', 'TDX'),
            ('input', 'tdx_id', region),
            ('input', 'tdx_basins_file', basin_name),
            ('input', 'tdx_streams_file', stream_name),
            ('input', 'tdx_basins_layer', basin_layer),
            ('input', 'tdx_streams_layer', stream_layer),
            ('outlet', 'selection_method', 'automatic'),
            ('outlet', 'manual_reach_id', 0),
            ('outlet', 'upstream_area_units', 'm2'),
            ('delineation', 'projected_crs', utm_crs(lat, lon)),
            ('delineation', 'outlet_subbasin_method', args.method),
            ('overlay', 'overlay_check', False),
            ('overlay', 'overlay_shapefile_directory', ''),
            ('overlay', 'overlay_shapefile', ''),
            ('overlay', 'calculate_area_difference', False),
            ('overlay', 'calculate_spatial_overlap', False),
            ('overlay', 'generate_comparison_map', False),
            # This is the desired mode, but NOT ACTIVE until Steps 03-07 support it.
            ('nested_gauges', 'mode', 'shared'),
            ('nested_gauges', 'enabled', False),
            ('nested_gauges', 'registry_file', str(args.registry.resolve())),
            ('nested_gauges', 'split_network_directory', str(args.registry.resolve().parent)),
            ('headwater', 'enabled', False),
            ('headwater', 'preserve_tdx_reach_id', True),
            ('headwater', 'method', 'copernicus_d8'),
            ('headwater', 'max_area_difference_pct', 20.0),
            # Read-only matching metadata for future integration.
            ('gauge_match', 'snap_status', snap_status),
            ('gauge_match', 'tdx_region', region),
            ('gauge_match', 'parent_tdx_reach_id', reach),
            ('gauge_match', 'shared_group_status', reg.get('group_status', '') if registry else ''),
            ('gauge_match', 'requires_review', bool(flags)),
        ]
        rendered = template
        for table, key, value in changes:
            rendered = put(rendered, table, key, value)
        parsed = tomllib.loads(rendered)
        if parsed['domain']['site_id'] != site or parsed['input']['tdx_id'] != region:
            raise RuntimeError(f'{site}: TOML roundtrip validation failed')
        if parsed['nested_gauges']['enabled'] or parsed['headwater']['enabled']:
            raise RuntimeError('Not-yet-implemented features must remain disabled')
        if parsed['overlay']['overlay_check']:
            raise RuntimeError('Sagehen overlay was not disabled')

        path = args.output_dir / f'control_{site}.toml'
        generated.append((path, rendered))
        records.append({
            'site_id': site,
            'control_path': str(path.resolve()),
            'snap_status': snap_status,
            'tdx_id': region,
            'parent_tdx_reach_id': reach,
            'region_basins_file': basin_name,
            'region_streams_file': stream_name,
            'region_metadata_ready': str(bool(entry and basin_layer and stream_layer)).lower(),
            'shared_group_status': reg.get('group_status', '') if registry else 'REGISTRY_NOT_AVAILABLE',
            'pending_flags': ';'.join(flags),
            'source_input_ready': str(bool(region and entry and basin_layer and stream_layer)).lower(),
            'nested_processing_enabled': 'false',
            'headwater_processing_enabled': 'false',
            'final_export_approved': 'false',
        })

    if args.write:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        collisions = [str(path) for path, _ in generated if path.exists()]
        manifest_path = args.output_dir / 'batch_controls_manifest.csv'
        if manifest_path.exists():
            collisions.append(str(manifest_path))
        if collisions and not args.overwrite:
            raise FileExistsError(f'{len(collisions)} output(s) already exist, e.g. {collisions[:3]}. '
                                  'Use --overwrite for generated controls.')
        for path, rendered in generated:
            path.write_text(rendered, encoding='utf-8')
        if records:
            with manifest_path.open('w', newline='', encoding='utf-8') as stream:
                writer = csv.DictWriter(stream, fieldnames=list(records[0]))
                writer.writeheader()
                writer.writerows(records)

    snap_counts = Counter(r['snap_status'] for r in records)
    issue_counts = Counter(flag for row in records for flag in row['pending_flags'].split(';') if flag)
    summary = {
        'mode': 'WRITE' if args.write else 'DRY_RUN',
        'target_controls': len(records),
        'snap_status_counts': dict(sorted(snap_counts.items())),
        'pending_flag_counts': dict(sorted(issue_counts.items())),
        'controls_with_region_metadata': sum(r['source_input_ready'] == 'true' for r in records),
        'all_generated_controls_export_approved': False,
        'nested_and_headwater_features_enabled': False,
        'Sagehen_template_modified': False,
        'output_directory': str(args.output_dir),
        'sample_sites': [r['site_id'] for r in records[:5]],
    }
    print(json.dumps(summary, indent=2))
    if args.write:
        print('MANIFEST:', args.output_dir / 'batch_controls_manifest.csv')


if __name__ == '__main__':
    main()
