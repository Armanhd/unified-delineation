#!/usr/bin/env python3
"""Resolve missing TDX source slopes from original connected TDX reaches.

Resolution order: native positive parent slope; nearest upstream positive native
slope (main branch by drainage area, if present); nearest downstream positive
native slope. Never fabricate a number or modify source GeoPackages.

Only filtered rows are read. Connectivity is taken from the original regional
TDX stream GeoPackage, not newly created fractional IDs.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import math
import sqlite3


def _name(cols, choices, required=True):
    options = {str(c).lower(): str(c) for c in cols}
    for candidate in choices:
        if candidate and candidate.lower() in options:
            return options[candidate.lower()]
    if required:
        raise KeyError(f'Cannot find {choices} in TDX columns {sorted(cols)}')
    return None


def _quote(identifier):
    return '"' + str(identifier).replace('"', '""') + '"'


def _positive(value):
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return result if math.isfinite(result) and result > 0 else None


class TDXSlopeResolver:
    """Read a bounded part of the original TDX network and resolve slopes."""

    def __init__(self, cfg, max_hops=20):
        inputs = cfg['input']
        region = str(inputs.get('tdx_id', '')).strip()
        self.path = Path(inputs['tdx_directory']) / (
            inputs.get('tdx_streams_file') or f'{region}-streamnet.gpkg')
        if not self.path.is_file():
            raise FileNotFoundError(f'TDX source GeoPackage missing: {self.path}')
        self.conn = sqlite3.connect(self.path.resolve().as_uri() + '?mode=ro', uri=True)
        tables = [r[0] for r in self.conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table'")]
        preferred = inputs.get('tdx_streams_layer')
        if preferred:
            self.table = _name(tables, [preferred])
        else:
            self.table = next((t for t in tables if 'streamnet' in t.lower()), None)
            if self.table is None:
                raise KeyError('No TDX streamnet table found')
        cols = [r[1] for r in self.conn.execute(f'PRAGMA table_info({_quote(self.table)})')]
        fields = cfg.get('tdx_fields', {})
        self.idcol = _name(cols, [fields.get('reach_id'), 'LINKNO'])
        self.slopecol = _name(cols, [fields.get('slope'), 'Slope'])
        self.downcol = _name(cols, [fields.get('downstream_id'), 'DSLINKNO', 'downstream_id', 'NextDownID'])
        self.areacol = _name(cols, ['USContArea', 'UpArea', 'upstream_area', 'DrainArea',
                                      'DrainageArea', 'Area', 'area'], required=False)
        self.cache = {}
        self.up_cache = {}
        self.max_hops = int(max_hops)
        self.table_sql = _quote(self.table)
        self.id_sql = _quote(self.idcol)
        self.down_sql = _quote(self.downcol)
        fields_sql = [self.id_sql, _quote(self.slopecol), self.down_sql]
        if self.areacol:
            fields_sql.append(_quote(self.areacol))
        self.fields_sql = ', '.join(fields_sql)

    def close(self):
        self.conn.close()

    def _consume(self, rows):
        result = []
        for r in rows:
            rid = int(r[0])
            down = int(r[2]) if r[2] is not None else 0
            area = _positive(r[3]) if self.areacol else None
            record = {'id': rid, 'slope': _positive(r[1]), 'down': down, 'area': area}
            self.cache[rid] = record
            result.append(record)
        return result

    def get_many(self, ids):
        missing = sorted(set(int(x) for x in ids) - self.cache.keys())
        for i in range(0, len(missing), 400):
            chunk = missing[i:i+400]
            params = ','.join('?' for _ in chunk)
            sql = (f'SELECT {self.fields_sql} FROM {self.table_sql} '
                   f'WHERE {self.id_sql} IN ({params})')
            self._consume(self.conn.execute(sql, chunk))
        absent = sorted(set(missing) - self.cache.keys())
        if absent:
            raise ValueError(f'TDX parent LINKNO missing from source: {absent[:12]}')
        return {i: self.cache[i] for i in ids}

    def upstream(self, rid):
        rid = int(rid)
        if rid not in self.up_cache:
            sql = (f'SELECT {self.fields_sql} FROM {self.table_sql} '
                   f'WHERE {self.down_sql} = ? AND {self.id_sql} != ?')
            rows = self._consume(self.conn.execute(sql, (rid, rid)))
            # Prefer a large main stem; in absence of area, deterministic LINKNO.
            rows.sort(key=lambda x: (-(x['area'] or 0.0), x['id']))
            self.up_cache[rid] = [r['id'] for r in rows]
        return self.up_cache[rid]

    def resolve(self, parents):
        """Return mapping parent LINKNO -> positive slope; and replacement metadata.

        All donors must have valid ORIGINAL TDX source slopes. Upstream is
        preferred across all distances, and the downstream search is used only
        if no suitable upstream donor is found within max_hops.
        """
        ids = sorted(set(int(p) for p in parents))
        original = self.get_many(ids)
        resolved = {}
        substitutions = []
        for rid in ids:
            value = original[rid]['slope']
            if value is not None:
                resolved[rid] = value
                continue
            donor, hops, method = self._find_upstream(rid)
            if donor is None:
                donor, hops, method = self._find_downstream(rid)
            if donor is None:
                raise ValueError(f'No valid original TDX slope upstream/downstream '
                                 f'within {self.max_hops} reaches for parent LINKNO {rid}')
            replacement = self.cache[donor]['slope']
            resolved[rid] = replacement
            substitutions.append({
                'parent_tdx_reach_id': rid,
                'original_tdx_slope': None,
                'assigned_slope': replacement,
                'donor_tdx_reach_id': donor,
                'method': method,
                'network_hops': hops,
                'production_usable': False,
            })
        metadata = {
            'source_file': str(self.path), 'source_layer': self.table,
            'original_id_field': self.idcol, 'original_slope_field': self.slopecol,
            'original_downstream_field': self.downcol,
            'upstream_branch_priority_field': self.areacol,
            'slope_policy': 'original_parent_then_upstream_then_downstream',
            'max_network_hops': self.max_hops,
            'missing_or_nonpositive_parent_slopes': [r['parent_tdx_reach_id'] for r in substitutions],
        }
        return resolved, substitutions, metadata

    def _find_upstream(self, origin):
        visited = {origin}
        frontier = [origin]
        for distance in range(1, self.max_hops + 1):
            next_frontier = []
            for current in frontier:
                for parent in self.upstream(current):
                    if parent not in visited:
                        visited.add(parent)
                        next_frontier.append(parent)
            if not next_frontier:
                break
            valid = [i for i in next_frontier if self.cache[i]['slope'] is not None]
            if valid:
                valid.sort(key=lambda x: (-(self.cache[x]['area'] or 0), x))
                return valid[0], distance, 'upstream_source_main_branch'
            frontier = next_frontier
        return None, None, None

    def _find_downstream(self, origin):
        visited = {origin}
        current = origin
        for distance in range(1, self.max_hops + 1):
            downstream_id = self.cache[current]['down']
            if downstream_id <= 0 or downstream_id in visited:
                break
            visited.add(downstream_id)
            try:
                self.get_many([downstream_id])
            except ValueError:
                break
            if self.cache[downstream_id]['slope'] is not None:
                return downstream_id, distance, 'downstream_source_main_branch'
            current = downstream_id
        return None, None, None
