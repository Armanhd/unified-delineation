#!/usr/bin/env python3
"""Safely patch current ARC Step06 exporter for an explicit preview-only slope fallback.

Backs up the original. Aborts if expected anchors have changed. NEVER edits TDX data.
"""
from pathlib import Path
import argparse
import shutil
import datetime

ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--script',type=Path,required=True)
ap.add_argument('--check',action='store_true',help='Verify patch anchors; do not write')
args=ap.parse_args()
p=args.script.resolve()
s=p.read_text()
if '--preview-slope-fallback' in s:
    print('Already has preview slope fallback option; no change')
    raise SystemExit(0)
anchors=[
    "    args = ap.parse_args()",
    "    reach_records = []",
    "        if slope is None:\n            raise ValueError(f'No valid positive TDX slope for reach {rid}; explicit DEM fallback required, not a made-up value')",
    "    reaches = gpd.GeoDataFrame(reach_records, crs=metric)",
]
notfound=[a[:60] for a in anchors if s.count(a)!=1]
if notfound:
    raise SystemExit('INCOMPATIBLE Step06 source; not modified. Unmatched/nonunique anchors: '+repr(notfound))
if args.check:
    print('Step06 source supports safe patch; no files changed')
    raise SystemExit(0)
s=s.replace("    args = ap.parse_args()", "    ap.add_argument('--preview-slope-fallback', choices=['none', 'neighbor-median'], default='none',\n                    help='Explicit preview-only network-neighbor slope substitution; NEVER production-approved')\n    args = ap.parse_args()")
s=s.replace("    reach_records = []", "    reach_records = []\n    missing_slope_ids = []\n    slope_preview_substitutions = []")
s=s.replace("        if slope is None:\n            raise ValueError(f'No valid positive TDX slope for reach {rid}; explicit DEM fallback required, not a made-up value')",
"        if slope is None:\n            if args.preview_slope_fallback != 'neighbor-median':\n                raise ValueError(f'No valid positive TDX slope for reach {rid}; explicit DEM fallback required, not a made-up value')\n            missing_slope_ids.append(rid)\n            slope = float('nan')")
replacement='''    # Preview-only fallback: interpolate from physically connected reaches, not arbitrary network mean.
    if missing_slope_ids:
        from collections import defaultdict, deque
        import statistics
        by_id = {int(rec['COMID']): rec for rec in reach_records}
        neighbors = defaultdict(set)
        for rec in reach_records:
            reach_id, down = int(rec['COMID']), int(rec['NextDownID'])
            if down in by_id:
                neighbors[reach_id].add(down)
                neighbors[down].add(reach_id)
        for rid in missing_slope_ids:
            explored = {rid}
            frontier = [rid]
            donors = []
            for distance in range(1, 4):
                next_frontier = sorted({j for i in frontier for j in neighbors[i] if j not in explored})
                explored.update(next_frontier)
                donors = [j for j in next_frontier
                          if pd.notna(by_id[j]['slope']) and np.isfinite(by_id[j]['slope'])
                          and float(by_id[j]['slope']) > 0]
                if donors:
                    val = float(statistics.median(float(by_id[j]['slope']) for j in donors))
                    by_id[rid]['slope'] = val
                    slope_preview_substitutions.append({'reach_id':rid,'preview_slope':val,
                        'donors':','.join(map(str,donors)),'network_distance':distance,
                        'production_usable':False})
                    break
                frontier=next_frontier
                if not frontier:
                    break
            if not donors:
                raise ValueError(f'No positive TDX slope donor within 3 connected reaches of {rid}; remain BLOCKED')
    # Explicit QA record; this flag must persist for review even if geometry is structurally valid.
    if slope_preview_substitutions:
        slope_qa_dir = root / 'diagnostics'
        slope_qa_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(slope_preview_substitutions).to_csv(
            slope_qa_dir / 'step06_preview_slope_substitutions.csv', index=False)
'''
s=s.replace("    reaches = gpd.GeoDataFrame(reach_records, crs=metric)",replacement+"    reaches = gpd.GeoDataFrame(reach_records, crs=metric)")
# Add audit flags to the Step06 metadata summary without altering downstream schema.
marker="               'slope_provenance':slope_provenance, 'exports':{}}"
if s.count(marker)!=1:
    raise SystemExit('Metadata summary anchor missing; refusing to write')
s=s.replace(marker,"               'slope_provenance':slope_provenance, 'exports':{},\n               'preview_only_slope_substitution':bool(slope_preview_substitutions),\n               'preview_slope_substitutions':slope_preview_substitutions,\n               'production_export_authorized':False}")
backup=p.with_name(p.name+'.before_preview_slope_'+datetime.datetime.now().strftime('%Y%m%d_%H%M%S'))
shutil.copy2(p,backup)
p.write_text(s)
print('PATCHED:',p)
print('BACKUP:',backup)
