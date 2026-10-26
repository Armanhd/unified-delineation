# Unified TDX-Hydro delineation: user manual (preview workflow)

**Version:** 2026-10-25. **Scope:** generic gauge catalogues; worked ARC example for CAMELS-SPAT. **Outputs are provisional and require GIS/hydrological review before CWARHM/SUMMA.** River-line headwater extensions are intentionally OUT of scope; perform them using your separate downstream Python workflow.

## 1. Concept and input requirements

A gauge catalogue must contain a **unique** `site_id`, `latitude`, `longitude` (WGS84, decimal degrees) and `expected_area_km2`. For CAMELS-SPAT the existing validated catalogue is `indexes/gauge_catalogs/camels-spat_gauge_catalog_validated.csv`. For another dataset, create a separate catalogue and a separate per-dataset controls directory and batch manifest. **Do not mix independent datasets in one batch manifest.**

Inputs:

- Original, read-only TDX `*-basins.gpkg` and `*-streamnet.gpkg` files.
- TDX regional index `indexes/tdx_regions.csv`; it must be regenerated after source files are fixed or replaced.
- Copernicus GLO-30 tiles for DEM-refined boundaries, or an explicit decision to use native TDX polygons.
- This repository's existing scripts under `1_delineation/tdx` and `1_delineation/common`; a functioning `unified_delineation` Python/GeoPandas environment; ARC Slurm.
- A representative existing per-basin **template TOML** (needed by Step 03a).

The stages are **02b global gauge snapping -> 02e consolidation -> 02f shared gauge registry -> 03a controls -> regional topology and split plan -> 02/03/04/05 per basin -> 05b/06b/07b (nested) OR 06/07 (ordinary) -> status collection -> ArcGIS review**. Use `--help` before running individual stages, especially when your installed code is newer than this guide.

## 2. Install the code changes (before restarting)

From the downloadable package, copy:

- `09_tdx_batch.py` -> `1_delineation/common/09_tdx_batch.py` (revised runner; Step03 `--allow-review` for **preview only**).
- `09_global_tdx_id_audit.py` -> `1_delineation/common/09_global_tdx_id_audit.py` (unchanged audit, included for convenience).
- `10_build_arcgis_review_layers.py` -> `1_delineation/common/10_build_arcgis_review_layers.py`.
- `patch_step06_preview_slope.py` -> `1_delineation/common/patch_step06_preview_slope.py`.

**Important Step06 caution.** The ordinary Step06 exporter has changed over your work. We have **not replaced it** with an older version. Instead the included patcher verifies exact anchors and makes a timestamped backup before injecting a **preview-only** neighbor-median slope fallback:

```bash
python 1_delineation/common/patch_step06_preview_slope.py --script 1_delineation/tdx/06_export_cwarhm_shapefiles.py --check
# Run the next command only if the check says it is compatible:
python 1_delineation/common/patch_step06_preview_slope.py --script 1_delineation/tdx/06_export_cwarhm_shapefiles.py
python 1_delineation/tdx/06_export_cwarhm_shapefiles.py --help
```

If anchor verification fails, **stop and use the unmodified exporter**. Do not blindly overwrite or downgrade your installed Step06. The fallback uses up to three graph edges of connected reaches; if no physically connected positive slope donor exists, it blocks rather than inventing a value. All replacements are written to `outputs/<SITE>/diagnostics/step06_preview_slope_substitutions.csv` and in the Step06 metadata summary. This is **not** production-authorized.

For preview slope support in the batch:

```bash
python 1_delineation/common/09_tdx_batch.py submit --label LABEL --preview-ordinary-slope-fallback --dry-run
```

Only use this flag **after** confirming that installed Step06 reports `--preview-slope-fallback` in its `--help` output. Without this opt-in, native zero slope still blocks rather than being filled silently.

## 3. Protect the old CAMELS-SPAT results before a fresh run

**Do not delete during installation.** First verify that no related jobs remain:

```bash
squeue -u "$USER"
```

The existing run `camels_tdx_full` contains useful evidence: 844 selected sites; 192 structurally validated; 652 failed/blocked; 582 not ready. The original shapefile scan found 194 unique sites with files; keep the old run if you want before/after comparisons.

Recommended: create a timestamped archive outside the live working directories. Check free disk and quotas (`du -sh outputs batch_tdx indexes`). Suggested archive targets: `batch_tdx/camels_tdx_full`, `outputs`, `0_control_files/tdx_camels_spat`, `indexes/gauge_catalogs/tdx_snapping`, `indexes/gauge_catalogs/tdx_shared_network`, and `indexes/tdx_regional_topology`. Exclude original TDX/Copernicus source files and the project source code. **Do not delete `indexes/gauge_catalogs/camels-spat_gauge_catalog_validated.csv` or `indexes/tdx_regions.csv` without first regenerating/replacing them.**

For the first rerun use a **new label** (e.g. `camels_tdx_clean_20261025`) rather than overwriting the old `camels_tdx_full` run. To regenerate all outputs, archive then explicitly remove or rename the generated directories only after checking their contents. Examples are in Section 8; do not run them before verifying source-path ownership.

## 4. CAMELS-SPAT: validate repaired TDX source files

```bash
cd /work/comphyd_lab/users/arman.haddadchi/unified_delineation_code
conda activate unified_delineation

# Identify required files and confirm file integrity for each TDX region.
find /work/comphyd_lab/data/geospatial-data/TDX_hydrofabrics/TDX_hydrofabrics \
  -maxdepth 1 -name '*-streamnet.gpkg' | wc -l
find /work/comphyd_lab/data/geospatial-data/TDX_hydrofabrics/TDX_hydrofabrics \
  -maxdepth 1 -name '*-basins.gpkg' | wc -l

python 1_delineation/tdx/01_build_tdx_index.py --help
```

Rebuild `indexes/tdx_regions.csv` using your installed Step01's documented arguments. **The exact invocation depends on the current Step01 version; confirm its `--help` first.** If a GeoPackage has been replaced, the old regional SQLite topology index and previous gauge-snaps must **not** be assumed valid. Rebuild them.

## 5. Re-create snapping and gauge registry from scratch

The existing CAMELS-SPAT Step02b Slurm array file is `1_delineation/tdx/run_02b_snap_gauge_array.sh` (or its actual installed location). Its `#SBATCH --array=0-71%8` assumes 1426 stations and 20 per task; update for a different catalogue. Verify all paths inside the script before submission.

```bash
mkdir -p slurm_logs/tdx_snapping
sbatch 1_delineation/tdx/run_02b_snap_gauge_array.sh
# After all snapping tasks finish:
python 1_delineation/tdx/02e_consolidate_gauge_snaps.py
python 1_delineation/tdx/02f_build_shared_reach_registry.py
```

Review outputs:

```
indexes/gauge_catalogs/tdx_snapping/global_500m/tdx_gauge_snaps_all.csv
indexes/gauge_catalogs/tdx_shared_network/tdx_gauge_reach_registry.csv
```

Classifications such as `PROVISIONAL`, `REVIEW`, `FAILED` are important. `FAILED` snaps are **not** silently accepted. `REVIEW` snapping is allowed only for preview (`--allow-review` in Step03).

## 6. Generate per-basin TOML controls and readiness manifest

```bash
python 1_delineation/common/03a_generate_batch_controls.py --help
```

Locate your **existing original template control file** and supply `--template`. The supported Step03a arguments include `--catalog`, `--snaps`, `--registry`, `--region-index`, `--output-dir`, `--write`, and `--overwrite`:

```bash
python 1_delineation/common/03a_generate_batch_controls.py \
  --template /ABSOLUTE/PATH/TO/YOUR_TEMPLATE.toml \
  --catalog indexes/gauge_catalogs/camels-spat_gauge_catalog_validated.csv \
  --snaps indexes/gauge_catalogs/tdx_snapping/global_500m/tdx_gauge_snaps_all.csv \
  --registry indexes/gauge_catalogs/tdx_shared_network/tdx_gauge_reach_registry.csv \
  --region-index indexes/tdx_regions.csv \
  --output-dir 0_control_files/tdx_camels_spat --write --overwrite
```

Inspect `0_control_files/tdx_camels_spat/batch_controls_manifest.csv`. The **number selected may change** now that TDX GeoPackages have been repaired. Do not assume 844 is still correct. Keep headwater river extension disabled; this is handled in a separate workflow.

## 7. Pilot, audit and full Slurm run

Start with three known nested gauges (if source-ready) and a few ordinary basins:

```bash
python 1_delineation/common/09_tdx_batch.py prepare \
  --label camels_pilot_clean --sites CAN_02HC025,CAN_02HC023,CAN_02HC047
python 1_delineation/common/09_tdx_batch.py submit --label camels_pilot_clean --dry-run
python 1_delineation/common/09_tdx_batch.py submit \
  --label camels_pilot_clean --region-parallel 2 --basin-parallel 4 --export-parallel 3 \
  --preview-ordinary-slope-fallback
```

Omit `--preview-ordinary-slope-fallback` if Step06 has not been patched. When finished:

```bash
python 1_delineation/common/09_tdx_batch.py collect --label camels_pilot_clean
cat batch_tdx/camels_pilot_clean/tdx_batch_results_summary.json
```

Then full CAMELS-SPAT:

```bash
python 1_delineation/common/09_tdx_batch.py prepare --label camels_tdx_clean_20261025 \
  --controls-dir 0_control_files/tdx_camels_spat
python 1_delineation/common/09_tdx_batch.py submit \
  --label camels_tdx_clean_20261025 --dry-run
python 1_delineation/common/09_tdx_batch.py submit \
  --label camels_tdx_clean_20261025 \
  --region-parallel 2 --basin-parallel 8 --export-parallel 6 \
  --preview-ordinary-slope-fallback
```

Use the opt-in fallback flag only after Step06 patch compatibility checks, as above. For regions with serious DEM NoData, the **default** still fails safely. A limited `--max-connected-void-cells N` can be used **only on exploratory preview runs** after studying gap locations and affected watershed boundaries; do not set a high universal threshold. The two cases over `dem_max_cells=30,000,000` require tiling/memory planning, not an arbitrary automatic raise.

Check jobs and failure logs:

```bash
squeue -u "$USER"
sacct -j JOBIDS --format=JobID,JobName,State,ExitCode,Elapsed -X
cat batch_tdx/camels_tdx_clean_20261025/tdx_batch_results_summary.json
```

The batch code uses Slurm `afterany` so other sites continue after failures; downstream tasks are then marked blocked if prerequisites are missing. Inspect `batch_tdx/LABEL/logs/<SITE>/` for detailed per-stage logs.

## 8. Starting truly clean (optional, destructive)

**Run only after archiving and verifying that you are at the intended project root.** The following resets CAMELS outputs/indices/controls; it does **not** delete the validated gauge catalogue, regional source TDX GeoPackages, DEM tiles, or Python source files.

```bash
cd /work/comphyd_lab/users/arman.haddadchi/unified_delineation_code
pwd
squeue -u "$USER"
# Check the targets before deleting:
du -sh outputs batch_tdx indexes/tdx_regional_topology \
  indexes/gauge_catalogs/tdx_snapping \
  indexes/gauge_catalogs/tdx_shared_network \
  0_control_files/tdx_camels_spat 2>/dev/null
# ONLY after you have archived them, intentionally execute:
rm -rf -- outputs batch_tdx indexes/tdx_regional_topology \
  indexes/gauge_catalogs/tdx_snapping \
  indexes/gauge_catalogs/tdx_shared_network \
  0_control_files/tdx_camels_spat
```

Rebuild Step01 region index as well if repaired sources change its content. Remove/rebuild `indexes/tdx_regions.csv` only if you know how to regenerate it. The `rm` block is a **proposed** user-executed action; no remote deletion has been performed.

## 9. Generate ArcGIS review layers

After collecting the batch results (Step 7.2), generate GIS layers for reviewing successful, failed, and pending watershed delineations.

The workflow retains the best available watershed polygon for each basin (`tdx_nested_preview` → `tdx_clipped` → `tdx_native` → partial Step05/Step04). Stations without valid watershed geometry are still included as gauge points.

### Run on Slurm

From the project root, using the batch label defined in Step 7.2:

```bash
mkdir -p slurm_logs

sbatch 1_delineation/tdx/run_10_arcgis_review.sbatch "$LABEL"
```

The Slurm script executes `1_delineation/tdx/10_build_arcgis_review_layers.py`.

### Outputs

All files are saved under `batch_tdx/<LABEL>/arcgis_review/`.

| File | Description |
|---|---|
| `basin_boundaries_for_review.shp` | Best available watershed polygon per station, including partial delineations |
| `all_station_points.shp` | Original gauge locations, including failed and pending stations |
| `all_basins_review_status.csv` | Processing status of all stations, including those without polygons |
| `basin_review.gpkg` | Watershed polygons in GeoPackage format |
| `station_review.gpkg` | Gauge locations and associated attributes |
| `review_layers_summary.json` | Summary of generated geometries and processing statuses |

Key attributes include `site_id`, `region_id`, `wf_status`, `snap_stat`, `geom_src`, `variant`, `has_poly`, and `problem`.

Check the results:

```bash
ls -lh batch_tdx/"$LABEL"/arcgis_review/
```

Open the GIS files in ArcGIS Pro or QGIS to inspect watershed boundaries, outlet positions, river networks, nesting, and failed delineations.

**Note:** These layers are intended for visual review, not production-approved CWARHM exports. The current Slurm wrapper requires its catalogue and environment paths to be configured for the target dataset and computing system.

## 10. Other datasets (colleague checklist)

1. Create a WGS84 station CSV with `site_id,latitude,longitude,expected_area_km2`. Validate IDs and areas.
2. Point Step02b/02e/02f to **dataset-specific** directories (rather than overwriting CAMELS-SPAT outputs).
3. Rebuild the region index after replacing TDX source GeoPackages.
4. Run Step03a with the other dataset's `--catalog` / `--snaps` / `--registry` / `--output-dir` and a suitable TOML template.
5. Run `09_tdx_batch.py prepare --label OTHER_DATASET_RUN --controls-dir 0_control_files/tdx_other_dataset`; it records the selected dataset's controls directory in `run_config.json` for all subsequent phases. Use a separate label and don't mix manifests.
6. Run small pilot; review outputs; then run the Slurm array with a unique batch label. Produce ArcGIS station/partial polygon files and inspect them manually.

## 11. Remaining known limitations

- `--allow-review` enables **diagnostic** snapping only; no automatically inferred gauge location becomes authoritative.
- Missing or corrupted TDX sources require a full region-index/resnap/control refresh.
- Nested plan/dependency blocks may persist if upstream stations were not completed.
- Some Step05 DEM NoData patches and two very large DEM windows remain unresolved, deliberately.
- Native TDX zero slopes can be replaced only with an explicit, flagged, preview-only connected-neighbor estimate; some failures will remain if no donor exists.
- The global ID audit checks **available** indexes, not necessarily all North American TDX regions; approve continental ID uniqueness separately.
- Structural validation PASS is not hydrological validation; always perform independent ArcGIS inspection.
