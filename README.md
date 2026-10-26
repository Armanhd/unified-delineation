# Unified Delineation

**Unified Delineation** is a research workflow for delineating watersheds, organizing river networks, and preparing geospatial inputs for distributed hydrological models. Its modular structure separates common data preparation and validation from delineation methods, allowing additional hydrofabric- and DEM-based methods to be incorporated over time.

The currently implemented and batch-tested workflow uses **TDX-Hydro** river networks and catchment polygons. It supports one gauge or a catalogue of gauges, including candidate gauge-to-river snapping, upstream-network tracing, nested/shared-reach analysis, per-basin CWARHM-style GIS exports, structural validation, and a combined ArcGIS review dataset.

> **Research / QA status:** Outputs are preliminary. A completed Slurm job, successful GIS export, or structural validation pass does **not** establish that a gauge location or watershed is hydrologically correct or approved for CWARHM/SUMMA production. Review status and any uncertainties must be preserved.

## What the workflow does

| Component | Main scripts | Function |
|---|---|---|
| **1. Gauge catalogue** | `1_delineation/common/00_prepare_gauge_catalog.py` | Check station IDs, WGS84 coordinates, and optional reference drainage area. |
| **2. TDX source inventory** | `1_delineation/tdx/01_build_tdx_index.py` | Index TDX-Hydro regional basin and stream-network GeoPackages; report metadata problems. |
| **3. Gauge snapping** | `1_delineation/tdx/02b_snap_gauge_catalog.py`, `run_02b_snap_gauge_array.sh` | Identify candidate TDX reaches near gauges. Default initial search is 500 m; failed sites can be retried at 2,000 m. |
| **4. Snapping QA and nested gauges** | `02e_consolidate_gauge_snaps.py`, `02f_build_shared_reach_registry.py` | Consolidate candidate matches, reconcile retries, and identify gauges on shared TDX parent reaches. |
| **5. Basin configuration** | `1_delineation/common/03a_generate_batch_controls.py` | Generate a basin-specific TOML control and source-readiness manifest. |
| **6. Regional/network preparation** | `1_delineation/tdx/04b_plan_shared_nested_splits.py`, `04c_validate_shared_split_geometry.py`, `1_delineation/common/09_global_tdx_id_audit.py` | Prepare nested-reach candidates, test geometry/connectivity, and check TDX identifier consistency. |
| **7. Basin delineation** | `02_select_tdx_region.py`, `03_select_tdx_outlet.py`, `04_trace_tdx_upstream.py`, `05_refine_tdx_outlet.py`, `05b_build_nested_hru_preview.py` | Select the outlet, trace upstream catchments/reaches, and prepare the selected full/native or applicable refinement/nested previews. |
| **8. GIS export and validation** | `06_export_cwarhm_shapefiles.py`, `06b_export_nested_cwarhm_preview.py`, `07_validate_and_visualise.py`, `07b_validate_nested_cwarhm_preview.py` | Export CWARHM-style catchments and river geometries, with QA/validation summaries. |
| **9. Batch management** | `1_delineation/common/09_tdx_batch.py` | Prepare tasks, submit dependency-controlled Slurm arrays, and collect results. |
| **10. ArcGIS review** | `1_delineation/tdx/10_build_arcgis_review_layers.py` | Assemble available results for *all* gauges, including unmatched/failed sites, into one review GeoPackage. |

The TDX pipeline uses native TDX hydrological connectivity as the basis for the upstream catchment. The tested `--method full` setup retains native TDX HRUs; it does **not** imply that DEM-based trimming of the outlet HRU has been applied. Other methods can be introduced under `1_delineation/` without changing the purpose of the common components.

## Repository layout

```text
unified-delineation/
├── .gitignore
├── README.md
├── 0_control_files/
│   └── control_USGS_10343500.toml   # Example/template settings
├── 1_delineation/
│   ├── common/                       # Catalogue, controls, Slurm, ID audit
│   ├── tdx/                          # TDX-Hydro stages and ArcGIS review
│   ├── merit/                        # Space for additional methods
│   └── dem/                          # Space for DEM-based workflows
├── environments/
│   ├── setup_conda_env.sh
│   ├── load_module_env.sh
│   └── unified_delineation.yml
└── docs/
    └── TDX_Hydro_Delineation_Manual.md
```

The `docs/` manual is a separate, detailed operating guide. The `merit/` and `dem/` paths do not necessarily contain complete, tested implementations. Generated inputs, temporary indexes, model-ready files, and batch outputs live outside version control.

## Inputs and prerequisites

- **Gauge catalogue (CSV):** `site_id`, `latitude`, `longitude` (EPSG:4326); optional `expected_area_km2` and `process_basin`.
- **TDX-Hydro inputs:** regional `*-basins.gpkg` and `*-streamnet.gpkg` files, supplied separately. Large source hydrofabrics are **not** included in this repository.
- **Control TOML:** a verified single-basin template containing source paths, delineation settings, and output settings.
- **Environment:** Python geospatial and raster libraries specified under `environments/`; the current parallel orchestration uses Slurm on an HPC system.

**First-time environment setup on the University of Calgary ARC system:**

```bash
cd /work/comphyd_lab/users/arman.haddadchi/unified_delineation_code
sbatch environments/setup_conda_env.sh
```

After the setup job completes successfully, activate in a new/current shell:

```bash
module load conda/base
conda activate unified_delineation
```

These absolute ARC paths and module names are installation-specific; other users must configure their own repository, TDX data, and environment locations.

## Running a TDX-Hydro batch: overview

The **[TDX-Hydro delineation manual](docs/TDX_Hydro_Delineation_Manual.md)** contains the complete commands, retry/reconciliation code, validation checks, and Slurm guidance. In brief:

1. Validate the station catalogue and index the TDX source data.
2. Snap stations to candidate river reaches at 500 m; optionally retry only failed sites at 2,000 m, preserving the first pass and the selected search radius.
3. Consolidate matches, build the shared-reach registry, and generate per-basin TOML controls.
4. Prepare the batch manifests using `09_tdx_batch.py prepare` **before** `submit`.
5. Dry-run the submission, then launch regional, ID audit, basin, export, and collection tasks with Slurm.
6. Review collected QA results; create the combined ArcGIS review GeoPackage as a separate final step.

For example, after controls are prepared:

```bash
python 1_delineation/common/09_tdx_batch.py prepare \
  --label my_tdx_run \
  --controls-dir 0_control_files/my_tdx_controls

python 1_delineation/common/09_tdx_batch.py submit \
  --label my_tdx_run \
  --region-parallel 2 \
  --basin-parallel 8 \
  --export-parallel 6 \
  --dry-run

# If the dry-run is correct, repeat without --dry-run.
# After batch processing:
python 1_delineation/common/09_tdx_batch.py collect --label my_tdx_run
```

The current 500 m snapping Slurm wrapper and ArcGIS Slurm wrapper may contain dataset- or ARC-specific paths. Check and adapt their catalogue, project root, environment, and batch label before independent use.

## Outputs: what is saved, and where?

All paths below are **relative to the project root** unless otherwise noted. Replace `<site_id>` with the gauge ID and `<label>` with the batch name supplied to `09_tdx_batch.py`. Output folders are created by the relevant stages; a failed or skipped stage may not produce all listed files.

### A. Dataset-wide preparation outputs

| Output | Path | Contents |
|---|---|---|
| Validated stations | `indexes/gauge_catalogs/<dataset>_gauge_catalog_validated.csv` *(name depends on catalogue)* | Validated gauge IDs and coordinates. |
| TDX regional index | `indexes/tdx_regions.csv` | Regional TDX source inventory. |
| TDX health report | `indexes/tdx_local_health.csv` | Metadata issues and source availability. |
| Initial snap batches | `indexes/gauge_catalogs/tdx_snapping/batch_500m/` | Per-batch `tdx_gauge_snaps_*.csv`, candidate CSVs, and JSON summaries. |
| Consolidated 500 m matches | `indexes/gauge_catalogs/tdx_snapping/consolidated_500m/tdx_gauge_snaps_all.csv` | One row per gauge, with matching status. |
| Optional 2,000 m retry | `indexes/gauge_catalogs/tdx_snapping/retry_2000m/` | Separate retry outputs for initially failed gauges. |
| Reconciled snaps | `indexes/gauge_catalogs/tdx_snapping/reconciled_500m_2000m/tdx_gauge_snaps_all.csv` | Selected matches, original status, retry flag and selected radius. |
| Shared-reach registry | `indexes/gauge_catalogs/tdx_shared_network/` | Candidate parent-reach groups, split plans, and nested topology information. |
| Basin-specific controls | `0_control_files/<dataset_controls>/control_<site_id>.toml` | Per-gauge configuration. |
| Control manifest | `0_control_files/<dataset_controls>/batch_controls_manifest.csv` | Source readiness and processing flags; includes gauges not ready for automated processing. |

The demonstrated directory names reflect the CAMELS-SPAT workflow; other datasets should use separate catalogue, snap, registry, and controls directories to prevent accidental mixing.

### B. Individual basin GIS and diagnostic outputs

**Main root:** `outputs/<site_id>/`

```text
outputs/<site_id>/
├── tdx_upstream_native.gpkg                 # Upstream TDX subbasins and reaches (Step 04)
├── tdx_upstream_refined.gpkg                # Step 05 refinement geometry, when available
├── tdx_nested_partition_preview.gpkg       # Nested/shared-gauge candidate geometry, if applicable
├── diagnostics/                             # Outlet selection, topology and process QA
├── cwarhm_exports/
│   ├── tdx_native/                          # Native TDX HRU variant, when exported
│   │   ├── shapefiles/
│   │   │   ├── catchment/<site_id>_catchment.shp
│   │   │   ├── river_basins/<site_id>_river_basins.shp
│   │   │   └── river_network/<site_id>_river_network.shp
│   │   └── metadata/                        # Export provenance and diagnostics
│   ├── tdx_clipped/                         # Clipped/refined preview, when exported
│   │   └── shapefiles/...
│   └── tdx_nested_preview/                  # Nested catchment preview, if applicable
│       └── shapefiles/...
└── validation/                              # Step 07 / Step 07b structural checks
```

**Meaning of the three GIS shapefile groups:**

| Folder | Feature | Typical use |
|---|---|---|
| `shapefiles/catchment/` | Individual HRU/subbasin polygons with reach/HRU attributes | Distributed catchment discretization for modelling. |
| `shapefiles/river_basins/` | One merged watershed-boundary polygon | View/check the entire delineated basin boundary. |
| `shapefiles/river_network/` | River-reach polylines and connectivity attributes | River-network and routing preparation. |

Shapefiles are multi-file datasets (`.shp`, `.shx`, `.dbf`, `.prj`, etc.); retain all components when copying GIS outputs. Export variants are method- and validation-dependent: **do not assume every basin has all three variants**. Under the `full` method, refined/intermediate geometry may still be produced for diagnostics, but the native-TDX-HRU configuration is the reference interpretation.

Individual basin validation and diagnostics typically include files under `outputs/<site_id>/validation/` and `outputs/<site_id>/diagnostics/`. Review these before accepting a shapefile as a hydrologically valid product.

### C. Slurm batch results and QA

**Main root:** `batch_tdx/<label>/`

| File or directory | Purpose |
|---|---|
| `sites.csv` | Selected source-ready basin tasks (not necessarily every catalogue gauge). |
| `regions.csv` | Regional task assignments. |
| `run_config.json` | Batch configuration and source-ready counts. |
| `slurm_logs/` | Per-phase Slurm standard output and error logs. |
| `status/` | Per-site and regional processing status JSONs. |
| `id_audit/` | Global identifier audit and conflict diagnostics. |
| `tdx_batch_results.csv` | Processing outcomes for selected basin tasks. |
| `tdx_batch_all_gauges_status.csv` | Status records for the complete input catalogue, including stations not selected for delineation. |
| `tdx_batch_results_summary.json` | Cross-site progress and status totals. |

Slurm uses dependencies between regional preparation, the ID audit, basin processing, export and final collection. Some dependencies are `afterany`, so later jobs may start even when upstream tasks failed. Inspect the collected status reports rather than relying only on a `COMPLETED` job state.

### D. Final combined ArcGIS review dataset

**Main root:** `batch_tdx/<label>/arcgis_review/`

The final review is generated by:

```text
1_delineation/tdx/10_build_arcgis_review_layers.py
```

It reads the completed batch status and available per-site GIS outputs. For the currently implemented review script, the files are:

```text
batch_tdx/<label>/arcgis_review/
├── camels_spat_tdx_review.gpkg      # Combined multi-layer GeoPackage
├── all_basins_review_status.csv    # All-station spatial review/status table
├── gis_read_errors.csv            # GIS read/geometry problems, if any
└── review_layers_summary.json     # Layer and feature counts
```

The name `camels_spat_tdx_review.gpkg` is **currently hard-coded and CAMELS-SPAT-specific**; change/parameterize it before calling the reviewer dataset-independent. The intended GeoPackage layers are:

| Layer | Contents |
|---|---|
| **`basin_boundaries`** | Available delineated boundaries, with site/status/source information. |
| **`river_network`** | Available river reaches, with identifiers and review metadata. |
| **`original_gauges`** | Original station point for **every input gauge**, including failed/unmatched gauges. |
| **`selected_outlets`** | Candidate snapped outlet points wherever a match was available. |

The review script preferentially uses available CWARHM-exported geometries, with Step 04/05 GIS data as labelled fallbacks when later exports are missing. It **does not manufacture a missing basin polygon or outlet**. Fields identify review status, selected source, snapping distance, and problems, so failed and incomplete basins can be inspected alongside successful ones in ArcGIS Pro or QGIS.

After `collect`, check the review script and Slurm wrapper for the correct project path and label:

```bash
python 1_delineation/tdx/10_build_arcgis_review_layers.py --help
sed -n '1,120p' 1_delineation/tdx/run_10_arcgis_review.sbatch
```

If the installed script supports the documented options, the direct command is:

```bash
python 1_delineation/tdx/10_build_arcgis_review_layers.py \
  --root "$PWD" \
  --label my_tdx_run \
  --catalog path/to/my_validated_catalog.csv \
  --output-dir batch_tdx/my_tdx_run/arcgis_review
```

For a large batch, use `1_delineation/tdx/run_10_arcgis_review.sbatch` **only after checking its currently hard-coded paths** (the earlier wrapper referenced the old `common/` script location).

## Validation and interpretation

- **`REVIEW`**: A candidate gauge/reach or output requires checking; it is **not** approval.
- **`FAILED`**: A stage did not produce a usable result. The original station can still appear in the final ArcGIS review layer.
- **Structural pass**: Geometry/topology checks passed the coded criteria, but manual hydrological review may still be required.
- **Production export authorization**: Kept separate from preview success; the current batch workflow marks outputs as **not production-approved**.

Recommended GIS checks include gauge-to-outlet location, delineated drainage area versus reference area, river connectivity, nested HRU partition consistency, missing/invalid geometries, and exceptionally large snapping distances.

## Reproducibility and contributions

The code is intended to support additional delineation methods under `1_delineation/merit/`, `1_delineation/dem/`, or later method directories, while sharing common catalogue, orchestration, QA, and export conventions. Current TDX scripts and HPC wrappers contain some ARC/CAMELS-SPAT assumptions and must be parameterized for general-purpose distribution. Do not represent them as fully portable until those assumptions have been removed and tested.

Generated datasets are intentionally excluded from Git. Keep raw TDX source files, station observations, geospatial outputs, and Slurm logs outside version control. Include example inputs and licensing information before a public release; do not redistribute TDX data unless its terms permit it.

For detailed, reproducible commands, use **[`docs/TDX_Hydro_Delineation_Manual.md`](docs/TDX_Hydro_Delineation_Manual.md)**.
