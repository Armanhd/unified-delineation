#!/bin/bash
#SBATCH --job-name=tdx_review
#SBATCH --time=08:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=2
#SBATCH --output=slurm_logs/tdx_review_%j.out
#SBATCH --error=slurm_logs/tdx_review_%j.err

set -euo pipefail
PROJECT=/work/comphyd_lab/users/arman.haddadchi/unified_delineation_code
LABEL="${1:-camels_tdx_full}"
cd "$PROJECT"
# Runs with the existing conda environment Python; no login-shell conda activation needed.
PYTHON=/home/arman.haddadchi/.conda/envs/unified_delineation/bin/python
"$PYTHON" -u 1_delineation/tdx/10_build_arcgis_review_layers.py \
    --root "$PROJECT" \
    --label "$LABEL" \
    --catalog indexes/gauge_catalogs/camels-spat_gauge_catalog_validated.csv \
    --chunk-size 10 \
    --progress-every 25
