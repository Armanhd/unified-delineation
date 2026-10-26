#!/bin/bash
#SBATCH --job-name=setup_delineation
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/setup_conda_%j.out
#SBATCH --error=logs/setup_conda_%j.err

set -e

ROOT=/work/comphyd_lab/users/arman.haddadchi/unified_delineation_code

cd "$ROOT"

# Load ARC Conda
module load conda/base

# Create environment if it does not exist
if ! conda env list | grep -q '^unified_delineation '; then
    conda env create \
        --solver libmamba \
        -f environments/unified_delineation.yml
fi

# Activate environment
eval "$(conda shell.bash hook)"
conda activate unified_delineation

# Verify
python --version

python -c "
import geopandas
import rasterio
import shapely
import pyproj
import fiona
import numpy
import pandas
import tomllib
print('Unified delineation environment ready.')
"