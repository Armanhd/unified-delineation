#!/bin/bash

# ARC Spack environment
source /work/comphyd_lab/local/modules/spack/2024v5/lmod-init-bash

module unuse "$MODULEPATH"
module use /work/comphyd_lab/local/modules/spack/2024v5/modules/linux-rocky8-x86_64/Core/

module purge

# Compiler and runtime
module load gcc/14.2.0

# Python and geospatial modules will be added
# after checking available ARC modules.

echo "Unified delineation HPC module environment initialized."
module list