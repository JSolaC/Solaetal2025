#!/bin/bash

#SBATCH --job-name=migration_model       # Job name
#SBATCH --cpus-per-task=12               # Use 8 CPUs for one task
#SBATCH --mem=200G                       # Total memory
#SBATCH --nodes=1                        # Request one node
#SBATCH --time=7-00:00:00                # Maximum runtime (1 day)
#SBATCH --output=log_%j.log              # Standard output and error combined
#SBATCH --mail-user=jordi.sola@kaust.edu.sa
#SBATCH --mail-type=ALL

# Go to your working directory
cd /home/codinajs/Models_new/MI_Model_k4_nospat_unscaled_filtered

# Ensure the .R directory exists and set Makevars
mkdir -p /home/codinajs/.R
echo -e "CXXFLAGS += -O2 -ftemplate-depth=1600\nCXX17FLAGS += -O2 -ftemplate-depth=1600" > /home/codinajs/.R/Makevars

# Initialize Conda
source ~/miniconda3/etc/profile.d/conda.sh

# Load environment and activate Conda
conda activate r_env

# Set STAN environment variables
export STAN_NUM_THREADS=12  # This controls threads within each chain
export OMP_NUM_THREADS=1   # Helps with BLAS parallelization if used

# Run the R script
Rscript MI_Model_IBEX_k4_nospat_unscaled_filtered.R

# Deactivate the Conda environment
conda deactivate