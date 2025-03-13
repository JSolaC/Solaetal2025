#!/bin/bash

#SBATCH --job-name=migration_model       # Job name
#SBATCH --cpus-per-task=12               # Use 12 CPUs for one task
#SBATCH --mem=300G                       # Total memory
#SBATCH --nodes=1                        # Request one node
#SBATCH --time=14-00:00:00                # Maximum runtime (1 day)
#SBATCH --output=log_%j.log              # Standard output and error combined
#SBATCH --mail-user=jordi.sola@kaust.edu.sa
#SBATCH --mail-type=ALL

# Go to your working directory
cd /home/codinajs/Models_new/MF_Model_nonspat_clustered_unscaled_filtered

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
Rscript MF_Model_IBEX_nonspat_clustered_unscaled_filtered.R

# Deactivate the Conda environment
conda deactivate