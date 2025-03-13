#!/bin/bash

#SBATCH --job-name=migration_model       # Job name
#SBATCH --cpus-per-task=4               # Use 8 CPUs for one task
#SBATCH --mem=400G                       # Total memory
#SBATCH --nodes=1                        # Request one node
#SBATCH --time=14-00:00:00                # Maximum runtime (14 days)
#SBATCH --output=log_%j.log              # Standard output and error combined
#SBATCH --mail-user=jordi.sola@kaust.edu.sa
#SBATCH --mail-type=ALL

# Increase stack size to avoid recursion issues
ulimit -s unlimited  # Set the stack size to unlimited to avoid Eigen template recursion issues

# Go to your working directory
cd /home/codinajs/Models_new/TB_Model_k4_nospat_unscaled

# Ensure the .R directory exists and set Makevars
mkdir -p /home/codinajs/.R
echo -e "CXXFLAGS += -O2 -ftemplate-depth=3200\nCXX17FLAGS += -O2 -ftemplate-depth=3200" > /home/codinajs/.R/Makevars

# Initialize Conda
source ~/miniconda3/etc/profile.d/conda.sh

# Load environment and activate Conda
conda activate r_env

# Set STAN environment variables
#export R_MAX_NUM_DLLS=999         # Increase the max number of dynamically loaded libraries

# Run the R script
Rscript TB_Model_IBEX_k4_nospat_unscaled.R

# Deactivate the Conda environment
conda deactivate