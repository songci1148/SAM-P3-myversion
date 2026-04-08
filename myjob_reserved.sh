#!/bin/bash
#
#SBATCH -A naiss2025-1-5
#SBATCH -J SAM
#SBATCH -t 03:30:00
#SBATCH -N 2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ci.song@misu.su.se
#SBARCH --reservation=bt-myres
#
# Run a single task in the foreground.
srun -n 64 ./SAM_ADV_UM5_SGS_TKE_RAD_RRTM_MICRO_P3MULTI_tke_frzcat_diag
#
# Script ends here
