#!/bin/bash
#SBATCH -e slurm_map_%A_%a.err
#SBATCH -o slurm_map_%A_%a.out
#SBATCH --array=0-900
#SBATCH -c 20
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmohn@mortonarb.org

cd mapcmd/
./map_S${SLURM_ARRAY_TASK_ID}
