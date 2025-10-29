#!/bin/bash
#SBATCH -e slurm_map_%A_%a.err
#SBATCH -o slurm_map_%A_%a.out
#SBATCH --array=0-15305
#SBATCH -c 1
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmohn@mortonarb.org

#cd phyCMDs
./phy${SLURM_ARRAY_TASK_ID}.cmd