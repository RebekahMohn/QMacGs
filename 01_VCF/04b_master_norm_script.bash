#!/bin/bash
#SBATCH -e slurm_map_%A_%a.err
#SBATCH -o slurm_map_%A_%a.out 
#SBATCH --array=1-721
#SBATCH -c 5
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmohn@mortonarb.org

#cd /cwork/ram163/OAK_GENOMES/14_MAPPING/
./norm_S${SLURM_ARRAY_TASK_ID}
