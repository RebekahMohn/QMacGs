#!/bin/bash
#SBATCH -e slurm_map_%A_%a.err
#SBATCH -o slurm_map_%A_%a.out 
#SBATCH --array=721
#SBATCH -c 20
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmohn@mortonarb.org

#cd /cwork/ram163/OAK_GENOMES/14_MAPPING/
./index2_S${SLURM_ARRAY_TASK_ID}
