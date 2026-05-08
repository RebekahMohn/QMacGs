#!/bin/bash
#SBATCH -e slurm_%A_%a.err
#SBATCH -o slurm_%A_%a.out 
#SBATCH --array=681-693
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmohn@mortonarb.org

cd /cwork/ram163/OAK_GENOMES/02_CLEAN_DATA/00_SCRIPTS/
./trim_filter_S${SLURM_ARRAY_TASK_ID}
