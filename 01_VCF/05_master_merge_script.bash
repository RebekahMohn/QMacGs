#!/bin/bash
#SBATCH -c 20
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmohn@mortonarb.org

cd /cwork/ram163/OAK_GENOMES/14_MAPPING/03_CALLED
bcftools merge -g /hpc/group/manoslab/QMACROCARPA_GENOMES/00_REFERENCES/Qmongolica_genome/Quercus_mongolica_genome.fasta -l samp_list -o ALL_LIST.gvcf.gz -Oz --threads 60
