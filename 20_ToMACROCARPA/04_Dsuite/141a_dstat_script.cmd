#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmohn@mortonarb.org

Dsuite Dinvestigate -w 1000,750 /cwork/ram163/OAK_GENOMES/24_MAPPING/MAC99SP5_MACREF.gvcf.gz /hpc/group/manoslab/QMACROCARPA_GENOMES/00a_SUBSETS/SetsGR_mac2.tsv /hpc/group/manoslab/QMACROCARPA_GENOMES/00a_SUBSETS/trios_alb.txt
Dsuite Dinvestigate -w 1000,750 /cwork/ram163/OAK_GENOMES/24_MAPPING/MAC99SP5_MACREF.gvcf.gz /hpc/group/manoslab/QMACROCARPA_GENOMES/00a_SUBSETS/SetsGR_mac2.tsv /hpc/group/manoslab/QMACROCARPA_GENOMES/00a_SUBSETS/trios_lob.txt
Dsuite Dinvestigate -w 1000,750 /cwork/ram163/OAK_GENOMES/24_MAPPING/MAC99SP5_MACREF.gvcf.gz /hpc/group/manoslab/QMACROCARPA_GENOMES/00a_SUBSETS/SetsGR_mac2.tsv /hpc/group/manoslab/QMACROCARPA_GENOMES/00a_SUBSETS/trios_mue.txt
Dsuite Dinvestigate -w 1000,750 /cwork/ram163/OAK_GENOMES/24_MAPPING/MAC99SP5_MACREF.gvcf.gz /hpc/group/manoslab/QMACROCARPA_GENOMES/00a_SUBSETS/SetsGR_mac2.tsv /hpc/group/manoslab/QMACROCARPA_GENOMES/00a_SUBSETS/trios_ste.txt
Dsuite Dinvestigate -w 1000,750 /cwork/ram163/OAK_GENOMES/24_MAPPING/MAC99SP5_MACREF.gvcf.gz /hpc/group/manoslab/QMACROCARPA_GENOMES/00a_SUBSETS/SetsGR_mac2.tsv /hpc/group/manoslab/QMACROCARPA_GENOMES/00a_SUBSETS/trios_lyr.txt
