#!/bin/bash
i=1
for file in /cwork/ram163/OAK_GENOMES/14_MAPPING/03_CALLED/*.vcf.gz;
do
        fq=`basename ${file} .vcf.gz`
	echo "cd /cwork/ram163/OAK_GENOMES/14_MAPPING/03_CALLED/" > /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/10_ToMONGOLICA/norm_S${i}
	echo "bcftools norm ${file} -f /hpc/group/manoslab/QMACROCARPA_GENOMES/00_REFERENCES/Qmongolica_genome/Quercus_mongolica_genome.fasta -Oz -o ${fq}_norm.vcf.gz --threads 5" >> /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/10_ToMONGOLICA/norm_S${i}
	echo "bcftools index /cwork/ram163/OAK_GENOMES/14_MAPPING/03_CALLED/${fq}_norm.vcf.gz --threads 5 " >> /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/10_ToMONGOLICA/norm_S${i}
        ((i=i+1))
done
