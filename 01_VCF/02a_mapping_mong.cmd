#!/bin/bash
i=1
for file in /cwork/ram163/OAK_GENOMES/02_CLEAN_DATA/02_CHLOROPLAST_FILTERED/S*_1.org_filtered.fq.gz;
do
        fq=`basename ${file} _1.org_filtered.fq.gz`
#	echo "#!/bin/bash" > /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/10_ToMONGOLICA/map_S${i}#
#	echo "#SBATCH --cpus-per-task=20" >> /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/10_ToMONGOLICA/map_S${i}
#	echo "#SBATCH --mem=65G"  >> /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/10_ToMONGOLICA/map_S${i}
#	echo "#SBATCH --mail-type=ALL"  >> /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/10_ToMONGOLICA/map_S${i}
#	echo "#SBATCH --mail-user=rmohn@mortonarb.org"  >> /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/10_ToMONGOLICA/map_S${i}
        echo "module load Bowtie2" > /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/10_ToMONGOLICA/map_S${i}
	echo "module load samtools" >> /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/10_ToMONGOLICA/map_S${i}
	echo "module load bcftools" >> /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/10_ToMONGOLICA/map_S${i}
	echo "bowtie2 -x /hpc/group/manoslab/QMACROCARPA_GENOMES/00_REFERENCES/Qmongolica_genome/Quercus_mongolica_genome -1 /cwork/ram163/OAK_GENOMES/02_CLEAN_DATA/02_CHLOROPLAST_FILTERED/${fq}_1.org_filtered.fq.gz -2 /cwork/ram163/OAK_GENOMES/02_CLEAN_DATA/02_CHLOROPLAST_FILTERED/${fq}_2.org_filtered.fq.gz --end-to-end --no-unal -p 20 | samtools  view -bS | samtools sort -@ 20 -m 3G -o /cwork/ram163/OAK_GENOMES/14_MAPPING/01_BAM_SORTED/${fq}.bam" >> /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/10_ToMONGOLICA/map_S${i}
	echo "cd /cwork/ram163/OAK_GENOMES/14_MAPPING/01_BAM_SORTED/" >> /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/10_ToMONGOLICA/map_S${i}
	echo "bcftools mpileup -f /hpc/group/manoslab/QMACROCARPA_GENOMES/00_REFERENCES/Qmongolica_genome/Quercus_mongolica_genome.fasta ${fq}.bam -g 5 -o /cwork/ram163/OAK_GENOMES/14_MAPPING/02_VARIANTS/${fq}.gvcf.gz -O z --threads 20" >> /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/10_ToMONGOLICA/map_S${i}
        ((i=i+1))
done
