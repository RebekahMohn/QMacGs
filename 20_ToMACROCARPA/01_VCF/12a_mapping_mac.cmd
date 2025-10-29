#!/bin/bash
i=1
for file in /cwork/ram163/OAK_GENOMES/02_CLEAN_DATA/02_CHLOROPLAST_FILTERED/S*_1.org_filtered.fq.gz;
do
    fq=`basename ${file} _1.org_filtered.fq.gz`
    echo "module load Bowtie2" > /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/mapcmd/map_S${i}
    echo "module load samtools" >> /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/mapcmd/map_S${i}
    echo "module load bcftools" >> /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/mapcmd/map_S${i}#
    echo "bowtie2 -x /hpc/group/manoslab/QMACROCARPA_GENOMES/00_REFERENCES/Qmacrocarpa_genome/quercus_macrocarpa_hap1 -1 /cwork/ram163/OAK_GENOMES/02_CLEAN_DATA/02_CHLOROPLAST_FILTERED/${fq}_1.org_filtered.fq.gz -2 /cwork/ram163/OAK_GENOMES/02_CLEAN_DATA/02_CHLOROPLAST_FILTERED/${fq}_2.org_filtered.fq.gz --end-to-end --no-unal -p 20 | samtools  view -bS | samtools sort -@ 20 -m 5G -o /cwork/ram163/OAK_GENOMES/24_MAPPING/01_BAM_SORTED/${fq}.bam" >> /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/mapcmd/map_S${i}
    echo "cd /cwork/ram163/OAK_GENOMES/24_MAPPING/01_BAM_SORTED/" >> /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/mapcmd/map_S${i}
    echo "bcftools mpileup -f /hpc/group/manoslab/QMACROCARPA_GENOMES/00_REFERENCES/Qmacrocarpa_genome/quercus_macrocarpa_hap1.fa ${fq}.bam -g 5 -Oz --threads 20 | bcftools call -Oz -m -g 5 | bcftools norm -f /hpc/group/manoslab/QMACROCARPA_GENOMES/00_REFERENCES/Qmacrocarpa_genome/quercus_macrocarpa_hap1.fa -Oz -o /cwork/ram163/OAK_GENOMES/24_MAPPING/03_CALLED/${fq}_norm.gvcf.gz --threads 20" >> /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/mapcmd/map_S${i}
	((i=i+1))
done

