module load bcftools
module load IQ-TREE
bcftools view -M2 -i 'TYPE!="indel" & QUAL>10 & FMT/DP>8 & FMT/DP<60' -o /cwork/ram163/GENOME_REGIONS/all/1633_Superscaffold4_24450001-24600001.vcf -S /hpc/group/manoslab/QMACROCARPA_GENOMES/00a_SUBSETS/Subset5.2Species99.tsv /hpc/group/manoslab/QMACROCARPA_GENOMES/14_VARIANTS_MONGOLICA/ALL_LIST.gvcf.gz -r Superscaffold4:24450001-24600001
cd /cwork/ram163/GENOME_REGIONS/all/
vcf2phylip.py -i /cwork/ram163/GENOME_REGIONS/all/1633_Superscaffold4_24450001-24600001.vcf --output-prefix 1633_Superscaffold4_24450001-24600001 -m 20 -f
iqtree2 -s 1633_Superscaffold4_24450001-24600001.min20.fasta -nt 5 -pre 1633_Superscaffold4_24450001-24600001 -B 1000 -mem 10G
 pxrr -t 1633_Superscaffold4_24450001-24600001.treefile -g SRR24768763.bam,SRR24768723.bam,SRR24768754.bam,SRR24768716.bam,SRR24768736.bam -o /cwork/ram163/GENOME_REGIONS/all/rooted/1633_Superscaffold4_24450001-24600001.rooted.tre
