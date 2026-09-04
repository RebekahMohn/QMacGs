module load bcftools
module load IQ-TREE
bcftools view -M2 -i 'TYPE!="indel" & QUAL>10 & FMT/DP>8 & FMT/DP<60' -o /cwork/ram163/GENOME_REGIONS/all/937_Superscaffold2_88500001-88650001.vcf -S /hpc/group/manoslab/QMACROCARPA_GENOMES/00a_SUBSETS/Subset5.2Species99.tsv /hpc/group/manoslab/QMACROCARPA_GENOMES/14_VARIANTS_MONGOLICA/ALL_LIST.gvcf.gz -r Superscaffold2:88500001-88650001
cd /cwork/ram163/GENOME_REGIONS/all/
vcf2phylip.py -i /cwork/ram163/GENOME_REGIONS/all/937_Superscaffold2_88500001-88650001.vcf --output-prefix 937_Superscaffold2_88500001-88650001 -m 20 -f
iqtree2 -s 937_Superscaffold2_88500001-88650001.min20.fasta -nt 5 -pre 937_Superscaffold2_88500001-88650001 -B 1000 -mem 10G
 pxrr -t 937_Superscaffold2_88500001-88650001.treefile -g SRR24768763.bam,SRR24768723.bam,SRR24768754.bam,SRR24768716.bam,SRR24768736.bam -o /cwork/ram163/GENOME_REGIONS/all/rooted/937_Superscaffold2_88500001-88650001.rooted.tre
