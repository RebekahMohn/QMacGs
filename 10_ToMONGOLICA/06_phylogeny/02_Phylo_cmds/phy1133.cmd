module load bcftools
module load IQ-TREE
bcftools view -M2 -i 'TYPE!="indel" & QUAL>10 & FMT/DP>8 & FMT/DP<60' -o /cwork/ram163/GENOME_REGIONS/all/1133_Superscaffold3_16650001-16800001.vcf -S /hpc/group/manoslab/QMACROCARPA_GENOMES/00a_SUBSETS/Subset5.2Species99.tsv /hpc/group/manoslab/QMACROCARPA_GENOMES/14_VARIANTS_MONGOLICA/ALL_LIST.gvcf.gz -r Superscaffold3:16650001-16800001
cd /cwork/ram163/GENOME_REGIONS/all/
vcf2phylip.py -i /cwork/ram163/GENOME_REGIONS/all/1133_Superscaffold3_16650001-16800001.vcf --output-prefix 1133_Superscaffold3_16650001-16800001 -m 20 -f
iqtree2 -s 1133_Superscaffold3_16650001-16800001.min20.fasta -nt 5 -pre 1133_Superscaffold3_16650001-16800001 -B 1000 -mem 10G
 pxrr -t 1133_Superscaffold3_16650001-16800001.treefile -g SRR24768763.bam,SRR24768723.bam,SRR24768754.bam,SRR24768716.bam,SRR24768736.bam -o /cwork/ram163/GENOME_REGIONS/all/rooted/1133_Superscaffold3_16650001-16800001.rooted.tre
