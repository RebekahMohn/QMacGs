module load bcftools
module load IQ-TREE
bcftools view -M2 -i 'TYPE!="indel" & QUAL>10 & FMT/DP>8 & FMT/DP<60' -o /cwork/ram163/GENOME_REGIONS/all/4617_Superscaffold11_19800001-19950001.vcf -S /hpc/group/manoslab/QMACROCARPA_GENOMES/00a_SUBSETS/Subset5.2Species99.tsv /hpc/group/manoslab/QMACROCARPA_GENOMES/14_VARIANTS_MONGOLICA/ALL_LIST.gvcf.gz -r Superscaffold11:19800001-19950001
cd /cwork/ram163/GENOME_REGIONS/all/
vcf2phylip.py -i /cwork/ram163/GENOME_REGIONS/all/4617_Superscaffold11_19800001-19950001.vcf --output-prefix 4617_Superscaffold11_19800001-19950001 -m 20 -f
iqtree2 -s 4617_Superscaffold11_19800001-19950001.min20.fasta -nt 5 -pre 4617_Superscaffold11_19800001-19950001 -B 1000 -mem 10G
 pxrr -t 4617_Superscaffold11_19800001-19950001.treefile -g SRR24768763.bam,SRR24768723.bam,SRR24768754.bam,SRR24768716.bam,SRR24768736.bam -o /cwork/ram163/GENOME_REGIONS/all/rooted/4617_Superscaffold11_19800001-19950001.rooted.tre
