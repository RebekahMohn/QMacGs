#!/bin/bash
#SBATCH --mem=607727
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmohn@mortonarb.org
#SBATCH --time=1-10:00:00
#Take Mac99 file and filter by depth and to only snps

module load bcftools
bcftools view -S /hpc/group/manoslab/QMACROCARPA_GENOMES/00a_SUBSETS/pure_Mac99a.tsv --threads 4 /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99SP5_MACREF.gvcf.gz | bcftools filter -Oz -i 'TYPE="snp" & QUAL>10 & FMT/DP>8 & FMT/DP<100' --threads 4 -o /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz

module load R/4.4.3
module load /opt/apps/modules-bak/R/4.5.1

#run rscript that will
##1. convert ped into lfmm
##2. read in the env. data
##3. run lfmm

##sample SYST-MOR-0007049 was missing a lot of data so was removed from this analysis
bcftools view -s SYST-MOR-0007049.bam --threads 4 /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP.vcf.gz -o /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz

module load /opt/apps/modules-bak/R/4.5.1

#Convert that file into ped
plink --vcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --recode12 --chr Chr01 --geno 0 --maf 0.01 --allow-extra-chr --out /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP_C1_maf01
Rscript gea_C1.R

plink --vcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --recode12 --chr Chr02  --geno 0 --maf 0.01 --allow-extra-chr --out /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP_C2_maf01
Rscript gea_C2.R

plink --vcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --recode12 --chr Chr03 --geno 0 --maf 0.01 --allow-extra-chr --out /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP_C3_maf01
Rscript gea_C3.R

plink --vcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --recode12 --chr Chr04 --geno 0 --maf 0.01 --allow-extra-chr --out /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP_C4_maf01
Rscript gea_C4.R

plink --vcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --recode12 --chr Chr05 --geno 0 --maf 0.01 --allow-extra-chr --out /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP_C5_maf01
Rscript gea_C5.R

plink --vcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --recode12 --chr Chr06 --geno 0 --maf 0.01 --allow-extra-chr --out /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP_C6_maf01
Rscript gea_C6.R

plink --vcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --recode12 --chr Chr07 --geno 0 --maf 0.01 --allow-extra-chr --out /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP_C7_maf01
Rscript gea_C7.R

plink --vcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --recode12 --chr Chr08 --geno 0 --maf 0.01 --allow-extra-chr --out /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP_C8_maf01
Rscript gea_C8.R

plink --vcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --recode12 --chr Chr09 --geno 0 --maf 0.01 --allow-extra-chr --out /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP_C9_maf01
Rscript gea_C9.R

plink --vcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --recode12 --chr Chr10 --geno 0 --maf 0.01 --allow-extra-chr --out /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP_C10_maf01
Rscript gea_C10.R

plink --vcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --recode12 --chr Chr11 --geno 0 --maf 0.01 --allow-extra-chr --out /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP_C11_maf01
Rscript gea_C11.R

plink --vcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --recode12 --chr Chr12 --geno 0 --maf 0.01 --allow-extra-chr --out /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP_C12_maf01
Rscript gea_C12.R
