module load VCFtools
module load R
cd /cwork/ram163/sweepfinder/
vcftools --gzvcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --keep mac_TN_BCP.pop.L --max-missing 1 --min-alleles 2 --max-alleles 2 --remove-indels --counts2 --chr Chr05 --from-bp 77248683 --out mac_TN_BCP.Chr05_77248683 --temp /scratch
Rscript /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/sweepfinder/sweepformat.R mac_TN_BCP.Chr05_77248683.frq.count
SweepFinder2 -lg 2000 mac_TN_BCP.Chr05_77248683.alfreq mac_TN_BCP.all.pcf mac_TN_BCP.Chr05_77248683.pcf.sf2
