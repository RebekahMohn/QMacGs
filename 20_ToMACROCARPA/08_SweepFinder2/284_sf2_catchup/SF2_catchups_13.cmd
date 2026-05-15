module load VCFtools
module load R
cd /cwork/ram163/sweepfinder/
vcftools --gzvcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --keep mac_SD_CSP.pop.L --max-missing 1 --min-alleles 2 --max-alleles 2 --remove-indels --counts2 --chr Chr02 --from-bp 75425996 --out mac_SD_CSP.Chr02_75425996 --temp /scratch
Rscript /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/sweepfinder/sweepformat.R mac_SD_CSP.Chr02_75425996.frq.count
SweepFinder2 -lg 2000 mac_SD_CSP.Chr02_75425996.alfreq mac_SD_CSP.all.pcf mac_SD_CSP.Chr02_75425996.pcf.sf2
