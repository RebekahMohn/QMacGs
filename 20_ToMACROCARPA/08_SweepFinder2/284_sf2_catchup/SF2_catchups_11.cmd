module load VCFtools
module load R
cd /cwork/ram163/sweepfinder/
vcftools --gzvcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --keep mac_OH_DPS.pop.L --max-missing 1 --min-alleles 2 --max-alleles 2 --remove-indels --counts2 --chr Chr02 --from-bp 67349750 --out mac_OH_DPS.Chr02_67349750 --temp /scratch
Rscript /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/sweepfinder/sweepformat.R mac_OH_DPS.Chr02_67349750.frq.count
SweepFinder2 -lg 2000 mac_OH_DPS.Chr02_67349750.alfreq mac_OH_DPS.all.pcf mac_OH_DPS.Chr02_67349750.pcf.sf2
