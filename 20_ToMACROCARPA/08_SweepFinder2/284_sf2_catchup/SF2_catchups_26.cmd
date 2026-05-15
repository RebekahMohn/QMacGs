module load VCFtools
module load R
cd /cwork/ram163/sweepfinder/
vcftools --gzvcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --keep mac_OH_DPS.pop.L --max-missing 1 --min-alleles 2 --max-alleles 2 --remove-indels --counts2 --chr Chr05 --from-bp 87074824 --out mac_OH_DPS.Chr05_87074824 --temp /scratch
Rscript /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/sweepfinder/sweepformat.R mac_OH_DPS.Chr05_87074824.frq.count
SweepFinder2 -lg 2000 mac_OH_DPS.Chr05_87074824.alfreq mac_OH_DPS.all.pcf mac_OH_DPS.Chr05_87074824.pcf.sf2
