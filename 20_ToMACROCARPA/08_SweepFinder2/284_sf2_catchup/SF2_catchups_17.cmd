module load VCFtools
module load R
cd /cwork/ram163/sweepfinder/
vcftools --gzvcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --keep mac_WI_RLF.pop.L --max-missing 1 --min-alleles 2 --max-alleles 2 --remove-indels --counts2 --chr Chr02 --from-bp 60735548 --out mac_WI_RLF.Chr02_60735548 --temp /scratch
Rscript /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/sweepfinder/sweepformat.R mac_WI_RLF.Chr02_60735548.frq.count
SweepFinder2 -lg 2000 mac_WI_RLF.Chr02_60735548.alfreq mac_WI_RLF.all.pcf mac_WI_RLF.Chr02_60735548.pcf.sf2
