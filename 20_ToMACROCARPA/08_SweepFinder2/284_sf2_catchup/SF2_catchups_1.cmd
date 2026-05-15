module load VCFtools
module load R
cd /cwork/ram163/sweepfinder/
vcftools --gzvcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --keep mac_IA_BSP.pop.L --max-missing 1 --min-alleles 2 --max-alleles 2 --remove-indels --counts2 --chr Chr02 --from-bp 73241929 --out mac_IA_BSP.Chr02_73241929 --temp /scratch
Rscript /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/sweepfinder/sweepformat.R mac_IA_BSP.Chr02_73241929.frq.count
SweepFinder2 -lg 2000 mac_IA_BSP.Chr02_73241929.alfreq mac_IA_BSP.all.pcf mac_IA_BSP.Chr02_73241929.pcf.sf2
