module load VCFtools
module load R
cd /cwork/ram163/sweepfinder/
vcftools --gzvcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --keep mac_IL_CHB.pop.L --max-missing 1 --min-alleles 2 --max-alleles 2 --remove-indels --counts2 --chr Chr02 --from-bp 64605666 --out mac_IL_CHB.Chr02_64605666 --temp /scratch
Rscript /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/sweepfinder/sweepformat.R mac_IL_CHB.Chr02_64605666.frq.count
SweepFinder2 -lg 2000 mac_IL_CHB.Chr02_64605666.alfreq mac_IL_CHB.all.pcf mac_IL_CHB.Chr02_64605666.pcf.sf2
