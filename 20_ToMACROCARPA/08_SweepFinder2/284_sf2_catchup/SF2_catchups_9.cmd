module load VCFtools
module load R
cd /cwork/ram163/sweepfinder/
vcftools --gzvcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --keep mac_MI_PCP.pop.L --max-missing 1 --min-alleles 2 --max-alleles 2 --remove-indels --counts2 --chr Chr02 --from-bp 59719518 --out mac_MI_PCP.Chr02_59719518 --temp /scratch
Rscript /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/sweepfinder/sweepformat.R mac_MI_PCP.Chr02_59719518.frq.count
SweepFinder2 -lg 2000 mac_MI_PCP.Chr02_59719518.alfreq mac_MI_PCP.all.pcf mac_MI_PCP.Chr02_59719518.pcf.sf2
