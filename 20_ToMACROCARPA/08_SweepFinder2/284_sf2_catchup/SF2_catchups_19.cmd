module load VCFtools
module load R
cd /cwork/ram163/sweepfinder/
vcftools --gzvcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --keep mac_MA_NMK.pop.L --max-missing 1 --min-alleles 2 --max-alleles 2 --remove-indels --counts2 --chr Chr04 --from-bp 72749539 --out mac_MA_NMK.Chr04_72749539 --temp /scratch
Rscript /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/sweepfinder/sweepformat.R mac_MA_NMK.Chr04_72749539.frq.count
SweepFinder2 -lg 2000 mac_MA_NMK.Chr04_72749539.alfreq mac_MA_NMK.all.pcf mac_MA_NMK.Chr04_72749539.pcf.sf2
