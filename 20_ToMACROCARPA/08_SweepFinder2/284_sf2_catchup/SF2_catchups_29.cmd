module load VCFtools
module load R
cd /cwork/ram163/sweepfinder/
vcftools --gzvcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --keep mac_TX_CHP.pop.L --max-missing 1 --min-alleles 2 --max-alleles 2 --remove-indels --counts2 --chr Chr05 --from-bp 85628803 --out mac_TX_CHP.Chr05_85628803 --temp /scratch
Rscript /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/sweepfinder/sweepformat.R mac_TX_CHP.Chr05_85628803.frq.count
SweepFinder2 -lg 2000 mac_TX_CHP.Chr05_85628803.alfreq mac_TX_CHP.all.pcf mac_TX_CHP.Chr05_85628803.pcf.sf2
