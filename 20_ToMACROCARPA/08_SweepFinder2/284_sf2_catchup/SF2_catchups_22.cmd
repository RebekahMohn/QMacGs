module load VCFtools
module load R
cd /cwork/ram163/sweepfinder/
vcftools --gzvcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --keep mac_TX_PLM.pop.L --max-missing 1 --min-alleles 2 --max-alleles 2 --remove-indels --counts2 --chr Chr04 --from-bp 79581705 --out mac_TX_PLM.Chr04_79581705 --temp /scratch
Rscript /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/sweepfinder/sweepformat.R mac_TX_PLM.Chr04_79581705.frq.count
SweepFinder2 -lg 2000 mac_TX_PLM.Chr04_79581705.alfreq mac_TX_PLM.all.pcf mac_TX_PLM.Chr04_79581705.pcf.sf2
