module load VCFtools
module load R
cd /cwork/ram163/sweepfinder/
vcftools --gzvcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --keep mac_TX_PLM.pop.L --max-missing 1 --min-alleles 2 --max-alleles 2 --remove-indels --counts2 --chr Chr02 --from-bp 60285535 --out mac_TX_PLM.Chr02_60285535 --temp /scratch
Rscript /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/sweepfinder/sweepformat.R mac_TX_PLM.Chr02_60285535.frq.count
SweepFinder2 -lg 2000 mac_TX_PLM.Chr02_60285535.alfreq mac_TX_PLM.all.pcf mac_TX_PLM.Chr02_60285535.pcf.sf2
