module load VCFtools
module load R
cd /cwork/ram163/sweepfinder/
vcftools --gzvcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --keep mac_OK_NOW.pop.L --max-missing 1 --min-alleles 2 --max-alleles 2 --remove-indels --counts2 --chr Chr02 --from-bp 66257717 --out mac_OK_NOW.Chr02_66257717 --temp /scratch
Rscript /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/sweepfinder/sweepformat.R mac_OK_NOW.Chr02_66257717.frq.count
SweepFinder2 -lg 2000 mac_OK_NOW.Chr02_66257717.alfreq mac_OK_NOW.all.pcf mac_OK_NOW.Chr02_66257717.pcf.sf2
