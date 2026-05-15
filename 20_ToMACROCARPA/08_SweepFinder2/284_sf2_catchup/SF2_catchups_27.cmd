module load VCFtools
module load R
cd /cwork/ram163/sweepfinder/
vcftools --gzvcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --keep mac_OK_NOW.pop.L --max-missing 1 --min-alleles 2 --max-alleles 2 --remove-indels --counts2 --chr Chr05 --from-bp 82118752 --out mac_OK_NOW.Chr05_82118752 --temp /scratch
Rscript /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/sweepfinder/sweepformat.R mac_OK_NOW.Chr05_82118752.frq.count
SweepFinder2 -lg 2000 mac_OK_NOW.Chr05_82118752.alfreq mac_OK_NOW.all.pcf mac_OK_NOW.Chr05_82118752.pcf.sf2
