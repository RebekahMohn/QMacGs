module load VCFtools
module load R
cd /cwork/ram163/sweepfinder/
vcftools --gzvcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --keep mac_IN_BOW.pop.L --max-missing 1 --min-alleles 2 --max-alleles 2 --remove-indels --counts2 --chr Chr05 --from-bp 82238754 --out mac_IN_BOW.Chr05_82238754 --temp /scratch
Rscript /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/sweepfinder/sweepformat.R mac_IN_BOW.Chr05_82238754.frq.count
SweepFinder2 -lg 2000 mac_IN_BOW.Chr05_82238754.alfreq mac_IN_BOW.all.pcf mac_IN_BOW.Chr05_82238754.pcf.sf2
