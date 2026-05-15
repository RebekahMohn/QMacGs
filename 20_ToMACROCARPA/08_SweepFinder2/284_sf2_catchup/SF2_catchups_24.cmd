module load VCFtools
module load R
cd /cwork/ram163/sweepfinder/
vcftools --gzvcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --keep mac_KY_BBF.pop.L --max-missing 1 --min-alleles 2 --max-alleles 2 --remove-indels --counts2 --chr Chr05 --from-bp 83852777 --out mac_KY_BBF.Chr05_83852777 --temp /scratch
Rscript /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/sweepfinder/sweepformat.R mac_KY_BBF.Chr05_83852777.frq.count
SweepFinder2 -lg 2000 mac_KY_BBF.Chr05_83852777.alfreq mac_KY_BBF.all.pcf mac_KY_BBF.Chr05_83852777.pcf.sf2
