module load VCFtools
module load R
cd /cwork/ram163/sweepfinder/
vcftools --gzvcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --keep mac_MB_SWP.pop.L --max-missing 1 --min-alleles 2 --max-alleles 2 --remove-indels --counts2 --chr Chr05 --from-bp 82808762 --out mac_MB_SWP.Chr05_82808762 --temp /scratch
Rscript /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/sweepfinder/sweepformat.R mac_MB_SWP.Chr05_82808762.frq.count
SweepFinder2 -lg 2000 mac_MB_SWP.Chr05_82808762.alfreq mac_MB_SWP.all.pcf mac_MB_SWP.Chr05_82808762.pcf.sf2
