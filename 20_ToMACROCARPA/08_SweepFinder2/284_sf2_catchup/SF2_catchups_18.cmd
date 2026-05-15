module load VCFtools
module load R
cd /cwork/ram163/sweepfinder/
vcftools --gzvcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --keep mac_KY_GRF.pop.L --max-missing 1 --min-alleles 2 --max-alleles 2 --remove-indels --counts2 --chr Chr04 --from-bp 76795638 --out mac_KY_GRF.Chr04_76795638 --temp /scratch
Rscript /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/sweepfinder/sweepformat.R mac_KY_GRF.Chr04_76795638.frq.count
SweepFinder2 -lg 2000 mac_KY_GRF.Chr04_76795638.alfreq mac_KY_GRF.all.pcf mac_KY_GRF.Chr04_76795638.pcf.sf2
