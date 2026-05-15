module load VCFtools
module load R
cd /cwork/ram163/sweepfinder/
vcftools --gzvcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --keep mac_KS_SMP.pop.L --max-missing 1 --min-alleles 2 --max-alleles 2 --remove-indels --counts2 --chr Chr02 --from-bp 60893553 --out mac_KS_SMP.Chr02_60893553 --temp /scratch
Rscript /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/sweepfinder/sweepformat.R mac_KS_SMP.Chr02_60893553.frq.count
SweepFinder2 -lg 2000 mac_KS_SMP.Chr02_60893553.alfreq mac_KS_SMP.all.pcf mac_KS_SMP.Chr02_60893553.pcf.sf2
