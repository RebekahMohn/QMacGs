#SweepFinder Code for running catchup
catchups<-read.delim(file="catchuppops",sep="\t",header=TRUE)
for(i in 1:length(catchups$pop)){
  ##need to add starting point##
  Rtext<-paste("module load VCFtools\n",
               "module load R\n",
               "cd /cwork/ram163/sweepfinder/\n",
               "vcftools --gzvcf /hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP2.vcf.gz --keep ",
               catchups$pop[i], ".pop.L --max-missing 1 --min-alleles 2 --max-alleles 2 --remove-indels --counts2 --chr ", catchups$chr[i],
               " --from-bp ",catchups$pos[i]," --out ",catchups$pop[i],".",catchups$chr[i],"_", catchups$pos[i]," --temp /scratch\n",
		"Rscript /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/20_ToMACROCARPA/sweepfinder/sweepformat.R ",catchups$pop[i],".",catchups$chr[i],"_", catchups$pos[i],".frq.count\n",
		"SweepFinder2 -lg 2000 ",catchups$pop[i],".",catchups$chr[i],"_", catchups$pos[i],".alfreq ",catchups$pop[i],".all.pcf ",catchups$pop[i],".",catchups$chr[i],"_", catchups$pos[i],".pcf.sf2",sep="")
	write(Rtext,paste("SF2_catchups_",i,".cmd",sep=""))      
      #                       
}
