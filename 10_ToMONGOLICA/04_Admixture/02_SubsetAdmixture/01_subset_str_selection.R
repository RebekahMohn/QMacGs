meta<-read.csv("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/metadata_sequenced_taxonomyfixed_outgroup.csv")
simp_meta<-data.frame(cbind(Seq=meta$Seq,colnum=meta$collectionNumber,species=meta$GenomicVoucherID_USE_THIS,state=meta$state.provIndIowanace,site=meta$site,lat=meta$latitude.orig,long=meta$longitude.orig))
simp_meta$Seq.bam<-paste(simp_meta$Seq,".bam",sep="")
samp_list<-c("SYST-MOR-0007212.bam", "SYST-MOR-0007185.bam",
             "SYST-MOR-0006994.bam", "SYST-MOR-0006970.bam","SYST-MOR-0007196.bam", "SYST-MOR-0006925.bam",
             "SYST-MOR-0007134.bam", "SYST-MOR-0006768.bam","SYST-MOR-0006754.bam", "SYST-MOR-0006760.bam",
             "SRR14632961.bam", "SRR14632950.bam","SRR14632932.bam", "SRR14632928.bam",
             "SRR14632968.bam", "SYST-MOR-0007046.bam","SYST-MOR-0007215.bam", "SYST-MOR-0006787.bam",
             "SYST-MOR-0007205.bam", "SYST-MOR-0007015.bam", "SYST-MOR-0006717.bam",
             "SYST-MOR-0007154.bam", "SYST-MOR-0006924.bam","SYST-MOR-0006718.bam", "SYST-MOR-0006773.bam",
             "SYST-MOR-0007080.bam", "SYST-MOR-0007061.bam","SYST-MOR-0006751.bam","SYST-MOR-0007120.bam",
             "SYST-MOR-0006984.bam","SYST-MOR-0007094.bam","SYST-MOR-0007240.bam","SYST-MOR-0007239.bam",
             "SYST-MOR-0007090.bam","SYST-MOR-0007092.bam","SYST-MOR-0006865.bam","SYST-MOR-0006742.bam",
             "SYST-MOR-0006886.bam","SYST-MOR-0006749.bam","SYST-MOR-0007173.bam","SYST-MOR-0007166.bam",
             "SYST-MOR-0007167.bam","SYST-MOR-0007160.bam","SYST-MOR-0007162.bam","SYST-MOR-0007178.bam")
simp_meta2<-subset(simp_meta,!Seq.bam %in% samp_list)
simp_meta3<-subset(simp_meta,!Seq.bam %in% samp_list)

for(i in 1:22){
  Seq4Subset<-sample(subset(simp_meta2,species=="Quercus macrocarpa Michx.")$Seq,15)
  simp_meta2<-subset(simp_meta2,!Seq %in% Seq4Subset)
  Seq4Subset<-c(Seq4Subset,sample(subset(simp_meta2,species=="Quercus alba L.")$Seq,15))
  Seq4Subset<-c(Seq4Subset,sample(subset(simp_meta2,species=="Quercus muehlenbergii Engelm.")$Seq,15))
  Seq4Subset<-c(Seq4Subset,sample(subset(simp_meta2,species=="Quercus bicolor Willd.")$Seq,15))
  Seq4Subset<-c(Seq4Subset,subset(simp_meta,species=="Quercus prinoides Willd.")$Seq)
  Seq4Subset<-c(Seq4Subset,sample(subset(simp_meta2,species=="Quercus stellata Wangenh.")$Seq,15))
  Seq4Subset<-c(Seq4Subset,subset(simp_meta2,species=="Quercus lyrata Walt.")$Seq)
  Seq4Subset<-c(Seq4Subset,subset(simp_meta2,species=="Quercus michauxii Nutt.")$Seq)
  Seq4Subset<-c(Seq4Subset,subset(simp_meta2,species=="Quercus montana Willd.")$Seq)
  Seq4Subset<-c(Seq4Subset,subset(simp_meta,species=="Quercus sinuata Walt.")$Seq)
  Seq4Subset<-c(Seq4Subset,subset(simp_meta2,species=="Quercus lobata  Nee.")$Seq)
  Seq4Subset<-c(Seq4Subset,subset(simp_meta2,species=="Quercus margaretiae Ashe ex Small")$Seq)
  Seq4Subset_samp<-paste(Seq4Subset,".bam",sep="")
  Seq4Subset_samp<-c(Seq4Subset_samp,samp_list)
  # simp_meta2<-subset(simp_meta2,!Seq.bam %in% Seq4Subset_samp)
  write.csv(Seq4Subset_samp,file=paste("C:/Users/rmohn/Desktop/10_Analysis/020_STR/Subset20perSpecies_",i,".csv",sep=""),row.names = FALSE,col.names = NA)
}

i=23
Seq4Subset<-subset(simp_meta2,species=="Quercus macrocarpa Michx.")$Seq
Seq4Subset<-c(Seq4Subset,sample(subset(simp_meta3,species=="Quercus macrocarpa Michx.")$Seq,8))
#simp_meta2<-subset(simp_meta2,!Seq %in% Seq4Subset)
Seq4Subset<-c(Seq4Subset,sample(subset(simp_meta2,species=="Quercus alba L.")$Seq,15))
Seq4Subset<-c(Seq4Subset,sample(subset(simp_meta2,species=="Quercus muehlenbergii Engelm.")$Seq,15))
Seq4Subset<-c(Seq4Subset,sample(subset(simp_meta2,species=="Quercus bicolor Willd.")$Seq,15))
Seq4Subset<-c(Seq4Subset,subset(simp_meta,species=="Quercus prinoides Willd.")$Seq)
Seq4Subset<-c(Seq4Subset,sample(subset(simp_meta2,species=="Quercus stellata Wangenh.")$Seq,15))
Seq4Subset<-c(Seq4Subset,subset(simp_meta2,species=="Quercus lyrata Walt.")$Seq)
Seq4Subset<-c(Seq4Subset,subset(simp_meta2,species=="Quercus michauxii Nutt.")$Seq)
Seq4Subset<-c(Seq4Subset,subset(simp_meta2,species=="Quercus montana Willd.")$Seq)
Seq4Subset<-c(Seq4Subset,subset(simp_meta,species=="Quercus sinuata Walt.")$Seq)
Seq4Subset<-c(Seq4Subset,subset(simp_meta2,species=="Quercus lobata  Nee.")$Seq)
Seq4Subset<-c(Seq4Subset,subset(simp_meta2,species=="Quercus margaretiae Ashe ex Small")$Seq)
Seq4Subset_samp<-paste(Seq4Subset,".bam",sep="")
Seq4Subset_samp<-c(Seq4Subset_samp,samp_list)
# simp_meta2<-subset(simp_meta2,!Seq.bam %in% Seq4Subset_samp)
write.csv(Seq4Subset_samp,file=paste("C:/Users/rmohn/Desktop/10_Analysis/020_STR/Subset20perSpecies_",i,".csv",sep=""),row.names = FALSE,col.names = NA)

for(i in 1:23){
  cat(paste("vcftools --gzvcf /hpc/group/manoslab/QMACROCARPA_GENOMES/14_VARIANTS_MONGOLICA/ALL_LIST.gvcf.gz --keep /hpc/group/manoslab/QMACROCARPA_GENOMES/00a_SUBSETS/Subset20Species_",i,".csv --max-alleles 2 --max-missing 0.50 --remove-indels --plink --out /cwork/ram163/OAK_GENOMES/15_STR/SS20_LIST_",i," --temp /scratch\nplink --file /cwork/ram163/OAK_GENOMES/15_STR/SS20_LIST_",i," --recode12 --out /cwork/ram163/OAK_GENOMES/15_STR/SS20_LIST_",i,"\nadmixture --cv /cwork/ram163/OAK_GENOMES/15_STR/SS20_LIST_",i,".ped 12 -j10 | tee log_SS20_",i,"_S12.out",sep=""),file=paste("C:/Users/rmohn/Desktop/10_Analysis/020_STR/STR_",i,".cmd",sep=""),append=FALSE)
}


