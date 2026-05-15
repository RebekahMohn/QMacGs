#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
#  args[2] = "out.txt"
}

#SweepFinder2 Configuring Data
filelist<-args[1] #list.files(".","*.count")
for(f in filelist){
  dat<-read.delim(file=f,header = T,sep = "\t",row.names = NULL)
  popname<-strsplit(f,"[.]")[[1]][1]
  csome<-strsplit(f,"[.]")[[1]][2]
  datrecon<-data.frame(position=dat$CHROM,x=dat$N_CHR,n=dat$N_ALLELES,folded=1)
  write.table(datrecon,paste(popname,".",csome,".alfreq",sep=""),quote=F, row.names = F,sep="\t")
}
