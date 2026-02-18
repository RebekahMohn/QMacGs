#Tajima's D
#Read in Tajima's D values for large populations
library(stringr)

TD<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/110_TajimasD/TD_ALL_POP.txt",sep=" ")
TDminus2<-subset(TD,!TajimaD>-2)
TDminus2wins<-aggregate(TajimaD~CHROM+BIN_START,data=TDminus2,FUN=mean)

gea<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_lfmmpval_cor.txt",sep=" ")

gea$CHR_Pos<-paste(gea$csome,gea$BP,sep="Q")


##Add the different speciees
i=1
genome_env<-read.table(paste("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_lfmmpval.txt",sep=""),sep=" ",header=TRUE)

genome_env$CHR_Pos<-paste(genome_env$V1,"Q",genome_env$V4,sep="")
genome_env$Chr2_Pos<-paste(genome_env$csome,"Q",genome_env$V4,sep="")
randsamples<-sample(1:length(genome_env$CHR_Pos),size = 1000,replace = F)
GEASAMPS<-data.frame(CHR=genome_env$V1[randsamples],Chr=genome_env$csome[randsamples],BP=genome_env$V4[randsamples])
GEASAMPS$TajimaD<-NA

for(j in 1:length(GEASAMPS$CHR)){
  chr1<-as.numeric(GEASAMPS$CHR[j])
  Chr<-GEASAMPS$Chr[j]
  pos<-as.numeric(GEASAMPS$BP[j])
  GEASAMPS$TajimaD[j]<-length(subset(TDminus2,CHROM==Chr&BIN_START<pos&BIN_START+5000>pos&!TajimaD>-2)$TajimaD)
}

ggplot(GEASAMPS)+
  geom_histogram(aes(x=TajimaD),binwidth = 1)

GEASAMPS_tajimaD<-length(subset(GEASAMPS,TajimaD>1)$TajimaD)

gea$TajimaD<-NA

for(j in 1:length(gea$CHR)){
  chr1<-gea$csome[j]
  pos<-gea$BP[j]
  gea$TajimaD[j]<-length(subset(TDminus2,CHROM==chr1&BIN_START<pos&BIN_START+5000>pos&!TajimaD>-2)$TajimaD)
}

gea_tajimaD<-subset(gea,TajimaD>1)
library(dplyr)
TD_ENV<-as.data.frame(bind_rows(TajimasD=table(gea_tajimaD$ENV),total=table(gea$ENV)))
rownames(TD_ENV)<-c("TajimasD", "total")

write.table(TD_ENV,"C:/Users/rmohn/Documents/GitHub/QMacGs/20_ToMACROCARPA/10_ComboAnalysis/Output/TajimasD_GEA.txt",sep="\t",quote=F)

#TajimasD
fisher.test(data.frame(c(19,1000),c(51,2617)))
fisher.test(data.frame(c(19,1000),c(6,540)))
fisher.test(data.frame(c(19,1000),c(16,148)))
fisher.test(data.frame(c(19,1000),c(0,10)))
fisher.test(data.frame(c(19,1000),c(2,31)))
fisher.test(data.frame(c(19,1000),c(1,73)))


