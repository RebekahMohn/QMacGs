
#Identifying whether adaptive are more likely to be introgressed than non-adaptive
##### 2000
gea<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_lfmmpval_cor.txt",sep=" ")
D_All_2000<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/107_Dsuite_MACREF/DPOPS/Dall_mac_2000.txt",header = TRUE)
D_SS25_2000<-subset(D_All_2000,d_f>0.25)

gea$CHR_Pos<-paste(gea$CHR,gea$BP,sep="Q")


##Add the different speciees
i=1
genome_env<-read.table(paste("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_lfmmpval.txt",sep=""),sep=" ",header=TRUE)
genome_env$CHR_Pos<-paste(genome_env$V1,"Q",genome_env$V4,sep="")
GEASAMPS<-data.frame(Chr_Pos=sample(genome_env$CHR_Pos,size = 1000,replace = F))
GEASAMPS$mueDstats<-0
GEASAMPS$steDstats<-0
GEASAMPS$lobDstats<-0
GEASAMPS$albDstats<-0
for(j in 1:length(GEASAMPS$Chr_Pos)){
  chr1<-as.numeric(str_split(GEASAMPS$Chr_Pos[j],"Q")[[1]][1])
  pos<-as.numeric(str_split(GEASAMPS$Chr_Pos[j],"Q")[[1]][2])
  GEASAMPS$mueDstats[j]<-length(subset(D_SS25_2000,chr3==chr1&windowStart<pos&windowEnd>pos&d_f>0.25&species=="mue")$d_f)
  GEASAMPS$steDstats[j]<-length(subset(D_SS25_2000,chr3==chr1&windowStart<pos&windowEnd>pos&d_f>0.25&species=="ste")$d_f)
  GEASAMPS$lobDstats[j]<-length(subset(D_SS25_2000,chr3==chr1&windowStart<pos&windowEnd>pos&d_f>0.25&species=="lob")$d_f)
  GEASAMPS$albDstats[j]<-length(subset(D_SS25_2000,chr3==chr1&windowStart<pos&windowEnd>pos&d_f>0.25&species=="alb")$d_f)
}

ggplot(GEASAMPS)+
  geom_histogram(aes(x=mueDstats),binwidth = 1)
ggplot(GEASAMPS)+
  geom_histogram(aes(x=lobDstats),binwidth = 1)
ggplot(GEASAMPS)+
  geom_histogram(aes(x=steDstats),binwidth = 1)
ggplot(GEASAMPS)+
  geom_histogram(aes(x=albDstats),binwidth = 1)


gea$mueDstats<-NA
gea$steDstats<-NA
gea$lobDstats<-NA
gea$albDstats<-NA

for(j in 1:length(gea$CHR)){
  chr1<-gea$CHR[j]
  pos<-gea$BP[j]
  gea$mueDstats[j]<-length(subset(D_SS25_2000,chr3==chr1&windowStart<pos&windowEnd>pos&d_f>0.25&species=="mue")$d_f)
  gea$steDstats[j]<-length(subset(D_SS25_2000,chr3==chr1&windowStart<pos&windowEnd>pos&d_f>0.25&species=="ste")$d_f)
  gea$lobDstats[j]<-length(subset(D_SS25_2000,chr3==chr1&windowStart<pos&windowEnd>pos&d_f>0.25&species=="lob")$d_f)
  gea$albDstats[j]<-length(subset(D_SS25_2000,chr3==chr1&windowStart<pos&windowEnd>pos&d_f>0.25&species=="alb")$d_f)
}

gea_dstat<-subset(gea,mueDstats>3|albDstats>3|steDstats>3|lobDstats>3)
gea_dstat_table<-bind_rows(mue=table(subset(gea_dstat,mueDstats>3)$ENV),ste=table(subset(gea_dstat,steDstats>3)$ENV),
          lob=table(subset(gea_dstat,lobDstats>3)$ENV),alb=table(subset(gea_dstat,albDstats>3)$ENV),total=table(gea$ENV))

write.table(gea_dstat_table,"C:/Users/rmohn/Documents/GitHub/QMacGs/20_ToMACROCARPA/10_ComboAnalysis/Output/gea_introgression.txt",sep="\t",quote=F,row.names = F)


length(subset(GEASAMPS,mueDstats>3)$mueDstats)
length(subset(GEASAMPS,steDstats>3)$steDstats)
length(subset(GEASAMPS,albDstats>3)$albDstats)
length(subset(GEASAMPS,lobDstats>3)$lobDstats)

table(subset(gea_dstat,mueDstats>3)$ENV)
length(subset(GEASAMPS,steDstats>3)$steDstats)
length(subset(GEASAMPS,albDstats>3)$albDstats)
length(subset(GEASAMPS,lobDstats>3)$lobDstats)

#mue
fisher.test(data.frame(c(45,1000),c(69,2617)))
fisher.test(data.frame(c(45,1000),c(179,540)))
fisher.test(data.frame(c(45,1000),c(25,148)))
fisher.test(data.frame(c(45,1000),c(0,10)))
fisher.test(data.frame(c(45,1000),c(14,31)))
fisher.test(data.frame(c(45,1000),c(18,73)))

#alb
fisher.test(data.frame(c(18,1000),c(19,2617)))
fisher.test(data.frame(c(18,1000),c(153,540)))
fisher.test(data.frame(c(18,1000),c(1,148)))
fisher.test(data.frame(c(18,1000),c(0,10)))
fisher.test(data.frame(c(18,1000),c(10,31)))
fisher.test(data.frame(c(18,1000),c(3,73)))


#alb
fisher.test(data.frame(c(62,1000),c(146,2617)))
fisher.test(data.frame(c(62,1000),c(123,540)))
fisher.test(data.frame(c(62,1000),c(3,148)))
fisher.test(data.frame(c(62,1000),c(0,10)))
fisher.test(data.frame(c(62,1000),c(11,31)))
fisher.test(data.frame(c(62,1000),c(3,73)))

#alb
fisher.test(data.frame(c(26,1000),c(38,2617)))
fisher.test(data.frame(c(26,1000),c(156,540)))
fisher.test(data.frame(c(26,1000),c(1,148)))
fisher.test(data.frame(c(26,1000),c(1,10)))
fisher.test(data.frame(c(26,1000),c(10,31)))
fisher.test(data.frame(c(26,1000),c(3,73)))