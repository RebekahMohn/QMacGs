
#Identifying whether adaptive are more likely to be introgressed than non-adaptive
##### 2000
gea<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_lfmmpval_cor.txt",sep=" ")
Dmac_All_2000<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/107_Dsuite_MACREF/DPOPS/Dall_mac_2000.txt",header = TRUE)
D_SS25_2000<-subset(Dmac_All_2000,d_f>0.25)

gea$CHR_Pos<-paste(gea$CHR,gea$BP,sep="Q")



lobPaint<-subset(read.delim("C:/Users/rmohn/Desktop/10_Analysis/113_painting/highpower/lob_mac_bic_ItrgByLat_samp.tsv",header = TRUE),(X.Chr01.=="Chr01"|X.Chr01.=="Chr02"|X.Chr01.=="Chr03"|
                                                                                                                                         X.Chr01.=="Chr04"|X.Chr01.=="Chr05"|X.Chr01.=="Chr06"|X.Chr01.=="Chr07"|X.Chr01.=="Chr08"|X.Chr01.=="Chr09"|X.Chr01.=="Chr10"|X.Chr01.=="Chr11"|X.Chr01.=="Chr12"))
colnames(lobPaint)<-c("chr","window","species","fixedBases","corgeneVsp","rsq","meanval","medianval","meanP","lobnumDev","sd")
lobPaint$csome_window<-paste(lobPaint$chr,"_",lobPaint$window,sep="")

stePaint<-subset(read.delim("C:/Users/rmohn/Desktop/10_Analysis/113_painting/highpower/ste_mac_bic_ItrgByLat_samp.tsv",header = TRUE),(X.Chr01.=="Chr01"|X.Chr01.=="Chr02"|X.Chr01.=="Chr03"|
                                                                                                                                         X.Chr01.=="Chr04"|X.Chr01.=="Chr05"|X.Chr01.=="Chr06"|X.Chr01.=="Chr07"|X.Chr01.=="Chr08"|X.Chr01.=="Chr09"|X.Chr01.=="Chr10"|X.Chr01.=="Chr11"|X.Chr01.=="Chr12"))
colnames(stePaint)<-c("chr","window","species","fixedBases","corgeneVsp","rsq","meanval","medianval","meanP","stenumDev","sd")
stePaint$csome_window<-paste(stePaint$chr,"_",stePaint$window,sep="")

muePaint<-subset(read.delim("C:/Users/rmohn/Desktop/10_Analysis/113_painting/highpower/mue_mac_bic_ItrgByLat_samp.tsv",header = TRUE),(X.Chr01.=="Chr01"|X.Chr01.=="Chr02"|X.Chr01.=="Chr03"|
                                                                                                                                         X.Chr01.=="Chr04"|X.Chr01.=="Chr05"|X.Chr01.=="Chr06"|X.Chr01.=="Chr07"|X.Chr01.=="Chr08"|X.Chr01.=="Chr09"|X.Chr01.=="Chr10"|X.Chr01.=="Chr11"|X.Chr01.=="Chr12"))
colnames(muePaint)<-c("chr","window","species","fixedBases","corgeneVsp","rsq","meanval","medianval","meanP","muenumDev","sd")
muePaint$csome_window<-paste(muePaint$chr,"_",muePaint$window,sep="")

albPaint<-subset(read.delim("C:/Users/rmohn/Desktop/10_Analysis/113_painting/highpower/alb_mac_bic_ItrgByLat_samp.tsv",header = TRUE),(X.Chr01.=="Chr01"|X.Chr01.=="Chr02"|X.Chr01.=="Chr03"|
                                                                                                                                         X.Chr01.=="Chr04"|X.Chr01.=="Chr05"|X.Chr01.=="Chr06"|X.Chr01.=="Chr07"|X.Chr01.=="Chr08"|X.Chr01.=="Chr09"|X.Chr01.=="Chr10"|X.Chr01.=="Chr11"|X.Chr01.=="Chr12"))
colnames(albPaint)<-c("achr","awindow","aspecies","afixedBases","acorgeneVsp","arsq","ameanval","amedianval","ameanP","albnumDev","asd")
albPaint$csome_window<-paste(albPaint$achr,"_",albPaint$awindow,sep="")

paintAll<-merge(unique(lobPaint),unique(stePaint),by="csome_window",all.x = TRUE)
paintAll<-merge(paintAll,unique(muePaint),by="csome_window",all.x = TRUE)
paintAll<-merge(paintAll,unique(albPaint),by="csome_window",all.x = TRUE)




##Add the different speciees
i=1
genome_env<-read.table(paste("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_lfmmpval.txt",sep=""),sep=" ",header=TRUE)
genome_env$CHR_Pos<-paste(genome_env$V1,"Q",genome_env$V4,sep="")
genome_env$csome_window<-paste(genome_env$V1,"_",round(gea1_c$BP/100000,0)*100000,sep="")

GEASAMPS<-data.frame(Chr_Pos=sample(genome_env$CHR_Pos,size = 1000,replace = F))
GEASAMPS$mueDstats<-0
GEASAMPS$steDstats<-0
GEASAMPS$lobDstats<-0
GEASAMPS$albDstats<-0
EASAMPS$paintStats<-0
for(j in 1:length(GEASAMPS$Chr_Pos)){
  chr1<-as.numeric(str_split(GEASAMPS$Chr_Pos[j],"Q")[[1]][1])
  pos<-as.numeric(str_split(GEASAMPS$Chr_Pos[j],"Q")[[1]][2])
  pos_round<-round(as.numeric(str_split(GEASAMPS$Chr_Pos[j],"Q")[[1]][2])/100000,0)*100000
  paintSS<-subset(paintAll,chr.x==chr1&window.x==pos_round)
  GEASAMPS$paintStats[j]<-length(subset(paintSS,(albnumDev>50|lobnumDev>50|muenumDev>50|stenumDev>50)&(fixedBases.x>=40&fixedBases.y>=40&fixedBases>=40&afixedBases>=40))$chr.x)
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

gea2<-gea
gea2$mueDstats<-NA
gea2$steDstats<-NA
gea2$lobDstats<-NA
gea2$albDstats<-NA

for(j in 1:length(gea$CHR)){
  chr1<-gea2$CHR[j]
  pos<-gea2$BP[j]
  gea2$mueDstats[j]<-length(subset(D_SS25_2000,chr3==chr1&windowStart<pos&windowEnd>pos&d_f>0.25&species=="mue")$d_f)
  gea2$steDstats[j]<-length(subset(D_SS25_2000,chr3==chr1&windowStart<pos&windowEnd>pos&d_f>0.25&species=="ste")$d_f)
  gea2$lobDstats[j]<-length(subset(D_SS25_2000,chr3==chr1&windowStart<pos&windowEnd>pos&d_f>0.25&species=="lob")$d_f)
  gea2$albDstats[j]<-length(subset(D_SS25_2000,chr3==chr1&windowStart<pos&windowEnd>pos&d_f>0.25&species=="alb")$d_f)
}

gea_dstat_2<-subset(gea2,mueDstats>3|albDstats>3|steDstats>3|lobDstats>3)
gea_dstat_2_table<-bind_rows(mue=table(subset(gea_dstat_2,mueDstats>3)$ENV),ste=table(subset(gea_dstat_2,steDstats>3)$ENV),
          lob=table(subset(gea_dstat_2,lobDstats>3)$ENV),alb=table(subset(gea_dstat_2,albDstats>3)$ENV),total=table(gea2$ENV))

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