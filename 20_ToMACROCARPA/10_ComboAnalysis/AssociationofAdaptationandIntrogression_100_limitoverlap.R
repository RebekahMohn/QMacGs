
#Identifying whether adaptive are more likely to be introgressed than non-adaptive
##### 100
library(dplyr)
library(stringr)

gea<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_lfmmpval_cor.txt",sep=" ")
D_All_100<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/107_Dsuite_MACREF/DPOPS/100_100/Dall_mac_100.txt",header = TRUE)
D_SS25_100<-subset(D_All_100,d_f>0.25)

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
  subset_D_100<-0
  subset_D_100<-subset(D_SS25_100,chr3==chr1&windowStart<pos&windowEnd>pos&d_f>0.25)
  subset_D_100b<-subset_D_100 %>% group_by(pop) %>% top_n(1, d_f)
  GEASAMPS$mueDstats[j]<-length(subset(subset_D_100b,species=="mue")$d_f)
  GEASAMPS$steDstats[j]<-length(subset(subset_D_100b,species=="ste")$d_f)
  GEASAMPS$lobDstats[j]<-length(subset(subset_D_100b,species=="lob")$d_f)
  GEASAMPS$albDstats[j]<-length(subset(subset_D_100b,species=="alb")$d_f)
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
  subset_D_100<-0
  subset_D_100<-subset(D_SS25_100,chr3==chr1&windowStart<pos&windowEnd>pos&d_f>0.25)
  subset_D_100b<-subset_D_100 %>% group_by(pop) %>% top_n(1, d_f)
  gea$mueDstats[j]<-length(subset(subset_D_100b,species=="mue")$d_f)
  gea$steDstats[j]<-length(subset(subset_D_100b,species=="ste")$d_f)
  gea$lobDstats[j]<-length(subset(subset_D_100b,species=="lob")$d_f)
  gea$albDstats[j]<-length(subset(subset_D_100b,species=="alb")$d_f)
}

gea_dstat<-subset(gea,mueDstats>3|albDstats>3|steDstats>3|lobDstats>3)
gea_dstat_table<-bind_rows(mue=table(subset(gea_dstat,mueDstats>3)$ENV),alb=table(subset(gea_dstat,albDstats>3)$ENV),
          lob=table(subset(gea_dstat,lobDstats>3)$ENV),ste=table(subset(gea_dstat,steDstats>3)$ENV),total=table(gea$ENV))

write.table(gea_dstat_table,"C:/Users/rmohn/Documents/GitHub/QMacGs/20_ToMACROCARPA/10_ComboAnalysis/Output/gea_introgression_100.txt",sep="\t",quote=F,row.names = F)


length(subset(GEASAMPS,mueDstats>3)$mueDstats)
length(subset(GEASAMPS,steDstats>3)$steDstats)
length(subset(GEASAMPS,albDstats>3)$albDstats)
length(subset(GEASAMPS,lobDstats>3)$lobDstats)

length(subset(GEASAMPS,mueDstats>3|steDstats>3|albDstats>3|lobDstats>3)$mueDstats)

table(subset(gea_dstat,mueDstats>3)$ENV)
length(subset(GEASAMPS,steDstats>3)$steDstats)
length(subset(GEASAMPS,albDstats>3)$albDstats)
length(subset(GEASAMPS,lobDstats>3)$lobDstats)

length(subset(gea_dstat,ENV=="e1")$mueDstats)
length(subset(gea_dstat,ENV=="e2")$mueDstats)
length(subset(gea_dstat,ENV=="e3")$mueDstats)
length(subset(gea_dstat,ENV=="e4")$mueDstats)
length(subset(gea_dstat,ENV=="e5")$mueDstats)
length(subset(gea_dstat,ENV=="e6")$mueDstats)

#mue
fisher.test(data.frame(c(467,2617),c(212,1000)))
fisher.test(data.frame(c(159,540),c(212,1000)))
fisher.test(data.frame(c(24,148),c(212,1000)))
fisher.test(data.frame(c(2,10),c(212,1000)))
fisher.test(data.frame(c(13,31),c(212,1000)))
fisher.test(data.frame(c(13,73),c(212,1000)))

#alb
fisher.test(data.frame(c(190,2617),c(102,1000)))
fisher.test(data.frame(c(62,540),c(102,1000)))
fisher.test(data.frame(c(12,148),c(102,1000)))
fisher.test(data.frame(c(1,10),c(102,1000)))
fisher.test(data.frame(c(2,31),c(102,1000)))
fisher.test(data.frame(c(6,73),c(102,1000)))


#lob
lobE1<-fisher.test(data.frame(c(732,2617),c(277,1000)))
lobE2<-fisher.test(data.frame(c(195,540),c(277,1000)))
lobE3<-fisher.test(data.frame(c(48,148),c(277,1000)))
lobE4<-fisher.test(data.frame(c(1,10),c(277,1000)))
lobE5<-fisher.test(data.frame(c(9,31),c(277,1000)))
lobE6<-fisher.test(data.frame(c(12,73),c(277,1000)))

lobE1$p.value
lobE2$p.value
lobE3$p.value
lobE4$p.value
lobE5$p.value
lobE6$p.value

lobE1$estimate
lobE2$estimate
lobE3$estimate
lobE4$estimate
lobE5$estimate
lobE6$estimate

#ste
steE1<-fisher.test(data.frame(c(348,2617),c(138,1000)))
steE2<-fisher.test(data.frame(c(114,540),c(138,1000)))
steE3<-fisher.test(data.frame(c(19,148),c(138,1000)))
steE4<-fisher.test(data.frame(c(4,10),c(138,1000)))
steE5<-fisher.test(data.frame(c(11,31),c(138,1000)))
steE6<-fisher.test(data.frame(c(6,73),c(138,1000)))


steE1$p.value
steE2$p.value
steE3$p.value
steE4$p.value
steE5$p.value
steE6$p.value

steE1$estimate
steE2$estimate
steE3$estimate
steE4$estimate
steE5$estimate
steE6$estimate



#any
anyE1<-fisher.test(data.frame(c(1222,2617),c(496,1000)))
anyE2<-fisher.test(data.frame(c(305,540),c(496,1000)))
anyE3<-fisher.test(data.frame(c(70,148),c(496,1000)))
anyE4<-fisher.test(data.frame(c(5,10),c(496,1000)))
anyE5<-fisher.test(data.frame(c(19,31),c(496,1000)))
anyE6<-fisher.test(data.frame(c(27,73),c(496,1000)))


anyE1$p.value
anyE2$p.value
anyE3$p.value
anyE4$p.value
anyE5$p.value
anyE6$p.value

anyE1$estimate
anyE2$estimate
anyE3$estimate
anyE4$estimate
anyE5$estimate
anyE6$estimate

