
#Identifying whether adaptive are more likely to be introgressed than non-adaptive
##### 2000
library(dplyr)
library(stringr)
library(ggplot2)

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
  subset_D_2000<-0
  subset_D_2000<-subset(D_SS25_2000,chr3==chr1&windowStart<pos&windowEnd>pos&d_f>0.25)
  subset_D_2000b<-subset_D_2000 %>% group_by(pop) %>% top_n(1, d_f)
  GEASAMPS$mueDstats[j]<-length(subset(subset_D_2000b,species=="mue")$d_f)
  GEASAMPS$steDstats[j]<-length(subset(subset_D_2000b,species=="ste")$d_f)
  GEASAMPS$lobDstats[j]<-length(subset(subset_D_2000b,species=="lob")$d_f)
  GEASAMPS$albDstats[j]<-length(subset(subset_D_2000b,species=="alb")$d_f)
}

# ggplot(GEASAMPS)+
#   geom_histogram(aes(x=mueDstats),binwidth = 1)
# ggplot(GEASAMPS)+
#   geom_histogram(aes(x=lobDstats),binwidth = 1)
# ggplot(GEASAMPS)+
#   geom_histogram(aes(x=steDstats),binwidth = 1)
# ggplot(GEASAMPS)+
#   geom_histogram(aes(x=albDstats),binwidth = 1)


gea$mueDstats<-NA
gea$steDstats<-NA
gea$lobDstats<-NA
gea$albDstats<-NA

for(j in 1:length(gea$CHR)){
  chr1<-gea$CHR[j]
  pos<-gea$BP[j]
  subset_D_2000<-0
  subset_D_2000<-subset(D_SS25_2000,chr3==chr1&windowStart<pos&windowEnd>pos&d_f>0.25)
  subset_D_2000b<-subset_D_2000 %>% group_by(pop) %>% top_n(1, d_f)
  gea$mueDstats[j]<-length(subset(subset_D_2000b,species=="mue")$d_f)
  gea$steDstats[j]<-length(subset(subset_D_2000b,species=="ste")$d_f)
  gea$lobDstats[j]<-length(subset(subset_D_2000b,species=="lob")$d_f)
  gea$albDstats[j]<-length(subset(subset_D_2000b,species=="alb")$d_f)
}

gea_dstat<-subset(gea,mueDstats>3|albDstats>3|steDstats>3|lobDstats>3)
gea_dstat_table<-bind_rows(mue=table(subset(gea_dstat,mueDstats>3)$ENV),ste=table(subset(gea_dstat,steDstats>3)$ENV),
          lob=table(subset(gea_dstat,lobDstats>3)$ENV),alb=table(subset(gea_dstat,albDstats>3)$ENV),total=table(gea$ENV))

write.table(gea_dstat_table,"C:/Users/rmohn/Documents/GitHub/QMacGs/20_ToMACROCARPA/10_ComboAnalysis/Output/gea_introgression_2000.txt",sep="\t",quote=F,row.names = F)


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
mueE1<-fisher.test(data.frame(c(66,2617),c(83,1000)))
mueE2<-fisher.test(data.frame(c(135,540),c(83,1000)))
mueE3<-fisher.test(data.frame(c(27,148),c(83,1000)))
mueE4<-fisher.test(data.frame(c(0,10),c(83,1000)))
mueE5<-fisher.test(data.frame(c(14,31),c(83,1000)))
mueE6<-fisher.test(data.frame(c(20,73),c(83,1000)))

mueE1$p.value
mueE2$p.value
mueE3$p.value
mueE4$p.value
mueE5$p.value
mueE6$p.value

mueE1$estimate
mueE2$estimate
mueE3$estimate
mueE4$estimate
mueE5$estimate
mueE6$estimate

#alb
albE1<-fisher.test(data.frame(c(20,2617),c(36,1000)))
albE2<-fisher.test(data.frame(c(117,540),c(36,1000)))
albE3<-fisher.test(data.frame(c(1,148),c(36,1000)))
albE4<-fisher.test(data.frame(c(0,10),c(36,1000)))
albE5<-fisher.test(data.frame(c(9,31),c(36,1000)))
albE6<-fisher.test(data.frame(c(1,73),c(36,1000)))

albE1$p.value
albE2$p.value
albE3$p.value
albE4$p.value
albE5$p.value
albE6$p.value

albE1$estimate
albE2$estimate
albE3$estimate
albE4$estimate
albE5$estimate
albE6$estimate


#lob
lobE1<-fisher.test(data.frame(c(159,2617),c(89,1000)))
lobE2<-fisher.test(data.frame(c(120,540),c(89,1000)))
lobE3<-fisher.test(data.frame(c(3,148),c(89,1000)))
lobE4<-fisher.test(data.frame(c(1,10),c(89,1000)))
lobE5<-fisher.test(data.frame(c(10,31),c(89,1000)))
lobE6<-fisher.test(data.frame(c(3,73),c(89,1000)))

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
steE1<-fisher.test(data.frame(c(38,2617),c(29,1000)))
steE2<-fisher.test(data.frame(c(114,540),c(29,1000)))
steE3<-fisher.test(data.frame(c(1,148),c(29,1000)))
steE4<-fisher.test(data.frame(c(1,10),c(29,1000)))
steE5<-fisher.test(data.frame(c(8,31),c(29,1000)))
steE6<-fisher.test(data.frame(c(2,73),c(29,1000)))

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
anyE1<-fisher.test(data.frame(c(236,2617),c(182,1000)))
anyE2<-fisher.test(data.frame(c(206,540),c(182,1000)))
anyE3<-fisher.test(data.frame(c(30,148),c(182,1000)))
anyE4<-fisher.test(data.frame(c(2,10),c(182,1000)))
anyE5<-fisher.test(data.frame(c(15,31),c(182,1000)))
anyE6<-fisher.test(data.frame(c(20,73),c(182,1000)))

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
