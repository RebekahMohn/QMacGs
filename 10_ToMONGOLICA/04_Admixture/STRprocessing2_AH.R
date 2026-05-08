# ADMIXTURE processing
# rmohn@mortonarb.org, 2025-10-30

# DESCRIPTION: This script takes a table of Q values, SE values 
# and calculates the confidence interval for each K value for each 
# sample. It then scores individuals based on these Confidence intervals 
# as pure (1 Q value with a maxCI > 0.999), uncertainPure (1 Q value with 
# a minCI > 0.001 and 0 Q value with a maxCI > 0.999 ), hybrid ( 2+ Q 
# value with a minCI > 0.001 and 0 Q value with a maxCI > 0.999), or 
# uncertain hybrids (  1+ Q value with a minCI > 0.001 and a different 1 
# Q value with a maxCI > 0.999).

# Assuming K=8

library(dplyr)
library(stringr)
library(reshape2)
library(readxl)
library(ggplot2)
library(tidyr)

library(rstudioapi)

# Getting the path of your current open file
current_path = 
  rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))



# read and process data
meta <- read_xlsx("../../00_METADATA/meta_pops.xlsx",sheet="Sheet1")
simp_meta <- data.frame(
  cbind(
    Seq=meta$Seq,
    colnum=meta$collectionNumber,
    species=meta$GenomicVoucherID_USE_THIS,
    state=meta$`state/provIndIowanace`,
    site=meta$siteABB,
    lat=meta$latitude.orig,
    long=meta$longitude.orig
    ))
STR_names <- read.delim("../../10_ToMONGOLICA/04_Admixture/SubsetFullSTR.txt", sep = "\t", header=FALSE)

All_STR8 <- read.delim("../../10_ToMONGOLICA/04_Admixture/ALLSTR_LIST_FILTER10000.8.Q",sep = " ", header=FALSE)

All_STR8$Sample <- STR_names$V1
All_STR8 <- All_STR8 %>%
  dplyr::mutate(across('Sample', str_replace, '.bam', ''))
All_S8_dat <- 
  merge(x=All_STR8, y=simp_meta, by.x="Sample", by.y="Seq",all.x=TRUE)

##############################
###Visualization Code #####
##############################

All_S8_long<-reshape2::melt(All_S8_dat,id.vars =  c("Sample","colnum","species","state","site","lat","long"))
All_S8_long$SamVar<-paste(All_S8_long$Sample,All_S8_long$variable,sep="_")
All_S8_long$species2<-paste("Q. ",sapply(strsplit(x=All_S8_long$species,split=c(" ")),"[",2),sep="")
ggplot(data=All_S8_long)+
  geom_bar(mapping=aes(x=paste(species,Sample),y=value,fill=variable),stat="identity",position="stack",width=1)+
  scale_fill_manual(values=c("#000000","#44AA99","#888888","#88CCEE","#CC6677","#661100","#DDCC77","#2C52CA","#999933","#882255","#AA4499","#117733"))+
  labs(x="samples",fill="population")+
  facet_grid(.~species2,scales = "free", space = "free")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),strip.text.y=element_text(angle=0),
        strip.text.x = element_text(angle = 90,face = "italic"))

######################
### STR_CI ###########
######################

STR_names <- read.delim("../../10_ToMONGOLICA/04_Admixture/SubsetFullSTR.txt",sep = "\t",header=FALSE)
All_STR8_se <- read.delim("../../10_ToMONGOLICA/04_Admixture/ALLSTR_LIST_FILTER10000.8.Q_se",sep = " ", header=FALSE)
All_STR8_se$Sample <- STR_names$V1
All_STR8_se <- All_STR8_se %>%
  dplyr::mutate(across('Sample',str_replace,'.bam',''))

## combine SE for macrocarpa V6 & V7


# SEM = sd / sqrt(n)
# so to calc the pooled SEM you do a weighted average:
# (1) divide each SEM by sample size
# (2) square the results to get the variance
# (3) average the variances
# (4) divide the averaged variance by the summed sample size

# but since the sample size here should be constant (it's a function of 
# the bootstrap replicates, not individuals or loci), you should be able 
# to just average the sem's. See my suggested edits below:


All_STR8_se$macComb <- mean(c(All_STR8_se$V6, All_STR8_se$V7))

All_S8_se_long<-reshape2::melt(All_STR8_se[,-c(6,7)],id.vars =  c("Sample"))
All_S8_se_long$SamVar <- paste(All_S8_se_long$Sample,All_S8_se_long$variable,sep="_")

##combine Q for macrocarpa
All_S8_dat$macComb <- All_S8_dat$V6+All_S8_dat$V7
All_S8_long <- reshape2::melt(All_S8_dat[,-c(7,8)],id.vars =  c("Sample","colnum","species","state","site","lat","long"))
All_S8_long$SamVar <- paste(All_S8_long$Sample,All_S8_long$variable,sep="_")


All_S8_dat_CI_long <- merge(All_S8_long,All_S8_se_long,by="SamVar")


##Calculate confidence intervals###
All_S8_dat_CI_long$minCI <- 
  All_S8_dat_CI_long$value.x - 1.96*All_S8_dat_CI_long$value.y
All_S8_dat_CI_long$maxCI <-
  All_S8_dat_CI_long$value.x + 1.96*All_S8_dat_CI_long$value.y


All_S8_dat_CI_long$strpop <- "sp"
All_S8_dat_CI_long[All_S8_dat_CI_long[,9]=="V1",16] <- "sinstemar"
All_S8_dat_CI_long[All_S8_dat_CI_long[,9]=="V2",16] <- "muepri"
All_S8_dat_CI_long[All_S8_dat_CI_long[,9]=="V3",16] <- "alb"
All_S8_dat_CI_long[All_S8_dat_CI_long[,9]=="V4",16] <- "biclyr"
All_S8_dat_CI_long[All_S8_dat_CI_long[,9]=="V5",16] <- "lob"
All_S8_dat_CI_long[All_S8_dat_CI_long[,9]=="macComb",16] <- "mac"
All_S8_dat_CI_long[All_S8_dat_CI_long[,9]=="V8",16] <- "micmon"

s8_maxCI_4wide<-as.data.frame(cbind(All_S8_dat_CI_long$strpop,All_S8_dat_CI_long$Sample.x,All_S8_dat_CI_long$maxCI))
s8_maxCI_wide<-spread(s8_maxCI_4wide,key = V1,value=as.numeric(V3))
write.csv(s8_maxCI_wide,"C:/Users/rmohn/Desktop/10_Analysis/02_STRUCTURE/FinalForCollabs/s8_max_CI.csv")

s8_minCI_4wide<-as.data.frame(cbind(All_S8_dat_CI_long$strpop,All_S8_dat_CI_long$Sample.x,All_S8_dat_CI_long$minCI))
s8_minCI_wide<-spread(s8_minCI_4wide,key = V1,value=as.numeric(V3))
write.csv(s8_minCI_wide,"C:/Users/rmohn/Desktop/10_Analysis/02_STRUCTURE/FinalForCollabs/s8_min_CI.csv")

############################################################################
### Scoring an individual contribution as likely pure or likely hybrid #####
############################################################################

All_S8_dat_CI_long$contribution<-"none"
All_S8_dat_CI_long$CIval<-0
for(i in 1:length(All_S8_dat_CI_long$minCI)){
  if(All_S8_dat_CI_long$maxCI[i]>0.999){
    All_S8_dat_CI_long$contribution[i]<-"pure?"
    All_S8_dat_CI_long$CIval[i]<-All_S8_dat_CI_long$maxCI[i]
  }else{
    if(All_S8_dat_CI_long$minCI[i]>0.001){
      All_S8_dat_CI_long$contribution[i]<-"hybrid?"
      All_S8_dat_CI_long$CIval[i]<-All_S8_dat_CI_long$minCI[i]
    }
  }
}

##################################################
###Calculating hybridization based on CI intervals
##################################################

# s8_hyb_4wide<-as.data.frame(cbind(All_S8_dat_CI_long$strpop,All_S8_dat_CI_long$Sample.x,All_S8_dat_CI_long$CIval))
# s8_hyb_wide<-spread(s8_hyb_4wide,key = V1,value=as.numeric(V3))
# 
# for(i in c("alb", "biclyr", "lob", "mac", "micmon", "muepri", "sinstemar")) {
#   s8_hyb_wide[[i]] <- as.numeric(s8_hyb_wide[[i]])
# }
# 
# s8_hyb_wide$pureCount <- s8_hyb_wide$hybCount <- 0
# s8_hyb_wide$hyb <- NA
# for(i in 1:length(s8_hyb_wide$V2)){
#   s8_hyb_wide$hybCount[i]<-sum(s8_hyb_wide[i,2:8]>.001&s8_hyb_wide[i,2:8]<.999,na.rm=TRUE)
#   s8_hyb_wide$pureCount[i]<-sum(s8_hyb_wide[i,2:8]>.999,na.rm=TRUE)
#   if(s8_hyb_wide$hybCount[i]> 1 &s8_hyb_wide$pureCount[i] < 1){
#     s8_hyb_wide$hyb[i]<-"hybrid"
#   }
#   if(s8_hyb_wide$hybCount[i]<2 & s8_hyb_wide$pureCount[i] < 1){
#     s8_hyb_wide$hyb[i]<-"uncertain pure"
#   }
#   if(s8_hyb_wide$hybCount[i]<1&s8_hyb_wide$pureCount[i]==1){
#     s8_hyb_wide$hyb[i]<-"pure"
#   }
#   if(s8_hyb_wide$pureCount[i]>1){
#     s8_hyb_wide$hyb[i]<-"?"
#   }
#   if(s8_hyb_wide$pureCount[i]>0&s8_hyb_wide$hybCount[i]>0){
#     s8_hyb_wide$hyb[i]<-"uncertain hybrid"
#   }
# }


##################################################
###Calculating hybridization based combined data
##################################################

s8_hyb_4wide<-as.data.frame(cbind(All_S8_dat_CI_long$strpop,All_S8_dat_CI_long$Sample.x,All_S8_dat_CI_long$contribution))
s8_hyb_wide<-spread(s8_hyb_4wide,key = V1,value=as.numeric(V3))

s8_hyb_wide$pureCount<-0
s8_hyb_wide$hybCount<-0
s8_hyb_wide$hyb<-NA
for(i in 1:length(s8_hyb_wide$V2)){
  s8_hyb_wide$hybCount[i]<-sum(s8_hyb_wide[i,2:8]=="hybrid?",na.rm=TRUE)
  s8_hyb_wide$pureCount[i]<-sum(s8_hyb_wide[i,2:8]=="pure?",na.rm=TRUE)
  if(s8_hyb_wide$hybCount[i]> 1 &s8_hyb_wide$pureCount[i] < 1){
    s8_hyb_wide$hyb[i]<-"hybrid"
  }
  if(s8_hyb_wide$hybCount[i]<2 & s8_hyb_wide$pureCount[i] < 1){
    s8_hyb_wide$hyb[i]<-"uncertain pure"
  }
  if(s8_hyb_wide$hybCount[i]<1&s8_hyb_wide$pureCount[i]==1){
    s8_hyb_wide$hyb[i]<-"pure"
  }
  if(s8_hyb_wide$pureCount[i]>1){
    s8_hyb_wide$hyb[i]<-"?"
  }
  if(s8_hyb_wide$pureCount[i]>0&s8_hyb_wide$hybCount[i]>0){
    s8_hyb_wide$hyb[i]<-"uncertain hybrid"
  }
}

# write results
All_S8_dat$V1<-round(All_S8_dat$V1,4)
All_S8_dat$V2<-round(All_S8_dat$V2,4)
All_S8_dat$V3<-round(All_S8_dat$V3,4)
All_S8_dat$V4<-round(All_S8_dat$V4,4)
All_S8_dat$V5<-round(All_S8_dat$V5,4)
All_S8_dat$V6<-round(All_S8_dat$V6,4)
All_S8_dat$V7<-round(All_S8_dat$V7,4)
All_S8_dat$V8<-round(All_S8_dat$V8,4)


s8_hyb_wide_v2 <- cbind(s8_hyb_wide, All_S8_dat[c(paste('V', 1:8, sep = ''), 'macComb')])
names(s8_hyb_wide_v2)[1] <- 'sample'
names(s8_hyb_wide_v2)[
  match(paste('V', 1:8, sep = ''), 
  names(s8_hyb_wide_v2))
  ] <- c(
    "sinstemar_1",
    "muepri_1",
    "alb_1",
    "biclyr_1",
    "lob_1",
    'mac1_1',
    'mac2_1',
    "micmon_1"
  )

NameConv<-read_excel("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/NameConv.xlsx")

s8_hyb_wide_v3<-merge(s8_hyb_wide_v2,NameConv, by.x="sample",by.y="Seq")

write.csv(s8_hyb_wide_v3, '../../10_ToMONGOLICA/04_Admixture/s8_hyb_wide_v3.csv')


########Determination
#make column for determination
#nameslist<-c("","","","","","","","","","",NA,"sinstemar","muepri","Q. alba","biclyr","Q. lobata",NA,NA,"micmon","Q. macrocarpa","")
nameslist<-c("sinstemar","muepri","Q. alba","biclyr","Q. lobata","micmon","Q. macrocarpa")


s8_hyb_wide_v2$determination<-"Q. sp"
#for hybrids 50x50
for(i in 1:length(s8_hyb_wide_v2$sample)){
  if(s8_hyb_wide_v2$hyb[i]=="hybrid"&sum(s8_hyb_wide_v2[i,c(12:16,19:20)]>0.4&s8_hyb_wide_v2[i,c(12:16,19:20)]<=0.6)>=1){
      s8_hyb_wide_v2$determination[i]<-paste(nameslist[s8_hyb_wide_v2[i,c(12:16,19:20)]>0.2&s8_hyb_wide_v2[i,c(12:16,19:20)]<=0.6],collapse =" x ")
  }else{
    s8_hyb_wide_v2$determination[i]<-paste(nameslist[s8_hyb_wide_v2[i,c(12:16,19:20)]>0.6],collapse="")
  }
}
#for 50x25x25 hybrids

#for pure
