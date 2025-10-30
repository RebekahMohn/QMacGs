

#######################################################
####### K=8   ##########################################
#########################################################
library(dplyr)
library(stringr)
library(reshape2)
library(readxl)
library(ggplot2)
meta<-read_xlsx("00_METADATA/meta_pops.xlsx",sheet="Sheet1")
simp_meta<-data.frame(cbind(Seq=meta$Seq,colnum=meta$collectionNumber,species=meta$GenomicVoucherID_USE_THIS,state=meta$`state/provIndIowanace`,site=meta$siteABB,lat=meta$latitude.orig,long=meta$longitude.orig))
STR_names<-read.delim("10_ToMONGOLICA/04_Admixture/SubsetFullSTR.txt",sep = "\t",header=FALSE)

All_STR8<-read.delim("10_ToMONGOLICA/04_Admixture/ALLSTR_LIST_FILTER10000.8.Q",sep = " ", header=FALSE)
All_STR8$Sample<-STR_names$V1
All_STR8<-All_STR8 %>%
  dplyr::mutate(across('Sample',str_replace,'.bam',''))
All_S8_dat<-merge(x=All_STR8,y=simp_meta,by.x="Sample",by.y="Seq",all.x=TRUE)

##############################
###Visualization Code #####
##############################

# All_S8_long<-reshape2::melt(All_S8_dat,id.vars =  c("Sample","colnum","species","state","site","lat","long"))
# All_S8_long$SamVar<-paste(All_S8_long$Sample,All_S8_long$variable,sep="_")
# ggplot(data=All_S8_long)+
#   geom_bar(mapping=aes(x=paste(species,Sample),y=value,fill=variable),stat="identity",position="stack",width=1)+
#   scale_fill_manual(values=c("#888888","#88CCEE","#CC6677","#117733","#661100","#DDCC77","#332288","#44AA99","#AA4499","#999933","#000000","#882255"))+
#   labs(title = "All Chromosomes; 8 populations")+
#   facet_grid(.~species,scales = "free", space = "free")+
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),strip.text.y=element_text(angle=0))

######################
### STR_CI ###########
######################
STR_names<-read.delim("10_ToMONGOLICA/04_Admixture/SubsetFullSTR.txt",sep = "\t",header=FALSE)
All_STR8_se<-read.delim("10_ToMONGOLICA/04_Admixture/ALLSTR_LIST_FILTER10000.8.Q_se",sep = " ", header=FALSE)
All_STR8_se$Sample<-STR_names$V1
All_STR8_se<-All_STR8_se %>%
  dplyr::mutate(across('Sample',str_replace,'.bam',''))

##combine SE for macrocarpa V6 & V7

All_STR8_se$macComb<-((sqrt(All_STR8_se$V6)+sqrt(All_STR8_se$V7))/2)^2
All_S8_se_long<-reshape2::melt(All_STR8_se[,-c(6,7)],id.vars =  c("Sample"))
All_S8_se_long$SamVar<-paste(All_S8_se_long$Sample,All_S8_se_long$variable,sep="_")

##combine Q for macrocarpa
All_S8_dat$macComb<-All_S8_dat$V6+All_S8_dat$V7
All_S8_long<-reshape2::melt(All_S8_dat[,-c(7,8)],id.vars =  c("Sample","colnum","species","state","site","lat","long"))
All_S8_long$SamVar<-paste(All_S8_long$Sample,All_S8_long$variable,sep="_")


All_S8_dat_CI_long<-merge(All_S8_long,All_S8_se_long,by="SamVar")


##Calculate confidence intervals###
All_S8_dat_CI_long$minCI<-All_S8_dat_CI_long$value.x-1.96*All_S8_dat_CI_long$value.y
All_S8_dat_CI_long$maxCI<-All_S8_dat_CI_long$value.x+1.96*All_S8_dat_CI_long$value.y


All_S8_dat_CI_long$strpop<-"sp"
All_S8_dat_CI_long[All_S8_dat_CI_long[,9]=="V1",16]<-"sinstemar"
All_S8_dat_CI_long[All_S8_dat_CI_long[,9]=="V2",16]<-"muepri"
All_S8_dat_CI_long[All_S8_dat_CI_long[,9]=="V3",16]<-"alb"
All_S8_dat_CI_long[All_S8_dat_CI_long[,9]=="V4",16]<-"biclyr"
All_S8_dat_CI_long[All_S8_dat_CI_long[,9]=="V5",16]<-"lob"
All_S8_dat_CI_long[All_S8_dat_CI_long[,9]=="macComb",16]<-"mac"
All_S8_dat_CI_long[All_S8_dat_CI_long[,9]=="V8",16]<-"micmon"





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
library(tidyr)

##################################################
###Calculating hybridization based on CI intervals
##################################################

s8_hyb_4wide<-as.data.frame(cbind(All_S8_dat_CI_long$strpop,All_S8_dat_CI_long$Sample.x,All_S8_dat_CI_long$CIval))
s8_hyb_wide<-spread(s8_hyb_4wide,key = V1,value=as.numeric(V3))

s8_hyb_wide$alb<-as.numeric(s8_hyb_wide$alb)
s8_hyb_wide$biclyr<-as.numeric(s8_hyb_wide$biclyr)
s8_hyb_wide$lob<-as.numeric(s8_hyb_wide$lob)
s8_hyb_wide$mac<-as.numeric(s8_hyb_wide$mac)
s8_hyb_wide$micmon<-as.numeric(s8_hyb_wide$micmon)
s8_hyb_wide$muepri<-as.numeric(s8_hyb_wide$muepri)
s8_hyb_wide$sinstemar<-as.numeric(s8_hyb_wide$sinstemar)

s8_hyb_wide$pureCount<-0
s8_hyb_wide$hybCount<-0
s8_hyb_wide$hyb<-NA
for(i in 1:length(s8_hyb_wide$V2)){
  s8_hyb_wide$hybCount[i]<-sum(s8_hyb_wide[i,2:8]>.001&s8_hyb_wide[i,2:8]<.999,na.rm=TRUE)
  s8_hyb_wide$pureCount[i]<-sum(s8_hyb_wide[i,2:8]>.999,na.rm=TRUE)
  if(s8_hyb_wide$hybCount[i]> 1 &s8_hyb_wide$pureCount[i] < 1){
    s8_hyb_wide$hyb[i]<-"hybrid"
  }
  if(s8_hyb_wide$hybCount[i]<2 & s8_hyb_wide$pureCount[i] < 1){
    s8_hyb_wide$hyb[i]<-"uncertainPure"
  }
  if(s8_hyb_wide$hybCount[i]<1&s8_hyb_wide$pureCount[i]==1){
    s8_hyb_wide$hyb[i]<-"pure"
  }
  if(s8_hyb_wide$pureCount[i]>1){
    s8_hyb_wide$hyb[i]<-"?"
  }
  if(s8_hyb_wide$pureCount[i]>0&s8_hyb_wide$hybCount[i]>0){
    s8_hyb_wide$hyb[i]<-"uncertainHybrid?"
  }
}


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
    s8_hyb_wide$hyb[i]<-"uncertainPure"
  }
  if(s8_hyb_wide$hybCount[i]<1&s8_hyb_wide$pureCount[i]==1){
    s8_hyb_wide$hyb[i]<-"pure"
  }
  if(s8_hyb_wide$pureCount[i]>1){
    s8_hyb_wide$hyb[i]<-"?"
  }
  if(s8_hyb_wide$pureCount[i]>0&s8_hyb_wide$hybCount[i]>0){
    s8_hyb_wide$hyb[i]<-"uncertainHybrid?"
  }
}