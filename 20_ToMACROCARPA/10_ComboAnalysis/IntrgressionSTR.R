###################################
## Dstat introgression ############
###################################

#################################################
## Macrocarpa Reference #########################
#################################################
#combining all species dstatistics+metadata into one spreadsheet for all species
 # Dlob_All<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/107_Dsuite_MACREF/DPOPS/bic_mac_lob_2000.txt",sep="\t")
 # Dlob_All$species<-rep("lob",length(Dlob_All$chr))
 # Dalb_All<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/107_Dsuite_MACREF/DPOPS/bic_mac_alb_2000.txt",sep="\t")
 # Dalb_All$species<-rep("alb",length(Dalb_All$chr))
 # Dste_All<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/107_Dsuite_MACREF/DPOPS/bic_mac_ste_2000.txt",sep="\t")
 # Dste_All$species<-rep("ste",length(Dste_All$chr))
 # Dmue_All<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/107_Dsuite_MACREF/DPOPS/bic_mac_mue_2000.txt",sep="\t")
 # Dmue_All$species<-rep("mue",length(Dmue_All$chr))
# 
# 
# D_All<-rbind(Dalb_All,Dlob_All,Dmue_All,Dste_All)
# D_All$windowMid<-rowMeans(D_All[,2:3])
# D_All$chr2<-D_All$chr
# D_All<-D_All %>%
#   dplyr::mutate(across('chr2',str_replace,'Chr',''))
# D_All$chr3<-as.numeric(D_All$chr2)
# write.table(D_All,file="C:/Users/rmohn/Desktop/10_Analysis/107_Dsuite_MACREF/DPOPS/Dall_mac_2000.txt",row.names=FALSE,quote=FALSE,sep="\t")

# Dbi_All$chr2<-Dbi_All$chr
# Dbi_All<-Dbi_All %>%
#   dplyr::mutate(across('chr2',str_replace,'Chr',''))
# Dbi_All$chr3<-as.numeric(Dbi_All$chr2)

##for 500,000 bp windows, saiving the max d_f for each population (this helps with visualizing)

D_All_2000<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/107_Dsuite_MACREF/DPOPS/Dall_mac_2000.txt",header = TRUE)

D25_mac_2000<-subset(D_All_2000,d_f>.25)

# Dsp_All_2000<-data.frame(pop=character(),chr=character(),windowStart=integer(),windowEnd=integer(),D=double(),f_d=double(),f_dM=double(),d_f=double(),lat=double(),long=double(),windowmid=integer(),species=character(),chr2=character(),chr3=integer())
# for(j in 1:12){
#   i=min(D_All_2000$windowmid)
#   while(i <= max(D_All_2000$windowmid)){
#     if(i %in% subset(Dalb_All,chr3==j)$windowmid){
#       Dwind<-rbind(subset(Dalb_All,windowmid==i&chr3==j),subset(Dlob_All,windowmid==i&chr3==j),subset(Dmue_All,windowmid==i&chr3==j),subset(Dste_All,windowmid==i&chr3==j))
#       Dsp_All_2000<-rbind(Dsp_All_2000,(Dwind %>% group_by(pop) %>% top_n(1, d_f)))
#     }
#     i=i+500000
#   }
# }

# Dsp_All_2000$g10<-Dsp_All_2000$d_f>0.25
# 
# Dsp_agg<-aggregate(g10~chr+windowmid+pop,data=Dsp_All_2000,FUN=max)
# Dsp_Sig<-aggregate(as.numeric(g10)~chr+windowmid,data=Dsp_agg,FUN=sum)
# 
# ggplot()+
#   geom_histogram(mapping=aes(x=Dsp_Sig$`as.numeric(g10)`),bins = 12)
# 
# length(subset(Dsp_Sig,`as.numeric(g10)`>=4)$chr)

###########################################################
# Introgression and genome structure ######################
###########################################################



library(readxl)
library(dplyr)
library(ggplot2)
#genomic regions
ChrStr<-read_excel("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/MacGenomeRegionsB.xlsx",sheet="Sheet3")
#assigning colors for use later
group.colors<-c(lyr="#000000",bic="#88CCEE",lob="#116644",mue="#AA5599",ste="#E69F00",alb="#332288",chromosome="black",none="white",CenArray="gray20",InterArray="gray90", IntraSpInv="gray70",InterSpInv="gray35",inversion="gray30")

#renaming column for ease of use
ChrStr$csome<-ChrStr$scaffold

#importing chromosome dimensions
C_Coords<-read_excel("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/mac_csome_dims.xlsx",sheet="Sheet1")

csomeSamp<-rep(C_Coords$chr,round(C_Coords$stop/100000,0))
invs<-read.delim(file="C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/Invs2b.txt")


InvITRGTest<-data.frame(invSize=integer(),samplenumber=integer(),invCsome=character(),invEnd=integer(),ITRGSite=integer())
for(j in unique(invs$span_bp.x)){
  randCsome<-sample(csomeSamp,200)
  for(i in 1:length(randCsome)){
    randSamps<-round(runif(min=j,max=C_Coords[C_Coords[,1]==randCsome[i],3][[1]],1),0)
    ITRGSites<-subset(D25_mac_2000,chr==randCsome[i]&windowEnd>randSamps-j&windowStart<randSamps&d_f>0.25)
    if(length(ITRGSites$windowmid>0)){
      albsigwins<-0
      muesigwins<-0
      stesigwins<-0
      lobsigwins<-0
      for(z in c((randSamps-j)+seq(from=5000,to=j,by=10000))){
        if(length(subset(ITRGSites, windowStart<z&windowEnd<z&species=="mue")$species)>3){muesigwins<-muesigwins+1}
        if(length(subset(ITRGSites, windowStart<z&windowEnd<z&species=="alb")$species)>3){albsigwins<-albsigwins+1}
        if(length(subset(ITRGSites, windowStart<z&windowEnd<z&species=="lob")$species)>3){lobsigwins<-lobsigwins+1}
        if(length(subset(ITRGSites, windowStart<z&windowEnd<z&species=="ste")$species)>3){stesigwins<-stesigwins+1}
      }
    }else{
      albsigwins<-0
      muesigwins<-0
      stesigwins<-0
      lobsigwins<-0
    }
    InvITRGTest<-rbind(cbind(invSize=j,samplenumber=i,invCsome=randCsome[i],invEnd=randSamps,albsigwins=as.numeric(albsigwins),muesigwins=as.numeric(muesigwins),stesigwins=as.numeric(stesigwins),lobsigwins=as.numeric(lobsigwins)),InvITRGTest) #sites=SNPsites,
  }
}

InvITRGTest[is.na(InvITRGTest[,6]),6]<-0
InvITRGTest[is.na(InvITRGTest[,7]),7]<-0
InvITRGTest[is.na(InvITRGTest[,8]),8]<-0
InvITRGTest[is.na(InvITRGTest[,9]),9]<-0
#InvITRGTest$palb<-InvITRGTest$alb/InvITRGTest$sites
#InvITRGTest$plob<-InvITRGTest$lob/InvITRGTest$sites
#InvITRGTest$pmue<-InvITRGTest$mue/InvITRGTest$sites
#InvITRGTest$pste<-InvITRGTest$ste/InvITRGTest$sites


####
InvITRGTest<-bind_rows(cbind(invSize=j,samplenumber=i,invCsome=randCsome[i],invEnd=randSamps,ITRGSite=Pops,INTS),InvITRGTest) #sites=SNPsites,


InvITRG<-data.frame(invSize=integer(),invnum=integer(),invCsome=character(),invEnd=integer(),albsigwins=numeric(),muesigwins=numeric(),lobsigwins=numeric(),stesigwins=numeric(),type=character())
for(j in 1:length(invs$chr_start)){
  #INTS<-data.frame(alb=0,lob=0,mue=0,ste=0)
  ITRGSites<-data.frame()
  ITRGSites<-subset(D25_mac_2000,chr==invs$chromosome[j]&windowEnd>invs$start[j]&windowStart<invs$end[j]&d_f>0.25)
  #SNPsites<-length(subset(Dsp_All_2000,csome==invs$chromosome[j]&V4>invs$start[j]&V4<invs$end[j])$V4)
  if(length(ITRGSites$windowmid>0)){
    albsigwins<-0
    muesigwins<-0
    stesigwins<-0
    lobsigwins<-0
    for(z in c((invs$start[j])+seq(from=5000,to=invs$span_bp.x[j],by=10000))){
      if(length(subset(ITRGSites, windowStart<z&windowEnd>z&species=="mue")$species)>3){muesigwins<-muesigwins+1}
      if(length(subset(ITRGSites, windowStart<z&windowEnd>z&species=="alb")$species)>3){albsigwins<-albsigwins+1}
      if(length(subset(ITRGSites, windowStart<z&windowEnd>z&species=="lob")$species)>3){lobsigwins<-lobsigwins+1}
      if(length(subset(ITRGSites, windowStart<z&windowEnd>z&species=="ste")$species)>3){stesigwins<-stesigwins+1}
    }
    
  }else{
    albsigwins<-0
    muesigwins<-0
    stesigwins<-0
    lobsigwins<-0
  }
  InvITRG<-rbind(cbind(invSize=as.numeric(invs$span_bp.x[j]),invnum=as.numeric(j),invCsome=invs$chromosome[j],invEnd=as.numeric(invs$end[j]),albsigwins=as.numeric(albsigwins),muesigwins=as.numeric(muesigwins),stesigwins=as.numeric(stesigwins),lobsigwins=as.numeric(lobsigwins),type=invs$INV[j]),InvITRG)
}



InvITRG[is.na(InvITRG[,6]),6]<-0
InvITRG[is.na(InvITRG[,7]),7]<-0
InvITRG[is.na(InvITRG[,8]),8]<-0
InvITRG[is.na(InvITRG[,9]),9]<-0
# InvITRG$pe1<-InvITRG$e1/InvITRG$sites
# InvITRG$pe2<-InvITRG$e2/InvITRG$sites
# InvITRG$pe3<-InvITRG$e3/InvITRG$sites
# InvITRG$pe4<-InvITRG$e4/InvITRG$sites
# InvITRG$pe5<-InvITRG$e5/InvITRG$sites
# InvITRG$pe6<-InvITRG$e6/InvITRG$sites


InvITRG$qalb<-NA
InvITRG$qlob<-NA
InvITRG$qmue<-NA
InvITRG$qste<-NA


for(k in 1:length(InvITRG$invSize)){
  InvITRG$qalb[k]<-length(subset(InvITRGTest,as.numeric(invSize)==as.numeric(InvITRG$invSize[k])&as.numeric(albsigwins)>=as.numeric(InvITRG$albsigwins[k]))$albsigwins)/200
  InvITRG$qmue[k]<-length(subset(InvITRGTest,as.numeric(invSize)==as.numeric(InvITRG$invSize[k])&as.numeric(muesigwins)>=as.numeric(InvITRG$muesigwins[k]))$muesigwins)/200
  InvITRG$qlob[k]<-length(subset(InvITRGTest,as.numeric(invSize)==as.numeric(InvITRG$invSize[k])&as.numeric(lobsigwins)>=as.numeric(InvITRG$lobsigwins[k]))$lobsigwins)/200
  InvITRG$qste[k]<-length(subset(InvITRGTest,as.numeric(invSize)==as.numeric(InvITRG$invSize[k])&as.numeric(stesigwins)>=as.numeric(InvITRG$stesigwins[k]))$stesigwins)/200
}


InvITRG$qlessalb<-NA
InvITRG$qlesslob<-NA
InvITRG$qlessmue<-NA
InvITRG$qlessste<-NA


for(k in 1:length(InvITRG$invSize)){
  InvITRG$qlessalb[k]<-length(subset(InvITRGTest,as.numeric(invSize)==as.numeric(InvITRG$invSize[k])&as.numeric(albsigwins)<=as.numeric(InvITRG$albsigwins[k]))$albsigwins)/200
  InvITRG$qlessmue[k]<-length(subset(InvITRGTest,as.numeric(invSize)==as.numeric(InvITRG$invSize[k])&as.numeric(muesigwins)<=as.numeric(InvITRG$muesigwins[k]))$muesigwins)/200
  InvITRG$qlesslob[k]<-length(subset(InvITRGTest,as.numeric(invSize)==as.numeric(InvITRG$invSize[k])&as.numeric(lobsigwins)<=as.numeric(InvITRG$lobsigwins[k]))$lobsigwins)/200
  InvITRG$qlessste[k]<-length(subset(InvITRGTest,as.numeric(invSize)==as.numeric(InvITRG$invSize[k])&as.numeric(stesigwins)<=as.numeric(InvITRG$stesigwins[k]))$stesigwins)/200
}


ggplot(data=InvITRG)+
  geom_violin(mapping=aes(x="alb",y=qalb))+
  geom_jitter(mapping=aes(x="alb",y=qalb,size=invSize),height=0,alpha=.5)+
  geom_violin(mapping=aes(x="lob",y=qlob))+
  geom_jitter(mapping=aes(x="lob",y=qlob,size=invSize),height=0,alpha=.5)+
  geom_violin(mapping=aes(x="mue",y=qmue))+
  geom_jitter(mapping=aes(x="mue",y=qmue,size=invSize),height=0,alpha=.5)+
  geom_violin(mapping=aes(x="ste",y=qste,size=invSize))+
  geom_jitter(mapping=aes(x="ste",y=qste,size=invSize),height=0,alpha=.5)+
  facet_grid(type~.)+
  labs(y="Quantile of introgression for inversion size")+
  theme_light()

#######
## Plot the values against the models to see if they make sense
#####



albinvSize<-ggplot(data=subset(InvITRGTest),aes(x=as.character(round(log(invSize),0)),y=alb))+
  geom_boxplot(alpha=.1)+
  #geom_quantile(quantiles=c(0.05,0.25,0.5,0.75,0.95))+
  #  ylim(0,0.1)+
  #geom_jitter(height=0,color="black",alpha=.1)+
  geom_jitter(data=subset(InvITRG),aes(x=as.character(round(log(invSize),0)),y=alb,color=qalb,shape=type),height=0,width=.4)+
  theme_light()+
  binned_scale(aesthetics = "color",
               scale_name = "stepsn", 
               palette = function(x) c("blue","orangered2","sienna2","orange","gold","yellow2"),
               breaks = c(0,0.025,0.25,0.5,0.75,0.975,1),
               limits = c(0, 1),
               show.limits = TRUE, 
               guide = "colorsteps",,name = "quantile"
  )+
  #scale_color_gradientn(colors=c("yellow","gold","gold2","orange","sienna2","tomato","red4","blue"),values = c(1,0.975,0.75,0.5,0.25,0.025,0))+
  labs(x="ln(inversion size)",y="pops x windows Q. alba")
mueinvSize<-ggplot(data=subset(InvITRGTest),aes(x=as.character(round(log(invSize),0)),y=mue))+
  geom_boxplot(alpha=.1)+
  #geom_quantile(quantiles=c(0.05,0.25,0.5,0.75,0.95))+
  #  ylim(0,0.1)+
  #geom_jitter(height=0,color="black",alpha=.1)+
  geom_jitter(data=subset(InvITRG),aes(x=as.character(round(log(invSize),0)),y=mue,color=qmue,shape=type),height=0,width=.4)+
  theme_light()+
  binned_scale(aesthetics = "color",
               scale_name = "stepsn", 
               palette = function(x) c("blue","orangered2","sienna2","orange","gold","yellow2"),
               breaks = c(0,0.025,0.25,0.5,0.75,0.975,1),
               limits = c(0, 1),
               show.limits = TRUE, 
               guide = "colorsteps",name = "quantile"
  )+
  #scale_color_gradientn(colors=c("yellow","gold","gold2","orange","sienna2","tomato","red4","blue"),values = c(1,0.975,0.75,0.5,0.25,0.025,0))+
  labs(x="ln(inversion size)",y="pops x windows Q. muehlenbergii")
steinvSize<-ggplot(data=subset(InvITRGTest),aes(x=as.character(round(log(invSize),0)),y=ste))+
  geom_boxplot(alpha=.1)+
  #geom_quantile(quantiles=c(0.05,0.25,0.5,0.75,0.95))+
  #  ylim(0,0.1)+
  #geom_jitter(height=0,color="black",alpha=.1)+
  geom_jitter(data=subset(InvITRG),aes(x=as.character(round(log(invSize),0)),y=ste,color=qste,shape=type),height=0,width=.4)+
  theme_light()+
  binned_scale(aesthetics = "color",
               scale_name = "stepsn", 
               palette = function(x) c("blue","orangered2","sienna2","orange","gold","yellow2"),
               breaks = c(0,0.025,0.25,0.5,0.75,0.975,1),
               limits = c(0, 1),
               show.limits = TRUE, 
               guide = "colorsteps",name = "quantile"
  )+
  #scale_color_gradientn(colors=c("yellow","gold","gold2","orange","sienna2","tomato","red4","blue"),values = c(1,0.975,0.75,0.5,0.25,0.025,0))+
  labs(x="ln(inversion size)",y="pops x windows Q. stellata")


lobinvSize<-ggplot(data=subset(InvITRGTest),aes(x=as.character(round(log(invSize),0)),y=lob))+
  geom_boxplot(alpha=.1)+
  #geom_quantile(quantiles=c(0.05,0.25,0.5,0.75,0.95))+
  #  ylim(0,0.1)+
  #geom_jitter(height=0,color="black",alpha=.1)+
  geom_jitter(data=subset(InvITRG),aes(x=as.character(round(log(invSize),0)),y=lob,color=qlob,shape=type),height=0,width=.4)+
  theme_light()+
  binned_scale(aesthetics = "color",
               scale_name = "stepsn", 
               palette = function(x) c("blue","orangered2","sienna2","orange","gold","yellow2"),
               breaks = c(0,0.025,0.25,0.5,0.75,0.975,1),
               limits = c(0, 1),
               show.limits = TRUE, 
               guide = "colorsteps",name = "quantile"
  )+
  #scale_color_gradientn(colors=c("yellow","gold","gold2","orange","sienna2","tomato","red4","blue"),values = c(1,0.975,0.75,0.5,0.25,0.025,0))+
  labs(x="ln(inversion size)",y="pops x windows Q. lobata")

library(ggpubr)
ggarrange(albinvSize,lobinvSize,mueinvSize,steinvSize,nrow=1,common.legend = TRUE,legend="right")

ggplot(data=InvITRGTest,aes(x=invSize,y=lob))+
  geom_bin_2d(binwidth=c(50000,10),alpha=0.5)+
  geom_quantile(quantiles=c(0.05,0.25,0.5,0.75,0.95))+
  #  ylim(0,0.1)+
  geom_point(data=InvITRG,aes(x=invSize,y=lob,color=type))+
  theme_light()+
  scale_fill_viridis_c()

ggplot(data=InvITRGTest,aes(x=invSize,y=ste))+
  geom_bin_2d(binwidth=c(500000,10),alpha=0.5)+
  geom_quantile(quantiles=c(0.05,0.25,0.5,0.75,0.95))+
  #  ylim(0,0.1)+
  geom_point(data=InvITRG,aes(x=invSize,y=ste,color=type))+
  theme_light()+
  scale_fill_viridis_c()

ggplot(data=InvITRGTest,aes(x=invSize,y=mue))+
  geom_bin_2d(binwidth=c(500000,10),alpha=0.5)+
  geom_quantile(quantiles=c(0.05,0.25,0.5,0.75,0.95))+
  #  ylim(0,0.1)+
  geom_point(data=InvITRG,aes(x=invSize,y=mue,color=type))+
  theme_light()+
  scale_fill_viridis_c()


table(subset(InvITRG,alb<=0.025|mue<=0.025|ste<=0.025|lob<=0.025)$type)


####Do for centromere interarray regions
library(readxl)
library(dplyr)
library(ggplot2)
#genomic regions


##############################
### Centromeres ##############
##############################

cents<-read_excel("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/MacGenomeRegionsB.xlsx",sheet="InterArrays")

###test centromeres
centITRGTest<-data.frame(centSize=integer(),samplenumber=integer(),centCsome=character(),centEnd=integer(),ITRGSite=integer())
for(j in unique(cents$Length)){
  randCsome<-sample(csomeSamp,200)
  for(i in 1:length(randCsome)){
    randSamps<-round(runif(min=j,max=C_Coords[C_Coords[,1]==randCsome[i],3][[1]],1),0)
    ITRGSites<-subset(D25_mac_2000,chr==randCsome[i]&windowEnd>randSamps-j&windowStart<randSamps&d_f>0.25)
    if(length(ITRGSites$windowmid>0)){
      albsigwins<-0
      muesigwins<-0
      stesigwins<-0
      lobsigwins<-0
      for(z in c((randSamps-j)+seq(from=5000,to=j,by=10000))){
        if(length(subset(ITRGSites, windowStart<z&windowEnd<z&species=="mue")$species)>3){muesigwins<-muesigwins+1}
        if(length(subset(ITRGSites, windowStart<z&windowEnd<z&species=="alb")$species)>3){albsigwins<-albsigwins+1}
        if(length(subset(ITRGSites, windowStart<z&windowEnd<z&species=="lob")$species)>3){lobsigwins<-lobsigwins+1}
        if(length(subset(ITRGSites, windowStart<z&windowEnd<z&species=="ste")$species)>3){stesigwins<-stesigwins+1}
      }
    }else{
      albsigwins<-0
      muesigwins<-0
      stesigwins<-0
      lobsigwins<-0
    }
    centITRGTest<-rbind(cbind(centSize=j,samplenumber=i,centCsome=randCsome[i],centEnd=randSamps,albsigwins=as.numeric(albsigwins),muesigwins=as.numeric(muesigwins),stesigwins=as.numeric(stesigwins),lobsigwins=as.numeric(lobsigwins)),centITRGTest) #sites=SNPsites,
  }
}

#centITRGTest$palb<-centITRGTest$alb/centITRGTest$sites
#centITRGTest$plob<-centITRGTest$lob/centITRGTest$sites
#centITRGTest$pmue<-centITRGTest$mue/centITRGTest$sites
#centITRGTest$pste<-centITRGTest$ste/centITRGTest$sites


####
#centITRGTest<-bind_rows(cbind(centSize=j,samplenumber=i,centCsome=randCsome[i],centEnd=randSamps,ITRGSite=Pops,INTS),centITRGTest) #sites=SNPsites,


centITRG<-data.frame(centSize=integer(),samplenumber=integer(),centCsome=character(),centEnd=integer(),ITRGSites=integer())
for(j in 1:length(cents$Start)){
  #INTS<-data.frame(alb=0,lob=0,mue=0,ste=0)
  ITRGSites<-data.frame()
  ITRGSites<-subset(D25_mac_2000,chr==cents$scaffold[j]&windowEnd>cents$Start[j]&windowStart<cents$Stop[j]&d_f>0.25)
  #SNPsites<-length(subset(Dsp_All_2000,csome==invs$chromosome[j]&V4>invs$start[j]&V4<invs$end[j])$V4)
  if(length(ITRGSites$windowmid>0)){
    albsigwins<-0
    muesigwins<-0
    stesigwins<-0
    lobsigwins<-0
    for(z in c((cents$Start[j])+seq(from=5000,to=cents$Length[j],by=10000))){
      if(length(subset(ITRGSites, windowStart<z&windowEnd>z&species=="mue")$species)>3){muesigwins<-muesigwins+1}
      if(length(subset(ITRGSites, windowStart<z&windowEnd>z&species=="alb")$species)>3){albsigwins<-albsigwins+1}
      if(length(subset(ITRGSites, windowStart<z&windowEnd>z&species=="lob")$species)>3){lobsigwins<-lobsigwins+1}
      if(length(subset(ITRGSites, windowStart<z&windowEnd>z&species=="ste")$species)>3){stesigwins<-stesigwins+1}
    }
    
  }else{
    albsigwins<-0
    muesigwins<-0
    stesigwins<-0
    lobsigwins<-0
  }
  centITRG<-rbind(cbind(centSize=as.numeric(cents$Length[j]),centnum=as.numeric(j),centCsome=cents$scaffold[j],centEnd=as.numeric(cents$Stop[j]),albsigwins=as.numeric(albsigwins),muesigwins=as.numeric(muesigwins),stesigwins=as.numeric(stesigwins),lobsigwins=as.numeric(lobsigwins),type="centInterArray"),centITRG)
}



centITRG[is.na(centITRG[,6]),6]<-0
centITRG[is.na(centITRG[,7]),7]<-0
centITRG[is.na(centITRG[,8]),8]<-0
centITRG[is.na(centITRG[,9]),9]<-0
# centITRG$pe1<-centITRG$e1/centITRG$sites
# centITRG$pe2<-centITRG$e2/centITRG$sites
# centITRG$pe3<-centITRG$e3/centITRG$sites
# centITRG$pe4<-centITRG$e4/centITRG$sites
# centITRG$pe5<-centITRG$e5/centITRG$sites
# centITRG$pe6<-centITRG$e6/centITRG$sites


centITRG$qalb<-NA
centITRG$qlob<-NA
centITRG$qmue<-NA
centITRG$qste<-NA


for(k in 1:length(centITRG$centSize)){
  centITRG$qalb[k]<-length(subset(centITRGTest,centSize==centITRG$centSize[k]&alb>=centITRG$alb[k])$alb)/200
  centITRG$qlob[k]<-length(subset(centITRGTest,centSize==centITRG$centSize[k]&lob>=centITRG$lob[k])$lob)/200
  centITRG$qmue[k]<-length(subset(centITRGTest,centSize==centITRG$centSize[k]&mue>=centITRG$mue[k])$mue)/200
  centITRG$qste[k]<-length(subset(centITRGTest,centSize==centITRG$centSize[k]&ste>=centITRG$ste[k])$ste)/200
}

for(k in 1:length(centITRG$centSize)){
  centITRG$qalb[k]<-length(subset(centITRGTest,as.numeric(centSize)==as.numeric(centITRG$centSize[k])&as.numeric(albsigwins)>=as.numeric(centITRG$albsigwins[k]))$albsigwins)/200
  centITRG$qmue[k]<-length(subset(centITRGTest,as.numeric(centSize)==as.numeric(centITRG$centSize[k])&as.numeric(muesigwins)>=as.numeric(centITRG$muesigwins[k]))$muesigwins)/200
  centITRG$qlob[k]<-length(subset(centITRGTest,as.numeric(centSize)==as.numeric(centITRG$centSize[k])&as.numeric(lobsigwins)>=as.numeric(centITRG$lobsigwins[k]))$lobsigwins)/200
  centITRG$qste[k]<-length(subset(centITRGTest,as.numeric(centSize)==as.numeric(centITRG$centSize[k])&as.numeric(stesigwins)>=as.numeric(centITRG$stesigwins[k]))$stesigwins)/200
}


albcentSize<-ggplot(data=subset(centITRGTest),aes(x=as.character(round(log(centSize),0)),y=alb))+
  geom_boxplot(alpha=.1)+
  #geom_quantile(quantiles=c(0.05,0.25,0.5,0.75,0.95))+
  #  ylim(0,0.1)+
  #geom_jitter(height=0,color="black",alpha=.1)+
  geom_jitter(data=subset(centITRG),aes(x=as.character(round(log(centSize),0)),y=alb,color=qalb,shape=type),height=0,width=.4)+
  theme_light()+
  binned_scale(aesthetics = "color",
               scale_name = "stepsn", 
               palette = function(x) c("blue","orangered2","sienna2","orange","gold","yellow2"),
               breaks = c(0,0.025,0.25,0.5,0.75,0.975,1),
               limits = c(0, 1),
               show.limits = TRUE, 
               guide = "colorsteps",,name = "quantile"
  )+
  #scale_color_gradientn(colors=c("yellow","gold","gold2","orange","sienna2","tomato","red4","blue"),values = c(1,0.975,0.75,0.5,0.25,0.025,0))+
  labs(x="ln(centromere inter-array size)",y="pops x windows Q. alba")

muecentSize<-ggplot(data=subset(centITRGTest),aes(x=as.character(round(log(centSize),0)),y=mue))+
  geom_boxplot(alpha=.1)+
  #geom_quantile(quantiles=c(0.05,0.25,0.5,0.75,0.95))+
  #  ylim(0,0.1)+
  #geom_jitter(height=0,color="black",alpha=.1)+
  geom_jitter(data=subset(centITRG),aes(x=as.character(round(log(centSize),0)),y=mue,color=qmue,shape=type),height=0,width=.4)+
  theme_light()+
  binned_scale(aesthetics = "color",
               scale_name = "stepsn", 
               palette = function(x) c("blue","orangered2","sienna2","orange","gold","yellow2"),
               breaks = c(0,0.025,0.25,0.5,0.75,0.975,1),
               limits = c(0, 1),
               show.limits = TRUE, 
               guide = "colorsteps",name = "quantile"
  )+
  #scale_color_gradientn(colors=c("yellow","gold","gold2","orange","sienna2","tomato","red4","blue"),values = c(1,0.975,0.75,0.5,0.25,0.025,0))+
  labs(x="ln(centromere inter-array size)",y="pops x windows Q. muehlenbergii")

stecentSize<-ggplot(data=subset(centITRGTest),aes(x=as.character(round(log(centSize),0)),y=ste))+
  geom_boxplot(alpha=.1)+
  #geom_quantile(quantiles=c(0.05,0.25,0.5,0.75,0.95))+
  #  ylim(0,0.1)+
  #geom_jitter(height=0,color="black",alpha=.1)+
  geom_jitter(data=subset(centITRG),aes(x=as.character(round(log(centSize),0)),y=ste,color=qste,shape=type),height=0,width=.4)+
  theme_light()+
  binned_scale(aesthetics = "color",
               scale_name = "stepsn", 
               palette = function(x) c("blue","orangered2","sienna2","orange","gold","yellow2"),
               breaks = c(0,0.025,0.25,0.5,0.75,0.975,1),
               limits = c(0, 1),
               show.limits = TRUE, 
               guide = "colorsteps",name = "quantile"
  )+
  #scale_color_gradientn(colors=c("yellow","gold","gold2","orange","sienna2","tomato","red4","blue"),values = c(1,0.975,0.75,0.5,0.25,0.025,0))+
  labs(x="ln(centromere inter-array)",y="pops x windows Q. stellata")


lobcentSize<-ggplot(data=subset(centITRGTest),aes(x=as.character(round(log(centSize),0)),y=lob))+
  geom_boxplot(alpha=.1)+
  #geom_quantile(quantiles=c(0.05,0.25,0.5,0.75,0.95))+
  #  ylim(0,0.1)+
  #geom_jitter(height=0,color="black",alpha=.1)+
  geom_jitter(data=subset(centITRG),aes(x=as.character(round(log(centSize),0)),y=lob,color=qlob,shape=type),height=0,width=.4)+
  theme_light()+
  binned_scale(aesthetics = "color",
               scale_name = "stepsn", 
               palette = function(x) c("blue","orangered2","sienna2","orange","gold","yellow2"),
               breaks = c(0,0.025,0.25,0.5,0.75,0.975,1),
               limits = c(0, 1),
               show.limits = TRUE, 
               guide = "colorsteps",name = "quantile"
  )+
  #scale_color_gradientn(colors=c("yellow","gold","gold2","orange","sienna2","tomato","red4","blue"),values = c(1,0.975,0.75,0.5,0.25,0.025,0))+
  labs(x="ln(centromere inter-array size)",y="pops x windows Q. lobata")

library(ggpubr)
ggarrange(albcentSize,lobcentSize,muecentSize,stecentSize,nrow=1,common.legend = TRUE,legend="right")

ggarrange(albinvSize,lobinvSize,mueinvSize,steinvSize,albcentSize,lobcentSize,muecentSize,stecentSize,nrow=2,ncol=4,common.legend = TRUE,legend="right",align = "v")










########################################
### Chromosome Figure #########################
#################################
D_All_2000<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/107_Dsuite_MACREF/DPOPS/Dall_mac_2000.txt",sep="\t")


Dsp_All_2000<-data.frame(pop=character(),chr=character(),windowStart=integer(),windowEnd=integer(),D=double(),f_d=double(),f_dM=double(),d_f=double(),lat=double(),long=double(),windowmid=integer(),species=character(),chr2=character(),chr3=integer())
for(j in 1:12){
  i=min(D_All$windowmid)
  while(i <= max(D_All$windowmid)){
    if(i %in% subset(Dalb_All,chr3==j)$windowmid){
      Dwind<-rbind(subset(Dalb_All,windowmid==i&chr3==j),subset(Dlob_All,windowmid==i&chr3==j),subset(Dmue_All,windowmid==i&chr3==j),subset(Dste_All,windowmid==i&chr3==j))
      Dsp_All_2000<-rbind(Dsp_All_2000,(Dwind %>% group_by(pop) %>% top_n(1, d_f)))
    }
    i=i+500000
  }
}


library(readxl)
group.colors<-c(lyr="#000000",lob="#88CCEE",bic="#116644",mue="#AA5599",ste="#E69F00",alb="#332288",chromosome="black",none="white",CenArray="gray20",centromere="darkcyan", IntraSpInv="blue4",InterSpInv="blue4",inversion="darkred")

popData<-read_excel("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/metadata_pops.xlsx",sheet="metadata_pops")
Dsp_NS_all<-merge(Dsp_All_2000,popData,by.x="pop",by.y="Pop")

cents<-read_excel("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/MacGenomeRegionsB.xlsx",sheet="CentRegion")
cents$chr<-cents$scaffold
invs<-read.delim(file="C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/Invs2b.txt")
invs$chr<-invs$chromosome
C_Coords<-read_excel("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/mac_csome_dims.xlsx",sheet="Sheet1")

Dsp_NS_all$count<-1
Dsp_Nall<-aggregate(count~chr+windowmid+species,data=subset(Dsp_NS_all,d_f>.25&MacPop=="N"),FUN=sum)
Dsp_Sall<-aggregate(count~chr+windowmid+species,data=subset(Dsp_NS_all,d_f>.25&MacPop=="S"),FUN=sum)


ggplot(data=Dsp_NS_all)+
  
  #geom_histogram(data=subset(Dsp_NS_all,d_f>.25&MacPop=="N"),mapping=aes(x=windowmid/1000000,fill=species),width=.5,binwidth = .5,position = "stack")+
  geom_col(data=Dsp_Sall,mapping=aes(x=windowmid/1000000,y=-count,fill=species),width=.5,position = "stack")+
  geom_col(data=Dsp_Nall,mapping=aes(x=windowmid/1000000,y=count,fill=species),width=.5,position = "stack")+
  
  geom_rect(data=C_Coords,aes(ymin=-.5,ymax=+.5,xmin=start/1000000,xmax=stop/1000000,fill="chromosome"),alpha=1)+
  
  geom_rect(data=cents,aes(ymin=-31,ymax=-28,xmin=Start/1000000,xmax=Stop/1000000,fill="centromere"),alpha=1)+
  geom_rect(data=invs,aes(ymin=28,ymax=31,xmin=start/1000000,xmax=end/1000000,fill="inversion"),alpha=1)+
  #geom_col(data=subset(highTD,chr=="Chr01"|chr=="Chr02"|chr=="Chr03"|chr=="Chr04"|chr=="Chr5"|chr=="Chr06"),mapping=aes(x=as.integer(START),y=-numpops),width=500000,alpha=.3)+
  #geom_vline(data=subset(Dsp_All_2000,chr3==j),mapping=aes(xintercept=i/1000000),color="red")+
  geom_blank(data=Dsp_All_2000,aes(x=windowmid/1000000))+
  scale_fill_manual(values=group.colors)+
  facet_grid(.~chr)+
  theme_light()+
  labs(x="Genomic window (Mbp)",y="Pops")+
  scale_x_continuous(breaks = seq(0,100,20)) +
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank(),panel.spacing.x=unit(0.2, "lines"),text = element_text(size = 8),panel.border = element_rect(color = "black", fill = NA, linewidth = .5))+
  #coord_cartesian(xlim = c(0, NA), expand = FALSE)+
  coord_flip()#xlim(c(0,NA))+
  
# ggplot(data=Dsp_All_2000)+
#   geom_histogram(data=subset(Dsp_NS_all,d_f>.25&MacPop=="N"),mapping=aes(y=windowmid/1000000,fill=species),width=.5,binwidth = .5,position = "stack")+
#   geom_col(data=Dsp_Sall,mapping=aes(y=windowmid/1000000,x=-count-12,fill=species),position = "stack")+
#   
#   geom_rect(data=C_Coords,aes(xmin=-6,xmax=0,ymin=start/1000000,ymax=stop/1000000,fill="chromosome"),alpha=1)+
#   
#   geom_rect(data=cents,aes(xmin=-8,xmax=-1,ymin=Start/1000000,ymax=Stop/1000000,fill="InterArray"),alpha=1)+
#   geom_rect(data=invs,aes(xmin=-15,xmax=-8,ymin=start/1000000,ymax=end/1000000,fill="inversion"),alpha=1)+
#   
#   #geom_histogram(data=subset(Dsp_All_2000,d_f>.25),mapping=aes(y=windowmid/1000000,fill=species),width=.5,binwidth = .5,position = "stack")+
#   #geom_col(data=subset(highTD,chr=="Chr01"|chr=="Chr02"|chr=="Chr03"|chr=="Chr04"|chr=="Chr5"|chr=="Chr06"),mapping=aes(y=as.integer(START),x=-numpops),width=500000,alpha=.3)+
#   #geom_vline(data=subset(Dsp_All_2000,chr3==j),mapping=aes(yintercept=i/1000000),color="red")+
#   geom_blank(data=Dsp_All_2000,aes(y=windowmid/1000000))+
#   scale_fill_manual(values=group.colors)+
#   facet_grid(.~chr, scales = "free_y",space="free_y")+
#   theme_light()+
#   labs(y="Genomic window (Mbp)",x="Pops")+
#   scale_y_continuous(breaks = seq(0,100,20)) +
#   theme(panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank(),panel.grid.minor.x=element_blank(),panel.spacing.y=unit(0, "lines"),text = element_text(size = 8),panel.border = element_rect(color = "black", fill = NA, linewidth = .5))+
#   coord_cartesian(ylim = c(0, NA), xlim=c(-15,55),expand = FALSE)#xlim(c(0,NA))+
# 


#################################################
## Mongolica Reference ##########################
#################################################

#combining all species dstatistics+metadata into one spreadsheet for all species
Dlob_AllMong<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_lobata.txt",sep="\t")
Dlob_AllMong$species<-rep("lob",length(Dlob_AllMong$chr))
Dalb_AllMong<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_alba.txt",sep="\t")
Dalb_AllMong$species<-rep("alb",length(Dalb_AllMong$chr))
Dste_AllMong<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_stellata.txt",sep="\t")
Dste_AllMong$species<-rep("ste",length(Dste_AllMong$chr))
Dmue_AllMong<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_muehlenbergii.txt",sep="\t")
Dmue_AllMong$species<-rep("mue",length(Dmue_AllMong$chr))

D_AllMong<-rbind(Dalb_AllMong,Dlob_AllMong,Dmue_AllMong,Dste_AllMong)
D_AllMong$windowMid<-rowMeans(D_AllMong[,2:3])


##for 500,000 bp windows, saiving the max d_f for each population (this helps with visualizing)
Dsp_All_2000Mong<-data.frame(pop=character(),chr=character(),windowStart=integer(),windowEnd=integer(),D=double(),f_d=double(),f_dM=double(),d_f=double(),lat=double(),long=double(),windowmid=integer(),species=character())
for(j in 1:12){
  i=min(D_AllMong$windowmid)
  while(i <= max(D_AllMong$windowmid)){
    if(i %in% subset(Dalb_AllMong,chr==paste("Superscaffold",j,sep=""))$windowmid){
      Dwind<-rbind(subset(Dalb_AllMong,windowmid==i&chr==paste("Superscaffold",j,sep="")),subset(Dlob_AllMong,windowmid==i&chr==paste("Superscaffold",j,sep="")),subset(Dmue_AllMong,windowmid==i&chr==paste("Superscaffold",j,sep="")),subset(Dste_AllMong,windowmid==i&chr==paste("Superscaffold",j,sep="")))
      Dsp_All_2000Mong<-rbind(Dsp_All_2000Mong,(Dwind %>% group_by(pop) %>% top_n(1, d_f)))
    }
    i=i+500000
  }
}


Dsp_All_2000Mong$g10<-Dsp_All_2000Mong$d_f>0.25

Dsp_aggMong<-aggregate(g10~chr+windowmid+pop,data=Dsp_All_2000Mong,FUN=max)
Dsp_SigMong<-aggregate(as.numeric(g10)~chr+windowmid,data=Dsp_agg,FUN=sum)

ggplot()+
  geom_histogram(mapping=aes(x=Dsp_SigMong$`as.numeric(g10)`),bins = 10)

length(subset(Dsp_SigMong,`as.numeric(g10)`>=4)$chr)

################################################
####### Mong PCA and dist range edge ###########
###############################################

library(readxl)
Oldpops<-read.csv("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/Pop_conversion.csv")
Dlob_Mong_All<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_lobata.txt",sep="\t")
Dlob_Mong_All$species<-rep("lob",length(Dlob_Mong_All$chr))



library(dplyr)
Dlob_Mong_All_win<-Dlob_Mong_All %>% group_by(pop,windowmid,chr) %>% top_n(1, d_f)
Dlob_Mong_All_sigwins<-subset(Dlob_Mong_All_win,d_f>0.25)
Dlob_Mong_All_sigwins$count<-1
Dlob_Mong_All_countsigwins<-aggregate(count~pop+lat,data=Dlob_Mong_All_sigwins,FUN=sum)
#merge distance to range edge
samps_edge_dist<-read.delim("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/dist_range_margin_centroid.txt",sep="\t")

Dlob_Mong_All_countsigwins<-merge(Dlob_Mong_All_countsigwins,Oldpops,by.x="pop",by.y="OLD_Pop")
Dlob_Mong_All_win_edge<-merge(Dlob_Mong_All_countsigwins,samps_edge_dist,by.x="NEW_Pop",by.y="Pop",all.x=T)
#do regression
lob_mong_distcent<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dlob_Mong_All_win_edge$count, x=Dlob_Mong_All_win_edge$dist_fromCenter_EM))+
  geom_point(mapping=aes(y = Dlob_Mong_All_win_edge$count, x=Dlob_Mong_All_win_edge$dist_fromCenter_EM,shape=Dlob_Mong_All_win_edge$MacPop.x, color=Dlob_Mong_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(x="distance from centroid",y="windows",subtitle="Q. lobata")

lob_mong_distedge<-ggplot()+
  geom_smooth(method=lm, mapping=aes(y = Dlob_Mong_All_win_edge$count, x=Dlob_Mong_All_win_edge$dist_fromEdge))+
  geom_jitter(mapping=aes(y = Dlob_Mong_All_win_edge$count, x=Dlob_Mong_All_win_edge$dist_fromEdge, shape=Dlob_Mong_All_win_edge$MacPop.x, color=Dlob_Mong_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(x="distance from edge",y="windows",subtitle="Q. lobata")

lob_mong_lat<-ggplot()+
  geom_smooth(method=lm, mapping=aes(y = Dlob_Mong_All_win_edge$count, x=Dlob_Mong_All_win_edge$latitude.orig))+
  geom_jitter(mapping=aes(y = Dlob_Mong_All_win_edge$count, x=Dlob_Mong_All_win_edge$latitude.orig, shape=Dlob_Mong_All_win_edge$MacPop.x, color=Dlob_Mong_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(x="latitude",y="windows",subtitle="Q. lobata")

##PCA
Dlob_Mong_All_win$chr_win<-paste(Dlob_Mong_All_win$chr,Dlob_Mong_All_win$windowmid,sep="_")
Dlob_Mong_All_win_small<-Dlob_Mong_All_win[,-c(1:6,9:12)]

library(tidyr)
Dlob_Mong_All_win_wide<-pivot_wider(Dlob_Mong_All_win_small,names_from=chr_win,values_from=d_f)
Dlob_Mong_All_win_wide<-Dlob_Mong_All_win_wide %>% 
  select(where(~!any(is.na(.))))


Dlob_Mong_PCA<-prcomp(Dlob_Mong_All_win_wide[,-1])
summary(Dlob_Mong_PCA)
Dlob_Mong_PCs<-as.data.frame(Dlob_Mong_PCA$x)
Dlob_Mong_PCs$Pop<-Dlob_Mong_All_win_wide$pop
Dlob_Mong_PCs<-merge(Dlob_Mong_PCs,Oldpops,by.x="Pop",by.y="OLD_Pop")
library(ggplot2)
lobMongPC1PC2<-ggplot(data=Dlob_Mong_PCs)+
  geom_point(aes(x=PC1, y=PC2,color=latitude.orig,size=longitude.orig,shape=MacPop),alpha=0.75)+
  scale_color_viridis_c()+
  labs(subtitle = "Q. lobata")

# library(ggplot2)
# ggplot(data=Dlob_Mong_PCs)+
#   geom_point(aes(x=PC1, y=latitude.orig))

#### Alba #####
library(readxl)


Dalb_Mong_All<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_alba.txt",sep="\t")
Dalb_Mong_All$species<-rep("alb",length(Dalb_Mong_All$chr))

library(dplyr)
Dalb_Mong_All_win<-Dalb_Mong_All %>% group_by(pop,windowmid,chr) %>% top_n(1, d_f)
Dalb_Mong_All_sigwins<-subset(Dalb_Mong_All_win,d_f>0.25)
Dalb_Mong_All_sigwins$count<-1
Dalb_Mong_All_countsigwins<-aggregate(count~pop+lat,data=Dalb_Mong_All_sigwins,FUN=sum)
#merge distance to range edge
samps_edge_dist<-read.delim("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/dist_range_margin_centroid.txt",sep="\t")

Dalb_Mong_All_countsigwins<-merge(Dalb_Mong_All_countsigwins,Oldpops,by.x="pop",by.y="OLD_Pop")
Dalb_Mong_All_win_edge<-merge(Dalb_Mong_All_countsigwins,samps_edge_dist,by.x="NEW_Pop",by.y="Pop",all.x=T)
#do regression
alb_mong_distcent<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dalb_Mong_All_win_edge$count, x=Dalb_Mong_All_win_edge$dist_fromCenter_EM))+
  geom_point(mapping=aes(y = Dalb_Mong_All_win_edge$count, x=Dalb_Mong_All_win_edge$dist_fromCenter_EM,shape=Dalb_Mong_All_win_edge$MacPop.x, color=Dalb_Mong_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(x="distance from centroid",y="windows",subtitle="Q. alba")

alb_mong_distedge<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dalb_Mong_All_win_edge$count, x=Dalb_Mong_All_win_edge$dist_fromEdge))+
  geom_jitter(mapping=aes(y = Dalb_Mong_All_win_edge$count, x=Dalb_Mong_All_win_edge$dist_fromEdge, shape=Dalb_Mong_All_win_edge$MacPop.x, color=Dalb_Mong_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(x="distance from edge",y="windows",subtitle="Q. alba")

alb_mong_lat<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dalb_Mong_All_win_edge$count, x=Dalb_Mong_All_win_edge$latitude.orig))+
  geom_jitter(mapping=aes(y = Dalb_Mong_All_win_edge$count, x=Dalb_Mong_All_win_edge$latitude.orig, shape=Dalb_Mong_All_win_edge$MacPop.x, color=Dalb_Mong_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(x="latitude",y="windows",subtitle="Q. alba")

##PCA
Dalb_Mong_All_win$chr_win<-paste(Dalb_Mong_All_win$chr,Dalb_Mong_All_win$windowmid,sep="_")
Dalb_Mong_All_win_small<-Dalb_Mong_All_win[,-c(1:6,9:12)]

library(tidyr)
Dalb_Mong_All_win_wide<-pivot_wider(Dalb_Mong_All_win_small,names_from=chr_win,values_from=d_f)
Dalb_Mong_All_win_wide<-Dalb_Mong_All_win_wide %>% 
  select(where(~!any(is.na(.))))


Dalb_Mong_PCA<-prcomp(Dalb_Mong_All_win_wide[,-1])
summary(Dalb_Mong_PCA)
Dalb_Mong_PCs<-as.data.frame(Dalb_Mong_PCA$x)
Dalb_Mong_PCs$Pop<-Dalb_Mong_All_win_wide$pop
Dalb_Mong_PCs<-merge(Dalb_Mong_PCs,Oldpops,by.x="Pop",by.y="OLD_Pop")
library(ggplot2)
albMongPC1PC2<-ggplot(data=Dalb_Mong_PCs)+
  geom_point(aes(x=PC1, y=PC2,color=latitude.orig,size=longitude.orig,shape=MacPop),alpha=0.75)+
  scale_color_viridis_c()+
  labs(subtitle = "Q. alba")

# library(ggplot2)
# ggplot(data=Dalb_Mong_PCs)+
#   geom_point(aes(x=PC1, y=latitude.orig))
# 
### Muehlenbergii ###
library(readxl)

Dmue_Mong_All<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_muehlenbergii.txt",sep="\t")
Dmue_Mong_All$species<-rep("mue",length(Dmue_Mong_All$chr))

library(dplyr)
Dmue_Mong_All_win<-Dmue_Mong_All %>% group_by(pop,windowmid,chr) %>% top_n(1, d_f)
Dmue_Mong_All_sigwins<-subset(Dmue_Mong_All_win,d_f>0.25)
Dmue_Mong_All_sigwins$count<-1
Dmue_Mong_All_countsigwins<-aggregate(count~pop+lat,data=Dmue_Mong_All_sigwins,FUN=sum)
#merge distance to range edge
samps_edge_dist<-read.delim("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/dist_range_margin_centroid.txt",sep="\t")

Dmue_Mong_All_countsigwins<-merge(Dmue_Mong_All_countsigwins,Oldpops,by.x="pop",by.y="OLD_Pop")
Dmue_Mong_All_win_edge<-merge(Dmue_Mong_All_countsigwins,samps_edge_dist,by.x="NEW_Pop",by.y="Pop",all.x=T)
#do regression
mue_mong_distcent<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dmue_Mong_All_win_edge$count, x=Dmue_Mong_All_win_edge$dist_fromCenter_EM))+
  geom_point(mapping=aes(y = Dmue_Mong_All_win_edge$count, x=Dmue_Mong_All_win_edge$dist_fromCenter_EM,shape=Dmue_Mong_All_win_edge$MacPop.x, color=Dmue_Mong_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(x="distance from centroid",y="windows",subtitle="Q. muehlenbergii")

mue_mong_distedge<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dmue_Mong_All_win_edge$count, x=Dmue_Mong_All_win_edge$dist_fromEdge))+
  geom_jitter(mapping=aes(y = Dmue_Mong_All_win_edge$count, x=Dmue_Mong_All_win_edge$dist_fromEdge, shape=Dmue_Mong_All_win_edge$MacPop.x, color=Dmue_Mong_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(x="distance from edge",y="windows",subtitle="Q. muehlenbergii")

mue_mong_lat<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dmue_Mong_All_win_edge$count, x=Dmue_Mong_All_win_edge$latitude.orig))+
  geom_jitter(mapping=aes(y = Dmue_Mong_All_win_edge$count, x=Dmue_Mong_All_win_edge$latitude.orig, shape=Dmue_Mong_All_win_edge$MacPop.x, color=Dmue_Mong_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(x="latitude",y="windows",subtitle="Q. muehlenbergii")

##PCA
Dmue_Mong_All_win$chr_win<-paste(Dmue_Mong_All_win$chr,Dmue_Mong_All_win$windowmid,sep="_")
Dmue_Mong_All_win_small<-Dmue_Mong_All_win[,-c(1:6,9:12)]

library(tidyr)
Dmue_Mong_All_win_wide<-pivot_wider(Dmue_Mong_All_win_small,names_from=chr_win,values_from=d_f)
Dmue_Mong_All_win_wide<-Dmue_Mong_All_win_wide %>% 
  select(where(~!any(is.na(.))))


Dmue_Mong_PCA<-prcomp(Dmue_Mong_All_win_wide[,-1])
summary(Dmue_Mong_PCA)
Dmue_Mong_PCs<-as.data.frame(Dmue_Mong_PCA$x)
Dmue_Mong_PCs$Pop<-Dmue_Mong_All_win_wide$pop
Dmue_Mong_PCs<-merge(Dmue_Mong_PCs,Oldpops,by.x="Pop",by.y="OLD_Pop")
library(ggplot2)
mueMongPC1PC2<-ggplot(data=Dmue_Mong_PCs)+
  geom_point(aes(x=PC1, y=PC2,color=latitude.orig,size=longitude.orig,shape=MacPop),alpha=0.75)+
  scale_color_viridis_c()+
  labs(subtitle = "Q. muehlenbergii")

# library(ggplot2)
# ggplot(data=Dmue_Mong_PCs)+
#   geom_point(aes(x=PC1, y=latitude.orig))

#### Stellata #####
library(readxl)

Dste_Mong_All<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_stellata.txt",sep="\t")
Dste_Mong_All$species<-rep("ste",length(Dste_Mong_All$chr))

library(dplyr)
Dste_Mong_All_win<-Dste_Mong_All %>% group_by(pop,windowmid,chr) %>% top_n(1, d_f)
Dste_Mong_All_sigwins<-subset(Dste_Mong_All_win,d_f>0.25)
Dste_Mong_All_sigwins$count<-1
Dste_Mong_All_countsigwins<-aggregate(count~pop+lat,data=Dste_Mong_All_sigwins,FUN=sum)
#merge distance to range edge
samps_edge_dist<-read.delim("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/dist_range_margin_centroid.txt",sep="\t")

Dste_Mong_All_countsigwins<-merge(Dste_Mong_All_countsigwins,Oldpops,by.x="pop",by.y="OLD_Pop")
Dste_Mong_All_win_edge<-merge(Dste_Mong_All_countsigwins,samps_edge_dist,by.x="NEW_Pop",by.y="Pop",all.x=T)
#do regression
ste_mong_distcent<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dste_Mong_All_win_edge$count, x=Dste_Mong_All_win_edge$dist_fromCenter_EM))+
  geom_point(mapping=aes(y = Dste_Mong_All_win_edge$count, x=Dste_Mong_All_win_edge$dist_fromCenter_EM,shape=Dste_Mong_All_win_edge$MacPop.x, color=Dste_Mong_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(x="distance from centroid",y="windows",subtitle="Q. stellata")

ste_mong_distedge<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dste_Mong_All_win_edge$count, x=Dste_Mong_All_win_edge$dist_fromEdge))+
  geom_jitter(mapping=aes(y = Dste_Mong_All_win_edge$count, x=Dste_Mong_All_win_edge$dist_fromEdge, shape=Dste_Mong_All_win_edge$MacPop.x, color=Dste_Mong_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(x="distance from edge",y="windows",subtitle="Q. stellata")

ste_mong_lat<-ggplot()+
  geom_smooth(method=lm,aes(y = Dste_Mong_All_win_edge$count,x=Dste_Mong_All_win_edge$latitude.orig))+
  geom_jitter(mapping=aes(y = Dste_Mong_All_win_edge$count, x=Dste_Mong_All_win_edge$latitude.orig, shape=Dste_Mong_All_win_edge$MacPop.x, color=Dste_Mong_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(x="latitude",y="windows",subtitle="Q. stellata")

##PCA
Dste_Mong_All_win$chr_win<-paste(Dste_Mong_All_win$chr,Dste_Mong_All_win$windowmid,sep="_")
Dste_Mong_All_win_small<-Dste_Mong_All_win[,-c(1:6,9:12)]

library(tidyr)
Dste_Mong_All_win_wide<-pivot_wider(Dste_Mong_All_win_small,names_from=chr_win,values_from=d_f)
Dste_Mong_All_win_wide<-Dste_Mong_All_win_wide %>% 
  select(where(~!any(is.na(.))))


Dste_Mong_PCA<-prcomp(Dste_Mong_All_win_wide[,-1])
summary(Dste_Mong_PCA)
Dste_Mong_PCs<-as.data.frame(Dste_Mong_PCA$x)
Dste_Mong_PCs$Pop<-Dste_Mong_All_win_wide$pop
Dste_Mong_PCs<-merge(Dste_Mong_PCs,Oldpops,by.x="Pop",by.y="OLD_Pop")
library(ggplot2)
steMongPC1PC2<-ggplot(data=Dste_Mong_PCs)+
  geom_point(aes(x=PC1, y=PC2,color=latitude.orig,size=longitude.orig,shape=MacPop),alpha=0.75)+
  scale_color_viridis_c()+
  labs(subtitle = "Q. stellata")

# library(ggplot2)
# ggplot(data=Dste_Mong_PCs)+
#   geom_point(aes(x=PC1, y=latitude.orig))


library(ggpubr)
ggarrange(albMongPC1PC2,lobMongPC1PC2,mueMongPC1PC2,steMongPC1PC2,ncol=2,nrow=2,common.legend = T,legend = "right")
ggarrange(alb_mong_distcent,lob_mong_distcent,mue_mong_distcent,ste_mong_distcent,common.legend = T,legend = "right")
ggarrange(alb_mong_distedge,lob_mong_distedge,mue_mong_distedge,ste_mong_distedge,common.legend = T,legend = "right")
ggarrange(alb_mong_lat,lob_mong_lat,mue_mong_lat,ste_mong_lat,common.legend = T,legend = "right")
