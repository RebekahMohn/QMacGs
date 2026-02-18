### look at the correlation between env. and df
library(biomartr)
library(ape)
library(phytools)
library(geiger)
library(readxl)
library(stringr)
library(ggpubr)

library(dplyr)
library(plyr)
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(adegenet)
library(data.table)
library(ggtree)
library(lme4)


states <- map_data("state")
counties <- map_data('county')
Canada<-map_data("world", "Canada")
group.colors<-c(mac="lightgreen",rub="gray40",lyr="#000000",bic="#88CCEE",
                lob="#116644",mue="#AA5599",pri="violet",ste="#E69F00",sin="gold",
                mar="gold3",alb="#332288",mic="steelblue",mon="steelblue4",e1="red4",
                e2="salmon",e3="gold4",e4="skyblue4",e5="blue4",e6="yellow2",`-`="white",
                CenArray="gray20",InterArray="gray80", intersp="gray70",hap="gray35",`+`='black')

gff_data<-read_gff("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/Hap1/quercus_macrocarpa_hap1_V1.0.gene.gff3")
adaptIntro<-read.delim(sep="\t","C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/adaptiveIntrogressed2_samp.tsv")
gene_annotation<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/Hap1/hap1_qmacrocarpa_EnTAP_annotated (1).tsv")
invs<-read.delim(file="C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/Invs2b.txt")
cents<-subset(read_excel("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/MacGenomeRegionsB.xlsx",sheet="Sheet1"),Type2=="Cent")
cents$csome<-cents$scaffold

BC6<-read.table("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/BioClim6vars.txt", sep="\t",header = T)
sampPop<-read.table("C:/Users/rmohn/Desktop/10_Analysis/107_DSuite_MACREF/SetsGR_mac2.tsv",sep="\t")
BC6_pops<-aggregate(cbind(wc2.1_30s_bio_2,wc2.1_30s_bio_6,wc2.1_30s_bio_8,wc2.1_30s_bio_13,wc2.1_30s_bio_14,wc2.1_30s_bio_18)~V2,data=merge(sampPop,BC6,by="V1"),FUN=mean)

gea_env<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_lfmmpval_cor.txt",sep=" ")


D_All_200<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/107_Dsuite_MACREF/DPOPS/Dall_mac_2000.txt",header = TRUE)


GEA_INT_DATABASE<-data.frame()

#

pdf("C:/Users/rmohn/Desktop/10_Analysis/120_Ad_Int/graphsOfRelationships2000b.pdf")
for(i in 1:length(adaptIntro$csome)){
  #Read the painting windows
  # mue_adInt<-read.delim(paste("C:/Users/rmohn/Desktop/10_Analysis/113_painting/AdIntWins/mue_",adaptIntro$V2[i],"_",adaptIntro$V1[i],".tsv",sep=""),sep="\t")
  # alb_adInt<-read.delim(paste("C:/Users/rmohn/Desktop/10_Analysis/113_painting/AdIntWins//alb_",adaptIntro$V2[i],"_",adaptIntro$V1[i],".tsv",sep=""),sep="\t")
  # lob_adInt<-read.delim(paste("C:/Users/rmohn/Desktop/10_Analysis/113_painting/AdIntWins/lob_",adaptIntro$V2[i],"_",adaptIntro$V1[i],".tsv",sep=""),sep="\t")
  # ste_adInt<-read.delim(paste("C:/Users/rmohn/Desktop/10_Analysis/113_painting/AdIntWins/ste_",adaptIntro$V2[i],"_",adaptIntro$V1[i],".tsv",sep=""),sep="\t")
  
  #Subset the GEA
  GeaSubset<-subset(gea_env,csome==adaptIntro$csome[i]&BP>adaptIntro$pos[i]-50000&BP<adaptIntro$pos[i]+50000)
  GeaSubset$gea_pos<-paste(GeaSubset$ENV,GeaSubset$BP,sep="_")
  
  #Subset the Dstats
  D_SS<-subset(D_All_2000,chr==adaptIntro$csome[i]&windowStart<adaptIntro$pos[i]+50000&windowEnd>adaptIntro$pos[i]-50000)
  D_SS$mac_pops<-paste("mac_",D_SS$pop,sep="")
  D_SS_GEA_MAX<-data.frame()
  D_SS_MAX_TAB<-data.frame()
  
  
  
  genesGEA<-data.frame()
  GEA_D<-data.frame(p_alb=numeric(),p_lob=numeric(),p_mue=numeric(),p_ste=numeric(),sigDGEA=character(),gea_pos=character())
  for(k in 1:length(GeaSubset$BP)){
    
    D_SS_GEA<-subset(D_SS,windowStart<GeaSubset$BP[k]&windowEnd>GeaSubset$BP[k])
    D_SS_GEA$distWind<-as.integer(D_SS_GEA$windowMid)-GeaSubset$BP[k]
    D_SS_GEA_2<-(D_SS_GEA %>% group_by(pop) %>% top_n(1, d_f))
    D_SS_GEA_MAX<-rbind(D_SS_GEA_2,D_SS_GEA_MAX)
    
    if(length(subset(D_SS_GEA_2,d_f>.25)$species)>0){
      D_SS_MAX_TAB<-as.data.frame(t(data.frame(table(subset(D_SS_GEA_2,d_f>.25)$species))[,2]))
      colnames(D_SS_MAX_TAB)<-t(data.frame(table(subset(D_SS_GEA_2,d_f>.25)$species))[,1])
      rownames(D_SS_MAX_TAB)<-NULL
    }

    
    if(GeaSubset$ENV[k]=="e1"){
      GEA_D_ASS<-merge(D_SS_GEA,BC6_pops[,1:2],by.x = "mac_pops",by.y="V2")
      
    } else if (GeaSubset$ENV[k]=="e2"){
      GEA_D_ASS<-merge(D_SS_GEA,BC6_pops[,c(1,3)],by.x = "mac_pops",by.y="V2")
      
    }  else if (GeaSubset$ENV[k]=="e3"){
      GEA_D_ASS<-merge(D_SS_GEA,BC6_pops[,c(1,4)],by.x = "mac_pops",by.y="V2")
      
    }  else if (GeaSubset$ENV[k]=="e4"){
      GEA_D_ASS<-merge(D_SS_GEA,BC6_pops[,c(1,5)],by.x = "mac_pops",by.y="V2")
      
    } else if (GeaSubset$ENV[k]=="e5"){
      GEA_D_ASS<-merge(D_SS_GEA,BC6_pops[,c(1,6)],by.x = "mac_pops",by.y="V2")
      
    } else if (GeaSubset$ENV[k]=="e6"){
      GEA_D_ASS<-merge(D_SS_GEA,BC6_pops[,c(1,7)],by.x = "mac_pops",by.y="V2")
      
    }
    ###The problem is that we need to differentiate which species it is
    if(length(D_SS_MAX_TAB)>0){
      pval_temp<-data.frame(row.names = c("GEA"))
      slope_GEA_temp<-data.frame(row.names = c("GEA"))
      for(l in colnames(D_SS_MAX_TAB)){
        SS_GEA_D<-subset(GEA_D_ASS,species==l)
        model1<-lm(SS_GEA_D$d_f~SS_GEA_D[,18])
        slope<-coefficients(model1)[2]
        pval_temp<-cbind(pval_temp,summary(model1)$coefficients[2,4])
        slope_GEA_temp<-cbind(slope_GEA_temp,slope)
        print(ggplot(data=SS_GEA_D,aes(x=d_f,y=SS_GEA_D[,18]))+ geom_point(aes(color=distWind))+geom_smooth(method="lm")+
               labs(caption=paste(slope,summary(model1)$coefficients[2,4],GeaSubset$ENV[k],GeaSubset$BP[k],GeaSubset$CHR, l)))
      }
      colnames(slope_GEA_temp)<-paste(colnames(D_SS_MAX_TAB),"_1",sep="")
      colnames(pval_temp)<-paste("p_",colnames(D_SS_MAX_TAB),sep="")
      slopePval<-cbind(slope_GEA_temp,pval_temp)
      GEA_D<-bind_rows(cbind(slopePval,D_SS_MAX_TAB,gea_pos=paste(GeaSubset$ENV[k],GeaSubset$BP[k],sep="_")),GEA_D)
    }
    
    
  }
  #GEA_INT_DATABASE<-bind_rows(GEA_INT_DATABASE,tempDF)
}
dev.off()
