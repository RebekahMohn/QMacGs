#Master GEA + INT windows


# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Biostrings")
# BiocManager::install("biomaRt")
# 
# install.packages("biomartr")
#install.packages("geiger")
#install.packages("phytools")
#install.packages("codetools")
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


# gea1<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllCe1_lfmmpval_cor.txt",sep="\t")
# gea2<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllCe2_lfmmpval_cor.txt",sep="\t")
# gea3<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllCe3_lfmmpval_cor.txt",sep="\t")
# gea4<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllCe4_lfmmpval_cor.txt",sep="\t")
# gea5<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllCe5_lfmmpval_cor.txt",sep="\t")
# gea6<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllCe6_lfmmpval_cor.txt",sep="\t")

gea_env<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_lfmmpval_cor.txt",sep=" ")
# DstatsMue<-read.delim()
# DstatsAlb<-read.delim()
# DstatsBic<-read.delim()
# DstatsLyr<-read.delim()
# DstatsSte<-read.delim()
# DstatsLob<-read.delim()
D_All_100<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/107_Dsuite_MACREF/DPOPS/100_100/Dall_mac.txt",header = TRUE)
#i=1
#D_All_100_2<-subset(D_All_100,chr==adaptIntro$csome[i]&windowStart<adaptIntro$stop[i]&windowEnd>adaptIntro$start[i])



taxa.list<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/111_GenomeRegionsPhy/taxa.list.txt")
convTable<-read.delim("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/oldMeta/metadata_figtree.txt")
taxa.vect<-taxa.list$samples
names(taxa.vect)<-taxa.list$taxon
macNsamp<-taxa.vect[13]
macNsamp<-unname(unlist(strsplit(macNsamp,",")))
macSsamp<-taxa.vect[14]
macSsamp<-unname(unlist(strsplit(macSsamp,",")))
macSamps<-c(macNsamp,macSsamp)


#leg<-cbind(percentile=c("<95%",">95%",">99%",">99.9%"),perc=c(94,96,99.5,99.99),long=c(-80,-85,-90,-95),lat=c(40,40,40,40))
leg<-cbind(percentile=c("<25%",">25%",">50%",">75%"),perc=c(94,96,99.5,99.99),long=c(-80,-85,-90,-95),lat=c(40,40,40,40))
leg2<-ggplot()+
  geom_point(data=leg,mapping=aes(x=long, y=lat,shape=percentile,size=percentile))+
  scale_shape_manual(values=c(17,16, 21, 21))+
  scale_size_manual(values=c(0.75,1,2.5,4))+
  theme(text = element_text(size = 8))+
  theme_light()
legmap<-as_ggplot(get_legend(leg2))

GEA_INT_DATABASE<-data.frame()
pdf("C:/Users/rmohn/Desktop/10_Analysis/120_Ad_Int/AdaptWins100_260205.pdf",height=6.5,width=9)

for(i in 1:length(adaptIntro$csome)){
  #Read the painting windows
  mue_adInt<-read.delim(paste("C:/Users/rmohn/Desktop/10_Analysis/113_painting/AdIntWins/mue_",adaptIntro$csome[i],"_",adaptIntro$pos[i],".tsv",sep=""),sep="\t")
  alb_adInt<-read.delim(paste("C:/Users/rmohn/Desktop/10_Analysis/113_painting/AdIntWins/alb_",adaptIntro$csome[i],"_",adaptIntro$pos[i],".tsv",sep=""),sep="\t")
  lob_adInt<-read.delim(paste("C:/Users/rmohn/Desktop/10_Analysis/113_painting/AdIntWins/lob_",adaptIntro$csome[i],"_",adaptIntro$pos[i],".tsv",sep=""),sep="\t")
  ste_adInt<-read.delim(paste("C:/Users/rmohn/Desktop/10_Analysis/113_painting/AdIntWins/ste_",adaptIntro$csome[i],"_",adaptIntro$pos[i],".tsv",sep=""),sep="\t")
  
  #Read the tree
  tree_top<-data.frame()
  tree_top2<-data.frame()
  treefile2<-list.files(path="C:/Users/rmohn/Desktop/10_Analysis/111_GenomeRegionsPhy/trees/",pattern=paste("*_",adaptIntro$csome[i],"_",adaptIntro$pos[i]+1,"-",adaptIntro$pos[i]+50001,".treefile",sep=""),full.names = TRUE)
  treefile1<-list.files(path="C:/Users/rmohn/Desktop/10_Analysis/111_GenomeRegionsPhy/trees/",pattern=paste("*_",adaptIntro$csome[i],"_",adaptIntro$pos[i]-50000+1,"-",adaptIntro$pos[i]+1,".treefile",sep=""),full.names = TRUE)
  
  if(length(treefile1)!=0){
    tree_top<-read.tree(treefile1)
    convTable2<-subset(convTable,Seq %in% tree_top$tip.label)
    
    tree_top_root<-midpoint_root(tree_top)
    tip_data<-data.frame(Seq=convTable2$Seq,Sp=factor(convTable2$Sp))
    tree_top_root_labs<-full_join(tree_top_root, tip_data, by = c("label"="Seq"))
    # tree_fig <- ggplot()+
    #   geom_blank()
    tree_fig <- ggplot(tree_top_root_labs,aes(x,y))+
      coord_cartesian(clip="off")+
      #xlim(-.05,0.05)+
      geom_tree()+
      #geom_tiplab(size=2)+
      theme_light()+
      geom_tippoint(aes(color=Sp),size=1)+
      scale_color_manual(values=group.colors)+
      theme(axis.title=element_blank(),axis.text = element_blank(),axis.line = element_blank(),
            axis.ticks = element_blank(),panel.grid = element_blank())
    
  }else{
    tree_fig <- ggplot()+
      geom_blank()
  }
  if(length(treefile2)!=0){
    tree_top2<-read.tree(treefile2)
    convTable2<-subset(convTable,Seq %in% tree_top2$tip.label)
    
    tree_top_root2<-data.frame()
    tree_top_root2<-midpoint_root(tree_top2)
    tree_top_root_labs2<-full_join(tree_top_root2, tip_data, by = c("label"="Seq"))
    tree_fig2 <- ggplot()+
      geom_blank()
    tree_fig2 <- ggplot(tree_top_root_labs2,aes(x,y))+
      coord_cartesian(clip="off")+
      #xlim(-.05,0.05)+
      geom_tree(size=.5)+
      #geom_tiplab(size=2)+
      theme_light()+
      geom_tippoint(aes(color=Sp),size=1)+
      scale_color_manual(values=group.colors)+
      theme(axis.title=element_blank(),axis.text = element_blank(),axis.line = element_blank(),
            axis.ticks = element_blank(),panel.grid = element_blank())
    
  }else{
    tree_fig2 <- ggplot()+
      geom_blank()
  }
  
  
  #Subset the GEA
  GeaSubset<-subset(gea_env,csome==adaptIntro$csome[i]&BP>adaptIntro$pos[i]-50000&BP<adaptIntro$pos[i]+50000)
  GeaSubset$gea_pos<-paste(GeaSubset$ENV,GeaSubset$BP,sep="_")
  
  #Subset the Dstats
  D_SS<-subset(D_All_100,chr==adaptIntro$csome[i]&windowStart<adaptIntro$pos[i]+50000&windowEnd>adaptIntro$pos[i]-50000)
  D_SS$mac_pops<-paste("mac_",D_SS$pop,sep="")
  D_SS_GEA_MAX<-data.frame()
  D_SS_MAX_TAB<-data.frame()
  
  
  
  #genes
  genes=subset(gff_data,seqid==adaptIntro$csome[i]&((start>adaptIntro$pos[i]-75000&start<adaptIntro$pos[i]+75000)|(end>adaptIntro$pos[i]-75000&end<adaptIntro$pos[i]+75000))&type=="mRNA")
  genes$ID<-str_split_i(str_split_i(genes$attribute,";",1),"=",2)
  regGenes<-merge(genes,gene_annotation,by.x="ID",by.y="Query_Sequence",all=FALSE)
  regGenes$mid<-(regGenes$start+regGenes$end)/2
  geneDescription<-paste(c(regGenes$SeqSearch_Description),collapse="; ")
  
  # genes=subset(gff_data,seqid==adaptIntro$csome[i]&((start>adaptIntro$pos[i]-75000&start<adaptIntro$pos[i]+75000)|(end>adaptIntro$pos[i]-75000&end<adaptIntro$pos[i]+75000))&type=="mRNA")
  # 
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
    # if(length(D_SS_MAX_TAB[1,]>0)){
    #   tempDF<-cbind(GeaSubset,D_SS_MAX_TAB)
    # }else{
    #   tempDF<-GeaSubset
    # }
    # 
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
        #print(ggplot(data=SS_GEA_D,aes(x=d_f,y=SS_GEA_D[,18]))+ geom_point(aes(color=distWind))+geom_smooth(method="lm")+
        #       labs(caption=paste(slope,summary(model1)$coefficients[2,4])))
      }
      colnames(slope_GEA_temp)<-paste(colnames(D_SS_MAX_TAB),"_1",sep="")
      colnames(pval_temp)<-paste("p_",colnames(D_SS_MAX_TAB),sep="")
      slopePval<-cbind(slope_GEA_temp,pval_temp)
      GEA_D<-bind_rows(cbind(slopePval,D_SS_MAX_TAB,gea_pos=paste(GeaSubset$ENV[k],GeaSubset$BP[k],sep="_")),GEA_D)
    }
    
    
    genesCloseTemp<-ddply(regGenes, .(seqid), function(z) {
      if(length(subset(z,(start<GeaSubset$BP[k]&end>GeaSubset$BP[k])|(end<GeaSubset$BP[k]&start>GeaSubset$BP[k]))$end)!=0){
        subset(z,(start<GeaSubset$BP[k]&end>GeaSubset$BP[k])|(end<GeaSubset$BP[k]&start>GeaSubset$BP[k]))
      }else{
        z$Start<-z$start
        z$Start[z$strand=="-"]<-z$end
       # print(abs(z$Start - GeaSubset$BP[k]) == min(abs(z$Start - GeaSubset$BP[k])))
        z[abs(z$Start - GeaSubset$BP[k]) == min(abs(z$Start - GeaSubset$BP[k])), ]
      }
    })
    
    
    # genesCloseTemp<-ddply(regGenes, .(seqid), function(z) {
    #   if(length(subset(z,(start<35181642&end>35181642)|(end<35181642&start>35181642))$end)==0){
    #     subset(z,(start<35181642&end>35181642)|(end<35181642&start>35181642))
    #   }else{
    #     StartBP<-z$start
    #     StartBP[z$strand=="-"]<-z$end
    #     z[abs(StartBP - 35181642) == min(abs(StartBP - 35181642)), ]
    #   }
    # })
    # 
    # genesCloseTemp<-ddply(regGenes, .(seqid), function(z) {
    #   if(length(subset(z,(start<GeaSubset$BP[k]&end>GeaSubset$BP[k])|(end<GeaSubset$BP[k]&start>GeaSubset$BP[k]))$end)!=0){
    #     subset(z,(start<GeaSubset$BP[k]&end>GeaSubset$BP[k])|(end<GeaSubset$BP[k]&start>GeaSubset$BP[k]))
    #   }else{
    #     z$Start<-z$start
    #     z$Start[z$strand=="-"]<-z$end
    #     print(abs(z$Start - GeaSubset$BP[k]) == min(abs(z$Start - GeaSubset$BP[k])))
    #     z[abs(z$Start - GeaSubset$BP[k]) == min(abs(z$Start - GeaSubset$BP[k])), ]
    #   }
    # })
    
    
    if(!is.na(genesCloseTemp[1,1])){
      genesCloseTemp$gea_pos<-GeaSubset$gea_pos[k]
      #genesCloseTemp$BP<-GeaSubset$BP[k]
      genesGEA<-bind_rows(genesGEA,genesCloseTemp)
      
      
    }
    
    
  }
  
  
  #genesGEA$gea_pos<-paste(genesGEA$ENV,genesGEA$BP,sep="_")
  
  #subset inversions
  invs_ss<-subset(invs,chromosome==adaptIntro$csome[i]&((start>adaptIntro$pos[i]-50000&start<adaptIntro$pos[i]+50000)|(end>adaptIntro$pos[i]-50000&end<adaptIntro$pos[i]+50000)|(start<adaptIntro$pos[i]-50000&end>adaptIntro$pos[i]+50000)))
  #subset centromeres
  cents_ss<-subset(cents,scaffold==adaptIntro$csome[i]&((Start>adaptIntro$pos[i]-50000&Start<adaptIntro$pos[i]+50000)|(Stop>adaptIntro$pos[i]-50000&Stop<adaptIntro$pos[i]+50000)|(Start<adaptIntro$pos[i]-50000&Stop>adaptIntro$pos[i]+50000)))
  
  # if(length(D_SS_MAX_TAB[1,]>0)){
  #   tempDF<-cbind(GeaSubset,D_SS_MAX_TAB)
  # }else{
  #   tempDF<-GeaSubset
  # }
  
  
  
  
  tempDF<-GeaSubset
  tempDF$invs<-invs_ss$INV[1]
  tempDF$cents<-cents_ss$Type[1]
  if(!is.na(genesCloseTemp[1,1])){
    tempDF<-merge(tempDF,genesGEA,by="gea_pos")
  }

  #tempDF<-merge(tempDF,genesGEA,by="gea_pos")
  if(length(GEA_D$p_alb)>0){
    GEA_D$lob_1[GEA_D$p_lob>0.05]<-0
    GEA_D$ste_1[GEA_D$p_ste>0.05]<-0
    GEA_D$mue_1[GEA_D$p_mue>0.05]<-0
    GEA_D$alb_1[GEA_D$p_alb>0.05]<-0
    GEA_D$sigDGEA<-"notSig"
    GEA_D$sigDGEA[GEA_D$p_lob<0.000005|GEA_D$p_ste<0.000005|GEA_D$p_mue<0.000005|GEA_D$p_alb<0.000005]<-"sig" #Is this correct for p_alba
    tempDF<-merge(tempDF,GEA_D,by="gea_pos",all = TRUE)
    #GEA_D$gea_pos<-paste(GeaSubset$ENV,GeaSubset$BP,sep="_")
    
  }else{
  tempDF$sigDGEA<-"-"
  }
  GEA_INT_DATABASE<-bind_rows(GEA_INT_DATABASE,tempDF)
  
  
  
  
  #visuals
  ggpaint<-ggplot()+
    labs(subtitle=paste(adaptIntro$csome[i]," ",adaptIntro$pos[i],sep=""),y="haplotypes sorted by latitude")+
    geom_tile(data=subset(lob_adInt,value==1), mapping=aes(x=as.numeric(pos),y=paste(latitude.orig,variable),fill="lob"))+
    geom_tile(data=subset(ste_adInt,value==1), mapping=aes(x=as.numeric(pos),y=paste(latitude.orig,variable),fill="ste"))+
    geom_tile(data=subset(mue_adInt,value==1), mapping=aes(x=as.numeric(pos),y=paste(latitude.orig,variable),fill="mue"))+
    geom_tile(data=subset(alb_adInt,value==1), mapping=aes(x=as.numeric(pos),y=paste(latitude.orig,variable),fill="alb"))+
    theme_light()+
    scale_fill_manual(values=group.colors)+
    theme(text = element_text(size=8),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),panel.spacing.x=unit(.2, "lines"),panel.spacing.y=unit(0, "lines"),
          axis.text.x = element_blank(),axis.title.x = element_blank())
  
  ggDstats<-ggplot()+
    geom_hline(yintercept=0.5,color="black",linewidth = 2)+
    geom_point(data=D_SS,mapping = aes(x=windowMid,y=d_f,group=pop,color=species))+
    geom_rect(data=D_SS_GEA_MAX,mapping = aes(xmin=windowStart,xmax=windowEnd,y=d_f,height=0.01, group=pop,fill=species))+
    theme_light()+
    scale_color_manual(values=group.colors)+
    scale_fill_manual(values=group.colors)+
    theme(text = element_text(size=8),panel.spacing.x=unit(.2, "lines"),panel.spacing.y=unit(0, "lines"),
          axis.text.x = element_blank(),axis.title.x = element_blank(),legend.position = "none")+
    coord_cartesian(ylim=c(-1,1),xlim = c(adaptIntro$pos[i]-50000,adaptIntro$pos[i]+50000))
  
  ggFeatures<-ggplot()+
    labs(x=NULL)+
    geom_rect(data=genes, aes(xmin=start,xmax=end,ymin=10,ymax=15,fill=strand),color="black")+
    geom_point(data=tempDF,mapping = aes(x=as.numeric(BP),y=-log10(pbonferroni),color=ENV,shape=sigDGEA),size=3)+
    geom_rect(data=invs_ss, aes(xmin=start,xmax=end,ymin=16,ymax=20,fill=INV))+
    geom_rect(data=cents_ss, aes(xmin=Start,xmax=Stop,ymin=21,ymax=25,fill=Type))+
    theme_light()+
    #scale_color_viridis_c()+
    scale_fill_manual(values=group.colors)+
    scale_color_manual(values=group.colors)+
    scale_shape_manual(values = c(notSig=4,sig=17,`-`=3))+
    theme(text = element_text(size=8),panel.grid.minor.y = element_blank(),
          panel.spacing.x=unit(.2, "lines"),panel.spacing.y=unit(0, "lines"),plot.caption = element_text(size=6))+
    #labs(caption=str_wrap(geneDescription, width=82))+
    coord_cartesian(xlim=c(adaptIntro$pos[i]-50000,adaptIntro$pos[i]+50000),ylim=c(1,25))
  
  LeftPanel<-ggpubr::ggarrange(ggpaint,ggDstats,ggFeatures,ncol=1,heights = c(4,2,2), align="v")
  
  sp_fig<-ggplot()+
    geom_polygon(data = counties, aes(x = long, y = lat, group = group),fill = 'white', color = "white", lwd = 0.1)+
    geom_polygon(data = states, aes(x = long, y = lat, group = group),
                 fill = NA, color = "gray85")+
    geom_polygon(data = Canada, aes(x = long, y = lat, group = group),
                 fill = "white", color = "gray85")+
    coord_map(xlim = c(-105,-68),ylim = c(29,50))+
    scale_size(limits=c(.1,1),range=c(2,8))+
    geom_point(data=subset(D_SS_GEA_MAX,d_f<.25),mapping=aes(x=long, y=lat), color = "black",size = .5,shape=17)+
    geom_point(data=subset(D_SS_GEA_MAX,d_f>=.25),mapping=aes(x=long, y=lat, color = species),size = 1, shape= 16)+
    geom_point(data=subset(D_SS_GEA_MAX,d_f>.50),mapping=aes(x=long, y=lat, color = species),size = 2.5, shape= 21)+
    geom_point(data=subset(D_SS_GEA_MAX,d_f>.75),mapping=aes(x=long, y=lat, color = species),size = 4, shape= 21)+
    #geom_point(data=subset(DwindMax,d_f>.75),mapping=aes(x=long, y=lat, color = species),size = 5.5, shape= 21)+
    scale_color_manual(values=group.colors)+
    theme(axis.title = element_blank(),text = element_text(size = 8))+
    guides(color="none",size="none",shape="none")
  
  map_fig<-ggpubr::ggarrange(sp_fig,legmap,ncol=2,nrow = 1,widths = c(4,1))
  trees<-ggpubr::ggarrange(tree_fig,tree_fig2,ncol=2,common.legend = T,legend = "right")
  
  #  add.scale.bar(cex=.8)
  #  nodelabels(text=rub_bt_l$node.label,cex=.75,frame = "none",adj=0)
  
  rightSide<-ggpubr::ggarrange(trees,map_fig,ncol=1)
  
  print(annotate_figure(ggpubr::ggarrange(LeftPanel,rightSide,ncol=2),bottom=text_grob(str_wrap(geneDescription, width=160),size=8)))
}

dev.off()

write.table(GEA_INT_DATABASE,file="C:/Users/rmohn/Desktop/10_Analysis/120_Ad_Int/Ad_Int_Database_samp100_260205.txt",row.names = FALSE, quote = FALSE, sep="\t")



GEA_INT_DATABASE<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/120_Ad_Int/Ad_Int_Database_samp.txt")
GEA_INT_DATABASE$lob_1[GEA_INT_DATABASE$p_lob>0.05]<-0
GEA_INT_DATABASE$ste_1[GEA_INT_DATABASE$p_ste>0.05]<-0
GEA_INT_DATABASE$mue_1[GEA_INT_DATABASE$p_mue>0.05]<-0
GEA_INT_DATABASE$alb_1[GEA_INT_DATABASE$p_alb>0.05]<-0


albGEAINT<-cbind(ENV=GEA_INT_DATABASE$ENV,Dsp=GEA_INT_DATABASE$alb,Dsp1=GEA_INT_DATABASE$alb_1,sp="alb")
mueGEAINT<-cbind(ENV=GEA_INT_DATABASE$ENV,Dsp=GEA_INT_DATABASE$mue,Dsp1=GEA_INT_DATABASE$mue_1,sp="mue")
steGEAINT<-cbind(ENV=GEA_INT_DATABASE$ENV,Dsp=GEA_INT_DATABASE$ste,Dsp1=GEA_INT_DATABASE$ste_1,sp="ste")
lobGEAINT<-cbind(ENV=GEA_INT_DATABASE$ENV,Dsp=GEA_INT_DATABASE$lob,Dsp1=GEA_INT_DATABASE$lob_1,sp="lob")


GEAINTDIR<-as.data.frame(rbind(albGEAINT,mueGEAINT,steGEAINT,lobGEAINT))
GEAINTDIR$neg<-NA
GEAINTDIR$neg[GEAINTDIR$Dsp1<0]<-"-"
GEAINTDIR$neg[GEAINTDIR$Dsp1>0]<-"+"
table(GEAINTDIR$sp,GEAINTDIR$neg,GEAINTDIR$ENV)






# #Calculating table
# 
# 
# 
# 
# for(i in 1:length(adaptIntro$csome)){
#   #Read the painting windows
#   mue_adInt<-read.delim(paste("C:/Users/rmohn/Desktop/10_Analysis/113_painting/AdIntWins/mue_",adaptIntro$csome[i],"_",adaptIntro$pos[i],".tsv",sep=""),sep="\t")
#   alb_adInt<-read.delim(paste("C:/Users/rmohn/Desktop/10_Analysis/113_painting/AdIntWins//alb_",adaptIntro$csome[i],"_",adaptIntro$pos[i],".tsv",sep=""),sep="\t")
#   lob_adInt<-read.delim(paste("C:/Users/rmohn/Desktop/10_Analysis/113_painting/AdIntWins/lob_",adaptIntro$csome[i],"_",adaptIntro$pos[i],".tsv",sep=""),sep="\t")
#   ste_adInt<-read.delim(paste("C:/Users/rmohn/Desktop/10_Analysis/113_painting/AdIntWins/ste_",adaptIntro$csome[i],"_",adaptIntro$pos[i],".tsv",sep=""),sep="\t")
#   
#   #Read the tree
#   tree_top<-data.frame()
#   tree_top2<-data.frame()
#   treefile2<-list.files(path="C:/Users/rmohn/Desktop/10_Analysis/111_GenomeRegionsPhy/trees/",pattern=paste("*_",adaptIntro$csome[i],"_",adaptIntro$pos[i]+1,"-",adaptIntro$pos[i]+50001,".treefile",sep=""),full.names = TRUE)
#   treefile1<-list.files(path="C:/Users/rmohn/Desktop/10_Analysis/111_GenomeRegionsPhy/trees/",pattern=paste("*_",adaptIntro$csome[i],"_",adaptIntro$pos[i]-50000+1,"-",adaptIntro$pos[i]+1,".treefile",sep=""),full.names = TRUE)
#   
#   if(length(treefile1)!=0){
#     tree_top<-read.tree(treefile1)
#     convTable2<-subset(convTable,Seq %in% tree_top$tip.label)
#     
#     tree_top_root<-midpoint_root(tree_top)
#     tip_data<-data.frame(Seq=convTable2$Seq,Sp=factor(convTable2$Sp))
#     tree_top_root_labs<-full_join(tree_top_root, tip_data, by = c("label"="Seq"))
#     tree_fig <- ggplot()+
#       geom_blank()
#     tree_fig <- ggplot(tree_top_root_labs,aes(x,y))+
#       coord_cartesian(clip="off")+
#       #xlim(-.05,0.05)+
#       geom_tree(size=.5)+
#       #geom_tiplab(size=2)+
#       theme_light()+
#       geom_tippoint(aes(color=Sp),size=1)+
#       scale_color_manual(values=group.colors)+
#       theme(axis.title=element_blank(),axis.text = element_blank(),axis.line = element_blank(),
#             axis.ticks = element_blank(),panel.grid = element_blank())
#     
#   }else{
#     tree_fig <- ggplot()+
#       geom_blank()
#     }
#   if(length(treefile2)!=0){
#     tree_top2<-read.tree(treefile2)
#     convTable2<-subset(convTable,Seq %in% tree_top2$tip.label)
#     
#     tree_top_root2<-data.frame()
#     tree_top_root2<-midpoint_root(tree_top2)
#     tree_top_root_labs2<-full_join(tree_top_root2, tip_data, by = c("label"="Seq"))
#     tree_fig2 <- ggplot()+
#       geom_blank()
#     tree_fig2 <- ggplot(tree_top_root_labs2,aes(x,y))+
#       coord_cartesian(clip="off")+
#       #xlim(-.05,0.05)+
#       geom_tree(size=.5)+
#       #geom_tiplab(size=2)+
#       theme_light()+
#       geom_tippoint(aes(color=Sp),size=1)+
#       scale_color_manual(values=group.colors)+
#       theme(axis.title=element_blank(),axis.text = element_blank(),axis.line = element_blank(),
#             axis.ticks = element_blank(),panel.grid = element_blank())
#     
#   }else{
#     tree_fig2 <- ggplot()+
#       geom_blank()
#   }
#   # tree_top_labs<-rename_taxa(tree_top,convTable2,key="Seq",value="Sp")
#   # tree_sisters<-data.frame()
#   # for(t in 1:length(macSamps)){
#   #   dropMac<-macSamps[-t]
#   #   #drop all except t
#   #   rtree_mac<-drop.tip(tree_top,dropMac)
#   #   #record sister taxa
#   #   sisters<-tryCatch(tips(rtree_mac, getSisters(rtree_mac,macSamps[t],mode=c("number"))),error = function(e) NA)
#   #   for(j in sisters){
#   #     names(taxa.vect)
#   #   }
#   #   tree_sisters<-rbind(tree_sisters, c(macSamps[t],sisters))
#   #   #If we figure it out, convert names to species
#   # }
#   
#   
#   #Subset the Dstats
#   D_SS<-subset(D_All_100,chr==adaptIntro$csome[i]&windowMid>adaptIntro$pos[i]-50000&windowMid<adaptIntro$pos[i]+50000)
#   D_SS_MAX<-(D_SS %>% group_by(pop) %>% top_n(1, d_f))
#   
#   #D_SS_MAX_TAB<-paste(names(table(subset(D_SS_MAX,d_f>.1)$species)),table(subset(D_SS_MAX,d_f>.1)$species), sep =": ")
#   D_SS_MAX_TAB<-data.frame()
#   if(length(subset(D_SS_MAX,d_f>.1)$species)>0){
#   
#   D_SS_MAX_TAB<-as.data.frame(t(data.frame(table(subset(D_SS_MAX,d_f>.1)$species))[,2]))
#   colnames(D_SS_MAX_TAB)<-t(data.frame(table(subset(D_SS_MAX,d_f>.1)$species))[,1])
#   rownames(D_SS_MAX_TAB)<-NULL
#   }
#   
#   #Subset the GEA
#   GeaSubset<-subset(gea_env,csome==adaptIntro$csome[i]&BP>adaptIntro$pos[i]-50000&BP<adaptIntro$pos[i]+50000)
#   GeaSubset$gea_pos<-paste(GeaSubset$ENV,GeaSubset$BP,sep="_")
# 
#   #genes
#   genes=subset(gff_data,seqid==adaptIntro$csome[i]&((start>adaptIntro$pos[i]-75000&start<adaptIntro$pos[i]+75000)|(end>adaptIntro$pos[i]-75000&end<adaptIntro$pos[i]+75000))&type=="mRNA")
#   genes$ID<-str_split_i(str_split_i(genes$attribute,";",1),"=",2)
#   regGenes<-merge(genes,gene_annotation,by.x="ID",by.y="Query_Sequence",all=FALSE)
#   regGenes$mid<-(regGenes$start+regGenes$end)/2
#   geneDescription<-paste(c(regGenes$SeqSearch_Description),collapse="; ")
#   
#   # genes=subset(gff_data,seqid==adaptIntro$csome[i]&((start>adaptIntro$pos[i]-75000&start<adaptIntro$pos[i]+75000)|(end>adaptIntro$pos[i]-75000&end<adaptIntro$pos[i]+75000))&type=="mRNA")
#   # 
#   genesGea<-data.frame()
#   for(k in 1:length(GeaSubset$BP)){
#     genesCloseTemp<-ddply(regGenes, .(seqid), function(z) {
#       z[abs(z$mid - GeaSubset$BP[k]) == min(abs(z$mid - GeaSubset$BP[k])), ]
#     })
#     if(!is.na(genesCloseTemp[1,1])){
#     genesCloseTemp$ENV<-GeaSubset$ENV[k]
#     genesCloseTemp$BP<-GeaSubset$BP[k]
#     genesGEA<- rbind(genesGea,genesCloseTemp)
#     genesGEA$gea_pos<-paste(genesGEA$ENV,genesGEA$BP,sep="_")
#     }
#   }
# 
#   
#   #subset inversions
#   invs_ss<-subset(invs,chromosome==adaptIntro$csome[i]&((start>adaptIntro$pos[i]-50000&start<adaptIntro$pos[i]+50000)|(end>adaptIntro$pos[i]-50000&end<adaptIntro$pos[i]+50000)|(start<adaptIntro$pos[i]-50000&end>adaptIntro$pos[i]+50000)))
#   #subset centromeres
#   cents_ss<-subset(cents,scaffold==adaptIntro$csome[i]&((Start>adaptIntro$pos[i]-50000&Start<adaptIntro$pos[i]+50000)|(Stop>adaptIntro$pos[i]-50000&Stop<adaptIntro$pos[i]+50000)|(Start<adaptIntro$pos[i]-50000&Stop>adaptIntro$pos[i]+50000)))
#   
#   if(length(D_SS_MAX_TAB[1,]>0)){
#     tempDF<-cbind(GeaSubset,D_SS_MAX_TAB)
#   }else{
#       tempDF<-GeaSubset
#     }
#     
#   tempDF$invs<-invs_ss$INV[1]
#   tempDF$cents<-cents_ss$Type[1]
#   tempDF<-merge(tempDF,genesGEA,by="gea_pos")
#   
#   GEA_INT_DATABASE<-bind_rows(GEA_INT_DATABASE,tempDF)
#   
#   
#   #visuals
#   # ggpaint<-ggplot()+
#   #   labs(subtitle=paste(adaptIntro$csome[i]," ",adaptIntro$pos[i],sep=""),y="haplotypes sorted by latitude")+
#   #   geom_tile(data=subset(lob_adInt,value==1), mapping=aes(x=as.numeric(pos),y=paste(latitude.orig,variable),fill="lob"))+
#   #   geom_tile(data=subset(ste_adInt,value==1), mapping=aes(x=as.numeric(pos),y=paste(latitude.orig,variable),fill="ste"))+
#   #   geom_tile(data=subset(mue_adInt,value==1), mapping=aes(x=as.numeric(pos),y=paste(latitude.orig,variable),fill="mue"))+
#   #   geom_tile(data=subset(alb_adInt,value==1), mapping=aes(x=as.numeric(pos),y=paste(latitude.orig,variable),fill="alb"))+
#   #   theme_light()+
#   #   scale_fill_manual(values=group.colors)+
#   #   theme(text = element_text(size=8),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),axis.text.y = element_blank(),
#   #         axis.ticks.y = element_blank(),panel.spacing.x=unit(.2, "lines"),panel.spacing.y=unit(0, "lines"),
#   #         axis.text.x = element_blank(),axis.title.x = element_blank())
#   # ggDstats<-ggplot()+
#   #   geom_point(data=D_SS,mapping = aes(x=windowMid,y=d_f,group=pop,color=species),size=1)+
#   #   theme_light()+
#   #   geom_hline(yintercept=0.1,color="red")+
#   #   scale_color_manual(values=group.colors)+
#   #   theme(text = element_text(size=8),panel.spacing.x=unit(.2, "lines"),panel.spacing.y=unit(0, "lines"),
#   #         axis.text.x = element_blank(),axis.title.x = element_blank())
#   # ggFeatures<-ggplot()+
#   #   labs(x=NULL)+
#   #   geom_rect(data=genes, aes(xmin=start,xmax=end,ymin=10,ymax=15,fill=strand))+
#   #   geom_point(data=GeaSubset,mapping = aes(x=as.numeric(BP),y=log10(pbonferroni),fill=ENV),size=3,shape=24)+
#   #   geom_rect(data=invs_ss, aes(xmin=start,xmax=end,ymin=16,ymax=20,fill=INV))+
#   #   geom_rect(data=cents_ss, aes(xmin=Start,xmax=Stop,ymin=21,ymax=25,fill=Type))+
#   #   theme_light()+
#   #   #scale_color_viridis_c()+
#   #   scale_fill_manual(values=group.colors)+
#   #   theme(text = element_text(size=8),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),axis.text.y = element_blank(),
#   #         axis.ticks.y = element_blank(),panel.spacing.x=unit(.2, "lines"),panel.spacing.y=unit(0, "lines"))+
#   #   labs(caption=str_wrap(geneDescription, width = 80))+
#   #   coord_cartesian(xlim=c(adaptIntro$pos[i]-50000,adaptIntro$pos[i]+50000))
#   # 
#   # LeftPanel<-ggarrange(ggpaint,ggDstats,ggFeatures,ncol=1,heights = c(4,2,2), align="v")
#   # 
#   # sp_fig<-ggplot()+
#   #   geom_polygon(data = counties, aes(x = long, y = lat, group = group),fill = 'white', color = "white", lwd = 0.1)+
#   #   geom_polygon(data = states, aes(x = long, y = lat, group = group),
#   #                fill = NA, color = "gray85")+
#   #   geom_polygon(data = Canada, aes(x = long, y = lat, group = group),
#   #                fill = "white", color = "gray85")+
#   #   coord_map(xlim = c(-105,-68),ylim = c(29,50))+
#   #   scale_size(limits=c(.1,1),range=c(2,8))+
#   #   geom_point(data=subset(D_SS_MAX,d_f<.10),mapping=aes(x=long, y=lat), color = "black",size = .5,shape=17)+
#   #   geom_point(data=subset(D_SS_MAX,d_f>=.10),mapping=aes(x=long, y=lat, color = species),size = 1, shape= 16)+
#   #   geom_point(data=subset(D_SS_MAX,d_f>.50),mapping=aes(x=long, y=lat, color = species),size = 2.5, shape= 21)+
#   #   geom_point(data=subset(D_SS_MAX,d_f>.75),mapping=aes(x=long, y=lat, color = species),size = 4, shape= 21)+
#   #   #geom_point(data=subset(DwindMax,d_f>.75),mapping=aes(x=long, y=lat, color = species),size = 5.5, shape= 21)+
#   #   scale_color_manual(values=group.colors)+
#   #   theme(axis.title = element_blank(),text = element_text(size = 8))+
#   #   guides(color="none",size="none",shape="none")
#   # 
#   # map_fig<-ggarrange(sp_fig,legmap,ncol=2,nrow = 1,widths = c(4,1))
#   # 
#   # 
#   # #  add.scale.bar(cex=.8)
#   # #  nodelabels(text=rub_bt_l$node.label,cex=.75,frame = "none",adj=0)
#   # tree_top_root<-midpoint_root(tree_top)
#   # 
#   # tip_data<-data.frame(Seq=convTable2$Seq,Sp=factor(convTable2$Sp))
#   # 
#   # tree_top_root_labs<-full_join(tree_top_root, tip_data, by = c("label"="Seq"))
#   # 
#   # tree_fig <- ggplot(tree_top_root_labs,aes(x,y))+
#   #   coord_cartesian(clip="off")+
#   #   #xlim(-.05,0.05)+
#   #   geom_tree(size=.5)+
#   #   #geom_tiplab(size=2)+
#   #   theme_light()+
#   #   geom_tippoint(aes(color=Sp),size=1)+
#   #   scale_color_manual(values=group.colors)+
#   #   theme(axis.title=element_blank(),axis.text = element_blank(),axis.line = element_blank(),
#   #         axis.ticks = element_blank(),panel.grid = element_blank())
#   # 
#   # tree_top_root2<-midpoint_root(tree_top2)
#   # 
#   # 
#   # tree_top_root_labs2<-full_join(tree_top_root2, tip_data, by = c("label"="Seq"))
#   # 
#   # tree_fig2 <- ggplot(tree_top_root_labs2,aes(x,y))+
#   #   coord_cartesian(clip="off")+
#   #   #xlim(-.05,0.05)+
#   #   geom_tree(size=.5)+
#   #   #geom_tiplab(size=2)+
#   #   theme_light()+
#   #   geom_tippoint(aes(color=Sp),size=1)+
#   #   scale_color_manual(values=group.colors)+
#   #   theme(axis.title=element_blank(),axis.text = element_blank(),axis.line = element_blank(),
#   #         axis.ticks = element_blank(),panel.grid = element_blank())
#   # 
#   # trees<-ggarrange(tree_fig,tree_fig2,ncol=2,common.legend = T,legend = "right")
#   # 
#   # rightSide<-ggarrange(trees,map_fig,ncol=1)
#   # 
#   # print(ggarrange(LeftPanel,rightSide,ncol=2))
#   
# }

#dev.off()

