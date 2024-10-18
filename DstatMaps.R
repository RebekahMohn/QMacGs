#DstatMaps
##This file is to visualize the results of Dsuite's Dinvestigate run on each population of Q. macrocarpa.
###Alb, lob, mue, ste: 
####P1=bic, P2=mac, P3=test
###Lyr
####P1=bic, P2=lyr, P3= mac
###Bic
####P1=lyr, P2=bic, P3=mac
####Note: Bic=-Lyr
##Alba
#combining all dstatistics and metadata for Q. alba
library(dplyr)
library(plyr)
meta_pop<-read.csv("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/metadata_pops.csv")
Dalb_All<-data.frame(pop=character(),chr=character(),windowStart=integer(),windowEnd=integer(),D=double(),f_d=double(),f_dM=double(),d_f=double(),lat=double(),long=double(),windowmid=integer())

for(i in 1:length(meta_pop$Pop)){
  pop<-meta_pop$Pop[i]
  lat<-meta_pop$latitude.orig[i]
  long<-meta_pop$longitude.orig[i]
  if(file.exists(paste("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_",pop,"_alba_localFstats__2000_2000.txt",sep=""))){
    Dalb_pop<-read.delim(paste("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_",pop,"_alba_localFstats__2000_2000.txt",sep=""))
    Dalb_pop$pop<-rep(pop,length(Dalb_pop$chr))
    Dalb_pop$lat<-rep(lat,length(Dalb_pop$chr))
    Dalb_pop$long<-rep(long,length(Dalb_pop$chr))
    Dalb_pop$windowmid<-round_any(rowMeans(Dalb_pop[,2:3]),500000)
    Dalb_All<-rbind(Dalb_All,Dalb_pop)
  }
  #Dalb_pop$Loc<-paste(Dalb_pop$CHROM,Dalb_pop$BIN_START,sep="_")
  
}
write.table(Dalb_All,file="C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_alba.txt",row.names=FALSE,quote=FALSE,sep="\t")
#I don't know that window number is the best method

##MUE
#combining all dstatistics and metadata for Q. muehlenbergii
library(plyr)
meta_pop<-read.csv("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/metadata_pops.csv")
macABF_mue<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_ABF_muehlenbergii_localFstats__2000_2000.txt")
Dmue_All<-data.frame(pop=character(),chr=character(),windowStart=integer(),windowEnd=integer(),D=double(),f_d=double(),f_dM=double(),d_f=double(),lat=double(),long=double(),windowmid=integer())

for(i in 1:length(meta_pop$Pop)){
  pop<-meta_pop$Pop[i]
  lat<-meta_pop$latitude.orig[i]
  long<-meta_pop$longitude.orig[i]
  if(file.exists(paste("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_",pop,"_muehlenbergii_localFstats__2000_2000.txt",sep=""))){
    Dmue_pop<-read.delim(paste("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_",pop,"_muehlenbergii_localFstats__2000_2000.txt",sep=""))
    Dmue_pop$pop<-rep(pop,length(Dmue_pop$chr))
    Dmue_pop$lat<-rep(lat,length(Dmue_pop$chr))
    Dmue_pop$long<-rep(long,length(Dmue_pop$chr))
    Dmue_pop$windowmid<-round_any(rowMeans(Dmue_pop[,2:3]),500000)
    Dmue_All<-rbind(Dmue_All,Dmue_pop)
  }
  #Dmue_pop$Loc<-paste(Dmue_pop$CHROM,Dmue_pop$BIN_START,sep="_")
  
}
write.table(Dmue_All,file="C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_muehlenbergii.txt",row.names=FALSE,quote=FALSE,sep="\t")


##STE
#combining all dstatistics and metadata for Q. stellata
library(plyr)
meta_pop<-read.csv("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/metadata_pops.csv")
macABF_ste<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_ABF_stellata_localFstats__2000_2000.txt")
Dste_All<-data.frame(pop=character(),chr=character(),windowStart=integer(),windowEnd=integer(),D=double(),f_d=double(),f_dM=double(),d_f=double(),lat=double(),long=double(),windowmid=integer())

for(i in 1:length(meta_pop$Pop)){
  pop<-meta_pop$Pop[i]
  lat<-meta_pop$latitude.orig[i]
  long<-meta_pop$longitude.orig[i]
  if(file.exists(paste("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_",pop,"_stellata_localFstats__2000_2000.txt",sep=""))){
    Dste_pop<-read.delim(paste("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_",pop,"_stellata_localFstats__2000_2000.txt",sep=""))
    Dste_pop$pop<-rep(pop,length(Dste_pop$chr))
    Dste_pop$lat<-rep(lat,length(Dste_pop$chr))
    Dste_pop$long<-rep(long,length(Dste_pop$chr))
    Dste_pop$windowmid<-round_any(rowMeans(Dste_pop[,2:3]),500000)
    Dste_All<-rbind(Dste_All,Dste_pop)
  }
  #Dste_pop$Loc<-paste(Dste_pop$CHROM,Dste_pop$BIN_START,sep="_")
  
}
write.table(Dste_All,file="C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_stellata.txt",row.names=FALSE,quote=FALSE,sep="\t")


##LOB
#combining all dstatistics and metadata for Q. lobata
library(plyr)
meta_pop<-read.csv("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/metadata_pops.csv")
macABF_lob<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_ABF_lobata_localFstats__2000_2000.txt")
Dlob_All<-data.frame(pop=character(),chr=character(),windowStart=integer(),windowEnd=integer(),D=double(),f_d=double(),f_dM=double(),d_f=double(),lat=double(),long=double(),windowmid=integer())

for(i in 1:length(meta_pop$Pop)){
  pop<-meta_pop$Pop[i]
  lat<-meta_pop$latitude.orig[i]
  long<-meta_pop$longitude.orig[i]
  if(file.exists(paste("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_",pop,"_lobata_localFstats__2000_2000.txt",sep=""))){
    Dlob_pop<-read.delim(paste("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_",pop,"_lobata_localFstats__2000_2000.txt",sep=""))
    Dlob_pop$pop<-rep(pop,length(Dlob_pop$chr))
    Dlob_pop$lat<-rep(lat,length(Dlob_pop$chr))
    Dlob_pop$long<-rep(long,length(Dlob_pop$chr))
    Dlob_pop$windowmid<-round_any(rowMeans(Dlob_pop[,2:3]),500000)
    Dlob_All<-rbind(Dlob_All,Dlob_pop)
  }
  #Dlob_pop$Loc<-paste(Dlob_pop$CHROM,Dlob_pop$BIN_START,sep="_")
  
}
write.table(Dlob_All,file="C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_lobata.txt",row.names=FALSE,quote=FALSE,sep="\t")

##bicolor/lyrata
###combining all dstatistics and metadata for Q. lyrata and Q. bicolor
library(plyr)
meta_pop<-read.csv("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/metadata_pops.csv")
macABF_biclyr<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_lyrata_macrocarpa_ABF_localFstats__2000_2000.txt")
Dbiclyr_All<-data.frame(pop=character(),chr=character(),windowStart=integer(),windowEnd=integer(),D=double(),f_d=double(),f_dM=double(),d_f=double(),lat=double(),long=double(),windowmid=integer())

for(i in 1:length(meta_pop$Pop)){
  pop<-meta_pop$Pop[i]
  lat<-meta_pop$latitude.orig[i]
  long<-meta_pop$longitude.orig[i]
  if(file.exists(paste("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_lyrata_macrocarpa_",pop,"_localFstats__2000_2000.txt",sep=""))){
    Dbiclyr_pop<-read.delim(paste("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_lyrata_macrocarpa_",pop,"_localFstats__2000_2000.txt",sep=""))
    Dbiclyr_pop$pop<-rep(pop,length(Dbiclyr_pop$chr))
    Dbiclyr_pop$lat<-rep(lat,length(Dbiclyr_pop$chr))
    Dbiclyr_pop$long<-rep(long,length(Dbiclyr_pop$chr))
    Dbiclyr_pop$windowmid<-round_any(rowMeans(Dbiclyr_pop[,2:3]),500000)
    Dbiclyr_All<-rbind(Dbiclyr_All,Dbiclyr_pop)
  }
}
write.table(Dbiclyr_All,file="C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_lyrata_macrocarpa.txt",row.names=FALSE,quote=FALSE,sep="\t")


##All
###combining all species dstatistics+metadata into one spreadsheet for all species
Dly_All<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_lyrata_macrocarpa.txt",sep="\t")
Dly_All$species<-rep("lyr",length(Dly_All$chr))
Dlob_All<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_lobata.txt",sep="\t")
Dlob_All$species<-rep("lob",length(Dlob_All$chr))
Dalb_All<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_alba.txt",sep="\t")
Dalb_All$species<-rep("alb",length(Dalb_All$chr))
Dste_All<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_stellata.txt",sep="\t")
Dste_All$species<-rep("ste",length(Dste_All$chr))
Dmue_All<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_muehlenbergii.txt",sep="\t")
Dmue_All$species<-rep("mue",length(Dmue_All$chr))
Dbi_All<-Dly_All
Dbi_All$d_f<-(-Dbi_All$d_f)
Dbi_All$species<-rep("bic",length(Dbi_All$chr))


D_All<-rbind(Dalb_All,Dbi_All,Dlob_All,Dly_All,Dmue_All,Dste_All)
D_All$windowMid<-rowMeans(D_All[,2:3])
write.table(D_All,file="C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/Dall_macrocarpa.txt",row.names=FALSE,quote=FALSE,sep="\t")

###for 500,000 bp windows, saiving the max d_f for each population (this helps with visualizing)

Dsp_All<-data.frame(pop=character(),chr=character(),windowStart=integer(),windowEnd=integer(),D=double(),f_d=double(),f_dM=double(),d_f=double(),lat=double(),long=double(),windowmid=integer(),species=character())
for(j in 1:12){
  i=min(D_All$windowmid)
  while(i <= max(D_All$windowmid)){
    if(i %in% subset(Dbi_All,chr==paste("Superscaffold",j,sep=""))$windowmid){
      Dwind<-rbind(subset(Dalb_All,windowmid==i&chr==paste("Superscaffold",j,sep="")),subset(Dbi_All,windowmid==i&chr==paste("Superscaffold",j,sep="")),subset(Dlob_All,windowmid==i&chr==paste("Superscaffold",j,sep="")),subset(Dly_All,windowmid==i&chr==paste("Superscaffold",j,sep="")),subset(Dmue_All,windowmid==i&chr==paste("Superscaffold",j,sep="")),subset(Dste_All,windowmid==i&chr==paste("Superscaffold",j,sep="")))
      Dsp_All<-rbind(Dsp_All,(Dwind %>% group_by(pop) %>% top_n(1, d_f)))
    }
    i=i+500000
  }
}

write.table(Dsp_All,file="C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/DspMax5E5_macrocarpa.txt",row.names=FALSE,quote=FALSE,sep="\t")

###Calculating DMSE for 150,000 bp windows
```{r}
library(dplyr)
library(data.table)
library(readxl)

D_All<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/Dall_macrocarpa.txt",header = TRUE)

DWindMSE<-data.frame(chr=integer(),window=integer(),alb=double(),bic=double(),lob=double(),lyr=double(),mue=double(),ste=double())

for(j in 1:12){
  i=75000
  while(i <= max(subset(D_All,chr==paste("Superscaffold",j,sep=""))$windowMid)){
    DallWind<-subset(D_All,windowMid>i-75000&windowMid<i+75000&chr==paste("Superscaffold",j,sep="")) %>% group_by(pop,species) %>% top_n(1, d_f)
    if(length(subset(DallWind,species=="alb")$d_f)>8|length(subset(DallWind,species=="bic")$d_f)>8|length(subset(DallWind,species=="lob")$d_f)>8|length(subset(DallWind,species=="mue")$d_f)>8|length(subset(DallWind,species=="mue")$d_f)>8){
      if(length(subset(DallWind,species=="alb")$d_f)>8){
        MSEalb<- mean(summary(lm(data=subset(DallWind,species=="alb"),d_f~windowMid))$residuals^2)
      }else{MSEalb<-NA}
      if(length(subset(DallWind,species=="bic")$d_f)>8){
        MSEbic<- mean(summary(lm(data=subset(DallWind,species=="bic"),d_f~windowMid))$residuals^2)
        MSElyr<- mean(summary(lm(data=subset(DallWind,species=="lyr"),d_f~windowMid))$residuals^2)
      }else{
        MSEbic<-NA
        MSElyr<-NA
      }
      if(length(subset(DallWind,species=="lob")$d_f)>8){
        MSElob<- mean(summary(lm(data=subset(DallWind,species=="lob"),d_f~windowMid))$residuals^2)
      }else{MSElob<-NA}
      if(length(subset(DallWind,species=="mue")$d_f)>8){
        MSEmue<- mean(summary(lm(data=subset(DallWind,species=="mue"),d_f~windowMid))$residuals^2)
      }else{MSEmue<-NA}
      if(length(subset(DallWind,species=="ste")$d_f)>8){
        MSEste<- mean(summary(lm(data=subset(DallWind,species=="ste"),d_f~windowMid))$residuals^2)
      }else{MSEste<-NA}
      MSE<-data.frame(chr=j,window=i,alb=MSEalb,bic=MSEbic,lob=MSElob,lyr=MSElyr,mue=MSEmue,ste=MSEste)
      DWindMSE<-rbind(DWindMSE,MSE)
    }
    i=i+75000
  }
}


write.table(DWindMSE,file="C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/Dmse_all_macrocarpa.txt",row.names=FALSE,quote=FALSE,sep="\t")


###Calculate windows of high DMSE
DMSE_All<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/Dmse_all_macrocarpa.txt",header = TRUE)
DMSE99<-quantile(c(DMSE_All$alb,DMSE_All$bic,DMSE_All$lob,DMSE_All$lyr,DMSE_All$mue,DMSE_All$ste),probs=c(0.99),na.rm=TRUE)

DMSE_highregions<-data.frame(chr=integer(),window=integer(),alb=double(),bic=double(),lob=double(),lyr=double(),mue=double(),ste=double(),windowstart=integer(),windowend=integer(),species=character())
for(j in 3:8){
  i<-1
  while(i < length(DMSE_All$window)){
    if(max(DMSE_All[i,j])>=DMSE99&!is.na(DMSE_All[i,j])){
      DMSE_highregions<-rbind(DMSE_highregions, cbind(DMSE_All[i,],windowstart=DMSE_All[i,2]-75000,windowend=DMSE_All[i,2]+75000,species=names(DMSE_All)[j]))
      #i=i+1
      k=0
      
      while(max(DMSE_All[i,j])>=DMSE99|is.na(DMSE_All[i,j])|k<3){
        i=i+1
        if(max(DMSE_All[i,j])>=DMSE99&!is.na(DMSE_All[i,j])){
          DMSE_highregions$windowend[length(DMSE_highregions$windowend)]<-DMSE_All[i,2]+75000
        }else{k=k+1}
      }
    }else{
      i=i+1
    }
  }
  
}
write.table(DMSE_highregions,"C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/DMSE_High_macrocarpa.txt",row.names=FALSE,quote=FALSE,sep="\t") 

###Visualizing the data
#install.packages("ggpubr")
library(ggpubr)
library(dplyr)
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(adegenet)
library(data.table)
library(readxl)
library(stringr)
states <- map_data("state")
counties <- map_data('county')
Canada<-map_data("world", "Canada")

D_All<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/Dall_macrocarpa.txt",header = TRUE)
D_All$chr2<-D_All$chr
D_All<-D_All %>%
  dplyr::mutate(across('chr2',str_replace,'Superscaffold',''))
Dsp_All<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/DspMax5E5_macrocarpa.txt",header = TRUE)
Dsp_All$chr2<-Dsp_All$chr
Dsp_All<-Dsp_All %>%
  dplyr::mutate(across('chr2',str_replace,'Superscaffold',''))
group.colors<-c(lyr="#000000",bic="#88CCEE",lob="#116644",mue="#AA5599",ste="#E69F00",alb="#332288")

highTD<-read_excel("C:/Users/rmohn/Desktop/10_Analysis/001_subsets/highTajima.D.xlsx",sheet = "Sheet1",col_names = TRUE,na = "NA")
highTD$chr<-highTD$CHROM
highTD$chr2<-highTD$CHROM
highTD<-highTD %>%
  dplyr::mutate(across('chr2',str_replace,'Superscaffold',''))
highTD$numpops<-rowSums(highTD[,7:50])

Dmse_highregions<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/DMSE_High_macrocarpa.txt",header=TRUE) 
Dmse_highregions_chr<-Dmse_highregions
Dmse_highregions_chr$chr2<-Dmse_highregions_chr$chr
Dmse_highregions_chr$chr<-paste("Superscaffold",Dmse_highregions_chr$chr2,sep="")

####making the figure legend
leg<-cbind(percentile=c("<95%",">95%",">99%",">99.9%"),perc=c(94,96,99.5,99.99),long=c(-80,-85,-90,-95),lat=c(40,40,40,40))
leg2<-ggplot()+
  geom_point(data=leg,mapping=aes(x=long, y=lat,shape=percentile,size=percentile))+
  scale_shape_manual(values=c(17,16, 21, 21))+
  scale_size_manual(values=c(0.75,1,2.5,4))+
  theme(text = element_text(size = 8))
  theme_light()
  
  
genes<-subset(read.delim(file = "C:/Users/rmohn/Downloads/Quercus_mongolica_gff3.gff3",sep="\t",header = FALSE),(V1=="Superscaffold1"|V1=="Superscaffold2"|V1=="Superscaffold3"|V1=="Superscaffold4"|V1=="Superscaffold5"|V1=="Superscaffold6"|V1=="Superscaffold7"|V1=="Superscaffold8"|V1=="Superscaffold9"|V1=="Superscaffold10"|V1=="Superscaffold11"|V1=="Superscaffold12")&V3=="CDS")
genes$chr<-genes$V1
genes$chr2<-genes$V1
genes<-genes %>%
  dplyr::mutate(across('chr2',str_replace,'Superscaffold',''))

cents<-read_excel("C:/Users/rmohn/Desktop/10_Analysis/11_GenomeRegions/DstatMSE_241011/Regions/Centromeres.xlsx")
cents$chr2<-cents$chr
cents$chr<-paste("Superscaffold",cents$chr,sep="")
  

for(j in 1:12){
  pdf(file=paste("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/FigsC",j,".pdf",sep=""),height=4.5,width=8)
  i=min(subset(Dsp_All,chr2==j)$windowmid)
  while(i <= max(subset(Dsp_All,chr2==j)$windowmid)){
    if(i %in% subset(Dsp_All,chr2==j)$windowmid){
      DwindMax<-subset(Dsp_All,windowmid==i&chr2==j)
      ####Map for each 500,000 bp genome region
      sp_fig<-ggplot()+
        geom_polygon(data = counties, aes(x = long, y = lat, group = group),fill = 'white', color = "white", lwd = 0.1)+
        geom_polygon(data = states, aes(x = long, y = lat, group = group),
                     fill = NA, color = "gray85")+
        geom_polygon(data = Canada, aes(x = long, y = lat, group = group),
                     fill = "white", color = "gray85")+
        coord_map(xlim = c(-105,-68),ylim = c(29,50))+
        scale_size(limits=c(.1,1),range=c(2,8))+
        geom_point(data=subset(DwindMax,d_f>=.24),mapping=aes(x=long, y=lat, color = species),size = 1, shape= 16)+
        geom_point(data=subset(DwindMax,d_f>.50),mapping=aes(x=long, y=lat, color = species),size = 2.5, shape= 21)+
        geom_point(data=subset(DwindMax,d_f>.79),mapping=aes(x=long, y=lat, color = species),size = 4, shape= 21)+
        #geom_point(data=subset(DwindMax,d_f>.75),mapping=aes(x=long, y=lat, color = species),size = 5.5, shape= 21)+
        geom_point(data=subset(DwindMax,d_f<.24),mapping=aes(x=long, y=lat), color = "black",size = .5,shape=17)+
        scale_color_manual(values=group.colors)+
        theme(axis.title = element_blank(),text = element_text(size = 8))+
        labs(caption=paste("Superscaffold ",j,": ",i,sep=""))+
        guides(color="none",size="none",shape="none")
      
      ####Graph of high D regions across the genome
      Csome_fig<-ggplot(data=Dsp_All)+
        #geom_rect(data=Dsd_highregions_chr,aes(ymin=-10,ymax=0,xmin=window,xmax=windowend),fill="#F0E442",alpha=.5)+
        geom_histogram(data=subset(Dsp_All,d_f>.5),mapping=aes(x=windowmid/1000000,fill=species),width=.5,binwidth = .5)+
        #geom_col(data=subset(highTD,chr=="Superscaffold1"|chr=="Superscaffold2"|chr=="Superscaffold3"|chr=="Superscaffold4"|chr=="Superscaffold5"|chr=="Superscaffold6"),mapping=aes(x=as.integer(START),y=-numpops),width=500000,alpha=.3)+
        geom_vline(data=subset(Dsp_All,chr2==j),mapping=aes(xintercept=i/1000000),color="red")+
        geom_blank(data=Dsp_All,aes(x=windowmid/1000000))+
        scale_fill_manual(values=group.colors)+
        facet_grid(.~as.numeric(chr2), scales = "free_x",space="free_x")+
        theme_light()+
        labs(x="Genomic window (Mbp)",y="Pops")+
        scale_x_continuous(breaks = seq(0,100,20)) +
        theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank(),panel.spacing.x=unit(0, "lines"),strip.text = element_blank(),text = element_text(size = 8),panel.border = element_rect(color = "black", fill = NA, linewidth = .5))+
        coord_cartesian(xlim = c(0, NA), ylim=c(0,55),expand = FALSE)#xlim(c(0,NA))
      
      ####Raw d_f in the window
      dfpoint<-ggplot(data=subset(D_All,chr2==j&windowMid>i-300000&windowMid<i+300000))+
        geom_rect(data=subset(Dmse_highregions_chr,chr2==j),aes(ymin=-1,ymax=1,xmin=windowstart/1000000,xmax=windowend/1000000),fill="darkgray",alpha=.25)+
        geom_rect(mapping=aes(ymin=-1,ymax=1,xmin=(i-250000)/1000000,xmax=(i+250000)/1000000),color="red",fill=NA)+
        coord_cartesian(xlim=c((i-450000)/1000000,(i+450000)/1000000), ylim=c(-1,1),expand = FALSE)+
        scale_color_manual(values=group.colors)+
        facet_grid(species~.)+
        geom_point(aes(x=windowMid/1000000,y=d_f,group=pop,color=species),size=.5)+
        #geom_vline(data=subset(D_All,chr==paste("Superscaffold",j,sep="")),mapping=aes(xintercept=i),color="red")+
      #geom_line(aes(x=windowMid,y=d_f,color=species,group=pop))+
        theme_light()+
        theme(panel.spacing.y=unit(0.5, "lines"),strip.text = element_blank(),text = element_text(size = 8))+
        guides(color="none")+
        ylab(bquote(d[f]))+
        xlab("window (Mbp)")

      ####Legend
      legmap<-as_ggplot(get_legend(leg2))

      ####Gene density in chromosome and centromere (estimated from Q. lobata)
      gen_dens<-ggplot(data=subset(genes,chr2==j))+
        geom_rect(data=subset(cents,chr2==j),mapping=aes(ymin=0,ymax=.035,xmin=centStart/1000000,xmax=centEnd/1000000),fill="lightblue",color="lightblue")+
        geom_vline(xintercept=i/1000000,color="red",linetype="dashed")+
        geom_density(data=subset(genes,chr2==j),aes(x=V4/1000000),bw=.5)+
        theme_light()+
        labs(x="window (Mbp)",y="gene density")+
        scale_x_continuous(breaks = seq(0,100,20)) +
          theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank(),text = element_text(size = 8))+
        coord_cartesian( ylim=c(0,.035),expand = FALSE)+
        theme(text = element_text(size = 8))+
        guides(color="none")

      ####Print out figure
      print(
        ggarrange(Csome_fig,
                  ggarrange(dfpoint,sp_fig,
                            ggarrange(legmap,gen_dens,ncol=1,heights=c(2,1)),
                            nrow=1,ncol=3,widths = c(1.5,3,1)),             
                  heights=c(1,3),ncol=1,nrow=2,common.legend = TRUE, legend = "right"))
      #ggsave(sp_fig,filename=paste("df_max_S",j,"_",i,".jpg",sep=""),device="jpeg",height=4,width=6,units="in")
    }
    i=i+500000
  }
  dev.off()
}
