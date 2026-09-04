#Visualing Sweepfinder Results

#Figure rows--samples, columns--chromosomes
# x axis position in chromosome
# y axis number
# color = p-value
library(ggplot2)
library(ggpubr)
library(readxl)

group.colors_e<-c(e1="skyblue", e2="salmon",e3="gold3",e4="skyblue4",e5="blue2",e6="red4",
                  Centromere="#008b8b",Inversion="#5b0000",Chromosome="gray40")
invs<-read.delim(file="C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/Invs2b.txt")
invs$csome<-invs$chromosome
cents<-subset(read_excel("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/MacGenomeRegionsB.xlsx",sheet="simple"),Type=="CenArray")
cents$csome<-cents$scaffold
C_Coords<-read_excel("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/mac_csome_dims.xlsx",sheet="Sheet1")
C_Coords$csome<-C_Coords$chr

sf2_pops<-read.csv("C:/Users/rmohn/Desktop/10_Analysis/115_SweepFinder/pops.txt")
# sampPop<-read.table("C:/Users/rmohn/Desktop/10_Analysis/107_DSuite_MACREF/SetsGR_mac2.tsv",sep="\t")
# popdata<-read_xlsx("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/meta_pops.xlsx",)

BC6<-read.table("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/BioClim6vars.txt", sep="\t",header = T)
sampPop<-read.table("C:/Users/rmohn/Desktop/10_Analysis/107_DSuite_MACREF/SetsGR_mac2.tsv",sep="\t")
BC6_pops<-aggregate(cbind(wc2.1_30s_bio_2,wc2.1_30s_bio_6,wc2.1_30s_bio_8,wc2.1_30s_bio_13,wc2.1_30s_bio_14,wc2.1_30s_bio_18)~V2,data=merge(sampPop,BC6,by="V1"),FUN=mean)


#all
all_sf1<-data.frame()
sf2_count<-0
for (i in 1:12){
  csome_tmp<-data.frame()
  if(i<10){
    chrom<-paste("Chr0",i,sep="")
    for (f in 1:length(sf2_pops[,1])) {
      tmppop<-read.delim(paste("C:/Users/rmohn/Desktop/10_Analysis/115_SweepFinder/",sf2_pops[f,1],".",chrom,".pcf.sf2",sep=""),header = T,sep="\t")
      tmppop$pop<-sf2_pops[f,1]
      tmppop$csome<-chrom
      csome_tmp<-rbind(csome_tmp,tmppop)
      # for(i in 1:12){
      #   if(i<10){
      #     
      #   }else{
      #     
      #   }
      # }
      sf2_count<-sum(sf2_count,length(tmppop$location))
    }
    
    csome_tmp<-subset(csome_tmp,LR>1)
    #filter out centromere regions
    csome_cents<-subset(cents,csome==chrom)
    for (k in 1:length(csome_cents$Type)){
      csome_tmp<-subset(csome_tmp,location>csome_cents$Stop[k]|location<csome_cents$Start[k])
    }
    all_sf1<-rbind(csome_tmp,all_sf1)
  }else{
    chrom<-paste("Chr",i,sep="")
    for (f in 1:length(sf2_pops[,1])) {
      tmppop<-read.delim(paste("C:/Users/rmohn/Desktop/10_Analysis/115_SweepFinder/",sf2_pops[f,1],".",chrom,".pcf.sf2",sep=""),header = T,sep="\t")
      tmppop$pop<-sf2_pops[f,1]
      tmppop$csome<-chrom
      csome_tmp<-rbind(csome_tmp,tmppop)
      sf2_count<-sum(sf2_count,length(tmppop$location))
      # for(i in 1:12){
      #   if(i<10){
      #     
      #   }else{
      #     
      #   }
      # }
    }
    csome_tmp<-subset(csome_tmp,LR>2)
    #filter out centromere regions
    csome_cents<-subset(cents,csome==chrom)
    for (k in 1:length(csome_cents$Type)){
      csome_tmp<-subset(csome_tmp,location>csome_cents$Stop[k]|location<csome_cents$Start[k])
    }
    all_sf1<-rbind(csome_tmp,all_sf1)
    
  }
}

#ForCatchup

catchup<-read.delim("C:/Users/rmohn/Documents/GitHub/QMacGs/20_ToMACROCARPA/08_SweepFinder2/284_sf2_catchup/catchuppops", sep="\t",header=T)


for(j in 1:length(catchup$chr)){
  chrom<-catchup$chr[j]
  tmppop<-read.delim(paste("C:/Users/rmohn/Documents/GitHub/QMacGs/20_ToMACROCARPA/08_SweepFinder2/285_SF2_Output/",catchup$pop[j],".",catchup$chr[j],"_",catchup$pos[j],".pcf.sf2",sep=""),header = T,sep="\t")
  tmppop<-tmppop[-1,]
  tmppop$pop<-catchup$pop[j]
  tmppop$csome<-catchup$chr[j]
  sf2_count<-sum(sf2_count,length(tmppop$location))
  #csome_tmp<-rbind(csome_tmp,tmppop)
  # for(i in 1:12){
  #   if(i<10){
  #     
  #   }else{
  #     
  #   }
  # }
#}
  csome_tmp<-subset(tmppop,LR>1)
#filter out centromere regions
  csome_cents<-subset(cents,csome==chrom)
  for (k in 1:length(csome_cents$Type)){
    csome_tmp<-subset(csome_tmp,location>csome_cents$Stop[k]|location<csome_cents$Start[k])
  }
  all_sf1<-rbind(csome_tmp,all_sf1)
}

#### All Samples

maxloc<-aggregate(location~pop+csome, data =all_sf1, FUN=max)

all_sf1_clim<-merge(all_sf1, BC6_pops,by.x = "pop",by.y="V2")
all_sf1_clim$pop_f<-factor(all_sf1_clim$pop,levels=c("mac_MB_SWP","mac_WI_RLF","mac_IA_BSP","mac_MA_NMK","mac_MI_PCP",
                                                                 "mac_SD_CSP","mac_IN_BOW","mac_OH_DPS","mac_IL_CHB","mac_KS_SMP",
                                                                 "mac_KY_BBF","mac_KY_GRF","mac_TN_BCP","mac_OK_NOW","mac_TX_CHP","mac_TX_PLM"))

all_sf1_clim$roundlocation<-round(all_sf1_clim$location/50000)*50000
all_sf1_clim$chrloc<-paste(all_sf1_clim$csome,all_sf1_clim$roundlocation)
all_sf1_clim_ag<-aggregate(list(alpha=all_sf1_clim$alpha,LR=all_sf1_clim$LR), by=list(roundlocation=all_sf1_clim$roundlocation,
                                                                                      pop=all_sf1_clim$pop, csome=all_sf1_clim$csome,wc2.1_30s_bio_6=all_sf1_clim$wc2.1_30s_bio_6),FUN=median)

all_sf1_clim_ag<-aggregate(list(alpha=all_sf1_clim$alpha,LR=all_sf1_clim$LR), by=list(roundlocation=all_sf1_clim$roundlocation,
                                                                                      pop=all_sf1_clim$pop, csome=all_sf1_clim$csome,wc2.1_30s_bio_6=all_sf1_clim$wc2.1_30s_bio_6),FUN=median)
all_sf1_clim_grp<-all_sf1_clim %>% group_by(chrloc) %>% top_n(1, LR)

write.table(all_sf1_clim, file="C:/Users/rmohn/Documents/GitHub/QMacGs/20_ToMACROCARPA/08_SweepFinder2/285_SF2_Output/all_sf2.pcf.sf2",sep="\t",quote=FALSE,row.names=FALSE )
write.table(all_sf1_clim_grp, file="C:/Users/rmohn/Documents/GitHub/QMacGs/20_ToMACROCARPA/08_SweepFinder2/285_SF2_Output/all_sf2_grp.pcf.sf2",sep="\t",quote=FALSE,row.names=FALSE )


pdf(file="C:/Users/rmohn/Desktop/10_Analysis/115_SweepFinder/sf2_run_all.pdf",height=8.5,width=11)
ggplot(data=all_sf1_clim)+
  geom_rect(data=subset(cents),mapping=aes(ymin=-2000.5,ymax=-1000.5,xmin=Start/1000000,xmax=Stop/1000000, fill="Centromere"))+
  geom_rect(data=subset(invs),mapping=aes(ymin=-1000.25,ymax=-0.25,xmin=start/1000000,xmax=end/1000000, fill="Inversion"))+
  geom_rect(data=subset(C_Coords),mapping=aes(ymin=-.25,ymax=2,xmin=start/1000000,xmax=stop/1000000, fill="Chromosome"))+
  geom_point(data=subset(all_sf1_clim),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="gray40")+
 # geom_point(data=subset(all_sf1_clim,pBH>=0.05),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="gray80")+
  facet_grid(pop_f~csome,scales="free_x",space="free_x")+
  scale_color_viridis_c(direction=-1)+
  scale_fill_manual(values=group.colors_e)+
  theme_light()+
  theme(panel.grid=element_blank(),panel.spacing.x=unit(0,"lines"),panel.spacing.y=unit(0.1,"lines"))+
  labs(x = "Mbps")

ggplot(data=subset(all_sf1_clim_grp))+
   geom_histogram(mapping=aes(x=LR),binwidth=100)


ggplot(data=subset(all_sf1_clim,csome=="Chr01"|csome=="Chr03"|csome=="Chr09"))+
  geom_rect(data=subset(cents,csome=="Chr01"|csome=="Chr03"|csome=="Chr09"),mapping=aes(ymin=-2000.5,ymax=-1000.5,xmin=Start/1000000,xmax=Stop/1000000, fill="Centromere"))+
  geom_rect(data=subset(invs,csome=="Chr01"|csome=="Chr03"|csome=="Chr09"),mapping=aes(ymin=-1000.25,ymax=-0.25,xmin=start/1000000,xmax=end/1000000, fill="Inversion"))+
  geom_rect(data=subset(C_Coords,csome=="Chr01"|csome=="Chr03"|csome=="Chr09"),mapping=aes(ymin=-.25,ymax=2,xmin=start/1000000,xmax=stop/1000000, fill="Chromosome"))+
  # geom_point(data=subset(all_sf1_clim,pBH<0.05&(csome=="Chr01"|csome=="Chr03"|csome=="Chr09")),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="red4")+
  # geom_point(data=subset(all_sf1_clim,pBH>=0.05&(csome=="Chr01"|csome=="Chr03"|csome=="Chr09")),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="gray80")+
  geom_point(data=subset(all_sf1_clim,csome=="Chr01"|csome=="Chr03"|csome=="Chr09"),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="gray40")+
  facet_grid(pop_f~csome,scales="free_x",space="free_x")+
  scale_color_viridis_c(direction=-1)+
  scale_fill_manual(values=group.colors_e)+
  theme_light()+
  theme(panel.grid=element_blank(),panel.spacing.x=unit(0,"lines"),panel.spacing.y=unit(0.1,"lines"))+
  labs(x = "Mbps")

ggplot(data=subset(all_sf1_clim_ag,csome=="Chr01"|csome=="Chr03"|csome=="Chr02"))+
  geom_rect(data=subset(cents,csome=="Chr01"|csome=="Chr03"|csome=="Chr02"),mapping=aes(ymin=-2000.5,ymax=-1000.5,xmin=Start/1000000,xmax=Stop/1000000, fill="Centromere"))+
  geom_rect(data=subset(invs,csome=="Chr01"|csome=="Chr03"|csome=="Chr02"),mapping=aes(ymin=-1000.25,ymax=-0.25,xmin=start/1000000,xmax=end/1000000, fill="Inversion"))+
  geom_rect(data=subset(C_Coords,csome=="Chr01"|csome=="Chr03"|csome=="Chr02"),mapping=aes(ymin=-.25,ymax=2,xmin=start/1000000,xmax=stop/1000000, fill="Chromosome"))+
  geom_point(data=subset(all_sf1_clim,pBH<0.05&(csome=="Chr01"|csome=="Chr03"|csome=="Chr02")),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="red4")+
  geom_point(data=subset(all_sf1_clim,pBH>=0.05&(csome=="Chr01"|csome=="Chr03"|csome=="Chr02")),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="gray80")+
  facet_grid(pop_f~csome,scales="free_x",space="free_x")+
  scale_color_viridis_c(direction=-1)+
  scale_fill_manual(values=group.colors_e)+
  theme_light()+
  theme(panel.grid=element_blank(),panel.spacing.x=unit(0,"lines"),panel.spacing.y=unit(0.1,"lines"))+
  labs(x = "Mbps")

ggplot(data=subset(all_sf1_clim_ag,csome=="Chr04"|csome=="Chr05"|csome=="Chr06"))+
  geom_rect(data=subset(cents,csome=="Chr04"|csome=="Chr05"|csome=="Chr06"),mapping=aes(ymin=-2000.5,ymax=-1000.5,xmin=Start/1000000,xmax=Stop/1000000, fill="Centromere"))+
  geom_rect(data=subset(invs,csome=="Chr04"|csome=="Chr05"|csome=="Chr06"),mapping=aes(ymin=-1000.25,ymax=-0.25,xmin=start/1000000,xmax=end/1000000, fill="Inversion"))+
  geom_rect(data=subset(C_Coords,csome=="Chr04"|csome=="Chr05"|csome=="Chr06"),mapping=aes(ymin=-.25,ymax=2,xmin=start/1000000,xmax=stop/1000000, fill="Chromosome"))+
  geom_point(data=subset(all_sf1_clim,pBH<0.05&(csome=="Chr04"|csome=="Chr05"|csome=="Chr06")),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="red4")+
  geom_point(data=subset(all_sf1_clim,pBH>=0.05&(csome=="Chr04"|csome=="Chr05"|csome=="Chr06")),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="gray80")+
  facet_grid(pop_f~csome,scales="free_x",space="free_x")+
  scale_color_viridis_c(direction=-1)+
  scale_fill_manual(values=group.colors_e)+
  theme_light()+
  theme(panel.grid=element_blank(),panel.spacing.x=unit(0,"lines"),panel.spacing.y=unit(0.1,"lines"))+
  labs(x = "Mbps")

ggplot(data=subset(all_sf1_clim,csome=="Chr07"|csome=="Chr08"|csome=="Chr09"))+
  geom_rect(data=subset(cents,csome=="Chr07"|csome=="Chr08"|csome=="Chr09"),mapping=aes(ymin=-2000.5,ymax=-1000.5,xmin=Start/1000000,xmax=Stop/1000000, fill="Centromere"))+
  geom_rect(data=subset(invs,csome=="Chr07"|csome=="Chr08"|csome=="Chr09"),mapping=aes(ymin=-1000.25,ymax=-0.25,xmin=start/1000000,xmax=end/1000000, fill="Inversion"))+
  geom_rect(data=subset(C_Coords,csome=="Chr07"|csome=="Chr08"|csome=="Chr09"),mapping=aes(ymin=-.25,ymax=2,xmin=start/1000000,xmax=stop/1000000, fill="Chromosome"))+
  geom_point(data=subset(all_sf1_clim,pBH<0.05&(csome=="Chr07"|csome=="Chr08"|csome=="Chr09")),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="red4")+
  geom_point(data=subset(all_sf1_clim,pBH>=0.05&(csome=="Chr07"|csome=="Chr08"|csome=="Chr09")),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="gray80")+
  facet_grid(pop_f~csome,scales="free_x",space="free_x")+
  scale_color_viridis_c(direction=-1)+
  scale_fill_manual(values=group.colors_e)+
  theme_light()+
  theme(panel.grid=element_blank(),panel.spacing.x=unit(0,"lines"),panel.spacing.y=unit(0.1,"lines"))+
  labs(x = "Mbps")

ggplot(data=subset(all_sf1_clim,csome=="Chr10"|csome=="Chr11"|csome=="Chr12"))+
  geom_rect(data=subset(cents,csome=="Chr10"|csome=="Chr11"|csome=="Chr12"),mapping=aes(ymin=-2000.5,ymax=-1000.5,xmin=Start/1000000,xmax=Stop/1000000, fill="Centromere"))+
  geom_rect(data=subset(invs,csome=="Chr10"|csome=="Chr11"|csome=="Chr12"),mapping=aes(ymin=-1000.25,ymax=-0.25,xmin=start/1000000,xmax=end/1000000, fill="Inversion"))+
  geom_rect(data=subset(C_Coords,csome=="Chr10"|csome=="Chr11"|csome=="Chr12"),mapping=aes(ymin=-.25,ymax=2,xmin=start/1000000,xmax=stop/1000000, fill="Chromosome"))+
  geom_point(data=subset(all_sf1_clim,pBH<0.05&(csome=="Chr10"|csome=="Chr11"|csome=="Chr12")),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="red4")+
  geom_point(data=subset(all_sf1_clim,pBH>=0.05&(csome=="Chr10"|csome=="Chr11"|csome=="Chr12")),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="gray80")+
  facet_grid(pop_f~csome,scales="free_x",space="free_x")+
  scale_color_viridis_c(direction=-1)+
  scale_fill_manual(values=group.colors_e)+
  theme_light()+
  theme(panel.grid=element_blank(),panel.spacing.x=unit(0,"lines"),panel.spacing.y=unit(0.1,"lines"))+
  labs(x = "Mbps")

# ggplot(data=all_sf1_clim)+
#   geom_rect(data=subset(invs),mapping=aes(ymin=-1.25,ymax=-0.25,xmin=start,xmax=end, fill="Inversion"))+
#   geom_rect(data=subset(cents),mapping=aes(ymin=-2.5,ymax=-1.5,xmin=Start,xmax=Stop, fill="Centromere"))+
#   geom_rect(data=subset(C_Coords),mapping=aes(ymin=-.25,ymax=1,xmin=start,xmax=stop, fill="Chromosome"))+
#   geom_point(mapping=aes(x=location,y=log(LR),color=-log(`alpha`)),size=.1)+
#   facet_grid(paste(wc2.1_30s_bio_6,pop)~csome,scales="free_x",space="free_x")+
#   scale_color_steps(low="gray80",high = "red3",midpoint = 2.995732)+
#   scale_fill_manual(values=group.colors_e)+
#   theme_light()+
#   theme(panel.grid=element_blank())


dev.off()




ggplot(data=subset(all_sf1_clim,csome=="Chr01"))+
  geom_rect(data=subset(cents,csome=="Chr01"),mapping=aes(ymin=-2000.5,ymax=-1000.5,xmin=Start/1000000,xmax=Stop/1000000, fill="Centromere"))+
  geom_rect(data=subset(invs,csome=="Chr01"),mapping=aes(ymin=-1000.25,ymax=-0.25,xmin=start/1000000,xmax=end/1000000, fill="Inversion"))+
  geom_rect(data=subset(C_Coords,csome=="Chr01"),mapping=aes(ymin=-.25,ymax=2,xmin=start/1000000,xmax=stop/1000000, fill="Chromosome"))+
  geom_point(data=subset(all_sf1_clim,(csome=="Chr01")),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="gray40")+
  #geom_point(data=subset(all_sf1_clim_ag,pBH>=0.05&(csome=="Chr01")),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="gray80")+
  facet_grid(bio6_f~csome,scales="free_x",space="free_x")+
  scale_color_viridis_c(direction=-1)+
  scale_fill_manual(values=group.colors_e)+
  theme_light()+
  theme(panel.grid=element_blank(),panel.spacing.x=unit(0,"lines"),panel.spacing.y=unit(0.1,"lines"))+
  labs(x = "Mbps")

ggplot(data=subset(all_sf1_clim,csome=="Chr03"))+
  geom_rect(data=subset(cents,csome=="Chr03"),mapping=aes(ymin=-2000.5,ymax=-1000.5,xmin=Start/1000000,xmax=Stop/1000000, fill="Centromere"))+
  geom_rect(data=subset(invs,csome=="Chr03"),mapping=aes(ymin=-1000.25,ymax=-0.25,xmin=start/1000000,xmax=end/1000000, fill="Inversion"))+
  geom_rect(data=subset(C_Coords,csome=="Chr03"),mapping=aes(ymin=-.25,ymax=2,xmin=start/1000000,xmax=stop/1000000, fill="Chromosome"))+
  geom_point(data=subset(all_sf1_clim,(csome=="Chr03")),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="gray40")+
  #geom_point(data=subset(all_sf1_clim_ag,pBH>=0.05&(csome=="Chr03")),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="gray80")+
  facet_grid(paste(wc2.1_30s_bio_6,pop)~csome,scales="free_x",space="free_x")+
  scale_color_viridis_c(direction=-1)+
  scale_fill_manual(values=group.colors_e)+
  theme_light()+
  theme(panel.grid=element_blank(),panel.spacing.x=unit(0,"lines"),panel.spacing.y=unit(0.1,"lines"))+
  labs(x = "Mbps")

ggplot(data=subset(all_sf1_clim,csome=="Chr09"&location<39000000&location>37500000))+
  #geom_rect(data=subset(cents,csome=="Chr09"&location<39000000&location>37500000),mapping=aes(ymin=-2000.5,ymax=-1000.5,xmin=Start/1000000,xmax=Stop/1000000, fill="Centromere"))+
  #geom_rect(data=subset(invs,csome=="Chr09"),mapping=aes(ymin=-1000.25,ymax=-0.25,xmin=start/1000000,xmax=end/1000000, fill="Inversion"))+
  #geom_rect(data=subset(C_Coords,csome=="Chr09"),mapping=aes(ymin=-.25,ymax=2,xmin=start/1000000,xmax=stop/1000000, fill="Chromosome"))+
  geom_point(data=subset(all_sf1_clim,(csome=="Chr09")&location<39000000&location>37500000),mapping=aes(x=location/1000000,y=LR),size=.01,color="gray40")+
  #geom_point(data=subset(all_sf1_clim_ag,pBH>=0.05&(csome=="Chr09")),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="gray80")+
  facet_grid(pop_f~csome,scales="free_x",space="free_x")+
  scale_color_viridis_c(direction=-1)+
  scale_fill_manual(values=group.colors_e)+
  theme_light()+
  theme(panel.grid=element_blank(),panel.spacing.x=unit(0,"lines"),panel.spacing.y=unit(0.1,"lines"))+
  labs(x = "Mbps")


ggplot(data=subset(all_sf1_clim_grp))+
  geom_rect(data=subset(cents),mapping=aes(ymin=-2000.5,ymax=-1000.5,xmin=Start/1000000,xmax=Stop/1000000, fill="Centromere"))+
  geom_rect(data=subset(invs),mapping=aes(ymin=-1000.25,ymax=-0.25,xmin=start/1000000,xmax=end/1000000, fill="Inversion"))+
  geom_rect(data=subset(C_Coords),mapping=aes(ymin=-.25,ymax=2,xmin=start/1000000,xmax=stop/1000000, fill="Chromosome"))+
  geom_point(data=subset(all_sf1_clim_grp),mapping=aes(x=location/1000000,y=LR,color=wc2.1_30s_bio_6),size=.01)+
  #geom_point(data=subset(all_sf1_clim_ag,pBH>=0.05&(csome=="Chr09")),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="gray80")+
  facet_grid(csome~.,scales="free_x",space="free_x")+
  scale_color_viridis_c()+
  scale_fill_manual(values=group.colors_e)+
  theme_light()+
  theme(panel.grid=element_blank(),panel.spacing.x=unit(0,"lines"),panel.spacing.y=unit(0.1,"lines"))+
  labs(x = "Mbps")


ggplot(data=subset(all_sf1_clim_grp,LR>4308.210060))+
  geom_rect(data=subset(cents),mapping=aes(ymin=-2000.5,ymax=-1000.5,xmin=Start/1000000,xmax=Stop/1000000, fill="Centromere"))+
  geom_rect(data=subset(invs),mapping=aes(ymin=-1000.25,ymax=-0.25,xmin=start/1000000,xmax=end/1000000, fill="Inversion"))+
  geom_rect(data=subset(C_Coords),mapping=aes(ymin=-.25,ymax=2,xmin=start/1000000,xmax=stop/1000000, fill="Chromosome"))+
  geom_point(data=subset(all_sf1_clim_grp,LR>4308.210060),mapping=aes(x=location/1000000,y=LR,color=wc2.1_30s_bio_6),size=.01)+
  #geom_point(data=subset(all_sf1_clim_ag,pBH>=0.05&(csome=="Chr09")),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="gray80")+
  facet_grid(csome~.,scales="free_x",space="free_x")+
  scale_color_viridis_c()+
  scale_fill_manual(values=group.colors_e)+
  theme_light()+
  theme(panel.grid=element_blank(),panel.spacing.x=unit(0,"lines"),panel.spacing.y=unit(0.1,"lines"))+
  labs(x = "Mbps")

### Centromeres ###
CentCoords<-read_excel()

CentCoords<-read_excel("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/MacGenomeRegionsB.xlsx",sheet="CentRegion")
CentCoords$csome<-cents$scaffold

sf2_cen<-data.frame()
for(j in 1:length(CentCoords$csome)){
  sf2_cen<-rbind(subset(all_sf1_clim,location<CentCoords$Stop[j]+1000000&location>CentCoords$Start[j]-1000000&csome==CentCoords$csome[j]),sf2_cen)  
}


ggplot(data=sf2_cen)+
  geom_rect(data=subset(cents),mapping=aes(ymin=-2000.5,ymax=-1000.5,xmin=Start/1000000,xmax=Stop/1000000, fill="Centromere"))+
  #geom_rect(data=subset(invs),mapping=aes(ymin=-1000.25,ymax=-0.25,xmin=start/1000000,xmax=end/1000000, fill="Inversion"))+
  #geom_rect(data=subset(C_Coords),mapping=aes(ymin=-.25,ymax=2,xmin=start/1000000,xmax=stop/1000000, fill="Chromosome"))+
  geom_point(data=subset(sf2_cen),mapping=aes(x=roundlocation/1000000,y=LR,color=log(LR)),size=.01)+
  # geom_point(data=subset(all_sf1_clim,pBH>=0.05),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="gray80")+
  facet_grid(pop_f~csome,scales="free_x",space="free_x")+
  scale_color_viridis_c(direction=-1)+
  scale_fill_manual(values=group.colors_e)+
  theme_light()+
  theme(panel.grid=element_blank(),panel.spacing.x=unit(0.2,"lines"),panel.spacing.y=unit(0.2,"lines"))+
  labs(x = "Mbps")


### No centromeres
CentCoords<-read_excel()

CentCoords<-read_excel("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/MacGenomeRegionsB.xlsx",sheet="CentRegion")
CentCoords$csome<-cents$scaffold

sf2_n_cen<-data.frame()
for(j in 1:length(CentCoords$csome)){
  sf2_n_cen<-rbind(subset(all_sf1_clim,(location>CentCoords$Stop[j]+1000000|location<CentCoords$Start[j]-1000000)&csome==CentCoords$csome[j]),sf2_n_cen)  
}


sf2_n_cen_grp<-data.frame()
for(j in 1:length(CentCoords$csome)){
  sf2_n_cen_grp<-rbind(subset(all_sf1_clim_grp,(location>CentCoords$Stop[j]+1000000|location<CentCoords$Start[j]-1000000)&csome==CentCoords$csome[j]),sf2_n_cen_grp)  
}


ggplot(data=sf2_n_cen)+
  geom_rect(data=subset(cents),mapping=aes(ymin=-2000.5,ymax=-1000.5,xmin=Start/1000000,xmax=Stop/1000000, fill="Centromere"))+
  #geom_rect(data=subset(invs),mapping=aes(ymin=-1000.25,ymax=-0.25,xmin=start/1000000,xmax=end/1000000, fill="Inversion"))+
  #geom_rect(data=subset(C_Coords),mapping=aes(ymin=-.25,ymax=2,xmin=start/1000000,xmax=stop/1000000, fill="Chromosome"))+
  geom_point(data=subset(sf2_n_cen),mapping=aes(x=roundlocation/1000000,y=LR,color=log(LR)),size=.01)+
  # geom_point(data=subset(all_sf1_clim,pBH>=0.05),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="gray80")+
  facet_grid(pop_f~csome,scales="free_x",space="free_x")+
  scale_color_viridis_c(direction=-1)+
  scale_fill_manual(values=group.colors_e)+
  theme_light()+
  theme(panel.grid=element_blank(),panel.spacing.x=unit(0.2,"lines"),panel.spacing.y=unit(0.2,"lines"))+
  labs(x = "Mbps")

quantile(subset(sf2_n_cen_grp)$LR,probs = c(0.01,0.05,0.025,0.5,0.75,0.95,0.99))


all_sf1_clim_grp<-read.delim(file="C:/Users/rmohn/Documents/GitHub/QMacGs/20_ToMACROCARPA/08_SweepFinder2/285_SF2_Output/all_sf2_grp.pcf.sf2",sep="\t")


ggplot(data=subset(all_sf1_clim_grp,LR>4158.719))+
  geom_rect(data=subset(cents),mapping=aes(ymin=-2000.5,ymax=-1000.5,xmin=Start/1000000,xmax=Stop/1000000, fill="Centromere"))+
  geom_rect(data=subset(invs),mapping=aes(ymin=-1000.25,ymax=-0.25,xmin=start/1000000,xmax=end/1000000, fill="Inversion"))+
  #geom_rect(data=subset(C_Coords),mapping=aes(ymin=-.25,ymax=2,xmin=start/1000000,xmax=stop/1000000, fill="Chromosome"))+
  geom_point(data=subset(all_sf1_clim_grp,LR>4158.719),mapping=aes(x=location/1000000,y=LR,color=wc2.1_30s_bio_6),size=.01,color="red2")+
  geom_point(data=subset(all_sf1_clim_grp,LR<=4158.719),mapping=aes(x=location/1000000,y=LR,color=wc2.1_30s_bio_6),size=.01,color="gray70")+
  #geom_point(data=subset(all_sf1_clim_ag,pBH>=0.05&(csome=="Chr09")),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="gray80")+
  #facet_wrap(csome~.,ncol=6)+
  facet_grid(csome~.,scales="free_x",space="free_x")+
  scale_color_viridis_c()+
  scale_fill_manual(values=group.colors_e)+
  theme_light()+
  theme(panel.grid=element_blank(),panel.spacing.x=unit(0,"lines"),panel.spacing.y=unit(0.1,"lines"))+
  labs(x = "Mbps")


all_sf1_clim<-read.delim( file="C:/Users/rmohn/Documents/GitHub/QMacGs/20_ToMACROCARPA/08_SweepFinder2/285_SF2_Output/all_sf2.pcf.sf2",sep="\t" )
all_sf1_clim$pop_f<-factor(all_sf1_clim$pop,levels=c("mac_MB_SWP","mac_WI_RLF","mac_IA_BSP","mac_MA_NMK","mac_MI_PCP",
                                                     "mac_SD_CSP","mac_IN_BOW","mac_OH_DPS","mac_IL_CHB","mac_KS_SMP",
                                                     "mac_KY_BBF","mac_KY_GRF","mac_TN_BCP","mac_OK_NOW","mac_TX_CHP","mac_TX_PLM"))

sf2_pop_fig<-ggplot(data=all_sf1_clim)+
  geom_rect(data=subset(cents),mapping=aes(ymin=-2000.5,ymax=-1000.5,xmin=Start/1000000,xmax=Stop/1000000, fill="Centromere"))+
  geom_rect(data=subset(invs),mapping=aes(ymin=-1000.25,ymax=-0.25,xmin=start/1000000,xmax=end/1000000, fill="Inversion"))+
  geom_rect(data=subset(C_Coords),mapping=aes(ymin=-.25,ymax=2,xmin=start/1000000,xmax=stop/1000000, fill="Chromosome"))+
  
  geom_point(data=subset(all_sf1_clim,LR<=4158.719),mapping=aes(x=location/1000000,y=LR,color=wc2.1_30s_bio_6),size=.01,color="gray85")+
  geom_point(data=subset(all_sf1_clim,LR>4158.719),mapping=aes(x=location/1000000,y=LR,color=wc2.1_30s_bio_6),size=.01,color="red2")+
  # geom_point(data=subset(all_sf1_clim,pBH>=0.05),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="gray80")+
  facet_grid(pop_f~csome,scales="free_x",space="free_x")+
  #scale_color_viridis_c()+
  scale_fill_manual(values=group.colors_e)+
  theme_light()+
  theme(panel.grid=element_blank(),panel.spacing.x=unit(0,"lines"),panel.spacing.y=unit(0.1,"lines"),text = element_text(size=9))+
  labs(x = "Mbps")


sf2_all_fig<-ggplot(data=subset(all_sf1_clim_grp,LR>4158.719))+
  geom_rect(data=subset(cents),mapping=aes(ymin=-2000.5,ymax=-1000.5,xmin=Start/1000000,xmax=Stop/1000000, fill="Centromere"))+
  geom_rect(data=subset(invs),mapping=aes(ymin=-1000.25,ymax=-0.25,xmin=start/1000000,xmax=end/1000000, fill="Inversion"))+
  #geom_rect(data=subset(C_Coords),mapping=aes(ymin=-.25,ymax=2,xmin=start/1000000,xmax=stop/1000000, fill="Chromosome"))+
  geom_point(data=subset(all_sf1_clim_grp,LR<=4158.719),mapping=aes(x=location/1000000,y=LR,color=wc2.1_30s_bio_6),size=.01,color="gray85")+
  geom_point(data=subset(all_sf1_clim_grp,LR>4158.719),mapping=aes(x=location/1000000,y=LR,color=wc2.1_30s_bio_6),size=.01)+
  #geom_point(data=subset(all_sf1_clim_ag,pBH>=0.05&(csome=="Chr09")),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="gray80")+
  #facet_wrap(csome~.,ncol=6)+
  facet_grid(.~csome,scales="free_x",space="free_x")+
  scale_colour_viridis_c()+#gradient(low = "skyblue", high = "red", na.value = NA)+
  scale_fill_manual(values=group.colors_e)+
  theme_light()+
  theme(strip.text.y=element_text(angle=0),panel.grid=element_blank(),panel.spacing.x=unit(0,"lines"),panel.spacing.y=unit(0.1,"lines"),text = element_text(size=9))+
  labs(x = "Mbps",color="min. temp.")

library(ggpubr)
pdf(file="C:/Users/rmohn/Desktop/10_Analysis/115_SweepFinder/sf2_all_pops_supp.pdf",height=8.5,width=11)
ggarrange(sf2_all_fig,sf2_pop_fig,ncol=1,heights = c(1,5),align = "h",labels = "auto",common.legend = T,legend = "right")
dev.off()

jpeg(file="C:/Users/rmohn/Desktop/10_Analysis/115_SweepFinder/sf2_all_pops_supp.jpg",height=8.5,width=11,units="in",res=600)
ggarrange(sf2_all_fig,sf2_pop_fig,ncol=1,heights = c(1,5),align = "h",labels = "auto",common.legend = T,legend = "right")
dev.off()

pdf(file="C:/Users/rmohn/Desktop/10_Analysis/115_SweepFinder/sf2_all_pops_intext.pdf",height=9,width=3.75)

ggplot(data=subset(all_sf1_clim_grp,LR>4158.719))+
  geom_rect(data=subset(cents),mapping=aes(ymin=-2000.5,ymax=-1000.5,xmin=Start/1000000,xmax=Stop/1000000, fill="Centromere"))+
  geom_rect(data=subset(invs),mapping=aes(ymin=-1000.25,ymax=-0.25,xmin=start/1000000,xmax=end/1000000, fill="Inversion"))+
  #geom_rect(data=subset(C_Coords),mapping=aes(ymin=-.25,ymax=2,xmin=start/1000000,xmax=stop/1000000, fill="Chromosome"))+
  geom_point(data=subset(all_sf1_clim_grp,LR>4158.719),mapping=aes(x=location/1000000,y=LR,color=">99%"),size=.01)+
  geom_point(data=subset(all_sf1_clim_grp,LR<=4158.719),mapping=aes(x=location/1000000,y=LR,color="<=99%"),size=.01)+
  #geom_point(data=subset(all_sf1_clim_ag,pBH>=0.05&(csome=="Chr09")),mapping=aes(x=roundlocation/1000000,y=LR),size=.01,color="gray80")+
  #facet_wrap(csome~.,ncol=6)+
  facet_grid(csome~.,scales="free_x",space="free_x")+
  scale_color_manual(values=c(">99%"="red2","<=99%"="gray70"))+
  scale_fill_manual(values=group.colors_e)+
  theme_light()+
  theme(panel.grid=element_blank(),panel.spacing.x=unit(0,"lines"),panel.spacing.y=unit(0.1,"lines"),legend.position = "bottom",text=element_text(size=10))+
  labs(x = "Mbps")

dev.off()