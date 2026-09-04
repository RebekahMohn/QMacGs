#Windows across the genome
#100,000 bp window
library(readxl)
library(ggplot2)


#genomic regions
ChrStr<-read_excel("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/MacGenomeRegionsB.xlsx",sheet="Sheet3")
#assigning colors for use later
group.colors<-c(lyr="#000000",bic="#88CCEE",lob="#116644",mue="#AA5599",ste="#E69F00",alb="#332288",none="white",CenArray="gray20",InterArray="gray80", IntraSpInv="gray70",InterSpInv="gray35")

#renaming column for ease of use
ChrStr$csome<-ChrStr$scaffold

#importing chromosome dimensions
C_Coords<-read_excel("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/mac_csome_dims.xlsx",sheet="Sheet1")

#make windows for exploratioin
windowsMaster<-data.frame(windows=numeric(),present=numeric())
for(i in 1:12){
  if(i<10){
    windowsMaster<-rbind(windowsMaster,cbind(seq(from=0,to=C_Coords$stop[i],by=100000),rep(paste("Chr0",i,sep=""),trunc(C_Coords$stop[i]/100000)+1)))
  }else{
    windowsMaster<-rbind(windowsMaster,cbind(seq(from=0,to=C_Coords$stop[i],by=100000),rep(paste("Chr",i,sep=""),trunc(C_Coords$stop[i]/100000)+1))) 
  }
}


#combine chromosome and windows into a column for use merging data sets
windowsMaster$csome_window<-paste(windowsMaster$V2,"_",windowsMaster$V1,sep="")

#GEA windows
gea1_c<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_e1_lfmmpval_cor.txt",sep=" ")
gea2_c<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_e2_lfmmpval_cor.txt",sep=" ")
gea3_c<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_e3_lfmmpval_cor.txt",sep=" ")
gea4_c<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_e4_lfmmpval_cor.txt",sep=" ")
gea5_c<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_e5_lfmmpval_cor.txt",sep=" ")
gea6_c<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_e6_lfmmpval_cor.txt",sep=" ")


gea1_c$ENV<-"e1"
gea2_c$ENV<-"e2"
gea3_c$ENV<-"e3"
gea4_c$ENV<-"e4"
gea5_c$ENV<-"e5"
gea6_c$ENV<-"e6"

gea_env<-rbind(gea1_c,gea2_c,gea3_c,gea4_c,gea5_c,gea6_c)

write.table(gea_env,"C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_lfmmpval_cor.txt",row.names = FALSE,quote = FALSE)

gea_env<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_lfmmpval_cor.txt",sep=" ")

#make a column for counting
gea1_c$e1<-1
gea2_c$e2<-1
gea3_c$e3<-1
gea4_c$e4<-1
gea5_c$e5<-1
gea6_c$e6<-1
#round the windows to 100,000
gea1_c$roundpos<-round(gea1_c$BP/100000,0)*100000
gea2_c$roundpos<-round(gea2_c$BP/100000,0)*100000
gea3_c$roundpos<-round(gea3_c$BP/100000,0)*100000
gea4_c$roundpos<-round(gea4_c$BP/100000,0)*100000
gea5_c$roundpos<-round(gea5_c$BP/100000,0)*100000
gea6_c$roundpos<-round(gea6_c$BP/100000,0)*100000

#combine chromosome and windows into a column for use merging data sets
gea1_c$csome_window<-paste(gea1_c$csome,"_",gea1_c$roundpos,sep="")
gea2_c$csome_window<-paste(gea2_c$csome,"_",gea2_c$roundpos,sep="")
gea3_c$csome_window<-paste(gea3_c$csome,"_",gea3_c$roundpos,sep="")
gea4_c$csome_window<-paste(gea4_c$csome,"_",gea4_c$roundpos,sep="")
gea5_c$csome_window<-paste(gea5_c$csome,"_",gea5_c$roundpos,sep="")
gea6_c$csome_window<-paste(gea6_c$csome,"_",gea6_c$roundpos,sep="")

#aggregate and count how many significant SNPs per 100,000 basepairs
gea1_c_agg<-aggregate(data=gea1_c, e1~csome_window,function(x) sum(x, na.rm=TRUE))
gea2_c_agg<-aggregate(data=gea2_c, e2~csome_window,function(x) sum(x, na.rm=TRUE))
gea3_c_agg<-aggregate(data=gea3_c, e3~csome_window,function(x) sum(x, na.rm=TRUE))
gea4_c_agg<-aggregate(data=gea4_c, e4~csome_window,function(x) sum(x, na.rm=TRUE))
gea5_c_agg<-aggregate(data=gea5_c, e5~csome_window,function(x) sum(x, na.rm=TRUE))
gea6_c_agg<-aggregate(data=gea6_c, e6~csome_window,function(x) sum(x, na.rm=TRUE))





gea_env_agg<-merge(gea1_c_agg,gea2_c_agg,by="csome_window",all = TRUE)
gea_env_agg<-merge(gea_env_agg,gea3_c_agg,by="csome_window",all = TRUE)
gea_env_agg<-merge(gea_env_agg,gea4_c_agg,by="csome_window",all = TRUE)
gea_env_agg<-merge(gea_env_agg,gea5_c_agg,by="csome_window",all = TRUE)
gea_env_agg<-merge(gea_env_agg,gea6_c_agg,by="csome_window",all = TRUE)


#merge full window data set with windows of significant gea results
windowsMaster2<-merge(windowsMaster,gea_env_agg,by="csome_window",all.x = TRUE)


#Marking windows associated with inversions or centromeres
windowsMaster2$hapinvs<-as.numeric(0)
windowsMaster2$albinvs<-0
windowsMaster2$cents<-0

# haps_invs<-subset(read.delim(file="C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/Inversions/haps_inv_characterization_large_inversions_simplified_coords.txt"),coordinate_system=="REF")
# intersp_invs<-subset(read.delim(file="C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/Inversions/qmac_v_qalba_INV_large_inversions_simplified_coords.txt"),coordinate_system=="REF")
# haps_invs$chr_start<-paste(haps_invs$chromosome,"_",haps_invs$start,sep="")
# intersp_invs$chr_start<-paste(intersp_invs$chromosome,"_",intersp_invs$start,sep="")
# invs<-merge(haps_invs,intersp_invs,by="chr_start",all = T)
# write.table(invs,"C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/Invs2.txt",sep="\t",quote = F,row.names = F)
invs<-read.delim(file="C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/Invs2b.txt")

cents<-read_excel("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/MacGenomeRegionsB.xlsx",sheet="CentRegion")

###adding the inversion types to the windows

for(i in 1:length(windowsMaster2$csome_window)){
  for(j in 1:length(invs$chromosome)){
    if(windowsMaster2$V2[i]==invs$chromosome[j]&as.numeric(windowsMaster2$V1[i])+50000>invs$start[j] & (as.numeric(windowsMaster2$V1[i])-50000)<invs$end[j]){
      if(invs$INV[j]=="hap"){
        windowsMaster2$hapinvs[i] <- 1
      }else{
        windowsMaster2$albinvs[i] <- 1
      }
    }
  }
  # for(j in 1:length(intersp_invs$chromosome)){
  #   if(windowsMaster2$V2[i]==intersp_invs$chromosome[j]&as.numeric(windowsMaster2$V1[i])+50000>intersp_invs$start[j] & (as.numeric(windowsMaster2$V1[i])-50000)<intersp_invs$end[j]){
  #     windowsMaster2$albinvs[i] <- 1
  #   }
  # }
  for(j in 1:length(cents$scaffold)){
    if(windowsMaster2$V2[i]==cents$scaffold[j] & as.numeric(windowsMaster2$V1[i])+50000>cents$Start[j] & (as.numeric(windowsMaster2$V1[i])-50000)<cents$Stop[j]){
      windowsMaster2$cents[i] <- 1
    }
  }
}


#windowsMaster2[windowsMaster2$hapinvs==1&windowsMaster2$albinvs==1,8]<-0
#windowsMaster2$gea<-0
windowsMaster2[is.na(windowsMaster2$e1),4]<-0
windowsMaster2[is.na(windowsMaster2$e2),5]<-0
windowsMaster2[is.na(windowsMaster2$e3),6]<-0
windowsMaster2[is.na(windowsMaster2$e4),7]<-0
windowsMaster2[is.na(windowsMaster2$e5),8]<-0
windowsMaster2[is.na(windowsMaster2$e6),9]<-0



#ADD PAINT WINDOWS to the Master table

lobPaint<-subset(read.delim("C:/Users/rmohn/Desktop/10_Analysis/113_painting/highpower/lob_mac_bic_ItrgByLat_samp.tsv",header = TRUE),(X.Chr01.=="Chr01"|X.Chr01.=="Chr02"|X.Chr01.=="Chr03"|
                   X.Chr01.=="Chr04"|X.Chr01.=="Chr05"|X.Chr01.=="Chr06"|X.Chr01.=="Chr07"|X.Chr01.=="Chr08"|X.Chr01.=="Chr09"|X.Chr01.=="Chr10"|X.Chr01.=="Chr11"|X.Chr01.=="Chr12"))
# ggplot(lobPaint)+
#   geom_point(aes(x=X.0.,y=NA_character_..2,color=NA_character_..2))+
#   scale_color_viridis_c(direction=-1)+
#   facet_grid(X.Chr01.~.)


colnames(lobPaint)<-c("chr","window","species","fixedBases","corgeneVsp","rsq","meanval","medianval","meanP","lobnumDev","sd")
lobPaint$csome_window<-paste(lobPaint$chr,"_",lobPaint$window,sep="")

stePaint<-subset(read.delim("C:/Users/rmohn/Desktop/10_Analysis/113_painting/highpower/ste_mac_bic_ItrgByLat_samp.tsv",header = TRUE),(X.Chr01.=="Chr01"|X.Chr01.=="Chr02"|X.Chr01.=="Chr03"|
                                                                                                                                   X.Chr01.=="Chr04"|X.Chr01.=="Chr05"|X.Chr01.=="Chr06"|X.Chr01.=="Chr07"|X.Chr01.=="Chr08"|X.Chr01.=="Chr09"|X.Chr01.=="Chr10"|X.Chr01.=="Chr11"|X.Chr01.=="Chr12"))
# ggplot(stePaint)+
#   geom_point(aes(x=X.0.,y=NA_character_..2,color=NA_character_..2))+
#   scale_color_viridis_c(direction=-1)+
#   facet_grid(X.Chr01.~.)
# 
colnames(stePaint)<-c("chr","window","species","fixedBases","corgeneVsp","rsq","meanval","medianval","meanP","stenumDev","sd")
stePaint$csome_window<-paste(stePaint$chr,"_",stePaint$window,sep="")
#colnames(stePaint)<-c("chr","window","species","fixedBases","corgeneVsp","rsq","meanval","meanP","stenumDev","csome_window")




muePaint<-subset(read.delim("C:/Users/rmohn/Desktop/10_Analysis/113_painting/highpower/mue_mac_bic_ItrgByLat_samp.tsv",header = TRUE),(X.Chr01.=="Chr01"|X.Chr01.=="Chr02"|X.Chr01.=="Chr03"|
                                                                                                                                   X.Chr01.=="Chr04"|X.Chr01.=="Chr05"|X.Chr01.=="Chr06"|X.Chr01.=="Chr07"|X.Chr01.=="Chr08"|X.Chr01.=="Chr09"|X.Chr01.=="Chr10"|X.Chr01.=="Chr11"|X.Chr01.=="Chr12"))
# ggplot(muePaint)+
#   geom_point(aes(x=X.0.,y=NA_character_..2,color=NA_character_..2))+
#   scale_color_viridis_c(direction=-1)+
#   facet_grid(X.Chr01.~.)

colnames(muePaint)<-c("chr","window","species","fixedBases","corgeneVsp","rsq","meanval","medianval","meanP","muenumDev","sd")
muePaint$csome_window<-paste(muePaint$chr,"_",muePaint$window,sep="")
#colnames(muePaint)<-c("chr","window","species","fixedBases","corgeneVsp","rsq","meanval","meanP","muenumDev","csome_window")



albPaint<-subset(read.delim("C:/Users/rmohn/Desktop/10_Analysis/113_painting/highpower/alb_mac_bic_ItrgByLat_samp.tsv",header = TRUE),(X.Chr01.=="Chr01"|X.Chr01.=="Chr02"|X.Chr01.=="Chr03"|
                                                                                                                                    X.Chr01.=="Chr04"|X.Chr01.=="Chr05"|X.Chr01.=="Chr06"|X.Chr01.=="Chr07"|X.Chr01.=="Chr08"|X.Chr01.=="Chr09"|X.Chr01.=="Chr10"|X.Chr01.=="Chr11"|X.Chr01.=="Chr12"))
# ggplot(albPaint)+
#   geom_point(aes(x=X.0.,y=NA_character_..2,color=NA_character_..2))+
#   scale_color_viridis_c(direction=-1)+
#   facet_grid(X.Chr01.~.)

#albPaint$csome_window<-paste(albPaint$X.Chr01.,"_",albPaint$X.0.,sep="")
colnames(albPaint)<-c("achr","awindow","aspecies","afixedBases","acorgeneVsp","arsq","ameanval","amedianval","ameanP","albnumDev","asd")
albPaint$csome_window<-paste(albPaint$achr,"_",albPaint$awindow,sep="")

windowsMaster3<-merge(windowsMaster2,unique(lobPaint),by="csome_window",all.x = TRUE)
windowsMaster3<-merge(windowsMaster3,unique(stePaint),by="csome_window",all.x = TRUE)
windowsMaster3<-merge(windowsMaster3,unique(muePaint),by="csome_window",all.x = TRUE)
windowsMaster3<-merge(windowsMaster3,unique(albPaint),by="csome_window",all.x = TRUE)



length(subset(windowsMaster3,(albnumDev>20|lobnumDev>20|muenumDev>20|stenumDev>20)&(fixedBases.x>=40&fixedBases.y>=40&fixedBases>=40&afixedBases>=40))$albnumDev)


#subsetting introgressed and gea positive loci
adaptiveIntrogressed<-subset(windowsMaster3,(e1>=1|e2>=1|e3>=1|e4>=1|e5>=1|e6>=1)&(albnumDev>20|lobnumDev>20|muenumDev>20|stenumDev>20)&(fixedBases.x>=40&fixedBases.y>=40&fixedBases>=40&afixedBases>=40))
write.table(sep="\t",adaptiveIntrogressed,"C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/adaptiveIntrogressed3_samp.tsv")


length(subset(windowsMaster3,(albnumDev>20|lobnumDev>20|muenumDev>20|stenumDev>20)&(fixedBases.x>=40&fixedBases.y>=40&fixedBases>=40&afixedBases>=40))$csome_window)
length(subset(windowsMaster3,(fixedBases.x>=40&fixedBases.y>=40&fixedBases>=40&afixedBases>=40))$csome_window)
length(subset(windowsMaster3,(fixedBases.x>=20&fixedBases.y>=20&fixedBases>=20&afixedBases>=20))$csome_window)


length(subset(windowsMaster3,(e1>=1|e2>=1|e3>=1|e4>=1|e5>=1|e6>=1)&fixedBases.x>=40&fixedBases.y>=40&fixedBases>=40&afixedBases>=40)$e1)


length(subset(windowsMaster3,(e1>=1|e2>=1|e3>=1|e4>=1|e5>=1|e6>=1))$e1)
length(subset(windowsMaster3,(e1>=1|e2>=1|e3>=1|e4>=1|e5>=1|e6>=1)&hapinvs==1)$e1)
length(subset(windowsMaster3,(e1>=1|e2>=1|e3>=1|e4>=1|e5>=1|e6>=1)&albinvs==1)$e1)
length(subset(windowsMaster3,(e1>=1|e2>=1|e3>=1|e4>=1|e5>=1|e6>=1)&cents==1)$e1)
length(subset(windowsMaster3,(e1>=1|e2>=1|e3>=1|e4>=1|e5>=1|e6>=1)&(hapinvs==1|albinvs==1))$e1)
length(subset(windowsMaster3,(e1>=1|e2>=1|e3>=1|e4>=1|e5>=1|e6>=1)&((hapinvs==1|albinvs==1)&cents==1))$e1)

length(subset(windowsMaster3,albinvs==1|hapinvs==1)$e1)
length(subset(windowsMaster3,cents==1)$e1)

length(subset(windowsMaster3,(e1>=1|e2>=1|e3>=1|e4>=1|e5>=1|e6>=1)&(albnumDev>50|lobnumDev>50|muenumDev>50|stenumDev>50)&(fixedBases.x>=40&fixedBases.y>=40&fixedBases>=40&afixedBases>=40))$e1)

###################

CsomeMapsAI<-ggplot(adaptiveIntrogressed)+
  geom_rect(mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=14.5,ymax = 15.5),fill="slategray4")+
  geom_rect(data=subset(adaptiveIntrogressed,hapinvs==1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=8,ymax = 10),fill="violet")+
  geom_rect(data=subset(adaptiveIntrogressed,albinvs==1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=10,ymax = 12), fill="darkgreen")+
  geom_rect(data=subset(adaptiveIntrogressed,cents==1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=12,ymax = 14), fill="blue4")+
  geom_rect(data=subset(adaptiveIntrogressed,e1>=1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=41,ymax = 45,fill=e1))+
  geom_rect(data=subset(adaptiveIntrogressed,e2>=1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=46,ymax = 50,fill=e2))+
  geom_rect(data=subset(adaptiveIntrogressed,e3>=1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=51,ymax = 55,fill=e3))+
  geom_rect(data=subset(adaptiveIntrogressed,e4>=1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=56,ymax = 60,fill=e4))+
  geom_rect(data=subset(adaptiveIntrogressed,e5>=1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=61,ymax = 65,fill=e5))+
  geom_rect(data=subset(adaptiveIntrogressed,e6>=1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=56,ymax = 70,fill=e6))+
  geom_rect(data=subset(adaptiveIntrogressed,lobnumDev>20),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=21,ymax = 25,color=lobnumDev))+
  geom_rect(data=subset(adaptiveIntrogressed,stenumDev>20),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=26,ymax = 30,color=stenumDev))+
  geom_rect(data=subset(adaptiveIntrogressed,muenumDev>20),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=31,ymax = 35,color=muenumDev))+
  geom_rect(data=subset(adaptiveIntrogressed,albnumDev>20),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=36,ymax = 40,color=albnumDev))+
  
  facet_wrap(V2~.,ncol=2,dir="v")+
  #scale_color_manual(values=c(hapInvs="coral4",macalbInvs="forestgreen",Cents="blue4"))+
  #scale_fill_viridis_c(direction=-1,)+
  scale_color_stepsn(n.breaks=8,colors=c("white","pink","violet","purple","purple4","black"))+
  theme_light()+
  theme(text = element_text(size=8),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        axis.ticks.y = element_blank(),panel.spacing.x=unit(.2, "lines"),panel.spacing.y=unit(0, "lines"),
        strip.text = element_text(margin = margin(t = 0.01, r = 0.01, b = 0.01, l = 0.01, unit = "in")))

jpeg("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/CsomeMap2AdapIntro2_samp.jpg",height=600,width=1200)
print(CsomeMapsAI)
dev.off()


#counting how windows are related


length(subset(windowsMaster2,albinvs==1)$V2)
length(subset(windowsMaster2,hapinvs==1)$V2)
length(subset(windowsMaster2,cents==1)$V2)
length(subset(windowsMaster2,x>=1)$V2)


length(subset(windowsMaster2,x>=1&cents==1)$V2)
length(subset(windowsMaster2,x>=1&albinvs==1)$V2)
length(subset(windowsMaster2,x>=1&hapinvs==1)$V2)


length(subset(windowsMaster2,cents==1&hapinvs==1)$V2)

length(subset(windowsMaster2,cents==1&albinvs==1)$V2)


length(subset(windowsMaster2,cents==1&hapinvs==1&x>=1)$V2)

length(subset(windowsMaster2,cents==1&albinvs==1&x>=1)$V2)



####probability of multiple individuals being positive in 7663 samples

ggplot()+
  geom_histogram(aes(x=rbinom(n=76630,size=378,prob=0.0001)))

quantile(rbinom(n=76630,size=378,prob=0.05),probs = c(0.999999,0.000001))



ggplot()+
  geom_histogram(aes(x=rbinom(n=76630,size=625,prob=0.00005)))

quantile(rbinom(n=76630,size=625,prob=0.00008),probs = c(0.999999,0.000001))



#visualization
CsomeMaps<-ggplot(windowsMaster3)+
  geom_rect(mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=4,ymax = 5),fill="slategray4")+
  geom_rect(data=subset(windowsMaster3,hapinvs==1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=5,ymax = 9),fill="violet")+
  geom_rect(data=subset(windowsMaster3,albinvs==1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=6,ymax = 10), fill="darkgreen")+
  geom_rect(data=subset(windowsMaster3,cents==1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=11,ymax = 15), fill="blue4")+
  #geom_rect(data=subset(windowsMaster3,val>=10),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=16,ymax = 20,fill=val))+
  geom_rect(data=subset(windowsMaster3,lobnumDev>40&lobnumDev<379),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=21,ymax = 25,color=lobnumDev))+
  geom_rect(data=subset(windowsMaster3,stenumDev>40&stenumDev<379),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=26,ymax = 30,color=stenumDev))+
  geom_rect(data=subset(windowsMaster3,muenumDev>40&muenumDev<379),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=31,ymax = 35,color=muenumDev))+
  geom_rect(data=subset(windowsMaster3,albnumDev>40&albnumDev<379),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=36,ymax = 40,color=albnumDev))+
  
  geom_rect(data=subset(adaptiveIntrogressed,e1>=1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=41,ymax = 45,fill=e1))+
  geom_rect(data=subset(adaptiveIntrogressed,e2>=1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=46,ymax = 50,fill=e2))+
  geom_rect(data=subset(adaptiveIntrogressed,e3>=1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=51,ymax = 55,fill=e3))+
  geom_rect(data=subset(adaptiveIntrogressed,e4>=1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=56,ymax = 60,fill=e4))+
  geom_rect(data=subset(adaptiveIntrogressed,e5>=1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=61,ymax = 65,fill=e5))+
  geom_rect(data=subset(adaptiveIntrogressed,e6>=1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=56,ymax = 70,fill=e6))+
  
  
  facet_wrap(V2~.,ncol=2)+
  #scale_fill_viridis_b(direction=-1,)+
  scale_color_stepsn(n.breaks=8,colors=c("white","white","pink","violet","purple","purple4","black"))+
  theme_light()+
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),axis.ticks.y = element_blank())

pdf("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/CsomeMapwide_samp.pdf",height=8.5,width=11)
print(CsomeMaps)
dev.off()

jpeg("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/CsomeMapwide_samp.jpg",height=8.5,width=11,units="in",res=100)
print(CsomeMaps)
dev.off()







CsomeMaps2<-ggplot(windowsMaster3)+
  geom_rect(mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,y="chromosome",height=1),fill="slategray4")+
  geom_rect(data=subset(windowsMaster3,hapinvs==1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,y="hapinvs",height = 1),fill="violet")+
  geom_rect(data=subset(windowsMaster3,albinvs==1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,y="albinvs",height = 1), fill="darkgreen")+
  geom_rect(data=subset(windowsMaster3,cents==1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,y="cents",height = 1), fill="blue4")+
  #geom_rect(data=subset(windowsMaster3,val>=3),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,y="gea",height = 1,fill=val))+
  geom_rect(data=subset(windowsMaster3,lobnumDev>50),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,y="lob",height = 1,color=lobnumDev))+
  geom_rect(data=subset(windowsMaster3,stenumDev>50),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,y="ste",height = 1,color=stenumDev))+
  geom_rect(data=subset(windowsMaster3,muenumDev>50),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,y="mue",height = 1,color=muenumDev))+
  geom_rect(data=subset(windowsMaster3,albnumDev>50),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,y="alb",height = 1,color=albnumDev))+
  
  facet_wrap(V2~.,ncol=2)+
  scale_fill_viridis_c(direction=-1,)+
  scale_color_stepsn(n.breaks=8,colors=c("white","white","pink","violet","purple","purple4","black"))+
  theme_light()+
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),axis.ticks.y = element_blank())

pdf("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/CsomeMapwide2_samp.pdf",height=8.5,width=11)
print(CsomeMaps2)
dev.off()

jpeg("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/CsomeMapwide2_samp.jpg",height=8.5,width=11,units="in",res=150)
print(CsomeMaps2)
dev.off()











CsomeMaps<-ggplot(windowsMaster3)+
  geom_rect(mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=14.5,ymax = 15.5),fill="slategray4")+
  geom_rect(data=subset(windowsMaster3,hapinvs==1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=8,ymax = 10),fill="violet")+
  geom_rect(data=subset(windowsMaster3,albinvs==1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=10,ymax = 12), fill="darkgreen")+
  geom_rect(data=subset(windowsMaster3,cents==1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=12,ymax = 14), fill="blue4")+
  geom_rect(data=subset(windowsMaster3,val>=10),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=16,ymax = 20,fill=val))+
  geom_rect(data=subset(windowsMaster3,lobnumDev>40&lobnumDev<379),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=21,ymax = 25,color=lobnumDev))+
  geom_rect(data=subset(windowsMaster3,stenumDev>40&stenumDev<379),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=26,ymax = 30,color=stenumDev))+
  geom_rect(data=subset(windowsMaster3,muenumDev>40&muenumDev<379),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=31,ymax = 35,color=muenumDev))+
  geom_rect(data=subset(windowsMaster3,albnumDev>40&albnumDev<379),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=36,ymax = 40,color=albnumDev))+
  
  facet_wrap(V2~.,ncol=2,dir="v")+
  scale_fill_viridis_c(direction=-1,)+
  scale_color_stepsn(n.breaks=8,colors=c("white","white","pink","violet","purple","purple4","black"))+
  theme_light()+
  theme(text = element_text(size=8),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        axis.ticks.y = element_blank(),panel.spacing.x=unit(.2, "lines"),panel.spacing.y=unit(0, "lines"),
        strip.text = element_text(margin = margin(t = 0.01, r = 0.01, b = 0.01, l = 0.01, unit = "in")))

pdf("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/CsomeMapwide_samp.pdf",height=8.5,width=11)
print(CsomeMaps)
dev.off()

jpeg("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/CsomeMapwide_samp.jpg",height=5,width=9,units="in",res=300)
print(CsomeMaps)
dev.off()

IntrMaps<-ggplot(windowsMaster3)+
  geom_rect(mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=14.5,ymax = 15.5),fill="slategray4")+
  geom_rect(data=subset(windowsMaster3,hapinvs==1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=8,ymax = 10,fill="hapInvs"))+
  geom_rect(data=subset(windowsMaster3,albinvs==1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=10,ymax = 12, fill="macalbInvs"))+
  geom_rect(data=subset(windowsMaster3,cents==1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=12,ymax = 14, fill="Cents"))+
  #geom_rect(data=subset(windowsMaster3,val>=10),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=16,ymax = 20,fill=val))+
  geom_rect(data=subset(windowsMaster3,lobnumDev>40&lobnumDev<379),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=21,ymax = 25,color=lobnumDev))+
  geom_rect(data=subset(windowsMaster3,stenumDev>40&stenumDev<379),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=26,ymax = 30,color=stenumDev))+
  geom_rect(data=subset(windowsMaster3,muenumDev>40&muenumDev<379),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=31,ymax = 35,color=muenumDev))+
  geom_rect(data=subset(windowsMaster3,albnumDev>40&albnumDev<379),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=36,ymax = 40,color=albnumDev))+
  
  facet_wrap(V2~.,ncol=2,dir="v")+
  scale_fill_manual(values=c(hapInvs="coral4",macalbInvs="forestgreen",Cents="blue4"))+
  #  scale_fill_viridis_c(direction=-1,)+
  scale_color_stepsn(n.breaks=8,colors=c("white","white","pink","violet","purple","purple4","black"))+
  theme_light()+
  theme(text = element_text(size=8),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        axis.ticks.y = element_blank(),panel.spacing.x=unit(.2, "lines"),panel.spacing.y=unit(0, "lines"),
        strip.text = element_text(margin = margin(t = 0.01, r = 0.01, b = 0.01, l = 0.01, unit = "in")))

# pdf("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/CsomeMapwide.pdf",height=8.5,width=11)
# print(CsomeMaps)
# dev.off()

jpeg("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/IntrMapwide_samp.jpg",height=5,width=9,units="in",res=300)
print(IntrMaps)
dev.off()


StrMaps<-ggplot(windowsMaster3)+
  geom_rect(mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=14.5,ymax = 15.5),fill="slategray4")+
  geom_rect(data=subset(windowsMaster3,hapinvs==1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=8,ymax = 10,fill="hapInvs"))+
  geom_rect(data=subset(windowsMaster3,albinvs==1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=10,ymax = 12, fill="macalbInvs"))+
  geom_rect(data=subset(windowsMaster3,cents==1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=12,ymax = 14, fill="Cents"))+  #geom_rect(data=subset(windowsMaster3,val>=10),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=16,ymax = 20,fill=val))+
  # geom_rect(data=subset(windowsMaster3,lobnumDev>40&lobnumDev<379),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=21,ymax = 25,color=lobnumDev))+
  # geom_rect(data=subset(windowsMaster3,stenumDev>40&stenumDev<379),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=26,ymax = 30,color=stenumDev))+
  # geom_rect(data=subset(windowsMaster3,muenumDev>40&muenumDev<379),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=31,ymax = 35,color=muenumDev))+
  # geom_rect(data=subset(windowsMaster3,albnumDev>40&albnumDev<379),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=36,ymax = 40,color=albnumDev))+
  
  facet_wrap(V2~.,ncol=2,dir="v")+
  scale_fill_manual(values=c(hapInvs="coral4",macalbInvs="forestgreen",Cents="blue4"))+
  #  scale_fill_viridis_c(direction=-1,)+
  # scale_color_stepsn(n.breaks=8,colors=c("white","white","pink","violet","purple","purple4","black"))+
  theme_light()+
  theme(text = element_text(size=8),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        axis.ticks.y = element_blank(),panel.spacing.x=unit(.2, "lines"),panel.spacing.y=unit(0, "lines"),
        strip.text = element_text(margin = margin(t = 0.01, r = 0.01, b = 0.01, l = 0.01, unit = "in")))

# pdf("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/CsomeMapwide.pdf",height=8.5,width=11)
# print(CsomeMaps)
# dev.off()

jpeg("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/StrMapwide_samp.jpg",height=5,width=9,units="in",res=300)
print(StrMaps)
dev.off()

GEAMaps<-ggplot(windowsMaster3)+
  geom_rect(mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=14.5,ymax = 15.5),fill="slategray4")+
  geom_rect(data=subset(windowsMaster3,hapinvs==1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=8,ymax = 10,color="hapInvs"))+
  geom_rect(data=subset(windowsMaster3,albinvs==1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=10,ymax = 12, color="macalbInvs"))+
  geom_rect(data=subset(windowsMaster3,cents==1),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=12,ymax = 14, color="Cents"))+  geom_rect(data=subset(windowsMaster3,val>=10),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=16,ymax = 20,fill=val))+
  #geom_rect(data=subset(windowsMaster3,lobnumDev>40&lobnumDev<379),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=21,ymax = 25,color=lobnumDev))+
  #geom_rect(data=subset(windowsMaster3,stenumDev>40&stenumDev<379),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=26,ymax = 30,color=stenumDev))+
  #geom_rect(data=subset(windowsMaster3,muenumDev>40&muenumDev<379),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=31,ymax = 35,color=muenumDev))+
  #geom_rect(data=subset(windowsMaster3,albnumDev>40&albnumDev<379),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=36,ymax = 40,color=albnumDev))+
  
  facet_wrap(V2~.,ncol=2,dir="v")+
  scale_color_manual(values=c(hapInvs="coral4",macalbInvs="forestgreen",Cents="blue4"))+
  scale_fill_viridis_c(direction=-1,)+
  #scale_color_stepsn(n.breaks=8,colors=c("white","white","pink","violet","purple","purple4","black"))+
  theme_light()+
  theme(text = element_text(size=8),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        axis.ticks.y = element_blank(),panel.spacing.x=unit(.2, "lines"),panel.spacing.y=unit(0, "lines"),
        strip.text = element_text(margin = margin(t = 0.01, r = 0.01, b = 0.01, l = 0.01, unit = "in")))

# pdf("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/CsomeMapwide.pdf",height=8.5,width=11)
# print(CsomeMaps)
# dev.off()

jpeg("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/GEAMapwide_samp.jpg",height=5,width=9,units="in",res=300)
print(GEAMaps)
dev.off()





ggplot(data = subset(windowsMaster3,x>=1))+
  geom_bin2d(mapping=aes(x=x,y=stenumDev),bins=10)+
  scale_fill_viridis_c(trans="log10")#+
#facet_wrap(albinvs+hapinvs~cents)
ggplot(data = subset(windowsMaster3,x>=1))+
  geom_bin2d(mapping=aes(x=x,y=muenumDev),bins=10)+
  scale_fill_viridis_c(trans="log10")#+
#facet_wrap(albinvs+hapinvs~cents)
ggplot(data = subset(windowsMaster3,x>=1))+
  geom_bin2d(mapping=aes(x=x,y=albnumDev),bins=10)+
  scale_fill_viridis_c(trans="log10")#+
#facet_wrap(albinvs+hapinvs~cents)
ggplot(data = subset(windowsMaster3,x>=1))+
  geom_bin2d(mapping=aes(x=x,y=lobnumDev),bins=10)+
  scale_fill_viridis_c(trans="log10")#+
#facet_wrap(albinvs+hapinvs~cents)



ggplot(subset(windowsMaster2,chr=="Chr07"))+
  geom_rect(mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=0,ymax = 1),fill="slategray4")+
  geom_rect(data=subset(windowsMaster2,hapinvs==1&chr=="Chr07"),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=1,ymax = 2),fill="purple4")+
  geom_rect(data=subset(windowsMaster2,albinvs==1&chr=="Chr07"),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=2,ymax = 3), fill="lightgreen")+
  geom_rect(data=subset(windowsMaster2,cents==1&chr=="Chr07"),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=3,ymax = 4), fill="blue4")+
  geom_rect(data=subset(windowsMaster2,x>=1&chr=="Chr07"),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=4,ymax = 5,fill=x))+
  geom_rect(data=subset(windowsMaster2,NA_character_..2.x>0&chr=="Chr07"),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=5,ymax = 6,color=NA_character_..2.x))+
  geom_rect(data=subset(windowsMaster2,NA_character_..2.y>0&chr=="Chr07"),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=6,ymax = 7,color=NA_character_..2.y))+
  geom_rect(data=subset(windowsMaster2,numDev>5&chr=="Chr07"),mapping = aes(xmin = as.numeric(V1)-50000,xmax = as.numeric(V1)+50000,ymin=7,ymax = 8,color=numDev))+
  
  facet_grid(V2~.)+
  scale_fill_viridis_b(direction=-1)+
  theme_light()+
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())+
  scale_color_viridis_b(direction=-1)


