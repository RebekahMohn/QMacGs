#tajimas and paint per region
library(readxl)
library(dplyr)
ChrStr<-read_excel("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/MacGenomeRegionsB.xlsx",sheet="Sheet3")
#ChrStr<-subset(ChrStr,Type2!="INVb")
#FstDxyA<-subset(read.csv("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/FstDxy_SP5_MACREF.csv",na.strings = "nan"),scaffold=="Chr01"|scaffold=="Chr02"|scaffold=="Chr03"|scaffold=="Chr04"|scaffold=="Chr05"|scaffold=="Chr06"|scaffold=="Chr07"|scaffold=="Chr08"|scaffold=="Chr09"|scaffold=="Chr10"|scaffold=="Chr11"|scaffold=="Chr12")
group.colors<-c(lyr="#000000",bic="#88CCEE",lob="#116644",mue="#AA5599",ste="#E69F00",alb="#332288",none="white",CenArray="gray20",InterArray="gray80", IntraSpInv="gray70",InterSpInv="gray35")


####paint slope
all_paint<-rbind(alb_mac_bic,lyr_mac_bic)
all_paint<-rbind(all_paint,lob_mac_bic)
all_paint<-rbind(all_paint,mue_mac_bic)
all_paint<-rbind(all_paint,ste_mac_bic)


all_paint<-rbind(alb_mac_bic,lob_mac_bic)
all_paint<-rbind(all_paint,mue_mac_bic)
all_paint<-rbind(all_paint,ste_mac_bic)

#all_paint_max<-aggregate(all_paint$slope,by=list(CHROM=all_paint$scaffold,window=all_paint$window),FUN=max)

all_paint_maxSlope<-all_paint %>%
  group_by(scaffold,window) %>%
  mutate(max = max(abs(slope),na.rm=T)) %>%
  ungroup()%>%
  filter(abs(slope)==max)
all_paint_maxSlope$CHROM<-all_paint_maxSlope$scaffold


all_paint_maxMean<-all_paint %>%
  group_by(scaffold,window) %>%
  mutate(max = max(mean,na.rm=T)) %>%
  ungroup()%>%
  filter(mean==max)
all_paint_maxMean$CHROM<-all_paint_maxMean$scaffold

all_paint_maxSlope$region<-rep("none",length(all_paint_maxSlope$scaffold))
for(i in 1:length(all_paint_maxSlope$scaffold)){
  for(j in 1:length(ChrStr$scaffold)){
    if(all_paint_maxSlope$scaffold[i]==ChrStr$scaffold[j]&all_paint_maxSlope$window[i]>ChrStr$Start[j] & (all_paint_maxSlope$window[i]+100000)<ChrStr$Stop[j]){
      all_paint_maxSlope$region[i] <- ChrStr$Type[j]
    }
  }
}


all_paint_maxMean$region<-rep("none",length(all_paint_maxMean$scaffold))
for(i in 1:length(all_paint_maxMean$scaffold)){
  for(j in 1:length(ChrStr$scaffold)){
    if(all_paint_maxMean$scaffold[i]==ChrStr$scaffold[j]&all_paint_maxMean$window[i]>ChrStr$Start[j] & (all_paint_maxMean$window[i]+100000)<ChrStr$Stop[j]){
      all_paint_maxMean$region[i]<-ChrStr$Type[j]
    }
  }
}

library(tidyverse)
library(egg)
library(multcompView)
library(reshape2)

all_paint_maxSlope$ChromRegion<-paste(all_paint_maxSlope$CHROM,all_paint_maxSlope$region,sep="_")
all_paint_maxSlope_ss<-subset(all_paint_maxSlope,region!="CenArray")

Letters.Slope<- data.frame(multcompLetters(TukeyHSD(aov(slope ~ ChromRegion, data = all_paint_maxSlope_ss))$ChromRegion[,4])$Letters)
colnames(Letters.Slope)[1]<-"Letter"
Letters.Slope$CHROM<-colsplit(rownames(Letters.Slope),names=c("CHROM","region"),pattern = "_")$CHROM
Letters.Slope$region<-colsplit(rownames(Letters.Slope),names=c("CHROM","region"),pattern = "_")$region


all_paint_maxSlope_C1<-subset(all_paint_maxSlope,CHROM=="Chr01")
all_paint_maxSlope_C1<-subset(all_paint_maxSlope,CHROM=="Chr02")
all_paint_maxSlope_C1<-subset(all_paint_maxSlope,CHROM=="Chr03")
all_paint_maxSlope_C1<-subset(all_paint_maxSlope,CHROM=="Chr04")
all_paint_maxSlope_C1<-subset(all_paint_maxSlope,CHROM=="Chr05")
all_paint_maxSlope_C1<-subset(all_paint_maxSlope,CHROM=="Chr06")
all_paint_maxSlope_C1<-subset(all_paint_maxSlope,CHROM=="Chr07")
all_paint_maxSlope_C1<-subset(all_paint_maxSlope,CHROM=="Chr08")
all_paint_maxSlope_C1<-subset(all_paint_maxSlope,CHROM=="Chr09")
all_paint_maxSlope_C1<-subset(all_paint_maxSlope,CHROM=="Chr10")
all_paint_maxSlope_C1<-subset(all_paint_maxSlope,CHROM=="Chr11")
all_paint_maxSlope_C1<-subset(all_paint_maxSlope,CHROM=="Chr12")



# TukeyHSD(aov(max ~ ChromRegion, data = all_paint_maxSlope_C1))$ChromRegion
# TukeyHSD(aov(max ~ ChromRegion, data = all_paint_maxSlope_C1))$ChromRegion
# TukeyHSD(aov(max ~ ChromRegion, data = all_paint_maxSlope_C1))$ChromRegion
# TukeyHSD(aov(max ~ ChromRegion, data = all_paint_maxSlope_C1))$ChromRegion
# TukeyHSD(aov(max ~ ChromRegion, data = all_paint_maxSlope_C1))$ChromRegion
# TukeyHSD(aov(max ~ ChromRegion, data = all_paint_maxSlope_C1))$ChromRegion
# TukeyHSD(aov(max ~ ChromRegion, data = all_paint_maxSlope_C1))$ChromRegion
# TukeyHSD(aov(max ~ ChromRegion, data = all_paint_maxSlope_C1))$ChromRegion
# TukeyHSD(aov(max ~ ChromRegion, data = all_paint_maxSlope_C1))$ChromRegion

all_paint_maxMean$ChromRegion<-paste(all_paint_maxMean$CHROM,all_paint_maxMean$region,sep="Q")


all_paint_maxMean_C1<-subset(all_paint_maxMean,CHROM=="Chr01")
all_paint_maxMean_C2<-subset(all_paint_maxMean,CHROM=="Chr02")
all_paint_maxMean_C3<-subset(all_paint_maxMean,CHROM=="Chr03")
all_paint_maxMean_C4<-subset(all_paint_maxMean,CHROM=="Chr04")
all_paint_maxMean_C5<-subset(all_paint_maxMean,CHROM=="Chr05")
all_paint_maxMean_C6<-subset(all_paint_maxMean,CHROM=="Chr06")
all_paint_maxMean_C7<-subset(all_paint_maxMean,CHROM=="Chr07")
all_paint_maxMean_C8<-subset(all_paint_maxMean,CHROM=="Chr08")
all_paint_maxMean_C9<-subset(all_paint_maxMean,CHROM=="Chr09")
all_paint_maxMean_C10<-subset(all_paint_maxMean,CHROM=="Chr10")
all_paint_maxMean_C11<-subset(all_paint_maxMean,CHROM=="Chr11")

all_paint_maxMean_C12<-subset(all_paint_maxMean,CHROM=="Chr12")





TukeyHSD(aov(max ~ ChromRegion, data = all_paint_maxMean_C1))$ChromRegion
TukeyHSD(aov(max ~ ChromRegion, data = all_paint_maxMean_C2))$ChromRegion
TukeyHSD(aov(max ~ ChromRegion, data = all_paint_maxMean_C3))$ChromRegion
TukeyHSD(aov(max ~ ChromRegion, data = all_paint_maxMean_C4))$ChromRegion
TukeyHSD(aov(max ~ ChromRegion, data = all_paint_maxMean_C5))$ChromRegion
TukeyHSD(aov(max ~ ChromRegion, data = all_paint_maxMean_C6))$ChromRegion
TukeyHSD(aov(max ~ ChromRegion, data = all_paint_maxMean_C7))$ChromRegion
TukeyHSD(aov(max ~ ChromRegion, data = all_paint_maxMean_C8))$ChromRegion
TukeyHSD(aov(max ~ ChromRegion, data = all_paint_maxMean_C9))$ChromRegion
TukeyHSD(aov(max ~ ChromRegion, data = all_paint_maxMean_C10))$ChromRegion
TukeyHSD(aov(max ~ ChromRegion, data = all_paint_maxMean_C11))$ChromRegion
TukeyHSD(aov(max ~ ChromRegion, data = all_paint_maxMean_C12))$ChromRegion
###Chr01
all_paint_maxSlope_C1<-subset(all_paint_maxSlope,CHROM=="Chr01")

Letters.SlopeC1<- data.frame(multcompLetters(TukeyHSD(aov(slope ~ ChromRegion, data = all_paint_maxSlope_C1))$ChromRegion[,4])$Letters)
colnames(Letters.SlopeC1)[1]<-"Letter"
Letters.SlopeC1$CHROM<-colsplit(rownames(Letters.SlopeC1),names=c("CHROM","region"),pattern = "_")$CHROM
Letters.SlopeC1$region<-colsplit(rownames(Letters.SlopeC1),names=c("CHROM","region"),pattern = "_")$region
# 
# 
# 
# ###Chr02
# all_paint_maxSlope_C2<-subset(all_paint_maxSlope,CHROM=="Chr02")
# 
# Letters.SlopeC2<- data.frame(multcompLetters(TukeyHSD(aov(max ~ ChromRegion, data = all_paint_maxSlope_C2))$ChromRegion[,4])$Letters)
# colnames(Letters.SlopeC2)[1]<-"Letter"
# Letters.SlopeC2$CHROM<-colsplit(rownames(Letters.SlopeC2),names=c("CHROM","region"),pattern = "_")$CHROM
# Letters.SlopeC2$region<-colsplit(rownames(Letters.SlopeC2),names=c("CHROM","region"),pattern = "_")$region
# 
# ###Chr03
# all_paint_maxSlope_C3<-subset(all_paint_maxSlope,CHROM=="Chr03")
# 
# Letters.SlopeC3<- data.frame(multcompLetters(TukeyHSD(aov(max ~ ChromRegion, data = all_paint_maxSlope_C3))$ChromRegion[,4])$Letters)
# colnames(Letters.SlopeC3)[1]<-"Letter"
# Letters.SlopeC3$CHROM<-colsplit(rownames(Letters.SlopeC3),names=c("CHROM","region"),pattern = "_")$CHROM
# Letters.SlopeC3$region<-colsplit(rownames(Letters.SlopeC3),names=c("CHROM","region"),pattern = "_")$region
# 
# ###Chr04
# all_paint_maxSlope_C4<-subset(all_paint_maxSlope,CHROM=="Chr04")
# 
# Letters.SlopeC4<- data.frame(multcompLetters(TukeyHSD(aov(slope ~ ChromRegion, data = all_paint_maxSlope_C4))$ChromRegion[,4])$Letters)
# colnames(Letters.SlopeC4)[1]<-"Letter"
# Letters.SlopeC4$CHROM<-colsplit(rownames(Letters.SlopeC4),names=c("CHROM","region"),pattern = "_")$CHROM
# Letters.SlopeC4$region<-colsplit(rownames(Letters.SlopeC4),names=c("CHROM","region"),pattern = "_")$region
# 
# ###Chr05
# all_paint_maxSlope_C5<-subset(all_paint_maxSlope,CHROM=="Chr05")
# 
# Letters.SlopeC5<- data.frame(multcompLetters(TukeyHSD(aov(slope ~ ChromRegion, data = all_paint_maxSlope_C5))$ChromRegion[,4])$Letters)
# colnames(Letters.SlopeC5)[1]<-"Letter"
# Letters.SlopeC5$CHROM<-colsplit(rownames(Letters.SlopeC5),names=c("CHROM","region"),pattern = "_")$CHROM
# Letters.SlopeC5$region<-colsplit(rownames(Letters.SlopeC5),names=c("CHROM","region"),pattern = "_")$region
# 
# ###Chr06
# all_paint_maxSlope_C6<-subset(all_paint_maxSlope,CHROM=="Chr06")
# 
# Letters.SlopeC6<- data.frame(multcompLetters(TukeyHSD(aov(slope ~ ChromRegion, data = all_paint_maxSlope_C6))$ChromRegion[,4])$Letters)
# colnames(Letters.SlopeC6)[1]<-"Letter"
# Letters.SlopeC6$CHROM<-colsplit(rownames(Letters.SlopeC6),names=c("CHROM","region"),pattern = "_")$CHROM
# Letters.SlopeC6$region<-colsplit(rownames(Letters.SlopeC6),names=c("CHROM","region"),pattern = "_")$region
# 
# ###Chr07
# all_paint_maxSlope_C7<-subset(all_paint_maxSlope,CHROM=="Chr07")
# 
# Letters.SlopeC7<- data.frame(multcompLetters(TukeyHSD(aov(slope ~ ChromRegion, data = all_paint_maxSlope_C7))$ChromRegion[,4])$Letters)
# colnames(Letters.SlopeC7)[1]<-"Letter"
# Letters.SlopeC7$CHROM<-colsplit(rownames(Letters.SlopeC7),names=c("CHROM","region"),pattern = "_")$CHROM
# Letters.SlopeC7$region<-colsplit(rownames(Letters.SlopeC7),names=c("CHROM","region"),pattern = "_")$region
# 
# ###C     hr08
# all_paint_maxSlope_C8<-subset(all_paint_maxSlope,CHROM=="Chr08")
# 
# Letters.SlopeC8<- data.frame(multcompLetters(TukeyHSD(aov(slope ~ ChromRegion, data = all_paint_maxSlope_C8))$ChromRegion[,4])$Letters)
# colnames(Letters.SlopeC8)[1]<-"Letter"
# Letters.SlopeC8$CHROM<-colsplit(rownames(Letters.SlopeC8),names=c("CHROM","region"),pattern = "_")$CHROM
# Letters.SlopeC8$region<-colsplit(rownames(Letters.SlopeC8),names=c("CHROM","region"),pattern = "_")$region
# 
###Chr09
all_paint_maxSlope_C9<-subset(all_paint_maxSlope,CHROM=="Chr09")

Letters.SlopeC9<- data.frame(multcompLetters(TukeyHSD(aov(slope ~ ChromRegion, data = all_paint_maxSlope_C9))$ChromRegion[,4])$Letters)
colnames(Letters.SlopeC9)[1]<-"Letter"
Letters.SlopeC9$CHROM<-colsplit(rownames(Letters.SlopeC9),names=c("CHROM","region"),pattern = "_")$CHROM
Letters.SlopeC9$region<-colsplit(rownames(Letters.SlopeC9),names=c("CHROM","region"),pattern = "_")$region

# ###Chr11
# all_paint_maxSlope_C11<-subset(all_paint_maxSlope,CHROM=="Chr11")
# 
# Letters.SlopeC11<- data.frame(multcompLetters(TukeyHSD(aov(slope ~ ChromRegion, data = all_paint_maxSlope_C11))$ChromRegion[,4])$Letters)
# colnames(Letters.SlopeC11)[1]<-"Letter"
# Letters.SlopeC11$CHROM<-colsplit(rownames(Letters.SlopeC11),names=c("CHROM","region"),pattern = "_")$CHROM
# Letters.SlopeC11$region<-colsplit(rownames(Letters.SlopeC11),names=c("CHROM","region"),pattern = "_")$region
# 
# ###Chr12
# all_paint_maxSlope_C12<-subset(all_paint_maxSlope,CHROM=="Chr12")
# 
# Letters.SlopeC12<- data.frame(multcompLetters(TukeyHSD(aov(slope ~ ChromRegion, data = all_paint_maxSlope_C12))$ChromRegion[,4])$Letters)
# colnames(Letters.SlopeC12)[1]<-"Letter"
# Letters.SlopeC12$CHROM<-colsplit(rownames(Letters.SlopeC12),names=c("CHROM","region"),pattern = "_")$CHROM
# Letters.SlopeC12$region<-colsplit(rownames(Letters.SlopeC12),names=c("CHROM","region"),pattern = "_")$region
# 
# ###Chr10
# all_paint_maxSlope_C10<-subset(all_paint_maxSlope,CHROM=="Chr10")
# 
# Letters.SlopeC10<- data.frame(multcompLetters(TukeyHSD(aov(slope ~ ChromRegion, data = all_paint_maxSlope_C10))$ChromRegion[,4])$Letters)
# colnames(Letters.SlopeC10)[1]<-"Letter"
# Letters.SlopeC10$CHROM<-colsplit(rownames(Letters.SlopeC10),names=c("CHROM","region"),pattern = "_")$CHROM
# Letters.SlopeC10$region<-colsplit(rownames(Letters.SlopeC10),names=c("CHROM","region"),pattern = "_")$region

# Letters.ConcUniq<- data.frame(multcompLetters(TukeyHSD(aov(`Reads concordant  1x`/`Total reads` ~ distMethod, data = Mapping1))$distMethod[,4])$Letters)
# colnames(Letters.ConcUniq)[1]<-"Letter"
# Letters.ConcUniq$Distance<-colsplit(rownames(Letters.ConcUniq),names=c("Distance","Method"),pattern = "Q")$Distance
# Letters.ConcUniq$Method<-colsplit(rownames(Letters.ConcUniq),names=c("Distance","Method"),pattern = "Q")$Method
# 

ggplot(data=subset(all_paint_maxSlope,region=="Cen"|region=="none"))+
  geom_violin(mapping = aes(y=abs(slope),x=region))+
  facet_grid(.~CHROM)+
  theme_light()

slope_bp<-ggplot(data=subset(all_paint_maxSlope,region!="CenArray"))+
  geom_hline(yintercept = 0.1,color="blue4")+
  geom_hline(yintercept = .5,color="blue4",alpha=.25)+
  geom_hline(yintercept = .9,color="blue4")+
  geom_boxplot(mapping = aes(y=abs(slope),x=region))+
  facet_grid(.~CHROM)+
  theme_light()+
  labs(x=NULL,y= "lat v. mean")+
  geom_text(data=subset(Letters.Slope,region!="CenArray"),aes(x=region, label=Letter,y=1.6),size=2, position = position_dodge(width = .75))+
  theme(axis.text.x = element_text(angle = 45, hjust =1),strip.background = element_blank(), strip.text.x = element_blank(),
        panel.spacing.x=unit(0, "lines"),text = element_text(size = 12),panel.border = element_rect(color = "black", fill = NA, linewidth = .5),panel.grid.major = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.minor.x = element_blank())




ggplot(data=subset(all_paint_maxMean,region=="Cen"|region=="none"))+
  geom_violin(mapping = aes(y=mean,x=region))+
  facet_grid(.~CHROM)+
  theme_light()

all_paint_maxMean$ChromRegion<-paste(all_paint_maxMean$CHROM,all_paint_maxMean$region,sep="_")


Letters.Mean<- data.frame(multcompLetters(TukeyHSD(aov(mean ~ ChromRegion, data = all_paint_maxMean))$ChromRegion[,4])$Letters)
colnames(Letters.Mean)[1]<-"Letter"
Letters.Mean$CHROM<-colsplit(rownames(Letters.Mean),names=c("CHROM","region"),pattern = "Q")$CHROM
Letters.Mean$region<-colsplit(rownames(Letters.Mean),names=c("CHROM","region"),pattern = "Q")$region


mean_bp<-ggplot(data=subset(all_paint_maxMean,region!="CenArray"))+
  geom_hline(yintercept = .1,color="blue4")+
  geom_hline(yintercept = .5,color="blue4",alpha=.25)+
  geom_hline(yintercept = .9,color="blue4")+
  geom_boxplot(mapping = aes(y=mean,x=region))+
  facet_grid(.~CHROM)+
  theme_light()+
  labs(x=NULL)+
  geom_text(data=subset(Letters.Mean,region!="CenArray"),aes(x=region, label=Letter,y=1.6),size=2, position = position_dodge(width = .75))+
  theme(axis.text.x = element_text(angle = 45, hjust =1),
        panel.spacing.x=unit(0, "lines"),text = element_text(size = 12),panel.border = element_rect(color = "black", fill = NA, linewidth = .5),panel.grid.major = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.minor.x = element_blank())


####Tajimas D
TD_2pops<-subset(TD_All_pop,(POP=="OK_NOW"|POP=="MN_QPA"))

TD_2pops$region<-rep("none",length(TD_2pops$CHROM))
for(i in 1:length(TD_2pops$CHROM)){
  for(j in 1:length(ChrStr$scaffold)){
    if(TD_2pops$CHROM[i]==ChrStr$scaffold[j]&as.numeric(TD_2pops$BIN_START[i])>ChrStr$Start[j] & (TD_2pops$BIN_START[i]+5000)<ChrStr$Stop[j]){
      TD_2pops$region[i] <- ChrStr$Type[j]
    }
  }
}

TD_2pops_ss<-subset(TD_2pops,region!="CenArray")

ggplot(data=subset(TD_2pops,region=="Cen"|region=="none"))+
  geom_boxplot(mapping = aes(y=TajimaD,x=region))+
  facet_grid(POP~CHROM)+
  theme_light()


TD_2pops_ss$ChromRegion<-paste(TD_2pops_ss$CHROM,TD_2pops_ss$region,sep="_")


Letters.TD_N<- data.frame(multcompLetters(TukeyHSD(aov(TajimaD ~ ChromRegion, data = subset(TD_2pops_ss,POP=="MN_QPA")))$ChromRegion[,4])$Letters)
colnames(Letters.TD_N)[1]<-"Letter"
Letters.TD_N$CHROM<-colsplit(rownames(Letters.TD_N),names=c("CHROM","region"),pattern = "_")$CHROM
Letters.TD_N$region<-colsplit(rownames(Letters.TD_N),names=c("CHROM","region"),pattern = "_")$region


TD_bp_N<-ggplot(data=subset(TD_2pops_ss,POP=="MN_QPA"))+
  geom_hline(yintercept = 2,color="blue4")+
  geom_hline(yintercept = 0,color="blue4",alpha=.25)+
  geom_hline(yintercept = -2,color="blue4")+
  geom_boxplot(mapping = aes(y=TajimaD,x=region))+
  facet_grid(.~CHROM)+
  labs(x = NULL,y="Tajima's D - North")+
  theme_light()+
  geom_text(data=subset(Letters.TD_N),aes(x=region, label=Letter,y=3.8),size=2, position = position_dodge(width = .75))+
  
  theme(axis.text.x = element_text(angle = 45, hjust =1),strip.background = element_blank(), strip.text.x = element_blank(),
        panel.spacing.x=unit(0, "lines"),text = element_text(size = 12),panel.border = element_rect(color = "black", fill = NA, linewidth = .5),panel.grid.major = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.minor.x = element_blank())


Letters.TD_S<- data.frame(multcompLetters(TukeyHSD(aov(TajimaD ~ ChromRegion, data = subset(TD_2pops_ss,POP=="OK_NOW")))$ChromRegion[,4])$Letters)
colnames(Letters.TD_S)[1]<-"Letter"
Letters.TD_S$CHROM<-colsplit(rownames(Letters.TD_S),names=c("CHROM","region"),pattern = "_")$CHROM
Letters.TD_S$region<-colsplit(rownames(Letters.TD_S),names=c("CHROM","region"),pattern = "_")$region


TD_bp_S<-ggplot(data=subset(TD_2pops_ss,POP=="OK_NOW"))+
  geom_hline(yintercept = 2,color="blue4")+
  geom_hline(yintercept = 0,color="blue4",alpha=.25)+
  geom_hline(yintercept = -2,color="blue4")+
  geom_boxplot(mapping = aes(y=TajimaD,x=region))+
  facet_grid(.~CHROM)+
  theme_light()+
  labs(x = "Genome Region",y="Tajima's D - South")+
  geom_text(data=subset(Letters.TD_S),aes(x=region, label=Letter,y=3.8),size=2, position = position_dodge(width = .75))+
  
  theme(axis.text.x = element_text(angle = 45, hjust =1),strip.background = element_blank(), strip.text.x = element_blank(),
        panel.spacing.x=unit(0, "lines"),text = element_text(size = 12),panel.border = element_rect(color = "black", fill = NA, linewidth = .5),panel.grid.major = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.minor.x = element_blank())


pdf("BoxPlot_PaintTD_Cen.pdf",height = 5, width = 8)
#ggarrange(mean_bp, slope_bp,TD_bp_N,TD_bp_S,ncol=1)
ggarrange(TD_bp_N,TD_bp_S,ncol=1)

dev.off()



###########INV


ggplot(data=subset(all_paint_maxSlope,region=="InterSpInv"|region=="IntraSpInv"|region=="none"))+
  geom_violin(mapping = aes(y=abs(slope),x=region))+
  facet_grid(.~CHROM)+
  theme_light()

slope_bp<-ggplot(data=subset(all_paint_maxSlope,region=="InterSpInv"|region=="IntraSpInv"|region=="none"))+
  geom_hline(yintercept = 0.1,color="blue4")+
  geom_hline(yintercept = .5,color="blue4",alpha=.25)+
  geom_hline(yintercept = .9,color="blue4")+
  geom_boxplot(mapping = aes(y=abs(slope),x=region))+
  facet_grid(.~CHROM)+
  theme_light()+
  labs(x=NULL,y= "lat v. mean")+
  #geom_text(data=subset(Letters.Slope,region=="InterSpInv"|region=="IntraSpInv"|region=="none"),aes(x=region, label=Letter,y=1.5),size=2, position = position_dodge(width = .75))+
  
  theme(axis.text.x = element_text(angle = 45, hjust =1),strip.background = element_blank(), strip.text.x = element_blank(),
        panel.spacing.x=unit(0, "lines"),text = element_text(size = 10),panel.border = element_rect(color = "black", fill = NA, linewidth = .5),panel.grid.major = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.minor.x = element_blank())




ggplot(data=subset(all_paint_maxMean,region=="Cen"|region=="none"))+
  geom_violin(mapping = aes(y=mean,x=region))+
  facet_grid(.~CHROM)+
  theme_light()

mean_bp<-ggplot(data=subset(all_paint_maxMean,region=="InterSpInv"|region=="IntraSpInv"|region=="none"))+
  geom_hline(yintercept = .1,color="blue4")+
  geom_hline(yintercept = .5,color="blue4",alpha=.25)+
  geom_hline(yintercept = .9,color="blue4")+
  geom_boxplot(mapping = aes(y=mean,x=region))+
  facet_grid(.~CHROM)+
  theme_light()+
  labs(x=NULL)+
  geom_text(data=subset(Letters.Mean,region=="InterSpInv"|region=="IntraSpInv"|region=="none"),aes(x=region, label=Letter,y=1),size=2, position = position_dodge(width = .75))+
  
  theme(axis.text.x = element_text(angle = 45, hjust =1),
        panel.spacing.x=unit(0, "lines"),text = element_text(size = 10),panel.border = element_rect(color = "black", fill = NA, linewidth = .5),panel.grid.major = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.minor.x = element_blank())


####Tajimas D


TD_bp_N<-ggplot(data=subset(TD_2pops,POP=="MN_QPA"&(region=="InterSpInv"|region=="IntraSpInv"|region=="none")))+
  geom_hline(yintercept = 2,color="blue4")+
  geom_hline(yintercept = 0,color="blue4",alpha=.25)+
  geom_hline(yintercept = -2,color="blue4")+
  geom_boxplot(mapping = aes(y=TajimaD,x=region))+
  facet_grid(.~CHROM)+
  labs(x = NULL,y="Tajima's D - North")+
  theme_light()+
  geom_text(data=subset(Letters.TD_N,region=="InterSpInv"|region=="IntraSpInv"|region=="none"),aes(x=region, label=Letter,y=3.8),size=2, position = position_dodge(width = .75))+
  
  theme(axis.text.x = element_text(angle = 45, hjust =1),strip.background = element_blank(), strip.text.x = element_blank(),
        panel.spacing.x=unit(0, "lines"),text = element_text(size = 10),panel.border = element_rect(color = "black", fill = NA, linewidth = .5),panel.grid.major = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.minor.x = element_blank())


TD_bp_S<-ggplot(data=subset(TD_2pops,POP=="OK_NOW"&(region=="InterSpInv"|region=="IntraSpInv"|region=="none")))+
  geom_hline(yintercept = 2,color="blue4")+
  geom_hline(yintercept = 0,color="blue4",alpha=.25)+
  geom_hline(yintercept = -2,color="blue4")+
  geom_boxplot(mapping = aes(y=TajimaD,x=region))+
  facet_grid(.~CHROM)+
  theme_light()+
  labs(x = "Genome Region",y="Tajima's D - South")+
  geom_text(data=subset(Letters.TD_S,region=="InterSpInv"|region=="IntraSpInv"|region=="none"),aes(x=region, label=Letter,y=3.8),size=2, position = position_dodge(width = .75))+
  
  theme(axis.text.x = element_text(angle = 45, hjust =1),strip.background = element_blank(), strip.text.x = element_blank(),
        panel.spacing.x=unit(0, "lines"),text = element_text(size = 10),panel.border = element_rect(color = "black", fill = NA, linewidth = .5),panel.grid.major = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.minor.x = element_blank())


pdf("BoxPlot_PaintTD_INV.pdf",height = 6, width = 10)
ggarrange(mean_bp, slope_bp,TD_bp_N,TD_bp_S,ncol=1)
dev.off()


###############
###ALL
################





ggplot(data=subset(all_paint_maxSlope,region=="InterSpInv"|region=="IntraSpInv"|region=="none"))+
  geom_violin(mapping = aes(y=abs(slope),x=region))+
  facet_grid(.~CHROM)+
  theme_light()

slope_bp<-ggplot(data=subset(all_paint_maxSlope))+
  geom_hline(yintercept = 0.1,color="blue4")+
  geom_hline(yintercept = .5,color="blue4",alpha=.25)+
  geom_hline(yintercept = .9,color="blue4")+
  geom_boxplot(mapping = aes(y=abs(slope),x=region))+
  theme_light()+
  labs(x=NULL,y= "lat v. mean")+
  theme(axis.text.x = element_text(angle = 45, hjust =1),strip.background = element_blank(), strip.text.x = element_blank(),
        panel.spacing.x=unit(0, "lines"),text = element_text(size = 12),panel.border = element_rect(color = "black", fill = NA, linewidth = .5),panel.grid.major = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.minor.x = element_blank())




ggplot(data=subset(all_paint_maxMean,region=="Cen"|region=="none"))+
  geom_violin(mapping = aes(y=mean,x=region))+
  facet_grid(.~CHROM)+
  theme_light()

mean_bp<-ggplot(data=subset(all_paint_maxMean))+
  geom_hline(yintercept = .1,color="blue4")+
  geom_hline(yintercept = .5,color="blue4",alpha=.25)+
  geom_hline(yintercept = .9,color="blue4")+
  geom_boxplot(mapping = aes(y=mean,x=region))+
  theme_light()+
  labs(x=NULL)+
  theme(axis.text.x = element_text(angle = 45, hjust =1),
        panel.spacing.x=unit(0, "lines"),text = element_text(size = 12),panel.border = element_rect(color = "black", fill = NA, linewidth = .5),panel.grid.major = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.minor.x = element_blank())


####Tajimas D


TD_bp_N<-ggplot(data=subset(TD_2pops,POP=="MN_QPA"))+
  geom_hline(yintercept = 2,color="blue4")+
  geom_hline(yintercept = 0,color="blue4",alpha=.25)+
  geom_hline(yintercept = -2,color="blue4")+
  geom_boxplot(mapping = aes(y=TajimaD,x=region))+
  labs(x = NULL,y="Tajima's D - South")+
  theme_light()+
  theme(axis.text.x = element_text(angle = 45, hjust =1),strip.background = element_blank(), strip.text.x = element_blank(),
        panel.spacing.x=unit(0, "lines"),text = element_text(size = 12),panel.border = element_rect(color = "black", fill = NA, linewidth = .5),panel.grid.major = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.minor.x = element_blank())


TD_bp_S<-ggplot(data=subset(TD_2pops,POP=="OK_NOW"))+
  geom_hline(yintercept = 2,color="blue4")+
  geom_hline(yintercept = 0,color="blue4",alpha=.25)+
  geom_hline(yintercept = -2,color="blue4")+
  geom_boxplot(mapping = aes(y=TajimaD,x=region))+
  theme_light()+
  labs(x = "Genome Region",y="Tajima's D - North")+
  theme(axis.text.x = element_text(angle = 45, hjust =1),strip.background = element_blank(), strip.text.x = element_blank(),
        panel.spacing.x=unit(0, "lines"),text = element_text(size = 12),panel.border = element_rect(color = "black", fill = NA, linewidth = .5),panel.grid.major = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.minor.x = element_blank())


pdf("BoxPlot_PaintTD_ALL.pdf",height = 7, width = 3)
ggarrange(mean_bp, slope_bp,TD_bp_N,TD_bp_S,ncol=1)
dev.off()



###########################
###TD vs. Max Slope#########
##############################
TD_2pops$Start1E6<-round(TD_2pops$BIN_START/100000,0)*100000
TD_2pops_mean<-aggregate(TD_2pops$TajimaD,by=list(TD_2pops$CHROM,TD_2pops$Start1E6,TD_2pops$POP),FUN=median,na.rm=TRUE)
TD_2pops_mean$loc<-paste(TD_2pops_mean$Group.1,TD_2pops_mean$Group.2,sep="_")

all_paint_maxSlope$loc<-paste(all_paint_maxSlope$scaffold,all_paint_maxSlope$window,sep="_")
TD_maxSlope<-merge(TD_2pops_mean,all_paint_maxSlope,by="loc",all = TRUE)

ggplot(TD_maxSlope,mapping=aes(x=x,y=abs(slope),color=Group.3))+
  geom_point(size=.1,alpha=.3)+
  geom_smooth(method = "lm")+
  facet_grid(Group.3~.)+
  theme_light()



all_paint_maxMean$loc<-paste(all_paint_maxMean$scaffold,all_paint_maxMean$window,sep="_")
TD_maxMean<-merge(TD_2pops_mean,all_paint_maxMean,by="loc",all = TRUE)

ggplot(TD_maxMean,mapping=aes(x=x,y=mean,color=Group.3))+
  geom_point(size=.1,alpha=.3)+
  geom_smooth(method = "lm")+
  facet_grid(Group.3~.)+
  theme_light()
