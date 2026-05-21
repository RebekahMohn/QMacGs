#SF2 + AI
library(readxl)

all_sf1_clim<-read.delim("C:/Users/rmohn/Documents/GitHub/QMacGs/20_ToMACROCARPA/08_SweepFinder2/285_SF2_Output/all_sf2.pcf.sf2",sep="\t")
all_sf1_clim_grp<-read.delim("C:/Users/rmohn/Documents/GitHub/QMacGs/20_ToMACROCARPA/08_SweepFinder2/285_SF2_Output/all_sf2_grp.pcf.sf2",sep="\t")

aisnps_100<-read_excel("C:/Users/rmohn/Desktop/10_Analysis/120_Ad_Int/Ad_Int_Database_samp100_260205.xlsx")

#
aisnps_100$round<-round(aisnps_100$BP/10000)*10000
aisnps_100$chr_round<-paste(aisnps_100$csome, aisnps_100$round)
aisnps_100_sf2<-merge(aisnps_100,all_sf1_clim_grp, by.x="chr_round",by.y="chrloc",all.x=T)

aisnps_100_sf2_sig<-subset(aisnps_100_sf2,sigDGEA=="sig")
quantile(aisnps_100_sf2$LR,na.rm=TRUE)
quantile(aisnps_100_sf2_sig$LR,na.rm=TRUE)
quantile(all_sf1_clim_grp$LR,na.rm=T)

qaisnps<-t(t(quantile(aisnps_100_sf2_sig$LR,probs=c(0.05,0.25,0.5,0.75,0.95))))

quantsmplr<-function(data,n){
  t(replicate(100,quantile(sample(data,n,replace = F),probs=c(0.05,0.25,0.5,0.75,0.95))))
}

#replicate(100,quantile(sample(all_sf1_clim_grp$LR,33,replace = F),probs=c(0.05,0.25,0.5,0.75,0.95)))

quantsample<-quantsmplr(all_sf1_clim_grp$LR,33)
#quantile(sample(all_sf1_clim_grp$LR,33,replace = F))
ggplot()+
  geom_boxplot(mapping=aes(x="5%",y=quantsample[,1]),color="gray40")+
  geom_boxplot(mapping=aes(x="25%",y=quantsample[,2]),color="gray40")+
  geom_boxplot(mapping=aes(x="50%",y=quantsample[,3]),color="gray40")+
  geom_boxplot(mapping=aes(x="75%",y=quantsample[,4]),color="gray40")+
  #geom_boxplot(mapping=aes(x="95%",y=quantsample[,5]),color="gray40")+
  geom_point(mapping=aes(x="5%",y=quantsample[,1]),color="gray40",alpha=0.1)+
  geom_point(mapping=aes(x="25%",y=quantsample[,2]),color="gray40",alpha=0.1)+
  geom_point(mapping=aes(x="50%",y=quantsample[,3]),color="gray40",alpha=0.1)+
  geom_point(mapping=aes(x="75%",y=quantsample[,4]),color="gray40",alpha=0.1)+
  #geom_point(mapping=aes(x="95%",y=quantsample[,5]),color="gray40",alpha=0.1)+
  
  geom_point(mapping=aes(x=row.names(qaisnps),y=qaisnps[,1]),color="red",size=1)+
  theme_light()

  
## conclusion at 25, 50, and 75 % quantiles the sweepfinder CLR is higher than randoms samples throughout the genome.


ggplot()+
  geom_point(data=subset(aisnps_100_sf2_sig,mue_1!="NA"),mapping=aes(x=wc2.1_30s_bio_6,y=as.numeric(mue_1),color=log(LR)))+
  scale_color_viridis_c()+
  theme_light()
  
ggplot()+
  geom_point(data=subset(aisnps_100_sf2_sig,ste_1!="NA"),mapping=aes(x=wc2.1_30s_bio_6,y=as.numeric(ste_1),color=log(LR))) +
  scale_color_viridis_c()+
  theme_light()

ggplot()+
  geom_point(data=subset(aisnps_100_sf2_sig,lob_1!="NA"),mapping=aes(x=wc2.1_30s_bio_6,y=as.numeric(lob_1),color=log(LR))) +
  scale_color_viridis_c()+
  theme_light()
