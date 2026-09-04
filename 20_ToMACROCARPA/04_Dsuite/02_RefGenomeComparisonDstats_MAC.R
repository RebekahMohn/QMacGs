#Reference genome and introgression comparison
library(dplyr)
library(ggplot2)
#################################################
## Macrocarpa Reference ##########################
#################################################

#combining all species dstatistics+metadata into one spreadsheet for all species
Dlob_AllMac<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/107_DSuite_MACREF/DPOPS/bic_mac_lob_2000.txt",sep="\t")
Dlob_AllMac$species<-rep("lob",length(Dlob_AllMac$chr))
Dalb_AllMac<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/107_DSuite_MACREF/DPOPS/bic_mac_alb_2000.txt",sep="\t")
Dalb_AllMac$species<-rep("alb",length(Dalb_AllMac$chr))
Dste_AllMac<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/107_DSuite_MACREF/DPOPS/bic_mac_ste_2000.txt",sep="\t")
Dste_AllMac$species<-rep("ste",length(Dste_AllMac$chr))
Dmue_AllMac<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/107_DSuite_MACREF/DPOPS/bic_mac_mue_2000.txt",sep="\t")
Dmue_AllMac$species<-rep("mue",length(Dmue_AllMac$chr))

D_AllMac<-rbind(Dalb_AllMac,Dlob_AllMac,Dmue_AllMac,Dste_AllMac)
D_AllMac$windowMid<-rowMeans(D_AllMac[,2:3])


##for 500,000 bp windows, saiving the max d_f for each population (this helps with visualizing)
Dsp_All_2000Mac<-data.frame(pop=character(),chr=character(),windowStart=integer(),windowEnd=integer(),D=double(),f_d=double(),f_dM=double(),d_f=double(),lat=double(),long=double(),windowmid=integer(),species=character())
for(j in 1:12){
  i=min(D_AllMac$windowmid)
  while(i <= max(D_AllMac$windowmid)){
    if(j<10){
      chro<-paste("Chr0",j,sep="")
    }else{
      chro<-paste("Chr",j,sep="")
    }
    
    if(i %in% subset(Dalb_AllMac,chr==chr)$windowmid){
      Dwind<-rbind(subset(Dalb_AllMac,windowmid==i&chr==chro),subset(Dlob_AllMac,windowmid==i&chr==chro),subset(Dmue_AllMac,windowmid==i&chr==chro),subset(Dste_AllMac,windowmid==i&chr==chro))
      Dsp_All_2000Mac<-rbind(Dsp_All_2000Mac,(Dwind %>% group_by(pop) %>% top_n(1, d_f)))
    }
    i=i+500000
  }
}


Dsp_All_2000Mac$g10<-Dsp_All_2000Mac$d_f>0.25

Dsp_aggMac<-aggregate(g10~chr+windowmid+pop,data=Dsp_All_2000Mac,FUN=max)
Dsp_SigMac<-aggregate(as.numeric(g10)~chr+windowmid,data=Dsp_aggMac,FUN=sum)

ggplot()+
  geom_histogram(mapping=aes(x=Dsp_SigMac$`as.numeric(g10)`),bins = 10)

length(subset(Dsp_SigMac,`as.numeric(g10)`>=4)$chr)

################################################
####### Mac PCA and dist range edge ###########
###############################################

library(readxl)
Pops<-read_xlsx("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/metadata_pops.xlsx",sheet="metadata_pops")
#Dlob_Mac_All<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_lobata.txt",sep="\t")
#Dlob_Mac_All$species<-rep("lob",length(Dlob_Mac_All$chr))



library(dplyr)
Dlob_Mac_All_win<-Dlob_AllMac %>% group_by(pop,windowmid,chr) %>% top_n(1, d_f)
Dlob_Mac_All_sigwins<-subset(Dlob_Mac_All_win,d_f>0.25)
Dlob_Mac_All_sigwins$count<-1
Dlob_Mac_All_countsigwins<-aggregate(count~pop+lat,data=Dlob_Mac_All_sigwins,FUN=sum)
#merge distance to range edge
#we need to recalculate this
#samps_edge_dist<-read.delim("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/dist_range_margin_centroid.txt",sep="\t")

Dlob_Mac_All_countsigwins<-merge(Dlob_Mac_All_countsigwins,Pops,by.x="pop",by.y="Pop")
Dlob_Mac_All_win_edge<-merge(Dlob_Mac_All_countsigwins,samps_edge_dist,by.x="pop",by.y="Pop",all.x=T)
#do regression
lob_Mac_distcent<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dlob_Mac_All_win_edge$count, x=Dlob_Mac_All_win_edge$dist_fromCenter_EM/100000))+
  geom_point(mapping=aes(y = Dlob_Mac_All_win_edge$count, x=Dlob_Mac_All_win_edge$dist_fromCenter_EM/100000, color=Dlob_Mac_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(color="latitude",x="distance from centroid",y="windows",subtitle="Q. lobata")

lob_Mac_distedge<-ggplot()+
  geom_smooth(method=lm, mapping=aes(y = Dlob_Mac_All_win_edge$count, x=Dlob_Mac_All_win_edge$dist_fromEdge/100000))+
  geom_jitter(mapping=aes(y = Dlob_Mac_All_win_edge$count, x=Dlob_Mac_All_win_edge$dist_fromEdge/100000, color=Dlob_Mac_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(color="latitude",x="distance from edge",y="windows",subtitle="Q. lobata")

lob_Mac_lat<-ggplot()+
  geom_smooth(method=lm, mapping=aes(y = Dlob_Mac_All_win_edge$count, x=Dlob_Mac_All_win_edge$lat_fixed))+
  geom_jitter(mapping=aes(y = Dlob_Mac_All_win_edge$count, x=Dlob_Mac_All_win_edge$lat_fixed, color=Dlob_Mac_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(color="latitude",x="latitude",y="windows",subtitle="Q. lobata")

##PCA
Dlob_Mac_All_win$chr_win<-paste(Dlob_Mac_All_win$chr,Dlob_Mac_All_win$windowmid,sep="_")
Dlob_Mac_All_win_small<-Dlob_Mac_All_win[,-c(1:6,9:14)]

library(tidyr)
Dlob_Mac_All_win_wide<-pivot_wider(Dlob_Mac_All_win_small,names_from=chr_win,values_from=d_f)
Dlob_Mac_All_win_wide<-Dlob_Mac_All_win_wide %>% 
  select(where(~!any(is.na(.))))


Dlob_Mac_PCA<-prcomp(Dlob_Mac_All_win_wide[,-1])
summary(Dlob_Mac_PCA)
Dlob_Mac_PCs<-as.data.frame(Dlob_Mac_PCA$x)
Dlob_Mac_PCs$Pop<-Dlob_Mac_All_win_wide$pop
Dlob_Mac_PCs<-merge(Dlob_Mac_PCs,Pops,by.x="Pop",by.y="Pop")
library(ggplot2)
lobMacPC1PC2<-ggplot(data=Dlob_Mac_PCs)+
  geom_point(aes(x=PC1, y=PC2,color=lat_fixed,size=long_fixed),alpha=0.75)+
  scale_color_viridis_c()+
  labs(color="latitude",subtitle = "Q. lobata")

# library(ggplot2)
# ggplot(data=Dlob_Mac_PCs)+
#   geom_point(aes(x=PC1, y=lat_fixed))

#### Alba #####
library(readxl)


#Dalb_Mac_All<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_alba.txt",sep="\t")
Dalb_AllMac$species<-rep("alb",length(Dalb_AllMac$chr))

library(dplyr)
Dalb_Mac_All_win<-Dalb_AllMac %>% group_by(pop,windowmid,chr) %>% top_n(1, d_f)
Dalb_Mac_All_sigwins<-subset(Dalb_Mac_All_win,d_f>0.25)
Dalb_Mac_All_sigwins$count<-1
Dalb_Mac_All_countsigwins<-aggregate(count~pop+lat,data=Dalb_Mac_All_sigwins,FUN=sum)
#merge distance to range edge
samps_edge_dist<-read.delim("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/dist_range_margin_centroid.txt",sep="\t")

Dalb_Mac_All_countsigwins<-merge(Dalb_Mac_All_countsigwins,Pops,by.x="pop",by.y="Pop")
Dalb_Mac_All_win_edge<-merge(Dalb_Mac_All_countsigwins,samps_edge_dist,by.x="pop",by.y="Pop",all.x=T)
#do regression
alb_Mac_distcent<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dalb_Mac_All_win_edge$count, x=Dalb_Mac_All_win_edge$dist_fromCenter_EM/100000))+
  geom_point(mapping=aes(y = Dalb_Mac_All_win_edge$count, x=Dalb_Mac_All_win_edge$dist_fromCenter_EM/100000, color=Dalb_Mac_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(color="latitude",x="distance from centroid",y="windows",subtitle="Q. alba")

alb_Mac_distedge<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dalb_Mac_All_win_edge$count, x=Dalb_Mac_All_win_edge$dist_fromEdge/100000))+
  geom_jitter(mapping=aes(y = Dalb_Mac_All_win_edge$count, x=Dalb_Mac_All_win_edge$dist_fromEdge/100000, color=Dalb_Mac_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(color="latitude",x="distance from edge",y="windows",subtitle="Q. alba")

alb_Mac_lat<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dalb_Mac_All_win_edge$count, x=Dalb_Mac_All_win_edge$lat_fixed))+
  geom_jitter(mapping=aes(y = Dalb_Mac_All_win_edge$count, x=Dalb_Mac_All_win_edge$lat_fixed, color=Dalb_Mac_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(color="latitude",x="latitude",y="windows",subtitle="Q. alba")

##PCA
Dalb_Mac_All_win$chr_win<-paste(Dalb_Mac_All_win$chr,Dalb_Mac_All_win$windowmid,sep="_")
Dalb_Mac_All_win_small<-Dalb_Mac_All_win[,-c(1:6,9:14)]

library(tidyr)
Dalb_Mac_All_win_wide<-pivot_wider(Dalb_Mac_All_win_small,names_from=chr_win,values_from=d_f)
Dalb_Mac_All_win_wide<-Dalb_Mac_All_win_wide %>% 
  select(where(~!any(is.na(.))))


# Dalb_Mac_PCA<-prcomp(Dalb_Mac_All_win_wide[,-1])
# summary(Dalb_Mac_PCA)
# Dalb_Mac_PCs<-as.data.frame(Dalb_Mac_PCA$x)
# Dalb_Mac_PCs$Pop<-Dalb_Mac_All_win_wide$pop
# Dalb_Mac_PCs<-merge(Dalb_Mac_PCs,Pops,by.x="Pop",by.y="Pop")
# library(ggplot2)
# albMacPC1PC2<-ggplot(data=Dalb_Mac_PCs)+
#   geom_point(aes(x=PC1, y=PC2,color=lat_fixed,size=long_fixed),alpha=0.75)+
#   scale_color_viridis_c()+
#   labs(color="latitude",subtitle = "Q. alba")
# 
# # library(ggplot2)
# # ggplot(data=Dalb_Mac_PCs)+
# #   geom_point(aes(x=PC1, y=lat_fixed))
# # 
### Muehlenbergii ###
library(readxl)

#Dmue_Mac_All<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_muehlenbergii.txt",sep="\t")
Dmue_AllMac$species<-rep("mue",length(Dmue_AllMac$chr))

library(dplyr)
Dmue_Mac_All_win<-Dmue_AllMac %>% group_by(pop,windowmid,chr) %>% top_n(1, d_f)
Dmue_Mac_All_sigwins<-subset(Dmue_Mac_All_win,d_f>0.25)
Dmue_Mac_All_sigwins$count<-1
Dmue_Mac_All_countsigwins<-aggregate(count~pop+lat,data=Dmue_Mac_All_sigwins,FUN=sum)
#merge distance to range edge
samps_edge_dist<-read.delim("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/dist_range_margin_centroid2.txt",sep="\t")

Dmue_Mac_All_countsigwins<-merge(Dmue_Mac_All_countsigwins,Pops,by.x="pop",by.y="Pop")
Dmue_Mac_All_win_edge<-merge(Dmue_Mac_All_countsigwins,samps_edge_dist,by.x="pop",by.y="Pop",all.x=T)
#do regression
mue_Mac_distcent<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dmue_Mac_All_win_edge$count, x=Dmue_Mac_All_win_edge$dist_fromCenter_EM/100000))+
  geom_point(mapping=aes(y = Dmue_Mac_All_win_edge$count, x=Dmue_Mac_All_win_edge$dist_fromCenter_EM/100000, color=Dmue_Mac_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(color="latitude",x="distance from centroid",y="windows",subtitle="Q. muehlenbergii")

mue_Mac_distedge<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dmue_Mac_All_win_edge$count, x=Dmue_Mac_All_win_edge$dist_fromEdge/100000))+
  geom_jitter(mapping=aes(y = Dmue_Mac_All_win_edge$count, x=Dmue_Mac_All_win_edge$dist_fromEdge/100000, color=Dmue_Mac_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(color="latitude",x="distance from edge",y="windows",subtitle="Q. muehlenbergii")

mue_Mac_lat<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dmue_Mac_All_win_edge$count, x=Dmue_Mac_All_win_edge$lat_fixed))+
  geom_jitter(mapping=aes(y = Dmue_Mac_All_win_edge$count, x=Dmue_Mac_All_win_edge$lat_fixed, color=Dmue_Mac_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(color="latitude",x="latitude",y="windows",subtitle="Q. muehlenbergii")

##PCA
Dmue_Mac_All_win$chr_win<-paste(Dmue_Mac_All_win$chr,Dmue_Mac_All_win$windowmid,sep="_")
Dmue_Mac_All_win_small<-Dmue_Mac_All_win[,-c(1:6,9:14)]

library(tidyr)
Dmue_Mac_All_win_wide<-pivot_wider(Dmue_Mac_All_win_small,names_from=chr_win,values_from=d_f)
Dmue_Mac_All_win_wide<-Dmue_Mac_All_win_wide %>% 
  select(where(~!any(is.na(.))))


Dmue_Mac_PCA<-prcomp(Dmue_Mac_All_win_wide[,-1])
summary(Dmue_Mac_PCA)
Dmue_Mac_PCs<-as.data.frame(Dmue_Mac_PCA$x)
Dmue_Mac_PCs$Pop<-Dmue_Mac_All_win_wide$pop
Dmue_Mac_PCs<-merge(Dmue_Mac_PCs,Pops,by.x="Pop",by.y="Pop")
library(ggplot2)
mueMacPC1PC2<-ggplot(data=Dmue_Mac_PCs)+
  geom_point(aes(x=PC1, y=PC2,color=lat_fixed,size=long_fixed),alpha=0.75)+
  scale_color_viridis_c()+
  labs(color="latitude",subtitle = "Q. muehlenbergii")

# library(ggplot2)
# ggplot(data=Dmue_Mac_PCs)+
#   geom_point(aes(x=PC1, y=lat_fixed))

#### Stellata #####
library(readxl)

#Dste_Mac_All<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/07_Dsuite/PopsRegions/bicolor_macrocarpa_stellata.txt",sep="\t")
Dste_AllMac$species<-rep("ste",length(Dste_AllMac$chr))

library(dplyr)
Dste_Mac_All_win<-Dste_AllMac %>% group_by(pop,windowmid,chr) %>% top_n(1, d_f)
Dste_Mac_All_sigwins<-subset(Dste_Mac_All_win,d_f>0.25)
Dste_Mac_All_sigwins$count<-1
Dste_Mac_All_countsigwins<-aggregate(count~pop+lat,data=Dste_Mac_All_sigwins,FUN=sum)
#merge distance to range edge
samps_edge_dist<-read.delim("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/dist_range_margin_centroid2.txt",sep="\t")

Dste_Mac_All_countsigwins<-merge(Dste_Mac_All_countsigwins,Pops,by.x="pop",by.y="Pop")
Dste_Mac_All_win_edge<-merge(Dste_Mac_All_countsigwins,samps_edge_dist,by.x="pop",by.y="Pop",all.x=T)
#do regression
ste_Mac_distcent<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dste_Mac_All_win_edge$count, x=Dste_Mac_All_win_edge$dist_fromCenter_EM/100000))+
  geom_point(mapping=aes(y = Dste_Mac_All_win_edge$count, x=Dste_Mac_All_win_edge$dist_fromCenter_EM/100000, color=Dste_Mac_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(color="latitude",x="distance from centroid",y="windows",subtitle="Q. stellata")

ste_Mac_distedge<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dste_Mac_All_win_edge$count, x=Dste_Mac_All_win_edge$dist_fromEdge/100000))+
  geom_jitter(mapping=aes(y = Dste_Mac_All_win_edge$count, x=Dste_Mac_All_win_edge$dist_fromEdge/100000, color=Dste_Mac_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(color="latitude",x="distance from edge",y="windows",subtitle="Q. stellata")

ste_Mac_lat<-ggplot()+
  geom_smooth(method=lm,aes(y = Dste_Mac_All_win_edge$count,x=Dste_Mac_All_win_edge$lat_fixed))+
  geom_jitter(mapping=aes(y = Dste_Mac_All_win_edge$count, x=Dste_Mac_All_win_edge$lat_fixed, color=Dste_Mac_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(color="latitude",x="latitude",y="windows",subtitle="Q. stellata")

##PCA
Dste_Mac_All_win$chr_win<-paste(Dste_Mac_All_win$chr,Dste_Mac_All_win$windowmid,sep="_")
Dste_Mac_All_win_small<-Dste_Mac_All_win[,-c(1:6,9:14)]

library(tidyr)
Dste_Mac_All_win_wide<-pivot_wider(Dste_Mac_All_win_small,names_from=chr_win,values_from=d_f)
Dste_Mac_All_win_wide<-Dste_Mac_All_win_wide %>% 
  select(where(~!any(is.na(.))))


Dste_Mac_PCA<-prcomp(Dste_Mac_All_win_wide[,-1])
summary(Dste_Mac_PCA)
Dste_Mac_PCs<-as.data.frame(Dste_Mac_PCA$x)
Dste_Mac_PCs$Pop<-Dste_Mac_All_win_wide$pop
Dste_Mac_PCs<-merge(Dste_Mac_PCs,Pops,by.x="Pop",by.y="Pop")
library(ggplot2)
steMacPC1PC2<-ggplot(data=Dste_Mac_PCs)+
  geom_point(aes(x=PC1, y=PC2,color=lat_fixed,size=long_fixed),alpha=0.75)+
  scale_color_viridis_c()+
  labs(color="latitude",subtitle = "Q. stellata")

# library(ggplot2)
# ggplot(data=Dste_Mac_PCs)+
#   geom_point(aes(x=PC1, y=lat_fixed))


library(ggpubr)
ggarrange(albMacPC1PC2,lobMacPC1PC2,mueMacPC1PC2,steMacPC1PC2,ncol=2,nrow=2,common.legend = T,legend = "right")
ggarrange(alb_Mac_distcent,lob_Mac_distcent,mue_Mac_distcent,ste_Mac_distcent,common.legend = T,legend = "right")
ggarrange(alb_Mac_distedge,lob_Mac_distedge,mue_Mac_distedge,ste_Mac_distedge,common.legend = T,legend = "right")
ggarrange(alb_Mac_lat,lob_Mac_lat,mue_Mac_lat,ste_Mac_lat,common.legend = T,legend = "right")
