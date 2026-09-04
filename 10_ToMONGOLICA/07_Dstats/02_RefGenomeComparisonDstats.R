#Reference genome and introgression comparison

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
Dsp_SigMong<-aggregate(as.numeric(g10)~chr+windowmid,data=Dsp_aggMong,FUN=sum)

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
samps_edge_dist<-read.delim("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/dist_range_margin_centroid2.txt",sep="\t")

Dlob_Mong_All_countsigwins<-merge(Dlob_Mong_All_countsigwins,Oldpops,by.x="pop",by.y="OLD_Pop")
Dlob_Mong_All_win_edge<-merge(Dlob_Mong_All_countsigwins,samps_edge_dist,by.x="NEW_Pop",by.y="Pop",all.x=T)
#do regression
lob_mong_distcent<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dlob_Mong_All_win_edge$count, x=Dlob_Mong_All_win_edge$dist_fromCenter_EM/100000))+
  geom_point(mapping=aes(y = Dlob_Mong_All_win_edge$count, x=Dlob_Mong_All_win_edge$dist_fromCenter_EM/100000, color=Dlob_Mong_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(color="latitude",x="distance from centroid",y="windows",subtitle="Q. lobata")

lob_mong_distedge<-ggplot()+
  geom_smooth(method=lm, mapping=aes(y = Dlob_Mong_All_win_edge$count, x=Dlob_Mong_All_win_edge$dist_fromEdge/100000))+
  geom_jitter(mapping=aes(y = Dlob_Mong_All_win_edge$count, x=Dlob_Mong_All_win_edge$dist_fromEdge/100000,  color=Dlob_Mong_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(color="latitude",x="distance from edge",y="windows",subtitle="Q. lobata")

lob_mong_lat<-ggplot()+
  geom_smooth(method=lm, mapping=aes(y = Dlob_Mong_All_win_edge$count, x=Dlob_Mong_All_win_edge$lat_fixed))+
  geom_jitter(mapping=aes(y = Dlob_Mong_All_win_edge$count, x=Dlob_Mong_All_win_edge$lat_fixed,  color=Dlob_Mong_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(color="latitude",x="latitude",y="windows",subtitle="Q. lobata")

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
  geom_point(aes(x=PC1, y=PC2,color=lat_fixed,size=long_fixed),alpha=0.75)+
  scale_color_viridis_c()+
  labs(color="latitude",subtitle = "Q. lobata")

# library(ggplot2)
# ggplot(data=Dlob_Mong_PCs)+
#   geom_point(aes(x=PC1, y=lat_fixed))

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
samps_edge_dist<-read.delim("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/dist_range_margin_centroid2.txt",sep="\t")

Dalb_Mong_All_countsigwins<-merge(Dalb_Mong_All_countsigwins,Oldpops,by.x="pop",by.y="OLD_Pop")
Dalb_Mong_All_win_edge<-merge(Dalb_Mong_All_countsigwins,samps_edge_dist,by.x="NEW_Pop",by.y="Pop",all.x=T)
#do regression
alb_mong_distcent<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dalb_Mong_All_win_edge$count, x=Dalb_Mong_All_win_edge$dist_fromCenter_EM/100000))+
  geom_point(mapping=aes(y = Dalb_Mong_All_win_edge$count, x=Dalb_Mong_All_win_edge$dist_fromCenter_EM/100000, color=Dalb_Mong_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(color="latitude",x="distance from centroid",y="windows",subtitle="Q. alba")

alb_mong_distedge<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dalb_Mong_All_win_edge$count, x=Dalb_Mong_All_win_edge$dist_fromEdge/100000))+
  geom_jitter(mapping=aes(y = Dalb_Mong_All_win_edge$count, x=Dalb_Mong_All_win_edge$dist_fromEdge/100000,  color=Dalb_Mong_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(color="latitude",x="distance from edge",y="windows",subtitle="Q. alba")

alb_mong_lat<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dalb_Mong_All_win_edge$count, x=Dalb_Mong_All_win_edge$lat_fixed))+
  geom_jitter(mapping=aes(y = Dalb_Mong_All_win_edge$count, x=Dalb_Mong_All_win_edge$lat_fixed,  color=Dalb_Mong_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(color="latitude",x="latitude",y="windows",subtitle="Q. alba")

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
  geom_point(aes(x=PC1, y=PC2,color=lat_fixed,size=long_fixed),alpha=0.75)+
  scale_color_viridis_c()+
  labs(color="latitude",subtitle = "Q. alba")

# library(ggplot2)
# ggplot(data=Dalb_Mong_PCs)+
#   geom_point(aes(x=PC1, y=lat_fixed))
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
samps_edge_dist<-read.delim("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/dist_range_margin_centroid2.txt",sep="\t")

Dmue_Mong_All_countsigwins<-merge(Dmue_Mong_All_countsigwins,Oldpops,by.x="pop",by.y="OLD_Pop")
Dmue_Mong_All_win_edge<-merge(Dmue_Mong_All_countsigwins,samps_edge_dist,by.x="NEW_Pop",by.y="Pop",all.x=T)
#do regression
mue_mong_distcent<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dmue_Mong_All_win_edge$count, x=Dmue_Mong_All_win_edge$dist_fromCenter_EM/100000))+
  geom_point(mapping=aes(y = Dmue_Mong_All_win_edge$count, x=Dmue_Mong_All_win_edge$dist_fromCenter_EM/100000, color=Dmue_Mong_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(color="latitude",x="distance from centroid",y="windows",subtitle="Q. muehlenbergii")

mue_mong_distedge<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dmue_Mong_All_win_edge$count, x=Dmue_Mong_All_win_edge$dist_fromEdge/100000))+
  geom_jitter(mapping=aes(y = Dmue_Mong_All_win_edge$count, x=Dmue_Mong_All_win_edge$dist_fromEdge/100000,  color=Dmue_Mong_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(color="latitude",x="distance from edge",y="windows",subtitle="Q. muehlenbergii")

mue_mong_lat<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dmue_Mong_All_win_edge$count, x=Dmue_Mong_All_win_edge$lat_fixed))+
  geom_jitter(mapping=aes(y = Dmue_Mong_All_win_edge$count, x=Dmue_Mong_All_win_edge$lat_fixed,  color=Dmue_Mong_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(color="latitude",x="latitude",y="windows",subtitle="Q. muehlenbergii")

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
  geom_point(aes(x=PC1, y=PC2,color=lat_fixed,size=long_fixed),alpha=0.75)+
  scale_color_viridis_c()+
  labs(color="latitude",subtitle = "Q. muehlenbergii")

# library(ggplot2)
# ggplot(data=Dmue_Mong_PCs)+
#   geom_point(aes(x=PC1, y=lat_fixed))

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
samps_edge_dist<-read.delim("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/dist_range_margin_centroid2.txt",sep="\t")

Dste_Mong_All_countsigwins<-merge(Dste_Mong_All_countsigwins,Oldpops,by.x="pop",by.y="OLD_Pop")
Dste_Mong_All_win_edge<-merge(Dste_Mong_All_countsigwins,samps_edge_dist,by.x="NEW_Pop",by.y="Pop",all.x=T)
#do regression
ste_mong_distcent<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dste_Mong_All_win_edge$count, x=Dste_Mong_All_win_edge$dist_fromCenter_EM/100000))+
  geom_point(mapping=aes(y = Dste_Mong_All_win_edge$count, x=Dste_Mong_All_win_edge$dist_fromCenter_EM/100000, color=Dste_Mong_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(color="latitude",x="distance from centroid",y="windows",subtitle="Q. stellata")

ste_mong_distedge<-ggplot()+
  geom_smooth(method=lm,mapping=aes(y = Dste_Mong_All_win_edge$count, x=Dste_Mong_All_win_edge$dist_fromEdge/100000))+
  geom_jitter(mapping=aes(y = Dste_Mong_All_win_edge$count, x=Dste_Mong_All_win_edge$dist_fromEdge/100000,  color=Dste_Mong_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(color="latitude",x="distance from edge",y="windows",subtitle="Q. stellata")

ste_mong_lat<-ggplot()+
  geom_smooth(method=lm,aes(y = Dste_Mong_All_win_edge$count,x=Dste_Mong_All_win_edge$lat_fixed))+
  geom_jitter(mapping=aes(y = Dste_Mong_All_win_edge$count, x=Dste_Mong_All_win_edge$lat_fixed,  color=Dste_Mong_All_win_edge$lat))+
  scale_color_viridis_b()+
  labs(color="latitude",x="latitude",y="windows",subtitle="Q. stellata")

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
  geom_point(aes(x=PC1, y=PC2,color=lat_fixed,size=long_fixed),alpha=0.75)+
  scale_color_viridis_c()+
  labs(color="latitude",color="latitude",subtitle = "Q. stellata")

# library(ggplot2)
# ggplot(data=Dste_Mong_PCs)+
#   geom_point(aes(x=PC1, y=lat_fixed))


library(ggpubr)
ggarrange(albMongPC1PC2,lobMongPC1PC2,mueMongPC1PC2,steMongPC1PC2,ncol=2,nrow=2,common.legend = T,legend = "right")
ggarrange(alb_mong_distcent,lob_mong_distcent,mue_mong_distcent,ste_mong_distcent,common.legend = T,legend = "right")
ggarrange(alb_mong_distedge,lob_mong_distedge,mue_mong_distedge,ste_mong_distedge,common.legend = T,legend = "right")
ggarrange(alb_mong_lat,lob_mong_lat,mue_mong_lat,ste_mong_lat,common.legend = T,legend = "right")


ggarrange(alb_Mac_distcent,lob_Mac_distcent,mue_Mac_distcent,ste_Mac_distcent,
          alb_mong_distcent,lob_mong_distcent,mue_mong_distcent,ste_mong_distcent,
          alb_Mac_distedge,lob_Mac_distedge,mue_Mac_distedge,ste_Mac_distedge,
          alb_mong_distedge,lob_mong_distedge,mue_mong_distedge,ste_mong_distedge,
          common.legend = T,legend = "right", labels = "auto")

