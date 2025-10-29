#install.packages("terra")
library(terra)
library(readxl)

sampleData<-read_excel("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/meta_pops.xlsx",sheet="Sheet1")
gpsPoints<-cbind(long=sampleData$longitude.orig,lat=sampleData$latitude.orig)
######################
######elevation#######
######################
elevRast<-terra::rast("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/wc2.1_30s_elev/wc2.1_30s_elev.tif")
elevPoints<-cbind(Seq=sampleData$Seq,long=sampleData$longitude.orig,lat=sampleData$latitude.orig,terra::extract(x=elevRast,y=gpsPoints))

library(ggplot2)
ggplot(data = subset(elevPoints,wc2.1_30s_elev < 1000))+
geom_point(mapping=aes(x=lat,y=wc2.1_30s_elev))
  

######################
######soil############
######################



#######################
######vapor############
#######################


######bioclim##########

library(car)

bio1_Rast<-terra::rast("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/wc2.1_30s_bio/wc2.1_30s_bio_1.tif")
bio2_Rast<-terra::rast("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/wc2.1_30s_bio/wc2.1_30s_bio_2.tif")
bio3_Rast<-terra::rast("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/wc2.1_30s_bio/wc2.1_30s_bio_3.tif")
bio4_Rast<-terra::rast("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/wc2.1_30s_bio/wc2.1_30s_bio_4.tif")
bio5_Rast<-terra::rast("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/wc2.1_30s_bio/wc2.1_30s_bio_5.tif")
bio6_Rast<-terra::rast("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/wc2.1_30s_bio/wc2.1_30s_bio_6.tif")
bio7_Rast<-terra::rast("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/wc2.1_30s_bio/wc2.1_30s_bio_7.tif")
bio8_Rast<-terra::rast("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/wc2.1_30s_bio/wc2.1_30s_bio_8.tif")
bio9_Rast<-terra::rast("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/wc2.1_30s_bio/wc2.1_30s_bio_9.tif")
bio10_Rast<-terra::rast("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/wc2.1_30s_bio/wc2.1_30s_bio_10.tif")
bio11_Rast<-terra::rast("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/wc2.1_30s_bio/wc2.1_30s_bio_11.tif")
bio12_Rast<-terra::rast("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/wc2.1_30s_bio/wc2.1_30s_bio_12.tif")
bio13_Rast<-terra::rast("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/wc2.1_30s_bio/wc2.1_30s_bio_13.tif")
bio14_Rast<-terra::rast("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/wc2.1_30s_bio/wc2.1_30s_bio_14.tif")
bio15_Rast<-terra::rast("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/wc2.1_30s_bio/wc2.1_30s_bio_15.tif")
bio16_Rast<-terra::rast("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/wc2.1_30s_bio/wc2.1_30s_bio_16.tif")
bio17_Rast<-terra::rast("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/wc2.1_30s_bio/wc2.1_30s_bio_17.tif")
bio18_Rast<-terra::rast("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/wc2.1_30s_bio/wc2.1_30s_bio_18.tif")
bio19_Rast<-terra::rast("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/wc2.1_30s_bio/wc2.1_30s_bio_19.tif")


bioclim_Points<-cbind(Seq=sampleData$Seq,long=sampleData$longitude.orig,lat=sampleData$latitude.orig,
                  bio1=terra::extract(x=bio1_Rast,y=gpsPoints),bio11=terra::extract(x=bio11_Rast,y=gpsPoints),
                  bio2=terra::extract(x=bio2_Rast,y=gpsPoints),bio12=terra::extract(x=bio12_Rast,y=gpsPoints),
                  bio3=terra::extract(x=bio3_Rast,y=gpsPoints),bio13=terra::extract(x=bio13_Rast,y=gpsPoints),
                  bio4=terra::extract(x=bio4_Rast,y=gpsPoints),bio14=terra::extract(x=bio14_Rast,y=gpsPoints),
                  bio5=terra::extract(x=bio5_Rast,y=gpsPoints),bio15=terra::extract(x=bio15_Rast,y=gpsPoints),
                  bio6=terra::extract(x=bio6_Rast,y=gpsPoints),bio16=terra::extract(x=bio16_Rast,y=gpsPoints),
                  bio7=terra::extract(x=bio7_Rast,y=gpsPoints),bio17=terra::extract(x=bio17_Rast,y=gpsPoints),
                  bio8=terra::extract(x=bio8_Rast,y=gpsPoints),bio18=terra::extract(x=bio18_Rast,y=gpsPoints),
                  bio9=terra::extract(x=bio9_Rast,y=gpsPoints),bio19=terra::extract(x=bio19_Rast,y=gpsPoints),
                  bio10=terra::extract(x=bio10_Rast,y=gpsPoints))


bioclim_Points<-na.omit(bioclim_Points)
cor_matrix<-cor(bioclim_Points[c("wc2.1_30s_bio_1","wc2.1_30s_bio_2","wc2.1_30s_bio_3", "wc2.1_30s_bio_4",
                                 "wc2.1_30s_bio_5","wc2.1_30s_bio_6","wc2.1_30s_bio_7","wc2.1_30s_bio_8",
                                 "wc2.1_30s_bio_9","wc2.1_30s_bio_10","wc2.1_30s_bio_11","wc2.1_30s_bio_12",
                                 "wc2.1_30s_bio_13","wc2.1_30s_bio_14","wc2.1_30s_bio_15","wc2.1_30s_bio_16",
                                 "wc2.1_30s_bio_17","wc2.1_30s_bio_18","wc2.1_30s_bio_19")])

cor_matrix<-cor(bioclim_Points[c("wc2.1_30s_bio_2","wc2.1_30s_bio_3",
                                 "wc2.1_30s_bio_5","wc2.1_30s_bio_6","wc2.1_30s_bio_8",
                                 "wc2.1_30s_bio_9","wc2.1_30s_bio_13","wc2.1_30s_bio_14",
                                "wc2.1_30s_bio_18")])


cor_matrix<-cor(bioclim_Points[c("wc2.1_30s_bio_2","wc2.1_30s_bio_6","wc2.1_30s_bio_8",
                                 "wc2.1_30s_bio_13","wc2.1_30s_bio_14",
                                 "wc2.1_30s_bio_18")])

#image(cor_matrix, main = "Correlation Matrix", col = colorRampPalette(c("blue", "white", "red"))(20))

library(corrplot)
corrplot(cor_matrix, method = 'ellipse', order = 'AOE')

###Only 9 Bioclim of interest
bioclim_ssPoints<-cbind(Seq=sampleData$Seq,long=sampleData$longitude.orig,lat=sampleData$latitude.orig,
                      bio2=terra::extract(x=bio2_Rast,y=gpsPoints),
                      bio3=terra::extract(x=bio3_Rast,y=gpsPoints),
                      bio5=terra::extract(x=bio5_Rast,y=gpsPoints),
                      bio6=terra::extract(x=bio6_Rast,y=gpsPoints),
                      bio8=terra::extract(x=bio8_Rast,y=gpsPoints),
                      bio9=terra::extract(x=bio9_Rast,y=gpsPoints),
                      bio13=terra::extract(x=bio13_Rast,y=gpsPoints),
                      bio14=terra::extract(x=bio14_Rast,y=gpsPoints),
                      bio18=terra::extract(x=bio18_Rast,y=gpsPoints))
                      
bioclim_ssPoints<-bioclim_ssPoints[-c(694:720),]
bioclim_ssPoints$Seq<-paste(bioclim_ssPoints$Seq,".bam",sep="")
bioclim_ssPoints<-bioclim_ssPoints[,-c(2,3)]

samples<-read.delim(file="C:/Users/rmohn/Desktop/10_Analysis/114_GEA/pure_Mac99a.tsv",sep="\t",header=FALSE)
bioclim_9ss<-merge(samples,bioclim_ssPoints,by.x="V1",by.y="Seq",all.x=TRUE)

write.table(bioclim_9ss,file = "C:/Users/rmohn/Desktop/10_Analysis/114_GEA/BioClim9vars.txt", sep="\t",quote=FALSE,row.names=FALSE)


###more strict Bioclim vars

bioclim_6ss<-bioclim_9ss[,-c(3,4,7)]
write.table(bioclim_6ss,file = "C:/Users/rmohn/Desktop/10_Analysis/114_GEA/BioClim6vars.txt", sep="\t",quote=FALSE,row.names=FALSE)



#devtools::install_github("bcm-uga/lfmm")

