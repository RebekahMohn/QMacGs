library(dplyr)
#library(reshape2)
library(tidyr)
library(ggpubr)
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(adegenet)
library(data.table)
library(readxl)
library(stringr)
library(scatterpie)
#install.packages("sf")
library(sf)

setwd("C:/Users/rmohn/Desktop/")

mac_range<- st_read("04_Mapping/quercus_shapefiles/quermacr.shp")
states <- map_data("state")
counties <- map_data('county')
Canada<-map_data("world", "Canada")
meta<-read_excel("00_Scripts_and_Labels/00_METADATA/meta_pops.xlsx",sheet="Sheet1")

simp_meta<-data.frame(cbind(Seq=meta$Seq,colnum=meta$collectionNumber,species=meta$GenomicVoucherID_USE_THIS,state=meta$State,site=meta$Pop,lat=meta$latitude.orig,long=meta$longitude.orig))

meta_pops<-read_excel("00_Scripts_and_Labels/00_METADATA/metadata_pops.xlsx",sheet="metadata_pops")


ggplot()+
  
  #geom_polygon(data = counties, aes(x = long, y = lat, group = group),fill = 'white', color = "gray", lwd = 0.1)+
  geom_polygon(data = states,aes(x = long, y = lat, group = group),
               fill = "white", color = "gray85")+
  geom_polygon(data = Canada,aes(x = long, y = lat, group = group),
               fill = "white", color = "gray85")+
  geom_sf(data=mac_range,fill="gray", color = "black",alpha=.2)+
  coord_sf(xlim = c(-103.5,-68),ylim = c(29,52))+
  geom_point(data=meta_pops,mapping = aes(x=longitude.orig,y=latitude.orig),color="gray65")+
  theme_light()+
  theme(panel.background=element_rect(fill="lightblue1"))

#############
### Hybrid Maps
##############

sp_colors<- c(alb="#888888",bic="#88CCEE",lob="#CC6677",lyr="#117733",mac="#661100",mic="#DDCC77",mue="#44AA99",mon="#332288",pri="#AA2299",sin="#999933",ste="#000000")

#alba x mac

alb_range<- st_read("04_Mapping/quercus_shapefiles/queralba.shp")
mac_range<- st_read("04_Mapping/quercus_shapefiles/quermacr.shp")
states <- map_data("state")
counties <- map_data('county')
Canada<-map_data("world", "Canada")
meta<-read_excel("00_Scripts_and_Labels/00_METADATA/meta_pops.xlsx",sheet="Sheet1")

simp_meta<-data.frame(cbind(Seq=meta$Seq,colnum=meta$collectionNumber,species=meta$GenomicVoucherID_USE_THIS,state=meta$State,site=meta$Pop,lat=meta$latitude.orig,long=meta$longitude.orig))

meta_pops<-read_excel("00_Scripts_and_Labels/00_METADATA/metadata_pops.xlsx",sheet="metadata_pops")


ggplot()+
  
  #geom_polygon(data = counties, aes(x = long, y = lat, group = group),fill = 'white', color = "gray", lwd = 0.1)+
  geom_polygon(data = states,aes(x = long, y = lat, group = group),
               fill = "white", color = "gray85")+
  geom_polygon(data = Canada,aes(x = long, y = lat, group = group),
               fill = "white", color = "gray85")+
  geom_sf(data=mac_range,aes(fill="mac", color = "mac"),alpha=.1)+
  geom_sf(data=alb_range,aes(fill="alb", color = "alb"),alpha=.2)+
  coord_sf(xlim = c(-103.5,-68),ylim = c(29,52))+
  geom_point(data=meta_pops,mapping = aes(x=longitude.orig,y=latitude.orig),color="gray65")+
  theme_light()+
  scale_color_manual(values=sp_colors)+
  scale_fill_manual(values=sp_colors)+
  theme(panel.background=element_rect(fill="gray90"))


#alb x mic
alb_range<- st_read("04_Mapping/quercus_shapefiles/queralba.shp")
mic_range<- st_read("04_Mapping/quercus_shapefiles/quermich.shp")
states <- map_data("state")
counties <- map_data('county')
Canada<-map_data("world", "Canada")
meta<-read_excel("00_Scripts_and_Labels/00_METADATA/meta_pops.xlsx",sheet="Sheet1")

simp_meta<-data.frame(cbind(Seq=meta$Seq,colnum=meta$collectionNumber,species=meta$GenomicVoucherID_USE_THIS,state=meta$State,site=meta$Pop,lat=meta$latitude.orig,long=meta$longitude.orig))

meta_pops<-read_excel("00_Scripts_and_Labels/00_METADATA/metadata_pops.xlsx",sheet="metadata_pops")


ggplot()+
  
  #geom_polygon(data = counties, aes(x = long, y = lat, group = group),fill = 'white', color = "gray", lwd = 0.1)+
  geom_polygon(data = states,aes(x = long, y = lat, group = group),
               fill = "white", color = "gray85")+
  geom_polygon(data = Canada,aes(x = long, y = lat, group = group),
               fill = "white", color = "gray85")+
  geom_sf(data=mic_range,aes(fill="mic", color = "mic"),alpha=.3)+
  geom_sf(data=alb_range,aes(fill="alb", color = "alb"),alpha=.2)+
  coord_sf(xlim = c(-103.5,-68),ylim = c(29,52))+
  geom_point(data=meta_pops,mapping = aes(x=longitude.orig,y=latitude.orig),color="gray65")+
  theme_light()+
  scale_color_manual(values=sp_colors)+
  scale_fill_manual(values=sp_colors)+
  theme(panel.background=element_rect(fill="gray90"))

#alb x mue
alb_range<- st_read("04_Mapping/quercus_shapefiles/queralba.shp")
mue_range<- st_read("04_Mapping/quercus_shapefiles/quermueh.shp")
states <- map_data("state")
counties <- map_data('county')
Canada<-map_data("world", "Canada")
meta<-read_excel("00_Scripts_and_Labels/00_METADATA/meta_pops.xlsx",sheet="Sheet1")

simp_meta<-data.frame(cbind(Seq=meta$Seq,colnum=meta$collectionNumber,species=meta$GenomicVoucherID_USE_THIS,state=meta$State,site=meta$Pop,lat=meta$latitude.orig,long=meta$longitude.orig))

meta_pops<-read_excel("00_Scripts_and_Labels/00_METADATA/metadata_pops.xlsx",sheet="metadata_pops")


ggplot()+
  
  #geom_polygon(data = counties, aes(x = long, y = lat, group = group),fill = 'white', color = "gray", lwd = 0.1)+
  geom_polygon(data = states,aes(x = long, y = lat, group = group),
               fill = "white", color = "gray85")+
  geom_polygon(data = Canada,aes(x = long, y = lat, group = group),
               fill = "white", color = "gray85")+
  geom_sf(data=mue_range,aes(fill="mue", color = "mue"),alpha=.2)+
  geom_sf(data=alb_range,aes(fill="alb", color = "alb"),alpha=.2)+
  coord_sf(xlim = c(-103.5,-68),ylim = c(29,52))+
  geom_point(data=meta_pops,mapping = aes(x=longitude.orig,y=latitude.orig),color="gray65")+
  theme_light()+
  scale_color_manual(values=sp_colors)+
  scale_fill_manual(values=sp_colors)+
  theme(panel.background=element_rect(fill="gray90"))

#bic x mac
bic_range<- st_read("04_Mapping/quercus_shapefiles/querbico.shp")
mac_range<- st_read("04_Mapping/quercus_shapefiles/quermacr.shp")
states <- map_data("state")
counties <- map_data('county')
Canada<-map_data("world", "Canada")
meta<-read_excel("00_Scripts_and_Labels/00_METADATA/meta_pops.xlsx",sheet="Sheet1")

simp_meta<-data.frame(cbind(Seq=meta$Seq,colnum=meta$collectionNumber,species=meta$GenomicVoucherID_USE_THIS,state=meta$State,site=meta$Pop,lat=meta$latitude.orig,long=meta$longitude.orig))

meta_pops<-read_excel("00_Scripts_and_Labels/00_METADATA/metadata_pops.xlsx",sheet="metadata_pops")


ggplot()+
  
  #geom_polygon(data = counties, aes(x = long, y = lat, group = group),fill = 'white', color = "gray", lwd = 0.1)+
  geom_polygon(data = states,aes(x = long, y = lat, group = group),
               fill = "white", color = "gray85")+
  geom_polygon(data = Canada,aes(x = long, y = lat, group = group),
               fill = "white", color = "gray85")+
  geom_sf(data=mac_range,aes(fill="mac", color = "mac"),alpha=.1)+
  geom_sf(data=bic_range,aes(fill="bic", color = "bic"),alpha=.2)+
  coord_sf(xlim = c(-103.5,-68),ylim = c(29,52))+
  geom_point(data=meta_pops,mapping = aes(x=longitude.orig,y=latitude.orig),color="gray65")+
  theme_light()+
  scale_color_manual(values=sp_colors)+
  scale_fill_manual(values=sp_colors)+
  theme(panel.background=element_rect(fill="gray90"))

#lyr x mac
lyr_range<- st_read("04_Mapping/quercus_shapefiles/querlyra.shp")
mac_range<- st_read("04_Mapping/quercus_shapefiles/quermacr.shp")
states <- map_data("state")
counties <- map_data('county')
Canada<-map_data("world", "Canada")
meta<-read_excel("00_Scripts_and_Labels/00_METADATA/meta_pops.xlsx",sheet="Sheet1")

simp_meta<-data.frame(cbind(Seq=meta$Seq,colnum=meta$collectionNumber,species=meta$GenomicVoucherID_USE_THIS,state=meta$State,site=meta$Pop,lat=meta$latitude.orig,long=meta$longitude.orig))

meta_pops<-read_excel("00_Scripts_and_Labels/00_METADATA/metadata_pops.xlsx",sheet="metadata_pops")


ggplot()+
  
  #geom_polygon(data = counties, aes(x = long, y = lat, group = group),fill = 'white', color = "gray", lwd = 0.1)+
  geom_polygon(data = states,aes(x = long, y = lat, group = group),
               fill = "white", color = "gray85")+
  geom_polygon(data = Canada,aes(x = long, y = lat, group = group),
               fill = "white", color = "gray85")+
  geom_sf(data=mac_range,aes(fill="mac", color = "mac"),alpha=.1)+
  geom_sf(data=lyr_range,aes(fill="lyr", color = "lyr"),alpha=.2)+
  coord_sf(xlim = c(-103.5,-68),ylim = c(29,52))+
  geom_point(data=meta_pops,mapping = aes(x=longitude.orig,y=latitude.orig),color="gray65")+
  theme_light()+
  scale_color_manual(values=sp_colors)+
  scale_fill_manual(values=sp_colors)+
  theme(panel.background=element_rect(fill="gray90"))

#lyr x mic
lyr_range<- st_read("04_Mapping/quercus_shapefiles/querlyra.shp")
mic_range<- st_read("04_Mapping/quercus_shapefiles/quermich.shp")
states <- map_data("state")
counties <- map_data('county')
Canada<-map_data("world", "Canada")
meta<-read_excel("00_Scripts_and_Labels/00_METADATA/meta_pops.xlsx",sheet="Sheet1")

simp_meta<-data.frame(cbind(Seq=meta$Seq,colnum=meta$collectionNumber,species=meta$GenomicVoucherID_USE_THIS,state=meta$State,site=meta$Pop,lat=meta$latitude.orig,long=meta$longitude.orig))

meta_pops<-read_excel("00_Scripts_and_Labels/00_METADATA/metadata_pops.xlsx",sheet="metadata_pops")


ggplot()+
  
  #geom_polygon(data = counties, aes(x = long, y = lat, group = group),fill = 'white', color = "gray", lwd = 0.1)+
  geom_polygon(data = states,aes(x = long, y = lat, group = group),
               fill = "white", color = "gray85")+
  geom_polygon(data = Canada,aes(x = long, y = lat, group = group),
               fill = "white", color = "gray85")+
  geom_sf(data=mic_range,aes(fill="mic", color = "mic"),alpha=.3)+
  geom_sf(data=lyr_range,aes(fill="lyr", color = "lyr"),alpha=.2)+
  coord_sf(xlim = c(-103.5,-68),ylim = c(29,52))+
  geom_point(data=meta_pops,mapping = aes(x=longitude.orig,y=latitude.orig),color="gray65")+
  theme_light()+
  scale_color_manual(values=sp_colors)+
  scale_fill_manual(values=sp_colors)+
  theme(panel.background=element_rect(fill="gray90"))

#mac x mue
mue_range<- st_read("04_Mapping/quercus_shapefiles/quermueh.shp")
mac_range<- st_read("04_Mapping/quercus_shapefiles/quermacr.shp")
states <- map_data("state")
counties <- map_data('county')
Canada<-map_data("world", "Canada")
meta<-read_excel("00_Scripts_and_Labels/00_METADATA/meta_pops.xlsx",sheet="Sheet1")

simp_meta<-data.frame(cbind(Seq=meta$Seq,colnum=meta$collectionNumber,species=meta$GenomicVoucherID_USE_THIS,state=meta$State,site=meta$Pop,lat=meta$latitude.orig,long=meta$longitude.orig))

meta_pops<-read_excel("00_Scripts_and_Labels/00_METADATA/metadata_pops.xlsx",sheet="metadata_pops")


ggplot()+
  
  #geom_polygon(data = counties, aes(x = long, y = lat, group = group),fill = 'white', color = "gray", lwd = 0.1)+
  geom_polygon(data = states,aes(x = long, y = lat, group = group),
               fill = "white", color = "gray85")+
  geom_polygon(data = Canada,aes(x = long, y = lat, group = group),
               fill = "white", color = "gray85")+
  geom_sf(data=mac_range,aes(fill="mac", color = "mac"),alpha=.1)+
  geom_sf(data=mue_range,aes(fill="mue", color = "mue"),alpha=.2)+
  coord_sf(xlim = c(-103.5,-68),ylim = c(29,52))+
  geom_point(data=meta_pops,mapping = aes(x=longitude.orig,y=latitude.orig),color="gray65")+
  theme_light()+
  scale_color_manual(values=sp_colors)+
  scale_fill_manual(values=sp_colors)+
  theme(panel.background=element_rect(fill="gray90"))
