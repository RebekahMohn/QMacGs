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
