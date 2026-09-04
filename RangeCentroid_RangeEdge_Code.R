#Set working directory
setwd("C:/Users/rmohn/Desktop/04_Mapping/shapefiles/quermacr/")

#Load in libraries
library(ggplot2)
library(terra)
library(sf)
library(httr)
library(readr)
#install.packages("rnaturalearth")
library(rnaturalearth)

#Additional packages that may be needed?
# library(readxl)
# library(dplyr)
# library(jtools)
# library(patchwork)
# library(ggeffects)
# library(here) 
# library(dplyr)
# library(phyloseq)
# library(tidyverse)
# library(stringr)



### ATTEMPT AT QUMA RANGE MODELING ########################

#Shapefile downloaded from: https://github.com/wpetry/USTreeAtlas/blob/main/shp/quermacr/quermacr.shp

# Load shapefile into R
bur_oak_range <- st_read("quermacr.shp")

# Plot the shapefile 
plot(st_geometry(bur_oak_range))

# Define coord ref system
if (is.na(st_crs(bur_oak_range))) {
  st_crs(bur_oak_range) <- 4326  # Assign CRS: WGS84 (EPSG: 4326)
}

# Get North American boundaries
nam <- ne_countries(continent = "North America", 
                    country = c("Canada", "United States of America", "Mexico"), 
                    returnclass = "sf")

# Transform the Quercus macrocarpa shapefile to match the CRS of the country boundaries
bur_oak_range <- st_transform(bur_oak_range, st_crs(nam))

# Create a plot with Quercus macrocarpa range and U.S. boundaries
#EM samples - all samples, all quercus species
ggplot() +
  geom_sf(data = bur_oak_range, fill = "green", color = "darkgreen") +  # oak range
  geom_sf(data = nam, fill = NA, color = "black") +  # North America 
  theme_minimal() +  
  labs(title = "EM Samples") + 
  coord_sf(xlim = c(-130, -60), ylim = c(25, 55)) #+ # Set long / lat range
#  geom_point(data = emWider_rootsITS, aes(x = long_fixed, y = lat_fixed), #add sampling points
 #            size = 1.5) 

#EM samples - QUMA only
#ggplot() +
  # geom_sf(data = bur_oak_range, fill = "green", color = "darkgreen") +  # oak range
  # geom_sf(data = nam, fill = NA, color = "black") + # North America 
  # theme_minimal() + 
  # labs(title = "EM Samples - QUMA only") + 
  # coord_sf(xlim = c(-130, -60), ylim = c(25, 55)) + # Set long / lat range
  # geom_point(data = emWider_rootsITS_QUMA, aes(x = long_fixed, y = lat_fixed), #add QUMA-only sampling points
  #            size = 1.5) 



#Calculating range centroid
#Code help from: https://gis.stackexchange.com/questions/43543/how-to-calculate-polygon-centroids-in-r-for-non-contiguous-shapes
#And help from chatGPT

# Transform the shapefile to a projected CRS (e.g., UTM)
bur_oak_range_utm <- st_transform(bur_oak_range, 32617)  # UTM zone 17N for example

# Calculate the centroids for each polygon in the range
centroids <- st_centroid(bur_oak_range_utm)

# Calculate the area of each polygon (these will be used as weights)
areas <- st_area(bur_oak_range_utm)

# Compute the weighted centroid (center-of-mass) by calculating the weighted average of the centroids
weighted_centroid <- centroids %>%
  st_cast("POINT") %>%
  st_as_sf() %>%
  st_set_geometry(NULL) %>%
  mutate(weight = as.numeric(areas)) %>%
  summarise(
    x = sum(weight * st_coordinates(centroids)[, 1]) / sum(weight),
    y = sum(weight * st_coordinates(centroids)[, 2]) / sum(weight))

# Convert weighted centroid back to sf for plotting or further analysis
weighted_centroid_sf <- st_sfc(st_point(c(weighted_centroid$x, weighted_centroid$y)), crs = st_crs(bur_oak_range_utm))

#EM samples with single centroid - all samples, all quercus
ggplot() +
  geom_sf(data = bur_oak_range, fill = "green", color = "darkgreen") +  # bur oak range
  geom_sf(data = nam, fill = NA, color = "black") +  # North America
  theme_minimal() + 
  labs(title = "EM Samples") + 
  geom_sf(data = weighted_centroid_sf, color = 'red', size = 3) + #add colored centroid
#  geom_point(data = emWider_rootsITS, aes(x = long_fixed, y = lat_fixed), #add sampling points
#             size = 1.5) +
  coord_sf(xlim = c(-130, -60), ylim = c(25, 55))  #Set lat, long range

#EM samples with single centroid - QUMA only
# ggplot() +
#   geom_sf(data = bur_oak_range, fill = "green", color = "darkgreen") +  # bur oak range
#   geom_sf(data = nam, fill = NA, color = "black") +  # North America
#   theme_minimal() +  
#   labs(title = "EM Samples - QUMA Only") + 
#   geom_sf(data = weighted_centroid_sf, color = 'red', size = 3) + #add colored centroid
#   geom_point(data = emWider_rootsITS_QUMA, aes(x = long_fixed, y = lat_fixed), #add sampling points
#              size = 1.5) +
#   coord_sf(xlim = c(-130, -60), ylim = c(25, 55))  #Set lat, long range
# 


#Calculating distance from centroids and range edge
# Convert the sampling points into an sf object and transform to the same CRS as the range
library(readxl)
samps<-subset(read_excel("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/meta_pops_updates05072026.xlsx",sheet="Sheet1"),!is.na(long_fixed))
sampling_points_EM <- st_as_sf(samps, coords = c("long_fixed", "lat_fixed"), crs = 4326)

# Transform sampling points to UTM (or the same CRS as the range)
sampling_points_utm_EM <- st_transform(sampling_points_EM, st_crs(bur_oak_range_utm))

# Calculate the distance from each sampling point to the weighted centroid
dist_fromCenter_EM <- st_distance(sampling_points_utm_EM, weighted_centroid_sf)


### Calculate the distance from each sampling point to the nearest range edge
# Get the boundary of the range (polygon)
range_boundary <- st_boundary(bur_oak_range_utm)

# Calculate the distance matrix from each sampling point to all components of the boundary
dist_matrix_EM <- st_distance(sampling_points_utm_EM, range_boundary)

# Apply the 'min' function across the columns (to find the nearest edge for each sampling point)
dist_fromEdge_EM <- apply(dist_matrix_EM, 1, min)

# Add the distances to the dataframe
sampling_points_EM$dist_fromCenter_EM <- as.numeric(dist_fromCenter_EM)  # Convert distance from centroid to numeric
sampling_points_EM$dist_fromEdge <- as.numeric(dist_fromEdge_EM)  # Convert distance from centroid to numeric


#colnames(sampling_points_EM[,3700:3769])

#colored by distance from edge - EM
ggplot() +
  geom_sf(data = bur_oak_range_utm, fill = 'white', color = 'black') +  # Range boundary
  geom_sf(data = weighted_centroid_sf, color = 'red', size = 3) +  # Centroid in red
  geom_sf(data = sampling_points_EM, aes(color = dist_fromEdge_EM), size = 2) +  # Sampling points colored by distance
  scale_color_viridis_c() +  # Color scale for distances
  theme_minimal() +
  labs(title = "EM sites - showing distance from edge")



#colored by distance from centroid - EM
ggplot() +
  geom_sf(data = bur_oak_range_utm, fill = 'white', color = 'black') +  # Range boundary
  geom_sf(data = weighted_centroid_sf, color = 'red', size = 3) +  # Centroid in red
  geom_sf(data = sampling_points_EM, aes(color = dist_fromCenter_EM), size = 2) +  # Sampling points colored by distance
  scale_color_viridis_c() +  # Color scale for distances
  theme_minimal() +
  labs(title = "EM sites - showing distance from centroid")

write.table(sampling_points_EM,"C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/dist_range_margin_centroid2.txt",sep="\t",quote = F,row.names = F)



library(readxl)
popls<-subset(read_excel("C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/metadata_pops.xlsx",sheet="metadata_pops"),!is.na(long_fixed))
pop_points_EM <- st_as_sf(popls, coords = c("long_fixed", "lat_fixed"), crs = 4326)

# Transform pop points to UTM (or the same CRS as the range)
pop_points_utm_EM <- st_transform(pop_points_EM, st_crs(bur_oak_range_utm))

# Calculate the distance from each pop point to the weighted centroid
dist_fromCenter_EM <- st_distance(pop_points_utm_EM, weighted_centroid_sf)


### Calculate the distance from each pop point to the nearest range edge
# Get the boundary of the range (polygon)
range_boundary <- st_boundary(bur_oak_range_utm)

# Calculate the distance matrix from each pop point to all components of the boundary
dist_matrix_EM <- st_distance(pop_points_utm_EM, range_boundary)

# Apply the 'min' function across the columns (to find the nearest edge for each pop point)
dist_fromEdge_EM <- apply(dist_matrix_EM, 1, min)

# Add the distances to the dataframe
pop_points_EM$dist_fromCenter_EM <- as.numeric(dist_fromCenter_EM)  # Convert distance from centroid to numeric
pop_points_EM$dist_fromEdge <- as.numeric(dist_fromEdge_EM)  # Convert distance from centroid to numeric


#colnames(pop_points_EM[,3700:3769])

#colored by distance from edge - EM
ggplot() +
  geom_sf(data = bur_oak_range_utm, fill = 'white', color = 'black') +  # Range boundary
  geom_sf(data = weighted_centroid_sf, color = 'red', size = 3) +  # Centroid in red
  geom_sf(data = pop_points_EM, aes(color = dist_fromEdge_EM), size = 2) +  # pop points colored by distance
  scale_color_viridis_c() +  # Color scale for distances
  theme_minimal() +
  labs(title = "EM sites - showing distance from edge")



#colored by distance from centroid - EM
ggplot() +
  geom_sf(data = bur_oak_range_utm, fill = 'white', color = 'black') +  # Range boundary
  geom_sf(data = weighted_centroid_sf, color = 'red', size = 3) +  # Centroid in red
  geom_sf(data = pop_points_EM, aes(color = dist_fromCenter_EM), size = 2) +  # pop points colored by distance
  scale_color_viridis_c() +  # Color scale for distances
  theme_minimal() +
  labs(title = "EM sites - showing distance from centroid")

write.table(pop_points_EM,"C:/Users/rmohn/Desktop/00_Scripts_and_Labels/00_METADATA/dist_range_margin_centroid2.txt",sep="\t",quote = F,row.names = F)
#Add columns back to main df (could also use new df but this keeps things consistent)
#emWider_rootsITS <- emWider_rootsITS %>%
#  left_join(pop_points_EM[, c("SampleID", "dist_fromCenter_EM", "dist_fromEdge_EM")], by = "SampleID")

#emWider_rootsITS_QUMA <- emWider_rootsITS_QUMA %>% #quercus only df
#  left_join(pop_points_EM[, c("SampleID", "dist_fromCenter_EM", "dist_fromEdge_EM")], by = "SampleID")
