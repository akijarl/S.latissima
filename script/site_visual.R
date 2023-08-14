setwd("C:/Users/aki/Desktop/S.latissima/analysis/")

library(sf)
library(terra)
library(dplyr)
library(spData)
#remotes::install_github("r-tmap/tmap@v4")
library(tmap)

isdat<-read_sf("C:/Users/aki/Downloads/data/ISL_adm0.shp")

bbox_new<-st_bbox(isdat)
bbox_new[1]<--24.5
bbox_new[2]<-63.5
bbox_new[3]<--20
bbox_new[4]<-65.5
bbox_new<-st_as_sfc(bbox_new)

sites<-read.csv(file = "sites.csv")
sites_sf <- st_as_sf(sites, coords = c("long", "lat"))

tm_shape(isdat,bbox=bbox_new)+
  tm_fill()+
  tm_borders()+
  tm_shape(sites_sf)+
  tm_dots(size=1.5,col="black")+
  tm_text("site",xmod=-.4,ymod=0.4)+
  tm_scale_bar(text.size = 0.7,position = c(0, -0.02))+
  #tm_scale_bar(text.size = 0.7,position = c("left", "bottom"))+
  tm_compass(type = "8star", size = 3, position = c(0.75, 0.7))



#library(maps)
#library(ggplot2)

#map_data_is <- map_data('world')[map_data('world')$region == "Iceland",]

#ggplot() +
  ## First layer: worldwide map
#  geom_polygon(data = map_data("world"),
#               aes(x=long, y=lat, group = group),
#               color = '#9c9c9c', fill = '#f3f3f3') +
  ## Second layer: Country map
#  geom_polygon(data = map_data_is,
#               aes(x=long, y=lat, group = group),
#               color = 'black', fill = 'grey') +
#  coord_map() +
#  coord_fixed(1.3,
#              xlim = c(-24.5, -20),
#              ylim = c(63.5, 65.5)) +
#  ggtitle("A map of Iceland") +
#  theme_classic()
