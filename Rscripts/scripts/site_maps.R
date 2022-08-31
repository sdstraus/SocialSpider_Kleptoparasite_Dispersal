library(tidyverse)
library(sf)
library(spData)
# (also try package - mapdata, maybe maptools)

# 
# ecuador <- read_sf("../data/paisaje/paisaje.shp")
coords <- read.csv("../data/eximius_coordinates.csv")

library(ggmap)
#devtools::install_github('oswaldosantos/ggsn')
library(ggsn)

coords <- coords %>% 
  mutate(locale = recode(locale,
                         sumaco = "Via Loreto",
                         archidona = "Archidona",
                         "jatun sacha" = "Jatun Sacha")) %>% 
  rename(Site = locale)

ecuador <- c(left = -77.85, bottom = -1.12, right = -77.5, top = -0.65)
all.sites <- get_map(ecuador) %>% ggmap()+
  geom_point(data = coords, aes(x = longitude, y = latitude, shape = Site, size = 0.5, stroke = 1))+
  scale_shape_manual(values=c(0,1,2))+
  xlab("Longitude")+
  ylab("Latitude")+
  scalebar(data = coords, dist = 5, dist_unit = "km", model = 'WGS84', transform = TRUE, 
           location = "bottomleft")+
  labs(size = NULL)
all.sites



sumaco <- c(left = -77.85, bottom = -0.75, right = -77.5, top = -0.65)
vl.sites <- get_map(sumaco) %>% ggmap()+
  geom_point(data = coords, aes(x = longitude, y = latitude, shape = Site, size = 0.5, stroke = 1))+
  scale_shape_manual(values=c(0,1,2))+
  xlab("Longitude")+
  ylab("Latitude")+
  scalebar(data = coords, dist = 5, dist_unit = "km", model = 'WGS84', transform = TRUE, 
           location = "bottomleft")+
  labs(size = NULL)
vl.sites

js <- c(left = -77.7, bottom = -1.1, right = -77.55, top = -1.0)
js.sites <- get_map(js) %>% ggmap()+
  geom_point(data = coords, aes(x = longitude, y = latitude, shape = Site, size = 0.5, stroke = 1))+
  scale_shape_manual(values=c(0,1,2))+
  xlab("Longitude")+
  ylab("Latitude")+
  scalebar(data = coords, dist = 5, dist_unit = "km", model = 'WGS84', transform = TRUE, 
           location = "bottomleft")+
  labs(size = NULL)
js.sites
