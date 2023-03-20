# qc and enhance alb tag data

library(tidyverse)
library(raster)
library(rnaturalearth)
library(rnaturalearthdata)
library(cmocean)
library(zoo)
library(here)
library(sf)
library(tidyverse)
library(viridis)

world <- ne_countries(scale = "medium", returnclass = "sf")
sf_use_s2(FALSE) # need to do this to remove spherical geometry

source(here("scripts", "functions", "remove_land.r"))

# loading in the processed alb tag data from 2003 - 2016
ALB_tag <- read.csv(here("data","ALB_tag", "AlbacoreHMMoce.csv")) %>%
    mutate(
        longitude = if_else(longitude >= 180, longitude - 360, longitude)
    ) %>% # convert 360 to -180
    mutate(
        year = lubridate::year(date),
        month = lubridate::month(date)
    ) %>%
    sf::st_as_sf(coords = c("longitude", "latitude"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0") %>%
    mutate(
        lon = sf::st_coordinates(.)[, 1],
        lat = sf::st_coordinates(.)[, 2]
    ) %>% 
    filter(
        lon <= 0 # just want it for the eastern pacific
    )

# seasonal distribution in by year plot
EEZ <- here("data", "EEZ", "eez.shp") %>%
    sf::st_read(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

ALB_tag %>% 
    ggplot() +
    geom_point(aes(lon, lat, fill = month), shape = 21) +
    scale_fill_viridis() +
    geom_sf(data = world, color = "black", fill = "grey") +
    geom_sf(data = EEZ, color = "black", fill = NA, size = 1) +
    #facet_wrap(~year) +
    coord_sf(xlim = c(-180, -117), ylim = c(10, 53.5), expand = TRUE) +
    labs(y = "Latitude", x = "Longitude") +
    theme_bw()

ggsave(here("Plots", "alb_tag_seasonal.png"),
    width = 10, height = 8, units = "in", dpi = 300
)

#########################
# removing points on land
#########################

# just in case
bathy <- here("data", "static", "NEP_bathy.nc")

ALB_tag<- remove_land(ALB_tag, bathy_file = bathy)

ALB_tag <- ALB_tag %>% 
            dplyr::select(lon, lat, date, year, month, z) %>% 
            mutate(date = as.Date(date),
                   id = row_number())

###################################################
# enhancing AIS with hycom & bathy data for the NEP
###################################################

HYCOM_NEP_dir<-"D:/HYCOM_NEP"

# filter for dates past 2012 because I dont have env data for before
ALB_tag  <- ALB_tag %>% 
            filter(year >= 2012) %>% 
            na.omit()

# function to enhance AIS data
ALB_tag2<-enhance_data(ALB_tag ,HYCOM_NEP_dir) #run it!

##############################
# enhancing distance from port
##############################
dis_port_path<-grep(list.files(path = "D:/HYCOM_NEP", full.names = TRUE, pattern = ".nc"), pattern = "dis_port", value = TRUE)
dis_port<-raster(dis_port_path)


ALB_tag2<-sf::st_as_sf(ALB_tag2, coords = c("lon","lat"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")


pts.ports<-raster::extract(dis_port,ALB_tag2)
fleet.ports<-data.frame(cbind(sf::st_coordinates(ALB_tag2),pts.ports))
variable<-"dis_port"
ALB_tag2[,variable]<-fleet.ports[,3]


###################################
# enhancing distance from seamounts
###################################
dis_seamount_path<-grep(list.files(path = "D:/HYCOM_NEP", full.names = TRUE, pattern = ".nc"), pattern = "seamount", value = TRUE)
dis_seamount<-raster(dis_seamount_path)


ALB_tag2<-sf::st_as_sf(ALB_tag2, coords = c("lon","lat"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")


pts.seamounts<-raster::extract(dis_seamount,ALB_tag2)
fleet.seamounts<-data.frame(cbind(sf::st_coordinates(ALB_tag2),pts.seamounts))
variable<-"dis_seamount"
ALB_tag2[,variable]<-fleet.seamounts[,3]

# clean data
ALB_tag3<-na.omit(ALB_tag2)
nrow(ALB_tag3)

# Save the data
saveRDS(ALB_tag3, here("data","ALB_tag","ALB_tag_2012to2016_enhanced.rds"))
