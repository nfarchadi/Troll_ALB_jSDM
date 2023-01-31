# qc alb tag data

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
source(here("scripts", "functions", "generate_abs.r"))


# loading in the processed logbooks from 2001 - 2019
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
  

# reducing spatial resolution to 0.25
# each row is the location for each day for that certain fish
res <- 0.25

ALB_tag <- ALB_tag %>%
    mutate(
        lat = floor(lat / res) * res + 0.5 * res,
        lon = floor(lon / res) * res + 0.5 * res
    )


# seasonal distribution in by year
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

# Since I reduced the resolution some points are now on land so we want to remove those. I made a function that uses bathymetry raster to determine if points are on land or not.

bathy <- here("data", "static", "NEP_bathy.nc")

ALB_tag<- remove_land(ALB_tag, bathy_file = bathy)

##############################
# generating background points
##############################

# Use one of your environmental rasters as a template. Be sure to have it as the same resolution as what your predicted outputs will be. Most likely this will be the raster with the lowest common resolution.

template<-raster(here("data","static","NEP_bathy.nc")) 
x<-template
res(x)<-0.25
template<-resample(template,x, method='bilinear')

#convex hull polygon which will be used to select background points from
hull<-terra::convHull(terra::vect(ALB_tag))
hull<-as(hull, "Spatial")


string=seq(1:ncell(template))
template[]=string #adds values to the raster

# need to use rast for template because 
template<-mask(template,hull) # only points within the convexhull

# extract unique values to presences
unique_pres<-raster::extract(template,
                             sf::st_as_sf(ALB_tag, coords = c("lon","lat"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
ALB_tag<-ALB_tag %>% mutate(unique=unique_pres,type=1)


# filter dates based on presence availability
udates <-unique(ALB_tag$date)
udates<-as.Date(udates)

# for loop below uses the generate_abs function that will generate the same number of pseudo-absences as the presences for each day from areas that there isn't a presence point

# ratio argument allows you to increase the number of pseudo-absences daily by a factor of whatever value you give it.

# Just in case I like to have more pseudo-absences so I sent the ratio argument to 3 to have a 1:3 presence:pseudo-absence ratio

absences<-data.frame()
for (i in 1:length(udates)) {
  ud<-udates[i]
  ab<-generate_abs(ALB_tag, ud, template, ratio = 3)
  absences<-rbind(absences,ab)
}

# need a z column to match with presence df
absences<- remove_land(absences, bathy_file = bathy)

absences <- absences %>% 
            sf::st_as_sf(coords = c("lon","lat"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0") %>% 
            mutate(
                lon = sf::st_coordinates(.)[, 1],
                lat = sf::st_coordinates(.)[, 2],
                date = as.Date(date)
                ) %>% 
            dplyr::select(lon, lat, type, date, year, month, z)
            
            
ALB_tag <- ALB_tag %>% 
            dplyr::select(lon, lat, type, date, year, month, z) %>% 
            mutate(date = as.Date(date))

#now we combine
Pres_Abs_ALB_tag_2003to2016<-rbind(ALB_tag,absences)


#make an id column. this will help when cbinding dataframes later
Pres_Abs_ALB_tag_2003to2016<-mutate(Pres_Abs_ALB_tag_2003to2016, id = row_number())


saveRDS(Pres_Abs_ALB_tag_2003to2016,here("data", "ALB_tag","Pres_Abs_ALB_tag_2003to2016_1to3ratio_absenceconstrained_convexhull.rds"))

