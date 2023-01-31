# Exploring the logbook data

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


# loading in the processed logbooks from 2001 - 2019
ALB_logbook <- readRDS(here("ALB_logbook", "Processed_ALBTrollLogbooks_2003_2016.rds")) %>%
    mutate(
        year = lubridate::year(date),
        month = lubridate::month(date)
    ) %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0") %>%
    mutate(
        lon = sf::st_coordinates(.)[, 1],
        lat = sf::st_coordinates(.)[, 2]
    )

alb_year_plot <- ALB_logbook %>%
    ggplot() +
    geom_point(aes(x = lon, y = lat)) +
    facet_wrap(~year) +
    geom_sf(data = world, color = "black", fill = "grey") +
    coord_sf(xlim = c(-180, -117), ylim = c(10, 53.5), expand = TRUE) +
    labs(y = "Latitude", x = "Longitude") +
    theme_bw()


ggsave(here("Plots", "alb_logbook_2003_2016.png"),
    width = 8, height = 5, units = "in", dpi = 300
)


# seasonal distribution patterns in the EEZ
EEZ <- here("data", "EEZ", "eez.shp") %>%
    sf::st_read(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

# ALB_logbook_EEZ <- st_join(ALB_logbook, EEZ) %>%
#     as_tibble() %>% 
#     filter(geoname == "United States Exclusive Economic Zone")

alb_logbook_seasonal<-ALB_logbook %>%
    ggplot() +
    geom_point(aes(lon, lat, fill = month), shape = 21) +
    scale_fill_viridis() +
    geom_sf(data = world, color = "black", fill = "grey") +
    geom_sf(data = EEZ, color = "black", fill = NA, size = 1) +
    facet_wrap(~year) +
    coord_sf(xlim = c(-180, -117), ylim = c(10, 53.5), expand = TRUE) +
    labs(y = "Latitude", x = "Longitude") +
    theme_bw()

ggsave(here("Plots", "alb_logbook_seasonal.png"),
    width = 10, height = 8, units = "in", dpi = 300
)
