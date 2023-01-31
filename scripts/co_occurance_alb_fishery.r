# co-occurance between alb and fishery

library(here)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

world <- ne_countries(scale = "medium", returnclass = "sf")
sf_use_s2(FALSE) # need to do this to remove spherical geometry


alb<-here("ALB_tag","Processed_ALB_tag_2003_2016.rds") %>% readRDS()
fishery<-here("ALB_logbook","Processed_ALBTrollLogbooks_2003_2016.rds") %>% readRDS() %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0") %>%
    mutate(
        lon = sf::st_coordinates(.)[, 1],
        lat = sf::st_coordinates(.)[, 2]
    )

# where do they overlap? - lets try daily and monthly
overlap <- left_join(alb %>% as_tibble(),
                     fishery %>% as_tibble(), 
                     by=c("date","lon","lat")
                     ) %>% 
            mutate(overlap = if_else(is.na(vesselID) == TRUE, 0, 1))

overlap_plot<-overlap %>%
    ggplot() +
    geom_point(aes(x = lon, y = lat,color = overlap)) +
    facet_wrap(~year) +
    viridis::scale_color_viridis() +
    geom_sf(data = world, color = "black", fill = "grey") +
    coord_sf(xlim = c(-180, -117), ylim = c(10, 53.5), expand = TRUE) +
    labs(y = "Latitude", x = "Longitude") +
    theme_bw()           

ggsave(here("Plots", "overlap.png"),
    width = 10, height = 8, units = "in", dpi = 300
)
