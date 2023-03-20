# co-occurance between alb and fishery

library(here)
library(tidyverse)
library(sf)
library(raster)
library(rnaturalearth)
library(rnaturalearthdata)

world <- ne_countries(scale = "medium", returnclass = "sf")
sf_use_s2(FALSE) # need to do this to remove spherical geometry

template<-raster(here("data","static","NEP_bathy.nc")) 
x<-template
res(x)<-0.25
template<-resample(template,x, method='bilinear')


alb<-here("data","ALB_tag","Processed_ALB_tag_2003_2016.rds") %>% readRDS()
alb <- alb %>% mutate(month_year = zoo::as.yearmon(date),
                      pres = 1)

# res <- 1

# alb <- alb %>%
#     mutate(
#         lat = floor(lat / res) * res + 0.5 * res,
#         lon = floor(lon / res) * res + 0.5 * res
#     )




fishery<-here("data","ALB_logbook","Processed_ALBTrollLogbooks_2003_2016.rds") %>% readRDS()

fishery <- fishery %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0") %>%
    mutate(
        lon = sf::st_coordinates(.)[, 1],
        lat = sf::st_coordinates(.)[, 2]
    ) %>% 
    mutate(month_year = zoo::as.yearmon(date),
           pres = 1)


# column for each month -- need to make it wide format
#alb %>% dplyr::select(lon,lat,pres,month_year) %>% sf::st_drop_geometry() %>% gather(month_year,pres,-lon,-lat) %>% view()
unique_month_year_alb <- alb %>% pull(month_year) %>% unique()
alb_stack <- stack()
for (i in 1:length(unique_month_year_alb)){
    alb_month_year <- alb %>% filter(month_year == unique_month_year_alb[i])
    alb_month_year <- alb_month_year %>% dplyr::select(lon,lat,pres) %>% sf::st_drop_geometry()
    alb_month_year_ras <- rasterize(alb_month_year[,1:2], template, field=1)
    alb_stack <- addLayer(alb_stack, alb_month_year_ras)
}

alb_stack<-setZ(alb_stack, unique_month_year_alb)
names(alb_stack)<-unique_month_year_alb


unique_month_year_fishery <- fishery %>% pull(month_year) %>% unique()
fishery_stack <- stack()
for (i in 1:length(unique_month_year_fishery)){
    fishery_month_year <- fishery %>% filter(month_year == unique_month_year_fishery[i])
    fishery_month_year <- fishery_month_year %>% dplyr::select(lon,lat,pres) %>% sf::st_drop_geometry()
    fishery_month_year_ras <- rasterize(fishery_month_year[,1:2], template, field=1)
    fishery_stack <- addLayer(fishery_stack, fishery_month_year_ras)
}

fishery_stack<-setZ(fishery_stack, unique_month_year_fishery)
names(fishery_stack)<-unique_month_year_fishery

same_month_year<-intersect(unique_month_year_alb,unique_month_year_fishery) %>% zoo::as.yearmon()

alb_stack_points<-alb_stack %>% rasterToPoints() %>% 
                    as.data.frame() %>%
                    rename_all(funs(c("x","y",as.character(unique_month_year_alb)))) %>% 
                    gather(date, pres, -x, -y) %>% 
                    mutate(date = zoo::as.yearmon(date),
                           month = lubridate::month(date),
                           year = lubridate::year(date)) %>% 
                    filter(pres > 0)

fishery_stack_points<-fishery_stack %>% rasterToPoints() %>% 
                    as.data.frame() %>%
                    rename_all(funs(c("x","y",as.character(unique_month_year_fishery)))) %>% 
                    gather(date, pres, -x, -y) %>% 
                    mutate(date = zoo::as.yearmon(date),
                           month = lubridate::month(date),
                           year = lubridate::year(date)) %>% 
                    filter(pres > 0)

# alb_fishery_points<-rbind(alb_stack_points,fishery_stack_points)

# alb_fishery_points %>% group_by(x,y,year)



overlap <- full_join(alb_stack_points %>% as_tibble(),
                                fishery_stack_points %>% as_tibble(),
                                by=c("x","y", "date", "month", "year")) %>%
                        group_by(x,y,month, year) %>% 
                        summarise(overlap = pres.x + pres.y, .groups = "drop") %>%
                        mutate(overlap = if_else(overlap > 0, 1, 0),
                               overlap = as.factor(overlap))

overlap_month_graph<-overlap %>%
    na.omit() %>% 
    ggplot() +
    geom_tile(aes(x = x, y = y, fill = overlap, width = 0.75, height = 0.75)) +
    facet_wrap(~month) +
    viridis::scale_color_viridis() +
    geom_sf(data = world, color = "black", fill = "grey") +
    coord_sf(xlim = c(-180, -117), ylim = c(10, 53.5), expand = TRUE) +
    labs(y = "Latitude", x = "Longitude") +
    theme_bw()

ggsave(here("Plots", "overlap_month_ras.png"),
    width = 10, height = 8, units = "in", dpi = 300
)

overlap_year_graph<-overlap %>%
    na.omit() %>% 
    ggplot() +
    geom_tile(aes(x = x, y = y, fill = overlap, width = 0.75, height = 0.75)) +
    facet_wrap(~year) +
    viridis::scale_color_viridis() +
    geom_sf(data = world, color = "black", fill = "grey") +
    coord_sf(xlim = c(-180, -117), ylim = c(10, 53.5), expand = TRUE) +
    labs(y = "Latitude", x = "Longitude") +
    theme_bw()

ggsave(here("Plots", "overlap_year_ras.png"),
    width = 10, height = 8, units = "in", dpi = 300
)













# where do they overlap? - lets try daily and monthly
overlap <- inner_join(alb %>% as_tibble(),
                     fishery %>% as_tibble(), 
                     by=c("date","lon","lat")
                     ) %>% 
            mutate(overlap = if_else(is.na(vesselID) == TRUE, 0, 1))

overlap %>%
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

