# Following aallyn's vignette but applying it to the ALB tag data
library(VAST)
library(tidyverse)
library(terra)
library(sf)
library(here)

source("C:/Users/nfarc/Desktop/RCodes_Shapefiles/TargetsSDM/R/vast_functions.r")

###########
# data prep
###########

#load in data
ALB_tag2<-here("data","ALB_tag","Pres_Abs_ALB_tag_20012to2016_1to3ratio_absenceconstrained_convexhull_enhanced.rds") %>% 
            readRDS() %>% 
            as.data.frame() %>% 
            na.omit()

#only need these select columns
ALB_tag2 <- ALB_tag2 %>% 
            dplyr::select(3:6,10:21,23,24) %>% 
            mutate(spec = "alb",
                   year_mon = lubridate::floor_date(date,"month")
                   )


## need to rescale data because VAST doesn't do it automatically
ALB_tag_scaled <- ALB_tag2 %>% 
                    mutate_at(.vars = c("ild", "n2", "sst_sd", "sst", "ssh_sd", "ssh", "eke", "salinity", "bathy", "rugosity", "dis_port", "dis_seamount"),
                              .funs = scale)
                    
# get rescaled params for later
rescale_params <- ALB_tag2 %>% 
                    summarize_at(.vars = c("ild", "n2", "sst_sd", "sst", "ssh_sd", "ssh", "eke", "salinity", "bathy", "rugosity", "dis_port", "dis_seamount"),
                              .funs = c("Mean"=mean, "SD" = sd), na.rm=TRUE)

#preping the df to include a "Pred_TF" column -- when Pred_TF == 1, the observation is only included in the predictions and NOT in the model fitting. This is a helpful way to do some quick model validation.
vast_ALB_data <- ALB_tag_scaled %>% 
                    mutate(Pred_TF = rep(0, nrow(ALB_tag_scaled)))
                    
# formating data for VAST  
vast_ALB_data <- vast_ALB_data  %>% 
  mutate(Year = lubridate::year(date),
         Count = if_else(Pres_abs == 1, as_units(1,"count"),as_units(0,"count") ),
         AreaSwept_km2 = 1) %>% 
  rename("Lat"="lat",
         "Lon"="lon") %>%
  dplyr::select(Year,Lat,Lon,Count,AreaSwept_km2,Pred_TF)


####################
# extrapolation grid
####################

# now we need to create our own extrapolation grid based on a region or someother shapefile
# this is following andrews 'vast_make_extrap_grd' function

# **if we have a shapefile like eez's it would start this way**

# region_shapefile <- st_read(here::here("data", "EEZ","eez.shp"))
# region_shapefile <- st_cast(region_shapefile, "POLYGON")
# region_shapefile <- region_shapefile[2,] # just was west coast eez

# for this example I am just going to use a convex hull

index_shapes <- c()
hull<-terra::convHull(terra::vect(vast_ALB_data %>% 
                              sf::st_as_sf(coords = c("Lon", "Lat"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
                              )
                  )

region_shapefile <- sf::st_as_sf(hull) 
# Transform crs of shapefile to common WGS84 lon/lat format.
region_wgs84 <- region_shapefile

# Get UTM zone
lon <- sum(st_bbox(region_wgs84)[c(1, 3)]) / 2
utm_zone <- floor((lon + 180) / 6) + 1

# Transform to the UTM zone
crs_utm <- st_crs(paste0("+proj=utm +zone=", utm_zone, " +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
region_utm <- st_transform(region_wgs84, crs = crs_utm)
  

# Make extrapolation grid with sf
cell_size <- 25000
region_grid <- st_as_sf(st_make_grid(region_utm, cellsize = cell_size, what = "centers"), crs = crs_utm)

# Now get only the points that fall within the shape polygon
points_keep <- data.frame("pt_row" = seq(from = 1, to = nrow(region_grid), by = 1), "in_out" = st_intersects(region_grid, region_utm, sparse = FALSE))
region_grid <- region_grid %>%
  mutate(., "in_poly" = st_intersects(region_grid, region_utm, sparse = FALSE)) %>%
  filter(., in_poly == TRUE)

# Convert back to WGS84 lon/lat, as that is what VAST expects.
extrap_grid <- region_grid %>%
  st_transform(., crs = 4326)

# Adding in the a strata/region component for stratified abundance. This will depend on index_shapes input. **I dont have an index_shapes input so I have a blank vector**
if (!is.null(index_shapes)) {
  extrap_grid <- extrap_grid %>%
    st_join(., index_shapes, join = st_within) %>%
    mutate(.,
      "Lon" = as.numeric(st_coordinates(.)[, 1]),
      "Lat" = as.numeric(st_coordinates(.)[, 2])
    ) %>%
    st_drop_geometry() %>%
    dplyr::select(., Lon, Lat, Region) %>%
    mutate(.,
      Area_km2 = ((cell_size / 1000)^2),
      STRATA = factor(Region, levels = index_shapes$Region, labels = index_shapes$Region)
    )
} else {
  extrap_grid <- extrap_grid %>%
    mutate(.,
      "Lon" = as.numeric(st_coordinates(.)[, 1]),
      "Lat" = as.numeric(st_coordinates(.)[, 2])
    ) %>%
    st_drop_geometry() %>%
    dplyr::select(., Lon, Lat) %>%
    mutate(., Area_km2 = ((cell_size / 1000)^2), row = 1:nrow(extrap_grid))
}

extrap_grid <- extrap_grid  %>% 
mutate(row = 1:nrow(extrap_grid))





 ################
 # model settings
 ################
 
# first setting the field and rho configuration settings - 
# Field config sets up the spatial/spatiotemporal components and how many factors should be estimate
# Rho config sets up the autoregressive structure on intercepts and spatiotemporal components
field_config <- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 0) # turn off last two for encounter only
rho_config <- c("Beta1" = 0, "Beta2" = 0, "Epsilon1" = 0, "Epsilon2" = 0) # need to turn off annual variation in intercept for positive catch model (Beta2 = 3)
obs_model <- c(1,0) # fitting only encounter only model -- first idex can be anything (we won't interper it), second index is 0 (indicates logit-link)
vast_settings<-FishStatsUtils::make_settings(n_x = 150, purpose = "index2", 
                                             Region = "USER", knot_method = "grid",
                                             bias.correct = FALSE, ObsModel = obs_model,
                                             Options = c("Calculate_Range" = TRUE),
                                             FieldConfig = field_config,
                                             RhoConfig = rho_config)

vast_settings


fit <- FishStatsUtils::fit_model(
 settings = vast_settings,
 Lat_i = vast_ALB_data$Lat,
 Lon_i = vast_ALB_data$Lon,
 t_i = vast_ALB_data$Year,
 b_i = vast_ALB_data$Count,
 a_i = vast_ALB_data$AreaSwept_km2,
 input_grid = extrap_grid,
 run_model = FALSE 
)

map_adjust<- fit$tmb_list$Map
map_adjust[['beta2_ft']] <- factor(rep(NA,5))


fit <- FishStatsUtils::fit_model(
 settings = vast_settings,
 Lat_i = vast_ALB_data$Lat,
 Lon_i = vast_ALB_data$Lon,
 t_i = vast_ALB_data$Year,
 b_i = vast_ALB_data$Count,
 a_i = vast_ALB_data$AreaSwept_km2,
 input_grid = extrap_grid,
 Map = map_adjust,
 run_model = TRUE 
)
summary(vast_ALB_data)
table(vast_ALB_data$Year)








vast0 <- vast_build_sdm(settings = vast_settings, extrap_grid = extrap_grid, 
                       sample_data = vast_ALB_data1, covariate_data = NULL, 
                       X1_formula = NULL, X2_formula = NULL, Q1_formula = NULL, 
                       Q2_formula = NULL, Xconfig_list = NULL, X_contrasts = NULL, 
                       index_shapes = NULL, spatial_info_dir = here::here(""))



#another method
## Specify model data and structure
?VAST::make_data()
  model_data <- VAST::make_data(
    b_i = vast_ALB_data$Count, #counts 
    a_i = vast_ALB_data$AreaSwept_km2, # no area is swept so treating it as standardized, 
    t_i = vast_ALB_data$Year,
    bias.correct=FALSE,
    ObsModel_ez = c(
      "PosDist" = 1, # positive catch rate ~ lognormal (but doesn't matter because only interested in encounter model)
      "Link" = 0 # encounter probability has logit link
    ),
    FieldConfig = c(
      "Omega1" = 1, # spatial variation in encounter probability
      "Epsilon1" = 1, # spatiotemporal variation in encounter probability
      "Omega2" = 0, # turning off spatial variation in positive catch rate
      "Epsilon2" = 0 # turning spatiotemporal variation in positive catch rate
    ), 
    spatial_list = EBS_knots[[i]], 
    Aniso = 1, # Geometric anisotropy,
    RhoConfig = c(
      "Beta1" = 0, # Fixed temporal effects in encounter probability
      "Beta2" = 3, # constant temporal effects in positive catch rate (need to turn off annual variation in intercept for positive catch model)
      "Epsilon1" = 0, # encounter probability, iid MVN
      "Epsilon2" = 0 # catch rate, iid MVN
    ))

