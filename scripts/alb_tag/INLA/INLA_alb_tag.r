## Attempting to use INLA for ALB tag data
## Using a log-Gaussian Cox process (point process model)

library(INLA)
library(tidyverse)
library(terra)
library(sf)
library(here)
library(inlabru)
library(raster)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(zoo)
sf_use_s2(FALSE)

#world <- ne_countries(scale = "medium", returnclass = "sf")
source(here("scripts", "functions","collinearity.r"))

###########
# data prep
###########

# load in data
ALB_tag<-here("data","ALB_tag","ALB_tag_2012to2016_enhanced.rds") %>% 
            readRDS() %>% 
            #sf::st_drop_geometry() %>% 
            #as.data.frame() %>% 
            na.omit()

# subseting the data to just the California EEZ -- just for exploration until I am ready to model the whole region
# EEZ<-here("data","shapefiles","eez.shp") %>% 
#     sf::st_read(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
# EEZ <- st_cast(EEZ, "POLYGON")
# EEZ <- EEZ[2,]

NEP<-here("data","shapefiles","NEP.shp") %>% 
    sf::st_read(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

# smaller domain so model runs quicker
NEP<-st_crop(NEP, c(xmin = -135, ymin = 32, xmax = -117, ymax = 55))


ALB_tag2 <- st_intersection(ALB_tag, NEP) %>% 
            sf::st_drop_geometry() %>% 
            dplyr::select(1:22)


ALB_tag2 %>% 
    ggplot() +
    geom_sf(data = NEP, color = "black", fill = "grey") +
    geom_point(aes(lon, lat)) +
    coord_sf(xlim = c(-135, -117), ylim = c(32, 55), expand = TRUE) +
    labs(y = "Latitude", x = "Longitude") +
    theme_bw()
    
################################
# make boundary from a shapefile
################################

coords <- ALB_tag2[c("lon","lat")] # other examples use only unique coords but doesn't seem to make a differece (a least when using a boundary)
# boundary <- inla.nonconvex.hull(as.matrix(coords), convex=-0.05, concave=-0.2)


# NEP<-here("data","shapefiles","NEP.shp") %>% 
#     sf::st_read(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

# # smaller domain so model runs quicker
# NEP<-st_crop(NEP, c(xmin = -130, ymin = 33, xmax = -120, ymax = 50))

bdry <- inla.sp2segment(NEP %>% as('Spatial'))
bdry$loc <- inla.mesh.map(bdry$loc)

###########
# make mesh
###########

#recommendations from https://rpubs.com/jafet089/886687

# Setting the max.edge
# Several studies, tutorials and group discussions have suggested the max.edge value to be between 1/3 to 1/10 times 
# smaller than the spatial range (https://haakonbakkagit.github.io/btopic104.html)
# since the range cannot be known until the model has been fitted we can make a somewhat conservative assumption that it
# is about 1/3 of the study area and thus specify the max.edge value to be 1/5 of the spatial range

mxedge <- diff(range(coords[,1]))/(3*5)

bound.outer <- diff(range(coords[,1]))/(3)


mesh1 <- inla.mesh.2d(boundary = bdry, #using the boundry agrument since we have a shapefile
                      max.edge=c(2,4)*mxedge, #specifies the max triangle edge length in inner domain and outer extension
                      offset = c(mxedge,bound.outer), #used to set the extension distance
                      cutoff = 0.3 # 1/5 of max.edge recomended. when using a boarder  value is no longer affected by the distance between points but the boundary polygon itself. Thus may be need to reduce cutoff value to achieve a higher the precision of the coastline
                      ) 

mesh1$crs <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

ggplot() +
     geom_sf(data=NEP,color='turquoise',fill='transparent')+  
     gg(mesh1) +
     geom_point(data=coords,aes(x = lon, y=lat), col='purple',size=1.7,alpha=0.5) 


########################################
# SPDE, integration points, and A matrix
########################################

# create SPDE object
spde <- inla.spde2.matern(mesh1, alpha=2)


# create integration stack -- following Suhaimi et al. 2021 method
max_x <- max(mesh1$loc[,1])
max_y <- max(mesh1$loc[,2])
loc.d <- t(matrix(c(0, 0, max_x, 0, max_x, max_y, 0, max_y, 0, 0), 2))

#make dual mesh
dd <- deldir::deldir(mesh1$loc[, 1], mesh1$loc[, 2])
tiles <- deldir::tile.list(dd)

#make domain into spatial polygon
domainSP <- SpatialPolygons(list(Polygons(
  list(Polygon(loc.d)), '0')))

#intersection between domain and dual mesh
poly.gpc <- as(domainSP@polygons[[1]]@Polygons[[1]]@coords, "gpc.poly")

# w now contains area of voronoi polygons
w <- sapply(tiles, function(p) rgeos::area.poly(rgeos::intersect(as(cbind(p$x, p$y), "gpc.poly"), poly.gpc)))

nv <- mesh1$n
n_alb_tag <- nrow(ALB_tag2)

# change data to include 0s for nodes and 1s for presences. 
# only necessary for "unstructured" data types (i.e. PO) - corresponds to y.pp in Suhaimi "unstructured" example
y.pp_alb_tag <- rep(0:1, c(nv, n_alb_tag)) ## nv - mesh nodes, n_alb_tag - alb presence

# add expectation vector (area for integration points/nodes and 0 for presences)
e.pp_alb_tag <- c(w, rep(0, n_alb_tag)) 

#  projection matrices for the integration points and the observation points
 
# diagonal matrix for integration point A matrix
imat <- Diagonal(nv, rep(1, nv))

# projection matrix for the integration points
data_alb_tag_A <- inla.spde.make.A(mesh1, as.matrix(coords))

# the entire projection matrix
A.pp_alb_tag <- rbind(imat, data_alb_tag_A)

#################################################################
# Extract environmental data at integration points and scale data
#################################################################

# scaling the tag data
ALB_tag2 <- ALB_tag2 %>% 
                    mutate_at(.vars = c("ild", "n2", "sst_sd", "sst", "ssh_sd", "ssh", "eke",
                                         "salinity", "bathy", "rugosity", "dis_port", "dis_seamount"),
                              .funs = scale)

# get covariates for integration points
# for these we will use a climotological average across the whole time period
clim_stack <- raster::stack("E:/HYCOM_NEP/clim_hycom/NEP_clim_scaled_hycom.grd")

env_covariates <- raster::extract(clim_stack, cbind(mesh1$loc[,1], mesh1$loc[,2])) %>%
                    as.data.frame()

# Checking for collinearity 
collinearity(na.omit(ALB_tag2 %>% dplyr::select(sst:rugosity)))
#############
# INLA stacks
#############

# Estimation stack
stk.pp <- inla.stack(
  data = list(y = y.pp_alb_tag, e = e.pp_alb_tag),
  A = list(1,A.pp_alb_tag),
  effects = list(
    list(intercept = rep(1, nv + n_alb_tag),
         eke = c(env_covariates$eke, ALB_tag2$eke),
         ild = c(env_covariates$ild, ALB_tag2$ild),
         n2 = c(env_covariates$n2, ALB_tag2$n2),
         ssh = c(env_covariates$ssh, ALB_tag2$ssh),
         ssh_sd = c(env_covariates$ssh_sd, ALB_tag2$ssh_sd),
         sst = c(env_covariates$sst, ALB_tag2$sst),
         sst_sd = c(env_covariates$sst_sd, ALB_tag2$sst_sd),
         bathy = c(env_covariates$bathy, ALB_tag2$bathy),
         rugosity = c(env_covariates$rugosity, ALB_tag2$rugosity)
         ),
    list(spatial_field = 1:spde$n.spde)
  ),
  tag = "alb_tag_pp"
)

# # Prediction stack
# df_pred <- data.frame(lon = mesh1$loc[,1],
#                       lat = mesh1$loc[,2])
# n_pred <- nrow(df_pred)
# A_pred <- Diagonal(n = n_pred)

# stk.pred <- inla.stack(
#   data = list(y = NA, e = NA),
#   A = list(1, A_pred),
#   effects = list(
#     list(intercept = rep(1, n_pred)),
#     list(spatial_field = 1:spde$n.spde))
#   )


# # Combine stacks
# stack <- inla.stack(stk.pp, stk.pred)


###############
# Model Formula
###############

form <- y ~ -1 + 
        intercept +
        f(inla.group(eke, n = 10, method = "quantile"), model = "rw2", constr=FALSE)+
        f(inla.group(ild, n = 10, method = "quantile"), model = "rw2", constr=FALSE)+
        f(inla.group(n2, n = 10, method = "quantile"), model = "rw2", constr=FALSE)+
        f(inla.group(ssh, n = 10, method = "quantile"), model = "rw2", constr=FALSE)+
        f(inla.group(ssh_sd, n = 10, method = "quantile"), model = "rw2", constr=FALSE)+
        f(inla.group(sst, n = 10, method = "quantile"), model = "rw2", constr=FALSE)+
        f(inla.group(sst_sd, n = 10, method = "quantile"), model = "rw2", constr=FALSE) +
        f(inla.group(bathy, n = 10, method = "quantile"), model = "rw2", constr=FALSE) +
        f(inla.group(rugosity, n = 10, method = "quantile"), model = "rw2", constr=FALSE)+
        f(spatial_field, model = spde)

###############
# Model Fitting
###############

t1 <- Sys.time()
output_all <- inla(form,
               family = "poisson",
               data = inla.stack.data(stk.pp),
               E = inla.stack.data(stk.pp)$e,
               control.predictor = list(A = inla.stack.A(stk.pp), compute = TRUE),
               control.compute = list(dic = TRUE,
                                      cpo = TRUE,
                                      waic = TRUE),
               verbose = TRUE, 
               safe = TRUE
               )
t2 <- Sys.time()
t2-t1
saveRDS(output_all, here("data","results","alb_tag_result_prelim_all_env.rds"))
output <- here("data","results","alb_tag_result_prelim.rds") %>% readRDS()
#########################
# Posterior Distributions
#########################

# posterior marginal distribution for the intercept 
plot(output_all$marginals.fixed$intercept,
     type = "l")

# tau (precision) parameter (controls the variance)
plot(output$marginals.hyperpar$'Theta1 for spatial_field', 
     type = "l",
     ylab = "density",
     xlab = "Precision")

# k (scale) parameter
plot(output$marginals.hyperpar$'Theta2 for spatial_field',
     type = "l")

# The precision and scale parameter in the SPDE forumla are hard to interpret
# Belw we transform them to more interpretable quantities (variance and range)
# the do.tranf = TRUE makes sure that marginals are calculated in the same scale as the data
output.field <- inla.spde2.result(inla = output_all,
                                  name = "spatial_field",
                                  spde = spde,
                                  do.transf=TRUE)


# these show the posterior for sigma^2 and the range
# confused how this sigma^2 differs from just sigma

# posterior distribution of the marginal variance of the latent field
plot(output.field$marginals.variance.nominal[[1]],
     type = "l",
     xlab = expression(sigma^2),
     ylab = expression(p(sigma^2*"|"*Z)))

# posterior distribution for marginal practical range parameter
# i.e. distance value above which spatial dependencies become negligiable
plot(output.field$marginals.range.nominal[[1]], 
     type = 'l',
     xlab = expression(l),
     ylab = expression(p(l*"|"*Z)))

# Checking if the range is smaller than the offset -- From Maria's example
range <- inla.emarginal(function(x) x, output.field$marginals.range.nominal[[1]]) #I think this is just get the mean value for the range parameter
range < max(diff(range(coords[,1])), diff(range(coords[,2])))*0.40 # YAY!

######################################
# Projection of Random Field on a Grid
######################################

# NEP_dim<-st_bbox(NEP) %>% matrix(ncol = 2)
# (dxy <- apply(NEP_dim,1, diff))
# (r <- dxy[1]/dxy[2])

# need to create a projector object
m<-150
proj.grid.mat <- 
  inla.mesh.projector(mesh1, 
                      xlim=c(st_bbox(NEP)$xmin,st_bbox(NEP)$xmax),
                      ylim=c(st_bbox(NEP)$ymin,st_bbox(NEP)$ymax),
                      dims=c(1, 1)*m)


# next clean the grid (i.e. set NA to the values outside boundary)
proj.grid.df<-data.frame(x = proj.grid.mat$lattice$loc[,1], y = proj.grid.mat$lattice$loc[,2])


proj.grid.sf<-st_as_sf(proj.grid.df, coords = c('x','y'),
          crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

# figuring out which points are in the NEP or not
ov <- sf::st_intersects(proj.grid.sf, NEP)
i.map <- lengths(ov) == 0 # 0 == not in the NEP

### Project the values of the mean, sd, of the spatial effect ###

# Posterior Mean
mean.g <- inla.mesh.project(proj.grid.mat, output_all$summary.random$spatial$mean)
mean.g <- cbind(proj.grid.df, mean=as.vector(mean.g))
mean.g[i.map,]<-NA

mean.g %>% 
    ggplot() +
    geom_sf(data = NEP, color = "black", fill = "grey") +
    geom_raster(aes(x, y, fill = mean)) +
    coord_sf(xlim = c(-135, -117), ylim = c(32, 55), expand = TRUE) +
    viridis::scale_fill_viridis(option = "turbo") +
    labs(y = "Latitude", x = "Longitude") +
    theme_bw()

# Posterior SD
sd.g <- inla.mesh.project(proj.grid.mat, output$summary.random$spatial$sd)
sd.g <- cbind(proj.grid.df, SD=as.vector(sd.g))
sd.g[i.map,]<-NA

sd.g %>% 
    ggplot() +
    geom_sf(data = NEP, color = "black", fill = "grey") +
    geom_raster(aes(x, y, fill = SD)) +
    coord_sf(xlim = c(-135, -117), ylim = c(32, 55), expand = TRUE) +
    viridis::scale_fill_viridis(option = "turbo") +
    labs(y = "Latitude", x = "Longitude") +
    theme_bw()

# Posterior q0.025
quantile_0.025 <- inla.mesh.project(proj.grid.mat, output$summary.random$spatial$`0.025quant`)
quantile_0.025[i.map]<-NA
plot(NEP %>% st_geometry())
fields::image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           quantile_0.025, add=TRUE)

# Posterior q0.975
quantile_0.975 <- inla.mesh.project(proj.grid.mat, output$summary.random$spatial$`0.975quant`)
quantile_0.975[i.map]<-NA
plot(NEP %>% st_geometry())
fields::image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           quantile_0.975, add=TRUE)


















############################################### TEMPLATE FOR BINOMIAL ################################################################





#########################
## Using a Binomial model
#########################

# load in data
ALB_tag<-here("data","ALB_tag","Pres_Abs_ALB_tag_20012to2016_1to3ratio_absenceconstrained_convexhull_enhanced.rds") %>% 
            readRDS() %>% 
            #sf::st_drop_geometry() %>% 
            #as.data.frame() %>% 
            na.omit()

# subseting the data to just the CCS -- just for exploration until I am ready to model the whole region
NEP<-here("data","shapefiles","NEP.shp") %>% 
    sf::st_read(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

# smaller domain so model runs quicker
NEP<-st_crop(NEP, c(xmin = -135, ymin = 32, xmax = -117, ymax = 55))


ALB_tag2 <- st_intersection(ALB_tag, NEP) %>% 
            sf::st_drop_geometry() %>% 
            dplyr::select(1:24)

# getting the data in a 1:1 ratio
udates <-unique(ALB_tag2$date)

# this code is finding the number of daily presences and randomly picking the same number of pseudo-absences
absences<-data.frame()
for (i in 1:length(udates)){
  subdate_presence<-filter(ALB_tag2, 
                           ALB_tag2$date== udates[i] & ALB_tag2$Pres_abs == 1)
  
  subdate_absences= ALB_tag2 %>%
    filter((ALB_tag2$Pres_abs == 0 & ALB_tag2$date == udates[i])) 
    
    if(nrow(subdate_presence) <= nrow(subdate_absences)){
         subdate_absences <-subdate_absences %>%
         .[sample(nrow(.),nrow(subdate_presence)),]
    } else {
         subdate_absences <-subdate_absences %>%
         .[sample(nrow(.),nrow(subdate_presence), replace = TRUE),]
    }
    
  absences<-rbind(absences,subdate_absences)
}

#subsetting the presences 
presence<-filter(ALB_tag2, ALB_tag2$Pres_abs == 1)

#making sure the number of presences and absences are the same
table(presence$year,presence$Pres_abs)
table(absences$year,absences$Pres_abs) 

ALB_tag2<-rbind(presence,absences)


################################
# make boundary from a shapefile
################################

coords <- ALB_tag2[c("lon","lat")] # other examples use only unique coords but doesn't seem to make a differece (a least when using a boundary)
# boundary <- inla.nonconvex.hull(as.matrix(coords), convex=-0.05, concave=-0.2)


# NEP<-here("data","shapefiles","NEP.shp") %>% 
#     sf::st_read(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

# # smaller domain so model runs quicker
# NEP<-st_crop(NEP, c(xmin = -130, ymin = 33, xmax = -120, ymax = 50))

bdry <- inla.sp2segment(NEP %>% as('Spatial'))
bdry$loc <- inla.mesh.map(bdry$loc)

###########
# make mesh
###########

#recommendations from https://rpubs.com/jafet089/886687

# Setting the max.edge
# Several studies, tutorials and group discussions have suggested the max.edge value to be between 1/3 to 1/10 times 
# smaller than the spatial range (https://haakonbakkagit.github.io/btopic104.html)
# since the range cannot be known until the model has been fitted we can make a somewhat conservative assumption that it
# is about 1/3 of the study area and thus specify the max.edge value to be 1/5 of the spatial range

mxedge <- diff(range(coords[,1]))/(3*5)

bound.outer <- diff(range(coords[,1]))/(3)


mesh1 <- inla.mesh.2d(boundary = bdry, #using the boundry agrument since we have a shapefile
                      max.edge=c(2,4)*mxedge, #specifies the max triangle edge length in inner domain and outer extension
                      offset = c(mxedge,bound.outer), #used to set the extension distance
                      cutoff = 0.3 # 1/5 of max.edge recomended. when using a boarder  value is no longer affected by the distance between points but the boundary polygon itself. Thus may be need to reduce cutoff value to achieve a higher the precision of the coastline
                      )

mesh1$crs <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

ggplot() +
     geom_sf(data=NEP,color='turquoise',fill='transparent')+  
     gg(mesh1) +
     geom_point(data=coords,aes(x = lon, y=lat), col='purple',size=1.7,alpha=0.5) 

########################################
# SPDE, spatial index, and A matrix
########################################

# create SPDE object
spde <- inla.spde2.matern(mesh1, alpha=2)

# add season term
yq <- as.yearqtr(as.yearmon(ALB_tag2$date, "%m/%d/%Y") + 1/12)
ALB_tag2$season <- format(yq, "%q") %>% as.numeric()
ALB_tag2 <- ALB_tag2 %>% 
               mutate(season = if_else(season == 3, 1, 2))

n_alb_tag <- length(ALB_tag2$Pres_abs)
n_seasons <- length(unique(ALB_tag2$season)) #only summar and fall which makes sense for when we just look at the cropped out region

# a spatial index that lists the spatial field index, which runs from 1 to mesh1$n for n_seasons 
# AND a spatial.field group which runs from 1 to n_years with each element replicated mesh1$n times
s_index <- inla.spde.make.index(name = "spatial_field",
                                n.spde = mesh1$n,
                                n.group = n_seasons)

# projection matrix for the integration points
data_alb_tag_A <- inla.spde.make.A(mesh = mesh1,
                                   loc = as.matrix(coords),
                                   
                                   group = ALB_tag2$season)

#################################################################
# Extract environmental data at integration points and scale data
#################################################################

# scaling the tag data
ALB_tag2 <- ALB_tag2 %>% 
                    mutate_at(.vars = c("ild", "n2", "sst_sd", "sst", "ssh_sd", "ssh", "eke",
                                         "salinity", "bathy", "rugosity", "dis_port", "dis_seamount"),
                              .funs = scale)

# Checking for collinearity 
collinearity(na.omit(ALB_tag2 %>% dplyr::select(sst:rugosity)))
#############
# INLA stacks
#############

# Estimation stack
stk.bi <- inla.stack(
  data = list(y = ALB_tag2$Pres_abs, Ntrials = rep(1, n_alb_tag)), #Ntrials could be just set to 1 also 
  A = list(1, data_alb_tag_A),
  effects = list(
       list(intercept = rep(1, n_alb_tag),
            eke = ALB_tag2$eke,
            ild = ALB_tag2$ild,
            n2 = ALB_tag2$n2,
            #ssh = inla.group(ALB_tag$ssh n = 25, method = "quantile"), # corellated with sst
            ssh_sd = ALB_tag2$ssh_sd,
            sst = ALB_tag2$sst,
            sst_sd = ALB_tag2$sst_sd,
            bathy = ALB_tag2$bathy, 
            rugosity = ALB_tag2$rugosity),
       #list(spatial_field = 1:spde$n.spde)
       s_index
    ),
  tag = "alb_tag_bi"
)


###############
# Model Formula
###############
# inla.group() bins data into groups according to the values of the covariate
# The option "quantile" uses equidistant quantiles in the probability space
form <- y ~ -1 + 
        intercept +
        f(inla.group(eke, n = 10, method = "quantile"), model = "rw2", constr=FALSE)+
        f(inla.group(ild, n = 10, method = "quantile"), model = "rw2", constr=FALSE)+
        f(inla.group(n2,, n = 10, method = "quantile"), model = "rw2", constr=FALSE)+
        #f(ssh, model = "rw2", constr=FALSE)+
        f(inla.group(ssh_sd, n = 10, method = "quantile"), model = "rw2", constr=FALSE)+
        f(inla.group(sst, n = 10, method = "quantile"), model = "rw2", constr=FALSE)+
        f(inla.group(sst_sd, n = 10, method = "quantile"), model = "rw2", constr=FALSE) +
        f(inla.group(bathy, n = 10, method = "quantile"), model = "rw2", constr=FALSE) +
        f(inla.group(rugosity, n = 10, method = "quantile"), model = "rw2", constr=FALSE)+
        f(spatial_field,
          model = spde,
          group = spatial_field.group,
          control.group = list(model = "exchangeable"))


###############
# Model Fitting
###############

t1 <- Sys.time()
output_bi <- inla(form,
                    family = "binomial",
                    data = inla.stack.data(stk.bi),
                    Ntrials = Ntrials,
                    control.predictor = list(A = inla.stack.A(stk.bi), compute = TRUE),
                    control.compute = list(dic = TRUE,
                                           cpo = TRUE,
                                           waic = TRUE),
                    # control.family = list(
                    #      list(link = "cloglog") #takes a lot longer using cloglog
                    #      ),
                    verbose = TRUE,
                    safe = TRUE
                    )
t2 <- Sys.time()
t2-t1
saveRDS(output_bi_2, here("data","results","alb_tag_result_prelim_all_env_sp_withpred.rds"))
output <- here("data","results","alb_tag_result_prelim_all_env_sp.rds") %>% readRDS()
output_bi <- output

#########################
# Posterior Distributions
#########################

# posterior marginal distribution for the intercept 
plot(output_bi$marginals.fixed$intercept,
     type = "l")

# tau (precision) parameter (controls the variance)
plot(output_bi$marginals.hyperpar$'Theta1 for spatial_field', 
     type = "l",
     ylab = "density",
     xlab = "Precision")

# k (scale) parameter
plot(output_bi$marginals.hyperpar$'Theta2 for spatial_field',
     type = "l")

# The precision and scale parameter in the SPDE forumla are hard to interpret
# Belw we transform them to more interpretable quantities (variance and range)
# the do.tranf = TRUE makes sure that marginals are calculated in the same scale as the data
output.field <- inla.spde2.result(inla = output_bi,
                                  name = "spatial_field",
                                  spde = spde,
                                  do.transf=TRUE)


# these show the posterior for sigma^2 and the range
# confused how this sigma^2 differs from just sigma

# posterior distribution of the marginal variance of the latent field
plot(output.field$marginals.variance.nominal[[1]],
     type = "l",
     xlab = expression(sigma^2),
     ylab = expression(p(sigma^2*"|"*Z)))

# posterior distribution for marginal practical range parameter
# i.e. distance value above which spatial dependencies become negligiable
plot(output.field$marginals.range.nominal[[1]], 
     type = 'l',
     xlab = expression(l),
     ylab = expression(p(l*"|"*Z)))

# Checking if the range is smaller than the offset -- From Maria's example
range <- inla.emarginal(function(x) x, output.field$marginals.range.nominal[[1]]) #I think this is just get the mean value for the range parameter
range < max(diff(range(coords[,1])), diff(range(coords[,2])))*0.40 # YAY!

######################################
# Projection of Random Field on a Grid
######################################

# NEP_dim<-st_bbox(NEP) %>% matrix(ncol = 2)
# (dxy <- apply(NEP_dim,1, diff))
# (r <- dxy[1]/dxy[2])

# need to create a projector object
m<-150
proj.grid.mat <- 
  inla.mesh.projector(mesh1, 
                      xlim=c(st_bbox(NEP)$xmin,st_bbox(NEP)$xmax),
                      ylim=c(st_bbox(NEP)$ymin,st_bbox(NEP)$ymax),
                      dims=c(1, 1)*m)


# next clean the grid (i.e. set NA to the values outside boundary)
proj.grid.df<-data.frame(x = proj.grid.mat$lattice$loc[,1], y = proj.grid.mat$lattice$loc[,2])


proj.grid.sf<-st_as_sf(proj.grid.df, coords = c('x','y'),
          crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

# figuring out which points are in the NEP or not
ov <- sf::st_intersects(proj.grid.sf, NEP)
i.map <- lengths(ov) == 0 # 0 == not in the NEP

# ### Project the values of the mean, sd, of the spatial effect - NOT SPATIOTEMPORAL ###

# # Posterior Mean
# mean.g <- inla.mesh.project(proj.grid.mat, output_bi$summary.random$spatial$mean)
# mean.g <- cbind(proj.grid.df, mean=as.vector(mean.g))
# mean.g[i.map,]<-NA

# mean.g %>% 
#     ggplot() +
#     geom_sf(data = NEP, color = "black", fill = "grey") +
#     geom_raster(aes(x, y, fill = mean)) +
#     coord_sf(xlim = c(-135, -117), ylim = c(32, 55), expand = TRUE) +
#     viridis::scale_fill_viridis(option = "turbo") +
#     labs(y = "Latitude", x = "Longitude") +
#     theme_bw()

# # Posterior SD
# sd.g <- inla.mesh.project(proj.grid.mat, output_bi$summary.random$spatial$sd)
# sd.g <- cbind(proj.grid.df, SD=as.vector(sd.g))
# sd.g[i.map,]<-NA

# sd.g %>% 
#     ggplot() +
#     geom_sf(data = NEP, color = "black", fill = "grey") +
#     geom_raster(aes(x, y, fill = SD)) +
#     coord_sf(xlim = c(-135, -117), ylim = c(32, 55), expand = TRUE) +
#     viridis::scale_fill_viridis(option = "turbo") +
#     labs(y = "Latitude", x = "Longitude") +
#     theme_bw()

# # Posterior q0.025
# quantile_0.025 <- inla.mesh.project(proj.grid.mat, output$summary.random$spatial$`0.025quant`)
# quantile_0.025[i.map]<-NA
# plot(NEP %>% st_geometry())
# fields::image.plot(proj.grid.mat$x, 
#            proj.grid.mat$y,
#            quantile_0.025, add=TRUE)

# # Posterior q0.975
# quantile_0.975 <- inla.mesh.project(proj.grid.mat, output$summary.random$spatial$`0.975quant`)
# quantile_0.975[i.map]<-NA
# plot(NEP %>% st_geometry())
# fields::image.plot(proj.grid.mat$x, 
#            proj.grid.mat$y,
#            quantile_0.975, add=TRUE)


#############################
# Visualize covariate effects
#############################
plot_data <- data.frame()
#cov_names <- output_bi$summary.random %>% names()
#cov_names <- cov_names[-length(cov_names)]
cov_names <- c("eke", "ild", "n2", "ssh_sd", "sst", "sst_sd", "bathy", "rugosity")



for (i in 1:length(cov_names)){
     Cov_summary.random <- output_bi$summary.random[[i]]
     cov_name<-cov_names[i]
     Mean_cov <- Cov_summary.random$mean
     Q_0.025_cov <- Cov_summary.random$`0.025quant`
     Q_0.975_cov <- Cov_summary.random$`0.975quant`
     Cov_value <- Cov_summary.random$ID *attr(ALB_tag2[,cov_name], 'scaled:scale') + attr(ALB_tag2[,cov_name], 'scaled:center')
     Cov <- rep(cov_name,length(Mean_cov))
     
     plot_cov_data <- data.frame(
          Cov = Cov, Mean_cov = Mean_cov,
          Q_0.025_cov = Q_0.025_cov, Q_0.975_cov = Q_0.975_cov,
          Cov_value = Cov_value)
     
     plot_data <- rbind(plot_data, plot_cov_data)
}


ggplot(plot_data) + 
  geom_line(aes(x=Cov_value,y=Mean_cov, color = Cov))+
  facet_wrap(~Cov, scales = "free")+
  #geom_line(aes(x=Cov_value,y=Q_0.025_cov),linetype ="dashed")+ 
  #geom_line(aes(x=Cov_value,y=Q_0.975_cov),linetype ="dashed")+
  theme_bw() +
  labs(y = "Marginal Effect", x = "Covariate Values") 


###############################################################
# Projection of Random Field on a Grid for Spatiotemporal Model
###############################################################

index_latent <- inla.stack.index(stk.bi, "alb_tag_bi")$data #index for the random field at the data locations
cor(ALB_tag2$Pres_abs, output_bi$summary.linear.predictor$mean[index_latent]) #correlation between the data response and the posterior mean of the latent values

lp_mean <- output_bi$summary.random$spatial$mean[index_latent]

# need to create a projector object
m<-150
proj.grid.mat <- inla.mesh.projector(mesh1, 
                      xlim=c(st_bbox(NEP)$xmin,st_bbox(NEP)$xmax),
                      ylim=c(st_bbox(NEP)$ymin,st_bbox(NEP)$ymax),
                      dims=c(1, 1)*m)

# next clean the grid (i.e. set NA to the values outside boundary)
proj.grid.df<-data.frame(x = proj.grid.mat$lattice$loc[,1], y = proj.grid.mat$lattice$loc[,2])

proj.grid.sf<-st_as_sf(proj.grid.df, coords = c('x','y'),
          crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

# figuring out which points are in the NEP or not
ov <- sf::st_intersects(proj.grid.sf, NEP)
i.map <- lengths(ov) == 0 # 0 == not in the NEP

### Project the values of the mean, sd, of the spatial effect -- SPATIOTEMPORAL MODEL ###

# Posterior Mean each time step (i.e. season)
xmean <- list()
for (j in 1:n_seasons){
     mean.g <- inla.mesh.project(proj.grid.mat, output_bi$summary.random$spatial_field$mean[s_index$spatial_field.group==j])
     mean.g <- cbind(proj.grid.df, mean=as.vector(mean.g))
     mean.g[i.map,]<-NA
     mean.g$season <- j
     xmean[[j]] <- mean.g
}
xmean<-bind_rows(xmean)

xmean %>% 
    ggplot() +
    geom_sf(data = NEP, color = "black") +
    geom_raster(aes(x, y, fill = mean)) +
    facet_wrap(~season)+
    coord_sf(xlim = c(-135, -117), ylim = c(32, 55), expand = TRUE) +
    viridis::scale_fill_viridis(option = "turbo") +
    labs(y = "Latitude", x = "Longitude") +
    theme_bw()


###############################################
# Prediction on a Grid for Spatiotemporal Model
###############################################
# read in env data - lets just look at one day (07/15/2013)
env_ras <- raster::stack("E:/HYCOM_NEP/hycom_combine_2013-07-15.grd")
env_ras <- crop(env_ras, NEP)
env_ras<-dropLayer(env_ras, c("water_u","water_v","ssh"))
env_ras<-scale(env_ras)

# reduce resolution to reduce computation -- FOR NOW
env_ras <- aggregate(env_ras, fact = 5, fun = mean)
dp <- rasterToPoints(env_ras)
dp<-dp %>% as.data.frame()
coords_pred<-dp[,1:2] # get the coordinates for a grid



# A matrix for prediction
data_alb_tag_Apred <-inla.spde.make.A(mesh = mesh1, 
                                      loc = as.matrix(coords_pred),
                                      group = 1, #selected season for prediction
                                      n.group = 2)

stk.pred<- inla.stack(
  data = list(y = NA),
  A = list(1, data_alb_tag_Apred), 
  effects = list(
       list(intercept = 1,
            eke = dp$eke,
            ild = dp$ild,
            n2 = dp$n2,
            #ssh = inla.group(ALB_tag$ssh n = 25, method = "quantile"), # correlated with sst
            ssh_sd = dp$ssh_sd,
            sst = dp$sst, 
            sst_sd = dp$sst_sd,
            bathy = dp$bathy,
            rugosity = dp$rugosity),
       #list(spatial_field = 1:spde$n.spde)
       s_index
    ),
  tag = "alb_tag_pred"
)

# Combine stacks
stk.full <- inla.stack(stk.bi, stk.pred)


# fit model
t1 <- Sys.time()
output_bi_pred <- inla(form,
                    family = "binomial",
                    data = inla.stack.data(stk.full, spde = spde),
                    Ntrials = Ntrials,
                    control.predictor = list(A = inla.stack.A(stk.full), compute = TRUE),
                    control.compute = list(dic = TRUE,
                                           cpo = TRUE,
                                           waic = TRUE),
                    # control.family = list(
                    #      list(link = "cloglog") #takes a lot longer using cloglog
                    #      ),
                    verbose = TRUE,
                    safe = TRUE
                    )
t2 <- Sys.time()
t2-t1


# predicting
index_pred <- inla.stack.index(stk.full, "alb_tag_pred")$data #index for the random field at the data locations
pred_mean_logit <- output_bi_pred$summary.fitted.values[index_pred, "mean"]
p.pred <- exp(pred_mean_logit)/(1 + exp(pred_mean_logit))

# x <- as.matrix(coords_pred)
# z <- as.matrix(p.pred)
pred_df <- data.frame(x = coords_pred$x, y = coords_pred$y,
           z = p.pred)
rasterFromXYZ(pred_df) %>% plot(zlim=c(0,1))
