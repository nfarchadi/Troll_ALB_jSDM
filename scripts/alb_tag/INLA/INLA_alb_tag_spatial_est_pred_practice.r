
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

# # add season term
# yq <- as.yearqtr(as.yearmon(ALB_tag2$date, "%m/%d/%Y") + 1/12)
# ALB_tag2$season <- format(yq, "%q") %>% as.numeric()
# ALB_tag2 <- ALB_tag2 %>% 
#                mutate(season = if_else(season == 3, 1, 2))

n_alb_tag <- length(ALB_tag2$Pres_abs)
# n_seasons <- length(unique(ALB_tag2$season)) #only summar and fall which makes sense for when we just look at the cropped out region

# a spatial index that lists the spatial field index, which runs from 1 to mesh1$n for n_seasons 
# AND a spatial.field group which runs from 1 to n_years with each element replicated mesh1$n times
s_index <- inla.spde.make.index(name = "spatial_field",
                                n.spde = mesh1$n)
lengths(s_index)
# projection matrix for the integration points
data_alb_tag_A <- inla.spde.make.A(mesh = mesh1,
                                   loc = as.matrix(coords))

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

##################
# Model Prediction
##################

# read in env data - lets just look at one day (07/15/2013)
env_ras <- raster::stack("E:/HYCOM_NEP/hycom_combine_2013-07-15.grd")
env_ras <- crop(env_ras, NEP)
env_ras<-raster::dropLayer(env_ras, c("water_u","water_v","ssh"))
env_ras<-scale(env_ras)
env_ras <- aggregate(env_ras, fact = 5, fun = mean)
dp <- rasterToPoints(env_ras)
dp<-dp %>% as.data.frame()
coords_pred<-dp[,1:2] # get the coordinates for a grid
Ap <- inla.spde.make.A(mesh = mesh1, loc = as.matrix(coords_pred))
            
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


# Prediction stack
stk.pred <- inla.stack(
  data = list(y = NA, Ntrials = NA),
  A = list(1, Ap),
  effects = list(
       list(intercept = 1,
            eke = dp$eke,
            ild = dp$ild,
            n2 = dp$n2, 
            #ssh = inla.group(ALB_tag$ssh n = 25, method = "quantile"), # corellated with sst
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

stk.full <- inla.stack(stk.bi, stk.pred)

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
        f(spatial_field, model = spde)


###############
# Model Fitting
###############

t1 <- Sys.time()
output_bi <- inla(form,
                    family = "binomial",
                    data = inla.stack.data(stk.full),
                    Ntrials = Ntrials,
                    control.predictor = list(link = 1, A = inla.stack.A(stk.full), compute = TRUE),
                    verbose = TRUE,
                    safe = TRUE
                    )
t2 <- Sys.time()
t2-t1




#remake A matrix for pred coords
Aprediction <- inla.spde.make.A(mesh = mesh1, loc = as.matrix(coords_pred))

# Prediction stack
stk.pred <- inla.stack(
  data = list(y = NA),
  A = list(1, Aprediction),
  effects = list(
       list(intercept = 1,
            eke = inla.group(dp$eke, n = 10, method = "quantile"),
            ild = inla.group(dp$ild, n = 10, method = "quantile"),
            n2 = inla.group(dp$n2, n = 10, method = "quantile"),
            #ssh = inla.group(ALB_tag$ssh n = 25, method = "quantile"), # corellated with sst
            ssh_sd = inla.group(dp$ssh_sd, n = 10, method = "quantile"),
            sst = inla.group(dp$sst, n = 10, method = "quantile"),
            sst_sd = inla.group(dp$sst_sd, n = 10, method = "quantile"),
            bathy = inla.group(dp$bathy, n = 10, method = "quantile"),
            rugosity = inla.group(dp$rugosity, n = 10, method = "quantile")),
       #list(spatial_field = 1:spde$n.spde)
       s_index
    ),
  tag = "alb_tag_pred"
)

stk.full <- inla.stack(stk.bi, stk.pred)


###############
# Model Fitting
###############

t1 <- Sys.time()
output_pred <- inla(form,
                    family = "binomial",
                    data = inla.stack.data(stk.full, spde = spde),
                    Ntrials = Ntrials,
                    control.predictor = list(A = inla.stack.A(stk.full), compute = FALSE),
                    control.compute = list(config = TRUE),
                    # control.family = list(
                    #      list(link = "cloglog") #takes a lot longer using cloglog
                    #      ),
                    verbose = TRUE,
                    safe = TRUE
                    )
t2 <- Sys.time()
t2-t1
