## Attempting to use INLA

library(INLA)
library(tidyverse)
library(terra)
library(sf)
library(here)

load("C:/Users/nfarc/Downloads/MDG_clean.Rdata")

###########
# data prep
###########

# load in data
ALB_tag2<-here("data","ALB_tag","Pres_Abs_ALB_tag_20012to2016_1to3ratio_absenceconstrained_convexhull_enhanced.rds") %>% 
            readRDS() %>% 
            as.data.frame() %>% 
            na.omit()

# # getting the data in a 1:1 ratio
# udates <-unique(ALB_tag2$date)

# absences<-data.frame()
# for (i in 1:length(udates)){
#   subdate_presence<-filter(ALB_tag2, 
#                            ALB_tag2$date== udates[i] & ALB_tag2$Pres_abs == 1)
  
#   subdate_absences= ALB_tag2 %>%
#     filter((ALB_tag2$Pres_abs == 0 & ALB_tag2$date == udates[i])) %>%
#     .[sample(nrow(.),nrow(subdate_presence)),]
  
#   absences<-rbind(absences,subdate_absences)
# }

# absences <- absences  %>% 
#             mutate(year = lubridate::year(date))

# #subsetting the presences 
# presence<-ALB_tag2  %>%  filter(Pres_abs == 1) %>% 
#             mutate(year = lubridate::year(date))

# #making sure the number of presences and absences are the same
# table(presence$year ,presence$Pres_abs)
# table(absences$year ,absences$Pres_abs) 

# ALB_tag3<-rbind(presence,absences)   

# ALB_tag3 %>% 
#     ggplot() +
#     geom_point(aes(lon, lat)) +
#     coord_sf(xlim = c(-180, -117), ylim = c(10, 53.5), expand = TRUE) +
#     labs(y = "Latitude", x = "Longitude") +
#     theme_bw()
 
#only need these select columns
ALB_tag3 <- ALB_tag3 %>% 
            filter(Pres_abs == 1) # only want presences
            dplyr::select(3:6,10:21,23,24) %>% 
            mutate(spec = "alb",
                   year_mon = lubridate::floor_date(date,"month")
                   )


###############
# make boundary
###############

coords <- ALB_tag3[c("lon","lat")] # other examples use only unique coords but doesn't seem to make a differece (a least when using a boundary)
boundary <- inla.nonconvex.hull(as.matrix(coords))

###########
# make mesh
###########

# try to see how mesh differs with different combinations or agruments


# mesh0 <- inla.mesh.2d(boundary = boundary, max.edge = c(3,5))
# plot(mesh0)
# points(coords[c("lon", "lat")], col = "red", type = "p")

# mesh1 <- inla.mesh.2d(loc = coords, max.edge = c(3))
# plot(mesh1)
# points(coords[c("lon", "lat")], col = "red", type = "p")

# mesh2 <- inla.mesh.2d(loc = coords, boundary = boundary,max.edge = c(3))
# plot(mesh2)
# points(coords[c("lon", "lat")], col = "red", type = "p")


# mesh3 <- inla.mesh.2d(loc = coords, boundary = boundary, max.edge = c(3, 4))
# plot(mesh3)
# points(coords[c("lon", "lat")], col = "red", type = "p")

# mesh4 <- inla.mesh.2d(loc = coords, boundary = boundary, max.edge = c(3, 4),
#                       cutoff = 1.5)
# plot(mesh4)
# points(coords[c("lon", "lat")], col = "red", type = "p")

# mesh5 <- inla.mesh.2d(boundary = boundary, max.edge = c(3, 4),
#                       cutoff = 1.5)
# plot(mesh5)
# points(coords[c("lon", "lat")], col = "red", type = "p")

mesh6 <- inla.mesh.2d(boundary = boundary,
                      max.edge = c(2, 4), #max edge length in interior and exterior domain
                      offset = c(1, 2), #how far I want to extend  
                      cutoff = 1.5) # min edge length
plot(mesh6)
###################
# SPDE and A matrix
###################

# create SPDE object
spde <- inla.spde2.matern(mesh6, alpha=2)

# need to construct a spatial (or spatiotemporal) index - NOT SURE IF THIS IS NEEDED
s_index <- inla.spde.make.index(
    name = "spatial.field",
    n.spde = mesh6$n # or spde$n.spde
    #n.group = n_year if spatiotemporal model
    )


# constructing the A matrix 
A <- inla.spde.make.A(
    mesh = mesh6,
    loc = as.matrix(coords)
    ) #if it was spatiotemporal it would also take time as a group

dim(A) #2237 observations x 789 nodes


############
# INLA stack
############

#first stack: estimation 
n_data <- nrow(ALB_tag3)
stack_est <- inla.stack(
    data = list(cnt = ALB_tag3$Pres_abs),
    A = list(A, 1),
    effects = list(s_index,
                   list(Intercept = rep(1, n_data))),
    tag = "est")

n=list(s_index,
     list(Intercept = rep(1, n_data))
)

n2=list(c(s_index,list(Intercept=1)))

n %>% str()
n2 %>% str()
