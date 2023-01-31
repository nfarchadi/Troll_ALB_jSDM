# enhancing alb tag data w/ HYCOM env data

#enhancing AIS with environmental data 

#loading in packages
library(rgdal)
library(rgeos)
library(tidyverse)
library(raster)
library(sf)
library(here)

source(here("scripts","functions","enhance_data.r"))

#load in the presence-absence data if its not already in your global environment 
ALB_tag<-here("data","ALB_tag",
               "Pres_Abs_ALB_tag_2003to2016_1to3ratio_absenceconstrained_convexhull.rds") %>% readRDS()


#changing the heading name
ALB_tag<-ALB_tag %>% rename(Pres_abs = type)


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

# save it!
saveRDS(ALB_tag2,here("data","ALB_tag","Pres_Abs_ALB_tag_20012to2016_1to3ratio_absenceconstrained_convexhull_enhanced.rds"))


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


# save data again just in case
saveRDS(ALB_tag2,here("data","ALB_tag","Pres_Abs_ALB_tag_20012to2016_1to3ratio_absenceconstrained_convexhull_enhanced.rds"))
table(ALB_tag2$Pres_abs,ALB_tag2$year)

#######################
# clean up the data set
#######################

ALB_tag3<-na.omit(ALB_tag2)

# getting the data in a 1:1 ratio
udates <-unique(ALB_tag3$date)

absences<-data.frame()
for (i in 1:length(udates)){
  subdate_presence<-filter(ALB_tag3, 
                           ALB_tag3$date== udates[i] & ALB_tag3$Pres_abs == 1)
  
  subdate_absences= ALB_tag3 %>%
    filter((ALB_tag3$Pres_abs == 0 & ALB_tag3$date == udates[i])) %>%
    .[sample(nrow(.),nrow(subdate_presence)),]
  
  absences<-rbind(absences,subdate_absences)
}

#subsetting the presences 
presence<-filter(ALB_tag3, ALB_tag3$Pres_abs == 1)    

#making sure the number of presences and absences are the same
table(presence$year,presence$Pres_abs)
table(absences$year,absences$Pres_abs) 

ALB_tag3_1to1ratio_w_env<-rbind(presence,absences)   

saveRDS(ALB_tag3_1to1ratio_w_env, here("data","AIS_processed","ALB_tag","Pres_Abs_2013to2020_NEP_USA_TROLL_onlyfishing_1to1ratio_absenceconstrained_convexhull_enhanced.rds"))