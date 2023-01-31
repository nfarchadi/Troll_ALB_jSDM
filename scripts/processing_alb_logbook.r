#Processing ALB troll logbook data for 2003 to 2016

library(tidyverse)
library(here)
library(raster)


#read ALB troll logbook data
ALB<-read.csv(here("ALB_logbook", "cleanedTrollLogbooks_1995_2019.csv"))

#filter from 2003 - 2016
ALB_03_16 <- ALB %>% 
  dplyr::select(-discarded) %>% #remove discard column
  na.omit() %>% #remove rows that have NA
  #filter(kept > 0) %>% #remove rows with 0 ALB kept
  mutate(year = lubridate::year(date),
         month = lubridate::month(date)) %>% #year and month column
  filter(year >= 2003 & year <=2016)

#remove records that have lat and lon as whole numbers - Barbs recommendation by assuming that these were approximate locations
ALB_03_16 <- ALB_03_16[-which(ALB_03_16$lat%%1==0 & ALB_03_16$lon%%1==0),]

#coarsen data to 0.25 res to align with expected accuracy in location reporting, and vessel movements while fishing (see Nieto et al. 2017)
res <- 0.25

ALB_03_16_0.25<-ALB_03_16 %>% 
  mutate(lat = floor(lat/res) * res + 0.5 * res,
          lon = floor(lon/res) * res + 0.5 * res) %>% 
  group_by(vesselID, date, lat, lon) %>% 
  summarise(kept = sum(kept, na.rm = T), .groups = "drop")

saveRDS(ALB_03_16_0.25, here("ALB_logbook","Processed_ALBTrollLogbooks_2003_2016.rds"))

ALB_logbook<- readRDS("C:/Users/nfarc/Desktop/Grad Work SDSU/Chapter2_JSDM/ALB_Logbook/Processed_ALBTrollLogbooks_2003_2016.rds")

ALB_logbook<-sf::st_as_sf(ALB_logbook, coords = c("lon","lat"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")



# ####------obtaining background points------####

# ## Use one of your environmental rasters as a template. Be sure to have it as the same resolution as what your predicted outputs will be. Most likely this will be the raster with the lowest common resolution. 
# template<-raster("E:/HYCOM_NEP/hycom_combine_2012-01-01.grd") 


# #convex hull polygon which will be used to select background points from
# hull<-terra::convHull(terra::vect(ALB_logbook))
# hull<-as(hull, "Spatial")


# #res(template)=5 ## change to fit your data
# string=seq(1:ncell(template))
# template[]=string #adds values to the raster
# #need to use rast for template because 
# template<-mask(template,hull)#only points within the convexhull

# # extract unique values to presences
# unique_pres<-raster::extract(template,
#                              ALB_logbook)
# ALB_logbook<-ALB_logbook %>% mutate(unique=unique_pres,type=1)


# ## filter dates based on presence availability
# udates <-unique(ALB_logbook$date)
# udates<-as.Date(udates)

# #function to generate absences daily. this will match the same number as your presences or you can make it sample more 
# generate_abs<-function(input_df, udates){    
#   subdate_presence<-filter(ALB_logbook, 
#                            ALB_logbook$date==udates)
  
#   subdate_absences=rasterToPoints(template) %>% as.data.frame()  %>% rename("unique"="water_u") %>% as.data.frame() %>%
#     filter(!(unique %in%subdate_presence$unique)) %>% mutate(type=0) %>%  ## only select pixels that don't have presences
#     .[sample(nrow(.),nrow(subdate_presence)*3),] %>% ## create 1:3 ratio. want more absences in case for later 
#     mutate(date=udates,geartype="background",flag="No_flag",
#            hours=NA,fishing_hours=0,mmsi_present=0,name=NA) %>% 
#     rename("lon"="x", "lat"="y")
  
#   subdate_absences$month<-lubridate::month(subdate_absences$date)
#   subdate_absences$year<-lubridate::year(subdate_absences$date)
#   subdate_absences<-subdate_absences[,c(5,1,2,7,6,8,9,10,11,12,13,3,4)]
  
#   return(subdate_absences)
# }





# ## now we have dates, get absences from areas that there isn't a presence point for each day
# absences<-data.frame()
# for (i in 1:length(udates)) {
#   ud<-udates[i]
#   ab<-generate_abs(ALB_logbook, ud)
#   absences<-rbind(absences,ab)
# }


# # # combining presence and background dfs
# # absences<-absences[,-c(4,5,6,8,9)]
# # #need a z column to match with presence df
# # absences<- remove_land(absences, bathy_file = bathy_file)
# # NEPtroll$date<-as.Date(NEPtroll$date)
# # NEPtroll<-NEPtroll[,c(1,3,2,4,6,5,8,9,7)]
# # #now we combine
# # Pres_Abs_NEPtroll_2012to2020<-rbind(NEPtroll,absences)


# #make an id column. this will help when cbinding dataframes later
# Pres_Abs_ALB_logbook_2009to2019<-mutate(Pres_Abs_ALB_logbook_2009to2019, id = row_number())


# saveRDS(Pres_Abs_ALB_logbook_2009to2019, "C:/Users/nfarc/Desktop/Grad Work SDSU/Chapter2_JSDM/Pres_Abs_2009to2019_ALBTROLL_logbook_1to3ratio_absenceconstrained_convexhull.rds")
