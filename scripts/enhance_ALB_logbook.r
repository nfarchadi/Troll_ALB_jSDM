#enhancing AIS with environmental data 

#loading in packages
library(rgdal)
library(rgeos)
library(tidyverse)
library(raster)
library(sf)
library(stars)
library(here)

ALB_logbook<- readRDS("C:/Users/nfarc/Desktop/Grad Work SDSU/Chapter2_JSDM/ALB_Logbook/Processed_ALBTrollLogbooks_2009_2019.rds")

ALB_logbook<-ALB_logbook %>% filter(date >= as.Date('2012-01-01') & 
                                      date < as.Date('2020-01-01'))

####----enhancing ALB logbook data with NEP env data----####

HYCOM_NEP_dir<-"E:/HYCOM_NEP/"


#function to enhance AIS data
enhance_ALB<-function(input_df, env_dir){
  dates <-unique(input_df$date) %>% as.Date() #get unique dates from df
  enhanced_df<-data.frame()#need a empty df
  
  #for loop that subsets input df daily 
  #then enhances that data specific to the same date of a raster file
  #rbinds it all at the end
  for (i in 1:length(dates)){ 
    day_subset<- filter(input_df, 
                        input_df$date==dates[i])
    
    day_subset<-sf::st_as_sf(day_subset, coords = c("lon","lat"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
    
    #bring in the same env date
    env_file_list <- grep(list.files(path = env_dir, full.names = TRUE, pattern = ".grd"), pattern = dates[i], value = TRUE)
    
    if(length(env_file_list)==0){
      next
    } 
    
    env_day_stack<-stack(env_file_list)
    
    pts.env<-raster::extract(env_day_stack,day_subset)
    
    day.pts.env<-cbind(sf::st_coordinates(day_subset),day_subset,pts.env)%>% dplyr::select(-geometry)
    
    enhanced_df<-rbind(enhanced_df,day.pts.env)
  }
  return(enhanced_df) #returns the fully enhanced df
}


ALB_enhanced<-enhance_ALB(ALB_logbook,HYCOM_NEP_dir) #run it!

saveRDS(ALB_enhanced,"C:/Users/nfarc/Desktop/Grad Work SDSU/Chapter2_JSDM/ALB_Logbook/ALBTrollLogbooks_2012_2019_enhanced_presencesonly.rds")
