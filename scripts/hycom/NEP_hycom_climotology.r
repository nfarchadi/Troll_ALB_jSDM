library(raster)
library(tidyverse)
library(foreach)
library(doParallel)

udates<-seq(as.Date("2012-01-01"),as.Date("2016-12-31"), by = "day")

HYCOM_NEP_dir<-"E:/HYCOM_NEP"

NEP_env_data<-list.files(HYCOM_NEP_dir, full.names = TRUE) %>%
  grep(pattern = ".grd", value = TRUE)
  
NEP_env_data<-NEP_env_data[1:(length(udates)-3)] # only need from 2012 - 2016

env_cov<-NEP_env_data[1] %>% stack() %>% names() 
env_cov<-env_cov[3:9] # only want from 3 to 10

for(i in 1:length(env_cov)){
    env_index <- i
    env_stack <- stack()
    for(ii in 1:length(NEP_env_data)){
        print(paste(env_cov[i],NEP_env_data[ii]))
        ras <- raster(NEP_env_data[ii], band = i)
        env_stack <- addLayer(env_stack, ras)
    }
    env_clim <- calc(env_stack, mean, na.rm=TRUE)
    names(env_clim)<-env_cov[i]
    writeRaster(env_clim, paste0(HYCOM_NEP_dir,"/clim_hycom/",env_cov[i],"_clim"),
                overwrite=TRUE)
}

# Now create a stack of the clim hycom data
NEP_clim_dir<-"E:/HYCOM_NEP/clim_hycom"
NEP_clim_env<-list.files(NEP_clim_dir, full.names = TRUE) %>%
  grep(pattern = ".grd", value = TRUE)

raster(NEP_env_data[1], band = 11)

clim_stack <- stack(NEP_clim_env)
clim_stack <- addLayer(clim_stack, 
                       raster(NEP_env_data[1], band = 11),
                       raster(NEP_env_data[1], band = 12)
                       )
writeRaster(clim_stack, paste0(NEP_clim_dir,"/NEP_clim_hycom"))
# scale the data
clim_stack_scaled<-scale(clim_stack)
writeRaster(clim_stack_scaled, paste0(NEP_clim_dir,"/NEP_clim_scaled_hycom"))
