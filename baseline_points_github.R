# preliminary
# compile baseline points corresponding to GBIF observational data
# for each species, all unique occurrence locations 2000-2004
# replicate the date for all years 2000-2020
# baseline data is a pseudo-observation dataset to simulate exposure to 
# climate change if a species did not move folliwng 2000-2004

library(tidyverse)
library(lubridate)
setwd("gbif/bg_points/")

# species list, change format of scientific names
list <- read_csv("species_list.csv")
list$Scientific <- gsub(" ","_",list$Scientific)
# loop seasons
for (season in c("breeding", "nonbreeding")){
  # Load pre-processed GBIF data, add month and year
  # The pre-processed GBIF dataset was assembled using SQL and annotated on 
  # Google Earth Engine, so there is no R code available
  alldata <- read_csv(paste0("data/",season,"_processed_20yr.csv"))
  alldata$month <- month(alldata$eventdate)
  alldata$year <- ifelse(alldata$month==12 & alldata$year>2009,alldata$year+1,alldata$year)
  # loop species
  for (sp in 1:nrow(list)){
    # cleanup environment following previous run
    gc(); if(sp==1){rm(sptab)}
    # try catch to keep running more species when one fails
    try({
      # species name formatting
      sp_name <- list$Scientific[sp]
      sp.name <- gsub("_"," ",sp_name)
      print(sp)
      # minimimum data requirement
      if(nrow(alldata[alldata$sci_name==sp.name,])>500){
        # grab 2000-2004 data, spatially filter to reduce spatial bias
        data <- alldata[alldata$sci_name==sp.name &
                          alldata$year<2005,] %>%
          dplyr::select(sci_name, eventdate, latitude, longitude) %>%
          group_by(latitude, longitude) %>% 
          sample_n(size = 1) %>%
          ungroup() 
        rm(sptab); gc()
        # for each point, replicate date for all years, add to running table
        # this generates 'baseline' data assuming pseudo-observations
        # corresponding to locations of 2000-2004 observations
        # we use this to estimate climate change exposure in absence of movement
        for (i in 1:nrow(data)){
        set <- seq(as.Date(paste0('2000-',str_sub(as.Date(data$eventdate[i]), start=-5))[[1]]),
                   length.out=21, by="years")
        tab <- merge(data.frame(cbind(sp.name, season, data$latitude[i], data$longitude[i])),
                     as.Date(set), all=T) 
        if(exists("sptab")){sptab <- rbind(sptab, tab)}else{sptab <- tab}
        }
        if(exists("fulltab")){fulltab <- rbind(fulltab, sptab)}else{fulltab <- sptab}
      }
    })
  }
  # name columns, create ID column, write table
  colnames(fulltab) <- c("species","season","latitude","longitude","eventdate")
  fulltab$ID <- 1:nrow(fulltab)
  write_csv(fulltab, paste0("~/gbif/bg_points/baseline_points_",season,".csv"))
}






