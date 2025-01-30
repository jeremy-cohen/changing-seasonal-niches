# building the bbs dataset
library(tidyverse)
library(raster)
library(rworldmap)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggpubr)
library(lubridate)
library(rgdal)


# species master list
list = read_csv("/vast/palmer/home.mccleary/jc3893/30x30/species_list_2024.csv")

### BBS compilation
setwd("/gpfs/gibbs/pi/jetz/data/species_datasets/occurrence/BBS/")
  lf = list.files("obs/", recursive=T, full.names=T)
  # load and combine observations
  # all existing rows signify presence
  for (i in 1:length(lf)){
    print(i)
    bbs1 = read_csv(lf[i]) 
    if(i==1){bbs=bbs1}else{bbs = rbind(bbs, bbs1)}
  }
  bbs$AOU = as.numeric(bbs$AOU)
  bbs$StateNum = as.numeric(bbs$StateNum)
  bbs$Route = as.numeric(bbs$Route)
  bbs$RPID = as.numeric(bbs$RPID)
  # load species key, attach names based on code
  splist = read_csv("SpeciesList.csv")[,c(2,3,7,8)]
  splist$scientificname = paste0(splist$Genus,"_",splist$Species)
  bbs = left_join(bbs, splist, "AOU")
  # load location key, attach coords based on code
  loclist = read_csv("Routes.csv")
  bbs = left_join(bbs, loclist, c("CountryNum", "StateNum", "Route"))
  # assign June as month
  bbs$month = 6
  # select columns
  bbs = bbs %>% 
    dplyr::select(RouteDataID, CountryNum, StateNum, Route, RPID,
                  Year, AOU, English_Common_Name, scientificname, 
                  RouteName, Latitude, Longitude, RouteTypeID, month)
  bbs = rename(bbs, "longitude"="Longitude", "latitude"="Latitude", "year"="Year")
  # restrict years
  bbs = filter(bbs, between(year, 1999, 2020))
  # save
  setwd("/gpfs/gibbs/pi/jetz/data/temporary_datasets/csn_jeremy/")
  write_csv(bbs, "bbs_breeding_experiment.csv")
  
  
  # annotate with ST relevant temp, evi, precip, and elev
  elev = raster("/vast/palmer/home.mccleary/jc3893/rasters/elevation_1KMmd_GMTEDmd.tif")
  # load bbs data
  bbs = read_csv("bbs_breeding_experiment.csv")
  # remove points over oceans by intersecting with continental polygons
  data(countriesLow)
  sf::sf_use_s2(FALSE)
  cc <- st_as_sf(countriesLow) %>%
    st_union()
  bbs = st_as_sf(bbs, coords = c("longitude", "latitude"),
                     crs = st_crs(cc))
  bbs.int = st_intersects(bbs, cc, sparse=F)
  bbs = bbs[as.vector(bbs.int),]
  bbs = data.frame(sf:::as_Spatial(bbs))
  bbs = rename(bbs, "lon"="coords.x1", "lat"="coords.x2")
  bbs$optional = NULL
 
  # annotate with monthly weather data by looping years and months
  for (year in 2000:2019){
    print(year)
    # year subset, add month column
    datayr = bbs[bbs$year==year,]
      dm = datayr[,c("lon","lat")]
      coordinates(dm) = c("lon","lat")
      
      # max temp, precip from CHELSA stored locally
      mc = "06" # define month as june
      pr = raster(paste0("/gpfs/gibbs/pi/jetz/data/environmental_datasets/CHELSA_v2-1/monthly/pr/CHELSA_pr_",
                           mc,"_",year,"_V.2.1.tif"))
      datayr$prec = raster::extract(pr, dm)
      tx = raster(paste0("/gpfs/gibbs/pi/jetz/data/environmental_datasets/CHELSA_v2-1/monthly/tasmax/CHELSA_tasmax_",
                         mc,"_",year,"_V.2.1.tif"))
      # adjust kelvin to Celsius
      datayr$tmax = (raster::extract(tx, dm)/10)-273
      
      # evi from modis stored locally
      # get year/date codes for all files
      list.evi = list.files("/gpfs/gibbs/pi/jetz/data/environmental_datasets/MODIS/MOD13C1/")
      evi.vals = gsub("A","",list.evi)
      evi.vals = as.numeric(gsub(".tif","",evi.vals))
      # get avg year/date code for this month
      day = 160 # early june date
      dayc = sprintf("%03d",day)
      dayf = mean(as.numeric(paste0(year, dayc)))
      # find closest file before that date
      wl = which(evi.vals < dayf)
      wm = which.min(dayf-wl)
      # get file, extract
      evi = raster(paste0("/gpfs/gibbs/pi/jetz/data/environmental_datasets/MODIS/MOD13C1/",
                          list.evi[wm]))
      datayr$modisevi = raster::extract(evi, dm)*.0001
      
      # elevation
      datayr$elev = raster::extract(elev, dm)
      
    # combine
    if(year==2000){datafinal=datayr}else{datafinal = rbind(datafinal, datayr)}
  }
  # write
  write_csv(datafinal, paste0("bbs_for_csn.csv"))



  # Baseline data generation from original data - see baseline_points_github.R for explanation
setwd("/gpfs/gibbs/pi/jetz/data/")
# import original seasonal data
  alldata = read_csv("temporary_datasets/csn_jeremy/bbs_for_csn.csv")
  # subset to 2000-2004 data
  data = alldata[alldata$year < 2005,]
  high = floor(nrow(data)/100000)
  rm(alldata)
  gc()
  # divide up by 100000 records at a time
  for (n in 0:high){
    if (file.exists(paste0("temporary_datasets/csn_jeremy/breeding_baseline_ebird_",n,".csv"))==F){
      rm(fulltab, datafinal)
      if(n==high){upper=nrow(data)}else{upper=((100000*n)+100000)}
      # loop through each batch
      for (i in ((100000*n)+1):upper){
        print(i)
        # generate series of dates corresponding to all 2000-2004 records
        set = 2000:2019
        tab = data.frame(data$scientificname[i],
                         data$lat[i], data$lon[i], set) 
        if(exists("fulltab")){fulltab = rbind(fulltab, tab)}else{fulltab = tab}
      }
      colnames(fulltab) = c("scientificname","lat","lon","year")
      
      ### Annotation of baseline data with environmental data as above
      for (year in sort(unique(fulltab$year))){ 
        print(year)
        rm(datacomb)
        datayr = fulltab[fulltab$year==year,]
          dm = datayr[,c("lon","lat")]
          coordinates(dm) = c("lon","lat")
          
          # tmax, prec stored locally
          mc = "06"
          pr = raster(paste0("environmental_datasets/CHELSA_v2-1/monthly/pr/CHELSA_pr_",
                             mc,"_",year,"_V.2.1.tif"))
          datayr$prec = raster::extract(pr, dm)
          tx = raster(paste0("environmental_datasets/CHELSA_v2-1/monthly/tasmax/CHELSA_tasmax_",
                             mc,"_",year,"_V.2.1.tif"))
          datayr$tmax = (raster::extract(tx, dm)/10)-273
          
          # evi
          # get year/date codes for all files
          list.evi = list.files("environmental_datasets/MODIS/MOD13C1/")
          evi.vals = gsub("A","",list.evi)
          evi.vals = as.numeric(gsub(".tif","",evi.vals))
          # get avg year/date code for this month
          day = 160 # early june date
          dayc = sprintf("%03d",day)
          dayf = mean(as.numeric(paste0(year, dayc)))
          # find closest file before that date
          wl = which(evi.vals < dayf)
          wm = which.min(dayf-wl)
          # get file, extract
          evi = raster(paste0("environmental_datasets/MODIS/MOD13C1/",
                              list.evi[wm]))
          datayr$modisevi = raster::extract(evi, dm)*.0001
          
        # combine
        if(exists("datafinal")){datafinal = rbind(datafinal, datayr)}else{datafinal=datayr}
      }
      write_csv(datafinal, paste0("temporary_datasets/csn_jeremy/breeding_baseline_bbs_",n,".csv"))
    }
  }
  
  # recombine all parts
  for (n in 0:high){ 
    part = read_csv(paste0("temporary_datasets/csn_jeremy/breeding_baseline_bbs_",n,".csv"))
    if(n==0){bl=part}else{bl=rbind(bl, part)}
  }
  rm(part)
  write_csv(bl, paste0("temporary_datasets/csn_jeremy/breeding_baseline_bbs_annotated.csv"))
  

