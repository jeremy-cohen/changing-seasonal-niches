# building the ebird dataset
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

# get ebird files from lab HPC
setwd("/gpfs/gibbs/pi/jetz/data/species_datasets/occurrence/ebird/ebird_Apr2024/global/raw/")
lf = list.files("all_files/")
# loop to pull 1000 files at a time and process
for (h in 0:4){
for (i in ((h*1000)+1):((h*1000)+1000)){
  # try catch to continue if it crashes
  try({
  print(i)
    # read both a data and index file
 eb2 = read_delim(paste0("all_files_newcol/", lf[i]))
 eb1 = read_delim(paste0("all_files/", lf[i])) %>%
   # combine
   left_join(eb2, by="molid") %>%
   # month and year columns
   mutate(month = month(eventDate),
          year = year(eventDate)) %>%
   # filter to our spatial extent, relevant seasons, species in our master list, relevant years
   # and filter out high-effort checklists incompatible with other data
   filter(between(latitude, -55, 87.11),
          between(longitude, -179.9, -42.68),
          month %in% c(1,2,6,7,8,12),
          scientificname %in% list$Scientific,
          duration_minutes<=300,
          effort_distance_km<=5,
          number_observers %in% 1:5,
          year %in% 2000:2020) %>%
   # keep only one checklist per sampling event to avoid duplicated lists
   group_by(sampling_event_identifier) %>%
   sample_n(size = 1) %>%
   ungroup() %>%
   # select only relevant columns
   dplyr::select(scientificname, longitude, latitude, eventDate, number_observers,
                 effort_distance_km, duration_minutes)
 # recombine
 if (i==((h*1000)+1)){eb = eb1}else{
   eb = rbind(eb, eb1)}
  })
}
  # write
  write_csv(eb, paste0("/gpfs/gibbs/pi/jetz/data/temporary_datasets/csn_jeremy/ebird_2024_",
            h,"_experiment.csv"))
  gc()
}

# rebuild each subset into master dataset
for (h in 0:4){
eb = read_csv(paste0("/gpfs/gibbs/pi/jetz/data/temporary_datasets/csn_jeremy/ebird_2024_",
            h,"_experiment.csv"))
if (h==0){ebird = eb}else{ebird = rbind(ebird, eb)}
}

# divide into seasonal datasets
setwd("/gpfs/gibbs/pi/jetz/data/temporary_datasets/csn_jeremy/")
ebird$month = month(ebird$eventDate)
breeding = ebird %>% filter(month %in% 6:8)
write_csv(breeding, "ebird_2024_breeding_experiment.csv")
nonbreeding = ebird %>% filter(month %in% c(1,2,12))
write_csv(nonbreeding, "ebird_2024_nonbreeding_experiment.csv")

# annotate with ST relevant temp, evi, precip, and elev
elev = raster("/vast/palmer/home.mccleary/jc3893/rasters/elevation_1KMmd_GMTEDmd.tif")
# loop seasons
for (season in c("breeding","nonbreeding")){
  # load seasonal data, add year
  ebird = read_csv(paste0("ebird_2024_",season,"_experiment.csv"))
  ebird$year = year(ebird$eventDate)
  # organize like original gbif dataset
  
  # remove points over oceans by intersecting with continental polygons
  data(countriesLow)
  sf::sf_use_s2(FALSE)
  cc <- st_as_sf(countriesLow) %>%
    st_union()
  ebird = st_as_sf(ebird, coords = c("longitude", "latitude"),
                     crs = st_crs(cc))
  ebird.int = st_intersects(ebird, cc, sparse=F)
  ebird = ebird[as.vector(ebird.int),]
  ebird = data.frame(sf:::as_Spatial(ebird))
  ebird = rename(ebird, "lon"="coords.x1", "lat"="coords.x2")
  ebird$optional = NULL
 
  # annotate with monthly weather data by looping years and months
  for (year in 2000:2019){
    print(year)
    # year subset, add month column
    datayr = ebird[ebird$year==year,]
    datayr$month = month(datayr$eventDate)
    # don't have data for 1999-2000 nonbreeding
    if(year==2000){datayr = datayr[datayr$month %in% c(6:12),]}
    # loop months
    mlist = sort(unique(datayr$month))
    for (month in mlist){
      print(month)
      # subset monthly data, separately subset out coordinates
      datam = datayr[datayr$month == month,]
      dm = datam[,c("lon","lat")]
      coordinates(dm) = c("lon","lat")
      
      # max temp, precip from CHELSA stored locally
      mc = sprintf("%02d",month)
      if(year==2019 & month > 6){ datam$prec = NA}else{
        pr = raster(paste0("/gpfs/gibbs/pi/jetz/data/environmental_datasets/CHELSA_v2-1/monthly/pr/CHELSA_pr_",
                           mc,"_",year,"_V.2.1.tif"))
        datam$prec = raster::extract(pr, dm)
      }
      tx = raster(paste0("/gpfs/gibbs/pi/jetz/data/environmental_datasets/CHELSA_v2-1/monthly/tasmax/CHELSA_tasmax_",
                         mc,"_",year,"_V.2.1.tif"))
      # adjust kelvin to Celsius
      datam$tmax = (raster::extract(tx, dm)/10)-273
      
      # evi from modis stored locally
      # get year/date codes for all files
      list.evi = list.files("/gpfs/gibbs/pi/jetz/data/environmental_datasets/MODIS/MOD13C1/")
      evi.vals = gsub("A","",list.evi)
      evi.vals = as.numeric(gsub(".tif","",evi.vals))
      # get avg year/date code for this month
      day = yday(datam$eventDate)
      dayc = sprintf("%03d",day)
      dayf = mean(as.numeric(paste0(year, dayc)))
      # find closest file before that date
      wl = which(evi.vals < dayf)
      wm = which.min(dayf-wl)
      # get file, extract
      evi = raster(paste0("/gpfs/gibbs/pi/jetz/data/environmental_datasets/MODIS/MOD13C1/",
                          list.evi[wm]))
      datam$modisevi = raster::extract(evi, dm)*.0001
      
      # elevation
      datam$elev = raster::extract(elev, dm)
      
      # combine
      if(month==min(mlist)){datacomb=datam}else{datacomb = rbind(datacomb, datam)}
    }
    # combine
    if(year==2000){datafinal=datacomb}else{datafinal = rbind(datafinal, datacomb)}
  }
  datafinal$month = NULL
  # write
  write_csv(datafinal, paste0("ebird_2024_",season,"_for_csn.csv"))
}


# Baseline data generation from original data - see baseline_points_github.R for explanation
setwd("/gpfs/gibbs/pi/jetz/data/")
for (season in c("breeding", "nonbreeding")){ 
  # import original seasonal data
  alldata = read_csv(paste0("temporary_datasets/csn_jeremy/ebird_2024_",season,"_for_csn.csv"))
  alldata$month = month(alldata$eventDate)
  alldata$adjyear = ifelse(alldata$month %in% c(1,2),alldata$year-1,alldata$year)
  # subset to 2000-2004 data
  data = alldata[alldata$adjyear < 2005,]
  high = floor(nrow(data)/100000)
  rm(alldata)
  # remove feb 29ths
  data = data[!(month(data$eventDate)==2 & day(data$eventDate)==29),]
  gc()
  # divide up by 100000 records at a time
  for (n in 0:high){
    if (file.exists(paste0("temporary_datasets/csn_jeremy/",season,
                           "_baseline_ebird_",n,".csv"))==F){
      rm(fulltab, datafinal)
      if(n==high){upper=nrow(data)}else{upper=((100000*n)+100000)}
      # loop through each batch
      for (i in ((100000*n)+1):upper){
        print(i)
        # generate series of dates corresponding to all 2000-2004 records
        date = data$eventDate[i]
        if(data$month[i] %in% c(1:2)){year(date) = 2001}else{year(date) = 2000}
        set = seq(date, (as.Date(date)+years(18)), by = 'years')
        tab = data.frame(data$scientificname[i],
                         data$lat[i], data$lon[i], set) 
        if(exists("fulltab")){fulltab = rbind(fulltab, tab)}else{fulltab = tab}
      }
      colnames(fulltab) = c("scientificname","lat","lon","eventDate")
      
      ### Annotation of baseline data with environmental data as above
      fulltab$year = year(fulltab$eventDate)
      for (year in sort(unique(fulltab$year))){ 
        print(year)
        rm(datacomb)
        datayr = fulltab[fulltab$year==year,]
        datayr$month = month(datayr$eventDate)
        # don't have data for 1999-2000 nonbreeding
        if(year==2000){datayr = datayr[datayr$month %in% c(6:12),]}
        # loop months
        mlist = sort(unique(datayr$month))
        for (month in mlist){
          datam = datayr[datayr$month == month,]
          dm = datam[,c("lon","lat")]
          coordinates(dm) = c("lon","lat")
          
          # tmax, prec stored locally
          mc = sprintf("%02d",month)
          pr = raster(paste0("environmental_datasets/CHELSA_v2-1/monthly/pr/CHELSA_pr_",
                             mc,"_",year,"_V.2.1.tif"))
          datam$prec = raster::extract(pr, dm)
          tx = raster(paste0("environmental_datasets/CHELSA_v2-1/monthly/tasmax/CHELSA_tasmax_",
                             mc,"_",year,"_V.2.1.tif"))
          datam$tmax = (raster::extract(tx, dm)/10)-273
          
          # evi
          # get year/date codes for all files
          list.evi = list.files("environmental_datasets/MODIS/MOD13C1/")
          evi.vals = gsub("A","",list.evi)
          evi.vals = as.numeric(gsub(".tif","",evi.vals))
          # get avg year/date code for this month
          day = yday(datam$eventDate)
          dayc = sprintf("%03d",day)
          dayf = mean(as.numeric(paste0(year, dayc)))
          # find closest file before that date
          wl = which(evi.vals < dayf)
          wm = which.min(dayf-wl)
          # get file, extract
          evi = raster(paste0("environmental_datasets/MODIS/MOD13C1/",
                              list.evi[wm]))
          datam$modisevi = raster::extract(evi, dm)*.0001
          
          # combine
          if(exists("datacomb")){datacomb = rbind(datacomb, datam)}else{datacomb=datam}
        }
        # combine
        if(exists("datafinal")){datafinal = rbind(datafinal, datacomb)}else{datafinal=datacomb}
      }
      write_csv(datafinal, paste0("temporary_datasets/csn_jeremy/",
                                  season,"_baseline_ebird_",n,".csv"))
    }
  }
  
  # recombine all parts
  for (n in 0:high){ 
    part = read_csv(paste0("temporary_datasets/csn_jeremy/",
                           season,"_baseline_ebird_",n,".csv"))
    if(n==0){bl=part}else{bl=rbind(bl, part)}
  }
  rm(part)
  write_csv(bl, paste0("temporary_datasets/csn_jeremy/",
                       season,"_baseline_ebird_annotated.csv"))
  
}
