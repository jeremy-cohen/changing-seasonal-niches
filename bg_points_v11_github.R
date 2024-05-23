# code sheet 1
# species-level trend workflow with baseline and observer points

  library(tidyverse)
  library(rgdal)
  library(sf)
  library(sp)
  library(MVNH)
  library(ggpubr)
  library(lubridate)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(lwgeom)
  library(raster)
  library(rgeos)
  library(broom)
  library(farver)
  library(Cairo)
  library(grid)
  library(gridExtra)
  library(foreach)
  library(parallel)
  library(doParallel)
  library(emmeans)
  library(interactions)
  library(fields)
  set.seed(1)
  
  # change computing instructions based on local or HPC
  loc <- "local" # local or hpc
  if(loc=="local"){
    setwd("~/stoat/bg_points/")
    memory.limit(300000)
  }else{ setwd("/gpfs/ysm/home/jc3893/gbif/bg_points/")
    n.cores <- 2
    my.cluster <- parallel::makeCluster(
      n.cores,
      type = "PSOCK")
    doParallel::registerDoParallel(cl = my.cluster)
  }
# fix ellipse plotting glitch in ggplot
source("fix_stat_ellipse.R")
# directories
dir.create("figs/")
dir.create("figs/trends/")
dir.create("outputs/")
dir.create("outputs/metrics/")
dir.create("outputs/models/")
dir.create("outputs/tables/")

 # load species list
  full_list <- read_csv("species_list.csv")
   full_list$Scientific <- gsub(" ","_",full_list$Scientific)
   list <- full_list

  # shapefiles
  countries <- readOGR("shapefiles", "countries")
  aba <- countries[c(9,38),] %>%
    st_as_sf()
  
 # seasonal loop
  foreach (season = c("breeding", "nonbreeding")) %dopar% { 
    library(tidyverse)
    library(rgdal)
    library(sf)
    library(sp)
    library(MVNH)
    library(ggpubr)
    library(lubridate)
    library(rnaturalearth)
    library(rnaturalearthdata)
    library(lwgeom)
    library(raster)
    library(rgeos)
    library(broom)
    library(farver)
    library(Cairo)
    library(grid)
    library(gridExtra)
    library(foreach)
    library(parallel)
    library(doParallel)
    library(emmeans)
    library(interactions)
    library(fields)
  set.seed(1)
    alldata <- read_csv(paste0("data/",season,"_processed_20yr.csv"))
    # push december events to next year so years are contiguous
    alldata$month <- month(alldata$eventdate)
    alldata$year <- ifelse(alldata$month==12 & alldata$year>2009,alldata$year+1,alldata$year)
    # fix resulting year issues in winter data
    if (season=="nonbreeding"){alldata$year <- alldata$year-1}
    
    # similarly load and modify baseline points for season
    bl <- read_csv(paste0("data/baseline_points_",season,"_bl.csv")) 
    bl$month <- month(bl$eventdate)
    bl$year <- year(bl$eventdate)
    bl$year <- ifelse(bl$month==12,bl$year+1,bl$year)
    if (season=="nonbreeding"){bl$year <- bl$year-1}
    bl <- bl[complete.cases(bl$tmax) &
               complete.cases(bl$modisevi),]
    
# loop through species to find niche and range shifts
   for (sp in 1:nrow(list)){
      gc(); rm(locdata, spdata)
      try({
        sp_name <- list$Scientific[sp]
        sp.name <- gsub("_"," ",sp_name)
        common <- list$Common[sp]
        code <- as.character(list$ebird_code[sp]) 
        if(nrow(alldata[alldata$sci_name==sp.name,])==0){zeroes <- c(zeroes,common)}
        # skip if already done
        if(file.exists(paste0("outputs/models/mod_",sp_name,"_",season,".csv"))==F & #CHANGE BACK TO F
           nrow(alldata[alldata$sci_name==sp.name,])>1000){
          # get species range (range polygons can be downloaded from ebird:
          #https://science.ebird.org/en/status-and-trends/download-data)
          if(loc=="local"){rangepath<-"~/ranges/range_polygons/"}else{
            rangepath<-"/gpfs/ysm/home/jc3893/30x30/range_polygons/"}
          total.range <- readOGR(unzip(
            paste0(rangepath,code,"-range-2020.gpkg.zip"),
            paste0(code,"-range-mr-2020.gpkg")), "range")
          
          # seasonal range, exclude range outside US/Canada
          range <- total.range[total.range$season_name==season |
                                 total.range$season_name=="resident",]
          if(nrow(range)==0){STOP}
          # Find overlap with US/Canada
          rangest <- st_as_sf(range) %>%
            st_intersection(aba) %>%
            st_union()
          box <- st_bbox(rangest)
          plot(rangest, col="bisque")
          print(common)
          # remove range file output
          try(file.remove(paste0(code,"-range-mr-2020.gpkg")))
          
          # all species occurrence points 
          spdata <- alldata[alldata$sci_name==sp.name,] 
          blsp <- bl[bl$species==sp.name,]
          
          # observer dataset creation
          # get unique locations per year for all data in species range
          locdata <- alldata %>% 
            mutate(lon = round(longitude, 2),
                   lat = round(latitude, 2),
                   year = year(eventdate)) %>% 
            group_by(lon, lat, year) %>% 
            sample_n(size = 1) %>% 
            ungroup()
          gc()
          
          # get seasonal pts for all species in range
          spdata <- st_as_sf(spdata, coords = c("longitude", "latitude"), 
                             crs = st_crs(rangest))
          spdata.int <- st_intersects(spdata, rangest, sparse=F)
          spdata <- spdata[as.vector(spdata.int),]
          spdata <- data.frame(sf:::as_Spatial(spdata))
          spdata$sci_name <- sp_name
          spdata$optional <- NULL
          spdata <- rename(spdata, "longitude"="coords.x1", "latitude"="coords.x2")
          
          if(nrow(spdata)<500){STOPLOWDATA}
          
          locdata <- st_as_sf(locdata, coords = c("longitude", "latitude"), 
                           crs = st_crs(rangest))
          locdata.int <- st_intersects(locdata, rangest, sparse=F)
          locdata <- locdata[as.vector(locdata.int),]
          locdata <- data.frame(sf:::as_Spatial(locdata))
          locdata$sci_name <- gsub(" ","_",locdata$sci_name)
          locdata$optional <- NULL
          locdata <- locdata[locdata$year<2021,]
          locdata <- rename(locdata, "longitude"="coords.x1", "latitude"="coords.x2")
          
          # setup baseline data for this species
          blsp <- st_as_sf(blsp, coords = c("longitude", "latitude"), 
                           crs = st_crs(rangest))
          bl.int <- st_intersects(blsp, rangest, sparse=F)
          blsp <- blsp[as.vector(bl.int),]
          blsp <- data.frame(sf:::as_Spatial(blsp))
          blsp$sci_name <- gsub(" ","_",blsp$species)
          blsp$optional <- NULL
          blsp <- blsp[blsp$year<2021,]
          blsp <- rename(blsp, "longitude"="coords.x1", "latitude"="coords.x2")
          blsp <- blsp[blsp$year %in% 2000:2020,]
          blsp <- blsp %>% 
            group_by(latitude, longitude) %>%
            dplyr::mutate(locID = cur_group_id()) %>%
            ungroup()
          # check for complete set of years for every location
          blsp2 <- blsp[blsp$year %in% 2000:2019,]
          blsp_locs <- blsp2 %>% 
            group_by(latitude, longitude, locID) %>%
            summarise(count= n())
          wlocs <- which(blsp_locs$count==20)
          blsp <- blsp[blsp$locID %in% wlocs,]
          
          # loop through years, compile centroids, create maps, niche plots
          for (year in 2000:2020){ # year=2018
            spdata_yr <- spdata[spdata$year==year,] 
            if(nrow(spdata_yr)<20){STOPLOWDATAYR}
            locdata_yr <- locdata[locdata$year==year,]  
            bl_yr <- blsp[blsp$year==year,]
            
            # similarity between species pts and baseline pts
            niche_vars <- c("tmax", "modisevi")
            mvnh.dis <- MVNH_dissimilarity(spdata_yr[,niche_vars], 
                                           bl_yr[,niche_vars],
                                           metric="Pianka")
            niche_pianka_bl <- mvnh.dis$`Pianka`[1]
            niche_mahal_bl <- mvnh.dis$`Mahalanobis_distance`[1]
            niche_dr_bl <- mvnh.dis$`Determinant_ratio`[1]
            spat_vars <- c("longitude", "latitude")
            mvnh.dis <- MVNH_dissimilarity(spdata_yr[,spat_vars], 
                                           bl_yr[,spat_vars],
                                           metric="Pianka")
            spat_pianka_bl <- mvnh.dis$`Pianka`[1]
            spat_mahal_bl <- mvnh.dis$`Mahalanobis_distance`[1]
            spat_dr_bl <- mvnh.dis$`Determinant_ratio`[1]
            metrics <- round(c(niche_pianka_bl, niche_mahal_bl, niche_dr_bl,
                               spat_pianka_bl, spat_mahal_bl, spat_dr_bl), 5)
            assign(paste0("met.",year), metrics)
            
            # figure comparing sp/bg within season
            if(season=="breeding"){col="darkorange1"}else{col="dodgerblue"}
            
            # centroids
            centroids <- rbind(c(median(locdata_yr$tmax), median(locdata_yr$modisevi),
                                 median(locdata_yr$prec)),
                               c(median(spdata_yr$tmax), median(spdata_yr$modisevi),
                               median(spdata_yr$prec)),
                                 c(median(bl_yr$tmax), median(bl_yr$modisevi),
                                   c(median(bl_yr$prec))))
            sds <- rbind(c(sd(locdata_yr$tmax)/sqrt(nrow(locdata_yr)), 
                           sd(locdata_yr$modisevi)/sqrt(nrow(locdata_yr)), 
                           sd(locdata_yr$prec)/sqrt(nrow(locdata_yr))),
                         c(sd(spdata_yr$tmax)/sqrt(nrow(spdata_yr)), 
                           sd(spdata_yr$modisevi)/sqrt(nrow(spdata_yr)), 
                           sd(spdata_yr$prec)/sqrt(nrow(spdata_yr))),
                         c(sd(bl_yr$tmax)/sqrt(nrow(bl_yr)), 
                           sd(bl_yr$modisevi)/sqrt(nrow(bl_yr)), 
                           sd(bl_yr$prec)/sqrt(nrow(bl_yr))))
            centroids <- round(data.frame(cbind(centroids, sds)), 5)
            colnames(centroids) <- c("x","y","z","sd.x","sd.y","sd.z")
            centroids$lab <- paste0(round(centroids$x,1),", ",round(centroids$y,2))
            centroids_text <- text_grob(c(paste0("Baseline: ",centroids$lab[1],
                                                 "\nSpecies: ",centroids$lab[3])))
            assign(paste0("c.",str_sub(year, 3, 4)), centroids_text)
            
            gg <- ggplot(bl_yr, aes(x = tmax, y = modisevi)) +
              geom_point(bl_yr, mapping=aes(x = tmax, y = modisevi),
                         col="gray70", alpha=0.2) +
              geom_point(spdata_yr, mapping=aes(x = tmax, y = modisevi),
                         col=col, alpha=0.2) +
              stat_clip_ellipse_niche(geom="polygon", level=.9, fill="gray70", color="black", alpha=0.5, size=.8) +
              stat_clip_ellipse_niche(spdata_yr, mapping=aes(x = tmax, y = modisevi),
                                geom="polygon", level=.9, fill=col, color="black", alpha=0.5, size=.8) +
              xlim(-20,50) +
              ylim(0,1) +
              ggtitle(year) +
              theme_bw() +
              theme(text = element_text(size=13),
                    panel.grid.major = element_line(),
                    panel.grid.minor = element_line(),
                    axis.line = element_line(color = "black"),
                    axis.text.x = element_text(color = "black", size=13),
                    axis.text.y = element_text(color = "black", size=13),
                    plot.title = element_text(size = 22, face = "bold", hjust=.5)) +
              geom_point(data=centroids[c(1,3),], aes(x, y), size=5, col=c("gray50", "black")) +
              xlab("MODIS - Max. Temp") +
              ylab("MODIS - EVI")
            assign(paste0("gg.",str_sub(year, 3, 4)), gg)
            
            # map figure comparing species/baseline points in space within season
            world <- ne_countries(scale='medium', returnclass = 'sf')
            lim <- st_bbox(rangest)
            
            # collate
            centroid <- rbind(c(median(locdata_yr$longitude), median(locdata_yr$latitude), median(locdata_yr$elev)),
                              c(median(spdata_yr$longitude), median(spdata_yr$latitude), median(spdata_yr$elev)),
                              c(median(bl_yr$longitude), median(bl_yr$latitude), median(bl_yr$elev)))
            ses <- rbind(c(sd(locdata_yr$longitude)/sqrt(nrow(locdata_yr)), 
                           sd(locdata_yr$latitude)/sqrt(nrow(locdata_yr)), 
                           sd(locdata_yr$elev)/sqrt(nrow(locdata_yr))),
                         c(sd(spdata_yr$longitude)/sqrt(nrow(spdata_yr)), 
                           sd(spdata_yr$latitude)/sqrt(nrow(spdata_yr)), 
                           sd(spdata_yr$elev)/sqrt(nrow(spdata_yr))),
                         c(sd(bl_yr$longitude)/sqrt(nrow(bl_yr)), 
                           sd(bl_yr$latitude)/sqrt(nrow(bl_yr)), 
                           sd(bl_yr$elev)/sqrt(nrow(bl_yr))))
            centroid <- round(data.frame(cbind(centroid, ses)), 5)
            colnames(centroid) <- c("x","y","elev","sd.x","sd.y","sd.elev")
            centroid$lab <- paste0(round(centroid$x,1),", ",round(centroid$y,2))
            centroid_text <- text_grob(c(paste0("Baseline: ",centroid$lab[1],
                                                "\nSpecies: ",centroid$lab[3])))
            assign(paste0("ct.",str_sub(year, 3, 4)), centroid_text)
            
            # sequence for pretty x axis plotting
            # need to round to smaller vals if small range species
            rnd <- ifelse((max(spdata$longitude)-min(spdata$longitude))<3,1,0)
            sequ <- seq(round(min(spdata$longitude),rnd), 
                round(max(spdata$longitude),rnd))
            sequ <- sequ[c(3,round(median(1:length(sequ),1)),length(sequ)-2)]
            
            map <- ggplot() +
              geom_sf(data=world, bg="gray95") + 
              geom_sf(data=rangest, bg="bisque1") +
              coord_sf(xlim=lim[c(1,3)], ylim=lim[c(2,4)],
                       expand=F) +
              theme_bw() +
              theme(text = element_text(size=13),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(color = "black"),
                    axis.text.x = element_text(color = "black", size=13),
                    axis.text.y = element_text(color = "black", size=13),
                    plot.title = element_text(size = 22, face = "bold", hjust=.5)) +
              geom_point(locdata_yr, mapping=aes(x=longitude, y=latitude),
                         col="gray20") +
              geom_point(bl_yr, mapping=aes(x=longitude, y=latitude)) +
              geom_point(spdata_yr, mapping=aes(x=longitude, y=latitude),
                         col=col) +
              scale_x_continuous(breaks=c(min(sequ),median(sequ),max(sequ))) +
              xlab("") +
              ylab("") +
              ggtitle(year)
            assign(paste0("map.",str_sub(year, 3, 4)), map)
            
            centroids <- c(as.numeric(centroids[1,1:6]), as.numeric(centroids[2,1:6]), as.numeric(centroids[3,1:6]),
                           as.numeric(centroid[1,1:6]), as.numeric(centroid[2,1:6]), as.numeric(centroid[3,1:6]))
            names(centroids) <- c("tmax.bg","modisevi.bg","prec.bg",
                                  "se.tmax.bg","se.modisevi.bg","se.prec.bg",
                                  "tmax.sp","modisevi.sp","prec.sp",
                                  "se.tmax.sp","se.modisevi.sp","se.prec.sp",
                                  "tmax.bl","modisevi.bl","prec.bl",
                                  "se.tmax.bl","se.modisevi.bl","se.prec.bl",
                                  "lon.bg","lat.bg","elev.bg","se.lon.bg","se.lat.bg","se.elev.bg",
                                  "lon.sp","lat.sp","elev.sp","se.lon.sp","se.lat.sp","se.elev.sp",
                                  "lon.bl","lat.bl","elev.bl","se.lon.bl","se.lat.bl","se.elev.bl")
            assign(paste0("c.",year), centroids)
            
          } # end year
          
          # save plots
          if(loc=="hpc"){
           Cairo(file=paste0("figs/",common,"_", season,".jpeg"), type="jpeg",
                 height=800, width=1200)
          }else{
           jpeg(file=paste0("figs/",common,"_", season,".jpeg"),
                height=800, width=1200)
          }
          grid.arrange(gg.00, gg.10, gg.20, c.00, c.10, c.20,
                       map.00, map.10, map.20, ct.00, ct.10, ct.20,
                       ncol=3, nrow=4, heights=c(5,1,7,1))

          dev.off()
         
          # annual output of centroids
          out <- data.frame(cbind(
            common, sp_name, season, 2000:2020, 
            rbind(met.2000, met.2001,met.2002,met.2003,met.2004,met.2005,
                  met.2006,met.2007,met.2008,met.2009,met.2010,
                  met.2011,met.2012,met.2013,met.2014,met.2015,
                  met.2016,met.2017,met.2018,met.2019,met.2020),
            rbind(c.2000, c.2001,c.2002,c.2003,c.2004,c.2005,
                  c.2006,c.2007,c.2008,c.2009,c.2010,
                  c.2011,c.2012,c.2013,c.2014,c.2015,
                  c.2016,c.2017,c.2018,c.2019,c.2020)))
          colnames(out)[1:10] <- c("common_name","scientific_name","season","year",
                                   "niche_pianka_bl","niche_mahal_bl","niche_dr_bl",
                                   "spat_pianka_bl","spat_mahal_bl","spat_dr_bl")
          # write out
          write_csv(out, paste0("outputs/metrics/out_",sp_name,"_",season,".csv"))
          
          # model influence of factors on decadal trend (both sp and bg)
          keepers <- c("year","tmax","modisevi","prec","latitude","longitude","elev")   
          spdata <- spdata[,keepers]
          spdata$status <- "focal"
          # cap at 50k due to memory issues
          if(nrow(spdata)>50000){spdata <- sample_n(spdata, 50000, replace=F)}
          # distances from 2000-2004 centroid for all points
          spdata.early <- spdata[spdata$year<2005,]
          acstart <- c(median(spdata.early$longitude),median(spdata.early$latitude))
          ac <- as.matrix(spdata[,c("longitude","latitude")])
          acstart <- do.call("rbind", replicate(nrow(ac), acstart, simplify = FALSE))
          spdata$dist <- rdist.earth(ac, acstart)[,1] # great circle distance
          # subset observer data to equalize with species data
          locdata <- locdata[,keepers]
          locdata$status <- "background"
          locdata <- sample_n(locdata, nrow(spdata), replace=T)
          # distances from 2000-2004 centroid for all points
          locdata.early <- locdata[locdata$year<2005,]
          acstart <- c(median(locdata.early$longitude),median(locdata.early$latitude))
          ac <- as.matrix(locdata[,c("longitude","latitude")])
          acstart <- do.call("rbind", replicate(nrow(ac), acstart, simplify = FALSE))
          locdata$dist <- rdist.earth(ac, acstart)[,1] # great circle distance
          # subset baseline data to equalize with species data
          blsp <- blsp[,keepers]
          blsp$status <- "baseline"
          blsp <- sample_n(blsp, nrow(spdata), replace=T)
          gc()
          # fuse
          spblbgdata <- rbind(subset(spdata, select=-dist), blsp, 
                              subset(locdata, select=-dist))
          spbgdata <- rbind(spdata, locdata)
          rm(spdata, blsp, locdata)
          gc()
          
          # Models interacting year with data type to estimate niche/range shifts
          size = 27
          # tmax (baseline and observer)
          tmax_lm <- lm(tmax~year*status, data=spblbgdata)
          em <- data.frame(emtrends(tmax_lm, "status", var="year"))
          tmax_lm_trends_bg <- em[1,c(2:3,5:6)]
          tmax_lm_trends_bl <- em[2,c(2:3,5:6)] 
          tmax_lm_trends_sp <- em[3,c(2:3,5:6)]   
          write_csv(em, paste0("outputs/tables/",sp_name,"_",season,"_tmax_lm_em.csv"))      
          tmax_tab <- tidy(tmax_lm)
          tmax_summary <- c(tmax_tab[6,c(2,3,5)], tmax_lm_trends_bg,
                            tmax_lm_trends_bl, tmax_lm_trends_sp)
          names(tmax_summary) <- c("tmax_coef","tmax_se","tmax_p",
                                   "tmax_coef_bg","tmax_se_bg",
                                   "tmax_lcl_bg","tmax_ucl_bg",
                                   "tmax_coef_bl","tmax_se_bl",
                                   "tmax_lcl_bl","tmax_ucl_bl",
                                   "tmax_coef_sp","tmax_se_sp",
                                   "tmax_lcl_sp","tmax_ucl_sp")      
          write_csv(tmax_tab, paste0("outputs/tables/",sp_name,"_",season,"_tmax_lm.csv"))   
          # model predictions
          spdata <- spbgdata[spbgdata$status=="focal",]
          spbgdata$pre <- predict(tmax_lm, spdata)
          plm <- lm(pre~year, data=spbgdata)
          # individual slopes plots
          tmax.ip <- interact_plot(tmax_lm, pred = year, modx = status, interval=T,
                                   colors=c("gray70","black",col), lty=c(2,2,2)) +
            xlab("Year") +
            ylab("Temperature (?C)") +
            theme_classic() +
            theme(legend.position="none",
                  axis.text = element_text(size = size),
                  axis.title = element_text(size = size),
                  axis.line = element_line(size=2),
                  axis.ticks = element_line(size=2)) +
            scale_x_continuous(breaks=c(2000,2010,2020))
          # tmax (observer only)
          tmax_lm <- lm(tmax~year*status, data=spbgdata)
          em <- data.frame(emtrends(tmax_lm, "status", var="year"))
          tmax_lm_trends_bg <- em[1,c(2:3,5:6)]
          tmax_lm_trends_sp <- em[2,c(2:3,5:6)]   
          write_csv(em, paste0("outputs/tables/",sp_name,"_",season,"_tmax_nobl_lm_em.csv"))      
          tmax_tab <- tidy(tmax_lm)
          tmax_summary_nobl <- c(tmax_tab[4,c(2,3,5)], tmax_lm_trends_bg,
                            tmax_lm_trends_sp)
          names(tmax_summary_nobl) <- c("tmax_coef_nobl","tmax_se_nobl","tmax_p_nobl",
                                   "tmax_coef_bg_nobl","tmax_se_bg_nobl",
                                   "tmax_lcl_bg_nobl","tmax_ucl_bg_nobl",
                                   "tmax_coef_sp_nobl","tmax_se_sp_nobl",
                                   "tmax_lcl_sp_nobl","tmax_ucl_sp_nobl")      
          write_csv(tmax_tab, paste0("outputs/tables/",sp_name,"_",season,"_tmax_nobl_lm.csv"))    
          
          # evi (baseline + observer)
          modisevi_lm <- lm(modisevi~year*status, data=spblbgdata)
          em <- data.frame(emtrends(modisevi_lm, "status", var="year"))
          modisevi_lm_trends_bg <- em[1,c(2:3,5:6)]
          modisevi_lm_trends_bl <- em[2,c(2:3,5:6)] 
          modisevi_lm_trends_sp <- em[3,c(2:3,5:6)]     
          write_csv(em, paste0("outputs/tables/",sp_name,"_",season,"_modisevi_lm_em.csv"))     
          modisevi_tab <- tidy(modisevi_lm)
          modisevi_summary <- c(modisevi_tab[6,c(2,3,5)], modisevi_lm_trends_bg,
                            modisevi_lm_trends_bl, modisevi_lm_trends_sp)
          names(modisevi_summary) <- c("modisevi_coef","modisevi_se","modisevi_p",
                                   "modisevi_coef_bg","modisevi_se_bg",
                                   "modisevi_lcl_bg","modisevi_ucl_bg",
                                   "modisevi_coef_bl","modisevi_se_bl",
                                   "modisevi_lcl_bl","modisevi_ucl_bl",
                                   "modisevi_coef_sp","modisevi_se_sp",
                                   "modisevi_lcl_sp","modisevi_ucl_sp")  
          write_csv(modisevi_tab, paste0("outputs/tables/",sp_name,"_",season,"_modisevi_lm.csv"))    
          # individual slopes plots
          modisevi.ip <- interact_plot(modisevi_lm, pred = year, modx = status, interval=T,
                                   colors=c("gray70","black",col), lty=c(2,2,2)) +
            xlab("Year") +
            ylab("Enhanced Vegetation Index") +
            theme_classic() +
            theme(legend.position="none",
                  axis.text = element_text(size = size),
                  axis.title = element_text(size = size))
          
          # evi (observer only)
          modisevi_lm <- lm(modisevi~year*status, data=spbgdata)
          em <- data.frame(emtrends(modisevi_lm, "status", var="year"))
          modisevi_lm_trends_bg <- em[1,c(2:3,5:6)]
          modisevi_lm_trends_sp <- em[2,c(2:3,5:6)]        
          write_csv(em, paste0("outputs/tables/",sp_name,"_",season,"_modisevi_nobl_lm_em.csv"))     
          modisevi_tab <- tidy(modisevi_lm)
          modisevi_summary_nobl <- c(modisevi_tab[4,c(2,3,5)], modisevi_lm_trends_bg,
                                modisevi_lm_trends_sp)
          names(modisevi_summary_nobl) <- c("modisevi_coef_nobl","modisevi_se_nobl","modisevi_p_nobl",
                                       "modisevi_coef_bg_nobl","modisevi_se_bg_nobl",
                                       "modisevi_lcl_bg_nobl","modisevi_ucl_bg_nobl",
                                       "modisevi_coef_sp_nobl","modisevi_se_sp_nobl",
                                       "modisevi_lcl_sp_nobl","modisevi_ucl_sp_nobl")  
          write_csv(modisevi_tab, paste0("outputs/tables/",sp_name,"_",season,"_modisevi_nobl_lm.csv")) 
          
          # prec (baseline + observer)
          prec_lm <- lm(prec~year*status, data=spblbgdata)
          em <- data.frame(emtrends(prec_lm, "status", var="year"))
          prec_lm_trends_bg <- em[1,c(2:3,5:6)]
          prec_lm_trends_bl <- em[2,c(2:3,5:6)] 
          prec_lm_trends_sp <- em[3,c(2:3,5:6)]     
          write_csv(em, paste0("outputs/tables/",sp_name,"_",season,"_prec_lm_em.csv"))     
          prec_tab <- tidy(prec_lm)
          prec_summary <- c(prec_tab[6,c(2,3,5)], prec_lm_trends_bg,
                                prec_lm_trends_bl, prec_lm_trends_sp)
          names(prec_summary) <- c("prec_coef","prec_se","prec_p",
                                       "prec_coef_bg","prec_se_bg",
                                       "prec_lcl_bg","prec_ucl_bg",
                                       "prec_coef_bl","prec_se_bl",
                                       "prec_lcl_bl","prec_ucl_bl",
                                       "prec_coef_sp","prec_se_sp",
                                       "prec_lcl_sp","prec_ucl_sp")  
          write_csv(prec_tab, paste0("outputs/tables/",sp_name,"_",season,"_prec_lm.csv"))    
          # individual slopes plots
          prec.ip <- interact_plot(prec_lm, pred = year, modx = status, interval=T,
                                       colors=c("gray70","black",col), lty=c(2,2,2)) +
            xlab("Year") +
            ylab("Mean Daily Precipitation (mm)") +
            theme_classic() +
            theme(legend.position="none",
                  axis.text = element_text(size = size),
                  axis.title = element_text(size = size))
          
          # prec (observer only)
          prec_lm <- lm(prec~year*status, data=spbgdata)
          em <- data.frame(emtrends(prec_lm, "status", var="year"))
          prec_lm_trends_bg <- em[1,c(2:3,5:6)]
          prec_lm_trends_sp <- em[2,c(2:3,5:6)]        
          write_csv(em, paste0("outputs/tables/",sp_name,"_",season,"_prec_nobl_lm_em.csv"))     
          prec_tab <- tidy(prec_lm)
          prec_summary_nobl <- c(prec_tab[4,c(2,3,5)], prec_lm_trends_bg,
                                     prec_lm_trends_sp)
          names(prec_summary_nobl) <- c("prec_coef_nobl","prec_se_nobl","prec_p_nobl",
                                            "prec_coef_bg_nobl","prec_se_bg_nobl",
                                            "prec_lcl_bg_nobl","prec_ucl_bg_nobl",
                                            "prec_coef_sp_nobl","prec_se_sp_nobl",
                                            "prec_lcl_sp_nobl","prec_ucl_sp_nobl")  
          write_csv(prec_tab, paste0("outputs/tables/",sp_name,"_",season,"_prec_nobl_lm.csv")) 
          
          # latitude (observer only)
          lat_lm <- lm(latitude~year*status, data=spbgdata)
          em <- data.frame(emtrends(lat_lm, "status", var="year"))
          lat_lm_trends_bg <- em[1,c(2:3,5:6)]
          lat_lm_trends_sp <- em[2,c(2:3,5:6)]
          write_csv(em, paste0("outputs/tables/",sp_name,"_",season,"_lat_lm_em.csv"))     
          lat_tab <- tidy(lat_lm)
          lat_summary <- c(lat_tab[4,c(2,3,5)], lat_lm_trends_bg, lat_lm_trends_sp)
          names(lat_summary) <- c("lat_coef","lat_se","lat_p",
                                   "lat_coef_bg","lat_se_bg",
                                   "lat_lcl_bg","lat_ucl_bg",
                                   "lat_coef_sp","lat_se_sp",
                                   "lat_lcl_sp","lat_ucl_sp")  
          write_csv(lat_tab, paste0("outputs/tables/",sp_name,"_",season,"_lat_lm.csv"))   
          # individual slopes plots
          lat.ip <- interact_plot(lat_lm, pred = year, modx = status, interval=T,
                                       colors=c("gray70",col), lty=c(2,2)) +
            xlab("Year") +
            ylab("Latitude") +
            theme_classic() +
            theme(legend.position="none",
                  axis.text = element_text(size = size),
                  axis.title = element_text(size = size),
                  axis.line = element_line(size=2),
                  axis.ticks = element_line(size=2)) +
            scale_x_continuous(breaks=c(2000,2010,2020))
          
          # elevation (observer only)
          elev_lm <- lm(elev~year*status, data=spbgdata)
          em <- data.frame(emtrends(elev_lm, "status", var="year"))
          elev_lm_trends_bg <- em[1,c(2:3,5:6)]
          elev_lm_trends_sp <- em[2,c(2:3,5:6)]
          write_csv(em, paste0("outputs/tables/",sp_name,"_",season,"_elev_lm_em.csv"))     
          elev_tab <- tidy(elev_lm)
          elev_summary <- c(elev_tab[4,c(2,3,5)], elev_lm_trends_bg, elev_lm_trends_sp)
          names(elev_summary) <- c("elev_coef","elev_se","elev_p",
                                  "elev_coef_bg","elev_se_bg",
                                  "elev_lcl_bg","elev_ucl_bg",
                                  "elev_coef_sp","elev_se_sp",
                                  "elev_lcl_sp","elev_ucl_sp")  
          write_csv(elev_tab, paste0("outputs/tables/",sp_name,"_",season,"_elev_lm.csv"))  
          # individual slopes plots
          elev.ip <- interact_plot(elev_lm, pred = year, modx = status, interval=T,
                                  colors=c("gray70",col), lty=c(2,2)) +
            xlab("Year") +
            ylab("Elevation (m)") +
            theme_classic() +
            theme(legend.position="none",
                  axis.text = element_text(size = size),
                  axis.title = element_text(size = size))
          
          # distance (observer only)
          dist_lm <- lm(dist~year*status, data=spbgdata)
          em <- data.frame(emtrends(dist_lm, "status", var="year"))
          dist_lm_trends_bg <- em[1,c(2:3,5:6)]
          dist_lm_trends_sp <- em[2,c(2:3,5:6)]
          write_csv(em, paste0("outputs/tables/",sp_name,"_",season,"_dist_lm_em.csv"))     
          dist_tab <- tidy(dist_lm)
          dist_summary <- c(dist_tab[4,c(2,3,5)], dist_lm_trends_bg, dist_lm_trends_sp)
          names(dist_summary) <- c("dist_coef","dist_se","dist_p",
                                   "dist_coef_bg","dist_se_bg",
                                   "dist_lcl_bg","dist_ucl_bg",
                                   "dist_coef_sp","dist_se_sp",
                                   "dist_lcl_sp","dist_ucl_sp")  
          write_csv(dist_tab, paste0("outputs/tables/",sp_name,"_",season,"_dist_lm.csv"))  
          # individual slopes plots
          dist.ip <- interact_plot(dist_lm, pred = year, modx = status, interval=T,
                                   colors=c("gray70",col), lty=c(2,2)) +
            xlab("Year") +
            ylab("GC distance (km)") +
            theme_classic() +
            theme(legend.position="none",
                  axis.text = element_text(size = size),
                  axis.title = element_text(size = size))
          
          # longitude (observer only) - need this only for arrow map
          lon_lm <- lm(longitude~year*status, data=spbgdata)
          em <- data.frame(emtrends(lon_lm, "status", var="year"))
          lon_lm_trends_bg <- em[1,c(2:3,5:6)]
          lon_lm_trends_sp <- em[2,c(2:3,5:6)]
          write_csv(em, paste0("outputs/tables/",sp_name,"_",season,"_lon_lm_em.csv"))     
          lon_tab <- tidy(lon_lm)
          lon_summary <- c(lon_tab[4,c(2,3,5)], lon_lm_trends_bg, lon_lm_trends_sp)
          names(lon_summary) <- c("lon_coef","lon_se","lon_p",
                                   "lon_coef_bg","lon_se_bg",
                                   "lon_lcl_bg","lon_ucl_bg",
                                   "lon_coef_sp","lon_se_sp",
                                   "lon_lcl_sp","lon_ucl_sp")  
          write_csv(lon_tab, paste0("outputs/tables/",sp_name,"_",season,"_lon_lm.csv"))  
          
          # save plots
          if(loc=="local"){
            gga <- ggarrange(tmax.ip, modisevi.ip, prec.ip, 
                             lat.ip, elev.ip, dist.ip, ncol=6, nrow=1)
            ggsave(paste0("figs/trends/",gsub(" ","_",common),"_",
                          season,"_trends_full.jpeg"), gga, "jpeg",
                   height=4, width=24, units="in")
            gga <- ggarrange(lat.ip, tmax.ip, ncol=2, nrow=1)
            ggsave(paste0("figs/trends/",gsub(" ","_",common),"_",
                          season,"_trends_two.jpeg"), gga, "jpeg",
                   height=4.5, width=8.5, units="in")
          }else{
            Cairo(paste0("figs/trends/",gsub(" ","_",common),"_",
                         season,"_trends_full.jpeg"), type="jpeg",
                  height=400, width=2200)
            grid.arrange(tmax.ip, modisevi.ip, prec.ip, 
                         lat.ip, elev.ip, dist.ip, ncol=6, nrow=1)
            dev.off()
            Cairo(paste0("figs/trends/",gsub(" ","_",common),"_",
                         season,"_trends_two.jpeg"), type="jpeg",
                  height=400, width=700)
            grid.arrange(lat.ip, tmax.ip, ncol=2, nrow=1)
            dev.off()
          }
          
          # comparing niche/range 2000-2004 vs 2016-2020
          spblbgdata <- spblbgdata[spblbgdata$elev>0,]
          early <- spblbgdata[spblbgdata$year<2005,]
          late <- spblbgdata[spblbgdata$year>2015,]
          
          for (type in 1:3){
            name <- c("background","baseline","focal")
            color <- c("gray70", "black", col)
            xmin <- quantile(spblbgdata$tmax,.001)-3; xmax <- quantile(spblbgdata$tmax,.999)+3
            ymin <- quantile(spblbgdata$modisevi,.001)-.03; ymax <- quantile(spblbgdata$modisevi,.999)+.03
           early_new <- early[early$status==name[type],]
           late_new <- late[late$status==name[type],]
           pts <- data.frame(rbind(c(median(early_new$tmax), median(early_new$modisevi)),
                          c(median(late_new$tmax), median(late_new$modisevi))))
          gg <- ggplot(early_new, aes(x = tmax, y = modisevi)) +
            stat_clip_ellipse_niche(geom="polygon", level=.9, fill=color[type], 
                              color="black", alpha=0.4, size=.8) +
            stat_clip_ellipse_niche(late_new, mapping=aes(x = tmax, y = modisevi),
                              geom="polygon", level=.9, fill=color[type], 
                              color="black", alpha=0.4, size=.8, lty=2) +
            theme_bw() +
            theme(text = element_text(size=18),
                  panel.grid.major = element_line(),
                  panel.grid.minor = element_line(),
                  axis.line = element_line(color = "black"),
                  axis.text.x = element_text(color = "black", size=18),
                  axis.text.y = element_text(color = "black", size=18)) +
            geom_point(data=pts, aes(X1, X2), size=5, col="gray50") +
            xlim(xmin, xmax) +
            ylim(ymin, ymax) +
            xlab("Max. daily temperature (C)") +
            ylab("EVI")
          assign(paste0("gg.niche.", type), gg)
          
          pts <- data.frame(rbind(c(median(early_new$latitude), median(early_new$elev)),
                                  c(median(late_new$latitude), median(late_new$elev))))
          xmin <- min(spblbgdata$latitude)-5; xmax <- max(spblbgdata$latitude)+5
          ymin <- 0; ymax <- quantile(spblbgdata$elev, .98)
          gg <- ggplot(early_new, aes(x = latitude, y = elev)) +
            stat_clip_ellipse_spat(geom="polygon", level=.9, fill=color[type], 
                              color="black", alpha=0.4, size=.8) +
            stat_clip_ellipse_spat(late_new, mapping=aes(x = latitude, y = elev),
                              geom="polygon", level=.9, fill=color[type], 
                              color="black", alpha=0.4, size=.8, lty=2) +
            theme_bw() +
            theme(text = element_text(size=18),
                  panel.grid.major = element_line(),
                  panel.grid.minor = element_line(),
                  axis.line = element_line(color = "black"),
                  axis.text.x = element_text(color = "black", size=18),
                  axis.text.y = element_text(color = "black", size=18)) +
            geom_point(data=pts, aes(X1, X2), size=5, col="gray50") +
            xlim(xmin, xmax) +
            ylim(ymin, ymax) +
            xlab("Latitude") +
            ylab("Elevation (m)")
          assign(paste0("gg.spat.", type), gg)
          }
          if(loc=="local"){
            gga <- ggarrange(gg.niche.3, gg.niche.2, gg.niche.1, 
                             gg.spat.3, gg.spat.2, gg.spat.1, ncol=3, nrow=2)
            ggsave(paste0("figs/trends/",gsub(" ","_",common),"_",season,"_comp.jpeg"), 
                   gga, "jpeg", height=8, width=12, units="in")
          }else{
            Cairo(paste0("figs/trends/",gsub(" ","_",common),"_",season,"_comp.jpeg"),
                   type="jpeg", height=600, width=900)
            grid.arrange(gg.niche.3, gg.niche.2, gg.niche.1, 
                         gg.spat.3, gg.spat.2, gg.spat.1, ncol=3, nrow=2)
            dev.off()
          }

          # model outputs table for species 
          model_out <- data.frame(c(common, sp_name, season, 
                                    tmax_summary, modisevi_summary, prec_summary,
                                    tmax_summary_nobl, modisevi_summary_nobl, prec_summary_nobl,
                                    lat_summary, elev_summary, dist_summary, lon_summary))
          rownames(model_out) <- NULL
          colnames(model_out) <- c("common_name","scientific_name","season",
                                   names(tmax_summary), names(modisevi_summary), 
                                   names(prec_summary), 
                                   names(tmax_summary_nobl), names(modisevi_summary_nobl), 
                                   names(prec_summary_nobl), 
                                   names(lat_summary), names(elev_summary), 
                                   names(dist_summary), names(lon_summary))
          write_csv(model_out, paste0("outputs/models/mod_",sp_name,"_",season,".csv"))
        }
      })
    } # end sp
  } # end season

