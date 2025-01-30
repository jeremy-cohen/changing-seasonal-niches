# Code sheet 2
# Relationships between spatial redistributions and niche shifts across species
# Code is identical across GBIF/ebird/BBS analyses

{library(tidyverse)
  library(sf)
  library(sp)
  library(MVNH)
  library(ggpubr)
  library(lubridate)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(lwgeom)
  library(raster)
  library(ggExtra)
  # code to fix ellipse glitch in ggplot
  source("~/stoat/fix_stat_ellipse.R")
  # working directory
  setwd("~/stoat/bg_points/")}

### PREP
# Trends across slopes (based on all data)
files = list.files("outputs/models/", full.names=T,  recursive=F)
for (i in 1:length(files)){
  file = read_csv(files[i])
  if(i==1){mod = file}else{mod = rbind(mod, file)}
}

# Species list - add column flagging interior west species to subset for elevation
splist = read_csv("~/ranges/species_list.csv") %>%
  dplyr::select(Common, IW)
mod = left_join(mod, splist, by=c("common_name"="Common"))

# Corrected values of niche shifts for observer/baseline
# Calculations correspond to those outlined in methods to adjust for observers/baseline
mod$tmax_obs_corr = mod$tmax_coef_bg - mod$tmax_coef_bl
mod$tmax_sp_corr = mod$tmax_coef_sp - mod$tmax_obs_corr
mod$tmax_sp_mit = mod$tmax_sp_corr - mod$tmax_coef_bl
mod$modisevi_obs_corr = mod$modisevi_coef_bg - mod$modisevi_coef_bl
mod$modisevi_sp_corr = mod$modisevi_coef_sp - mod$modisevi_obs_corr
mod$modisevi_sp_mit = mod$modisevi_sp_corr - mod$modisevi_coef_bl
mod$prec_obs_corr = mod$prec_coef_bg - mod$prec_coef_bl
mod$prec_sp_corr = mod$prec_coef_sp - mod$prec_obs_corr
mod$prec_sp_mit = mod$prec_sp_corr - mod$prec_coef_bl
# write csv
write_csv(mod, "outputs/mod_slopes.csv")

###### Trends across slopes table (based on annual avgs)
mod = read_csv("outputs/mod_slopes.csv")
# change distance to abs value
mod$dist_coef = abs(mod$dist_coef)

### Species-level niche/range shift tables
t1 = dplyr::select(mod, c(scientific_name, common_name, season,
                           tmax_sp_corr, tmax_sp_mit,
                           modisevi_sp_corr, modisevi_sp_mit,
                           prec_sp_corr, prec_sp_mit,
                           lat_coef, elev_coef, dist_coef, mahal))
# absolute value of distance moved (negative values not relevant)
t1$dist_coef = abs(t1$dist_coef)
# multiply per-year slopes by 20 to get change over 20 years
t1 = cbind(t1[,1:3],round(t1[4:13]*20,2))
# column names, compile seasonal tables
colnames(t1) = c("Species","Common name","Season","Temperature",
                  "Mitigated Temperature","EVI","Mitigated EVI",
                  "Precipitation","Mitigated Precipitation",
                  "Latitude","Elevation","Distance","Mahalanobis distance")
br = t1[t1$Season=="breeding",]
br$Season=NULL
write_csv(br, "outputs/sp_table_breeding.csv")
nbr = t1[t1$Season=="nonbreeding",]
nbr$Season=NULL
write_csv(nbr, "outputs/sp_table_nonbreeding.csv")


### Density plots and means/ses for niche/range shifts
# loop seasons
for (season in c("breeding","nonbreeding")){
  mods = mod[mod$season==season,]
# What's the distribution of background trends?
  # Loop 3 niche dimensions and 3 spatial dimensions
for (i in c("lat","elev","dist","tmax","modisevi","prec")){
  # per-year change to 20 year change adjustment
  if(i %in% c("lat","elev","dist")){
mod$dens = mod[,paste0(i,"_coef_bg")][[1]] *20
  }else{
    mod$dens = mod[,paste0(i,"_obs_corr")][[1]] *20}
if(i=="lat"){xlab="Δ Latitude (°)"}
if(i=="elev"){xlab="Δ Elevation (m)"}
if(i=="dist"){xlab="Δ Distance (km)";
mod$dens = abs(mod$dens)}
if(i=="tmax"){xlab="Δ Max. temperature (°C)"}
if(i=="modisevi"){xlab="Δ EVI"}
  if(i=="prec"){xlab="Δ Daily precipitation (mm)"}
  # get medians
  summer_med = median(mod[mod$season=="breeding",]$dens)
  winter_med = median(mod[mod$season=="nonbreeding",]$dens)
  # density plots
gg = ggplot(mod, aes(x=dens, fill=season)) +
  geom_density(alpha=.5, size=1) +
  geom_vline(xintercept=0, color="gray", linetype="dashed", cex=.7) +
  geom_vline(xintercept=summer_med, color="darkorange1", linetype="dashed", cex=1) +
  geom_vline(xintercept=winter_med, color="dodgerblue", linetype="dashed", cex=1) +
  scale_fill_manual(values = c("darkorange1","dodgerblue")) +
  theme(legend.position="none",
        axis.text.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black", size=16),
        axis.ticks.y=element_blank()) + 
  theme(text = element_text(size=20)) +
  xlab(xlab) +
  ylab("Density")
assign(paste0("gg.",i), gg)
}
  # arrange and save plots
gga = ggarrange(gg.lat, gg.elev, gg.dist, 
                 gg.tmax, gg.modisevi, gg.prec, ncol=6, nrow=1)
ggsave("figs/trends/spatial_sampling_trends.jpeg", 
       gga, "jpeg", height=4, width=24, units="in")

# What's the distribution of baseline trends?
for (i in c("tmax","modisevi","prec")){
  mod$dens = mod[,paste0(i,"_coef_bl")][[1]] *20
  if(i=="tmax"){xlab="Δ Max. temperature (°C)"}
  if(i=="modisevi"){xlab="Δ EVI"}
  if(i=="prec"){xlab="Δ Daily precipitation (mm)"}
  summer_med = median(mod[mod$season=="breeding",]$dens)
  winter_med = median(mod[mod$season=="nonbreeding",]$dens)
  gg = ggplot(mod, aes(x=dens, fill=season)) +
    geom_density(alpha=.5, size=1) +
    geom_vline(xintercept=0, color="gray", linetype="dashed", cex=.7) +
    geom_vline(xintercept=summer_med, color="darkorange1", linetype="dashed", cex=1) +
    geom_vline(xintercept=winter_med, color="dodgerblue", linetype="dashed", cex=1) +
    scale_fill_manual(values = c("darkorange1","dodgerblue")) +
    theme(legend.position="none",
          axis.text.y=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(color = "black"),
          axis.text.x = element_text(color = "black", size=16),
          axis.ticks.y=element_blank()) + 
    theme(text = element_text(size=20)) +
    xlab(xlab) +
    ylab("Density")
  assign(paste0("gg.",i), gg)
}
gga = ggarrange(gg.tmax, gg.modisevi, gg.prec, ncol=3, nrow=1)
ggsave("figs/trends/baseline_trends.jpeg", 
       gga, "jpeg", height=4, width=12, units="in")

# What's the distribution of adjusted species shifts?
for (i in c("lat","elev","dist","tmax","modisevi","prec")){
  if(i %in% c("lat","elev","dist")){
    mod$dens = mod[,paste0(i,"_coef")][[1]] *20
  }else{
    mod$dens = mod[,paste0(i,"_sp_corr")][[1]] *20}
  if(i=="lat"){xlab="Δ Latitude (°)"}
  if(i=="elev"){xlab="Δ Elevation (m)"}
  if(i=="dist"){xlab="Δ Distance (km)";
  mod$dens = abs(mod$dens)}
  if(i=="tmax"){xlab="Δ Max. temperature (°C)"}
  if(i=="modisevi"){xlab="Δ EVI"}
  if(i=="prec"){xlab="Δ Daily precipitation (mm)"}
  summer_med = median(mod[mod$season=="breeding",]$dens)
  winter_med = median(mod[mod$season=="nonbreeding",]$dens)
  gg = ggplot(mod, aes(x=dens, fill=season)) +
    geom_density(alpha=.5, size=1) +
    geom_vline(xintercept=0, color="gray", linetype="dashed", cex=.7) +
    geom_vline(xintercept=summer_med, color="darkorange1", linetype="dashed", cex=1) +
    geom_vline(xintercept=winter_med, color="dodgerblue", linetype="dashed", cex=1) +
    scale_fill_manual(values = c("darkorange1","dodgerblue")) +
    theme(legend.position="none",
          axis.text.y=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(color = "black"),
          axis.text.x = element_text(color = "black", size=16),
          axis.ticks.y=element_blank()) + 
    theme(text = element_text(size=20)) +
    xlab(xlab) +
    ylab("Density")
  assign(paste0("gg.",i), gg)
}
gga = ggarrange(gg.lat, gg.elev, gg.dist, gg.tmax, gg.modisevi, gg.prec, ncol=6, nrow=1)
ggsave("figs/trends/species_trends.jpeg", 
       gga, "jpeg", height=4, width=24, units="in")

# What's the distribution of mitigated species niche shifts?
for (i in c("tmax","modisevi","prec")){
    mod$dens = mod[,paste0(i,"_sp_mit")][[1]] *20
    if(i=="tmax"){xlab="Δ Max. temperature (°C)"}
    if(i=="modisevi"){xlab="Δ EVI"}
    if(i=="prec"){xlab="Δ Daily precipitation (mm)"}
    summer_med = median(mod[mod$season=="breeding",]$dens)
    winter_med = median(mod[mod$season=="nonbreeding",]$dens)
  gg = ggplot(mod, aes(x=dens, fill=season)) +
    geom_density(alpha=.5, size=1) +
    geom_vline(xintercept=0, color="gray", linetype="dashed", cex=.7) +
    geom_vline(xintercept=summer_med, color="darkorange1", linetype="dashed", cex=1) +
    geom_vline(xintercept=winter_med, color="dodgerblue", linetype="dashed", cex=1) +
    scale_fill_manual(values = c("darkorange1","dodgerblue")) +
    theme(legend.position="none",
          axis.text.y=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(color = "black"),
          axis.text.x = element_text(color = "black", size=16),
          axis.ticks.y=element_blank()) + 
    theme(text = element_text(size=20)) +
    xlab(xlab) +
    ylab("Density")
  assign(paste0("gg.",i), gg)
}
gga = ggarrange(gg.tmax, gg.modisevi, gg.prec, ncol=3, nrow=1)
ggsave("figs/trends/species_mit_trends.jpeg", 
       gga, "jpeg", height=4, width=12, units="in")

} # end season

### create table of mean/se niche and range shift across species
output <- matrix(nrow=0, ncol=5)
# loop seasons
for (season in c("breeding","nonbreeding")){
  # subset season
  mods <- mod[mod$season==season,]
  # What's the distribution of backgrounds?
  for (i in c("lat","elev","dist","tmax","modisevi","prec")){
    if(i %in% c("lat","elev","dist")){
      # per-year to 20 year trend adjustment
      mods$dens <- mods[,paste0(i,"_coef_bg")][[1]] *20
      if(i=="dist"){mods$dens <- abs(mods$dens)}
    }else{
      mods$dens <- mods[,paste0(i,"_obs_corr")][[1]] *20}
    out <- c(season,i,"observer", median(mods$dens),
             sd(mods$dens)/sqrt(nrow(mods)))
    output <- rbind(output, out)
  }
  # What's the distribution of baselines?
  for (i in c("tmax","modisevi","prec")){
    mods$dens <- mods[,paste0(i,"_coef_bl")][[1]] *20
    out <- c(season,i,"baseline", median(mods$dens),
             sd(mods$dens)/sqrt(nrow(mods)))
    output <- rbind(output, out)
  }
  # What's the distribution of adjusted species shifts?
  for (i in c("lat","elev","dist","tmax","modisevi","prec")){
    if(i %in% c("lat","elev","dist")){
      mods$dens <- mods[,paste0(i,"_coef")][[1]] *20
      if(i=="dist"){mods$dens <- abs(mods$dens)}
    }else{
      mods$dens <- mods[,paste0(i,"_sp_corr")][[1]] *20}
    out <- c(season,i,"realized", median(mods$dens), 
             sd(mods$dens)/sqrt(nrow(mods)))
    output <- rbind(output, out)
  }
  # What's the distribution of mitigated species niche shifts?
  for (i in c("tmax","modisevi","prec")){
    mods$dens <- mods[,paste0(i,"_sp_mit")][[1]] *20
    out <- c(season,i,"mitigated", median(mods$dens), 
             sd(mods$dens)/sqrt(nrow(mods)))
    output <- rbind(output, out)
  }
  # Add in mahalanobis distance in n-dimensional space ('realized' only)
  out = c(season,"mahal","realized", median(mods$mahal, na.rm=T),
          sd(mods$mahal, na.rm=T)/sqrt(nrow(mods)))
  output <- rbind(output, out)
}
# name columns and compile table
output = data.frame(output)
rownames(output) = NULL
colnames(output) = c("season","variable","dataset","mean","SE")
output$variable = plyr::revalue(
  output$variable, c("tmax"="Max. Temperature", "lat"="Latitude",
                     "elev"="Elevation","prec"="Precipitation",
                     "modisevi"="EVI","dist"="Distance",
                     "mahal"="Mahalanobis Distance"))
output = cbind(output[,1:3],round(as.numeric(output[,4]),3),
                round(as.numeric(output[,5]),3))
output = merge(output[output$season=="breeding",],
                output[output$season=="nonbreeding",],
                by=c("variable","dataset"))
output$season.x = NULL; output$season.y = NULL
colnames(output) = c("Variable","Trend","Mean","SE","Mean","SE")
write_csv(output, "outputs/means_table.csv")


##############################################

##### Niche vs space plots (Fig 3; supplementary figures)

### space vs mitigated niche loss
# loop seasons
  for (season in c("breeding","nonbreeding")){
    # subset by season
    sub_data = mod[mod$season==season,]
    # colors
    if(season=="breeding"){col="darkorange2"}else{col="dodgerblue"}
    # loop combos of space and niche shifts
    for (spvar in c("lat","dist","elev")){
    for (var in c("tmax","modisevi","prec")){
      # adjustments and labels
        if(spvar=="elev"){sub_data = sub_data[sub_data$IW=="yes",]}
      sub_data$spat = sub_data[,paste0(spvar,"_coef")][[1]] *20
      sub_data$env = sub_data[,paste0(var,"_sp_mit")][[1]] *20
      if(var=="tmax"){ylab="Mitigated Δ Max. temperature (°C)"}
      if(var=="modisevi"){ylab="Mitigated Δ EVI"}
      if(var=="prec"){ylab="Mitigated Δ Daily precipitation (mm)"}
      if(season=="breeding"){xlab=""}else{
        if(spvar=="lat"){xlab="Δ Latitude (°)"}
        if(spvar=="elev"){xlab="Δ Elevation (m)"}
        if(spvar=="dist"){xlab="Δ Distance (km)"}
      }
      # color gradient for points based on closeness to no niche shift
      low=col; mid="white"; high="gray30"
      # plot across species
      gg = ggplot(sub_data, aes(spat, env)) +
        geom_point(aes(col=env), size=4, pch=21, stroke=1.5, fill="black") +
        theme(text = element_text(size=15),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(color = "black"),
              axis.text.y = element_text(color = "black", size=15),
              axis.text.x = element_text(color = "black", size=15),
              legend.position="none") +
        geom_smooth(color="black", method="lm") +
        geom_hline(yintercept=0, lty=2) +
        geom_vline(xintercept=0, lty=2) +
        scale_colour_gradient2(low=low, mid=mid, high=high, midpoint=0) +
        xlab(xlab) +
        ylab(ylab)
      gg = ggMarginal(gg, type="density", size=10, fill=col, alpha=.3)
      assign(paste0("dif.",var,".",spvar,".",season), gg)
    }
  }
  }

# space vs realized niche loss
for (season in c("breeding","nonbreeding")){
  sub_data = mod[mod$season==season,]
  if(season=="breeding"){col="darkorange2"}else{col="dodgerblue"}
  for (spvar in c("lat","dist","elev")){
  for (var in c("tmax","modisevi","prec")){
      if(spvar=="elev"){sub_data = sub_data[sub_data$IW=="yes",]}
        sub_data$spat = sub_data[,paste0(spvar,"_coef")][[1]] *20
      sub_data$env = sub_data[,paste0(var,"_sp_corr")][[1]] *20 
if(var=="tmax"){ylab="Realized Δ Max. temperature (°C)"}
      if(var=="modisevi"){ylab="Realized Δ EVI"}
      if(var=="prec"){ylab="Realized Δ Daily precipitation (mm)"}
      if(season=="breeding"){xlab=""}else{
        if(spvar=="lat"){xlab="Δ Latitude (°)"}
        if(spvar=="elev"){xlab="Δ Elevation (m)"}
        if(spvar=="dist"){xlab="Δ Distance (km)"}
      }
      low="white"; mid=col; high="white"
    gg = ggplot(sub_data, aes(spat,env)) +
      geom_point(aes(fill=env), size=4, pch=21, col="black") +
      geom_smooth(color="black", method="lm") +
      theme(text = element_text(size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(color = "black"),
            axis.text.y = element_text(color = "black", size=15),
            axis.text.x = element_text(color = "black", size=15),
            legend.position="none") +
      geom_hline(yintercept=0, lty=2) +
      geom_vline(xintercept=0, lty=2) +
      scale_fill_gradient2(low=low, mid=mid, high=high, 
                           midpoint=0) +
      xlab(xlab) +
      ylab(ylab)
    gg = ggMarginal(gg, type="density", size=10, fill=col, alpha=.3)
    assign(paste0("sp.",var,".",spvar,".",season), gg)
  }
}
}

# arrange, save all plots for all combos of niche/space shifts
{gga = ggarrange(sp.tmax.lat.breeding, dif.tmax.lat.breeding, 
                  sp.tmax.lat.nonbreeding, dif.tmax.lat.nonbreeding)
ggsave("figs/trends/cross_sp_trends_tmax_lat.jpeg", 
       gga, "jpeg", height=8, width=10, units="in")
gga = ggarrange(sp.tmax.elev.breeding, dif.tmax.elev.breeding, 
                 sp.tmax.elev.nonbreeding, dif.tmax.elev.nonbreeding,
                 labels=c("a)","b)","c)","d)"))
ggsave("figs/trends/cross_sp_trends_tmax_elev.jpeg", 
       gga, "jpeg", height=8, width=10, units="in")
gga = ggarrange(sp.tmax.dist.breeding, dif.tmax.dist.breeding,
                 sp.tmax.dist.nonbreeding, dif.tmax.dist.nonbreeding,
                 labels=c("a)","b)","c)","d)"))
ggsave("figs/trends/cross_sp_trends_tmax_dist.jpeg", 
       gga, "jpeg", height=8, width=10, units="in")

gga = ggarrange(sp.modisevi.lat.breeding, dif.modisevi.lat.breeding, 
                 sp.modisevi.lat.nonbreeding, dif.modisevi.lat.nonbreeding,
                 labels=c("a)","b)","c)","d)"))
ggsave("figs/trends/cross_sp_trends_modisevi_lat.jpeg", 
       gga, "jpeg", height=8, width=10, units="in")
gga = ggarrange(sp.modisevi.elev.breeding, dif.modisevi.elev.breeding, 
                 sp.modisevi.elev.nonbreeding, dif.modisevi.elev.nonbreeding,
                 labels=c("a)","b)","c)","d)"))
ggsave("figs/trends/cross_sp_trends_modisevi_elev.jpeg", 
       gga, "jpeg", height=8, width=10, units="in")
gga = ggarrange(sp.modisevi.dist.breeding, dif.modisevi.dist.breeding, 
                 sp.modisevi.dist.nonbreeding, dif.modisevi.dist.nonbreeding,
                 labels=c("a)","b)","c)","d)"))
ggsave("figs/trends/cross_sp_trends_modisevi_dist.jpeg", 
       gga, "jpeg", height=8, width=10, units="in")

gga = ggarrange(sp.prec.lat.breeding, dif.prec.lat.breeding,
                 sp.prec.lat.nonbreeding, dif.prec.lat.nonbreeding,
                 labels=c("a)","b)","c)","d)"))
ggsave("figs/trends/cross_sp_trends_prec_lat.jpeg", 
       gga, "jpeg", height=8, width=10, units="in")
gga = ggarrange(sp.prec.elev.breeding, dif.prec.elev.breeding, 
                 sp.prec.elev.nonbreeding, dif.prec.elev.nonbreeding,
                 labels=c("a)","b)","c)","d)"))
ggsave("figs/trends/cross_sp_trends_prec_elev.jpeg", 
       gga, "jpeg", height=8, width=10, units="in")
gga = ggarrange(sp.prec.dist.breeding, dif.prec.dist.breeding, 
                 sp.prec.dist.nonbreeding, dif.prec.dist.nonbreeding,
                 labels=c("a)","b)","c)","d)"))
ggsave("figs/trends/cross_sp_trends_prec_dist.jpeg", 
       gga, "jpeg", height=8, width=10, units="in")
}


### space vs n-dimensional niche loss
mod2 = mod[mod$mahal<1,] # cut 2 species with astronomically high values (all others between 0-1)
# seasonal loop
for (season in c("breeding","nonbreeding")){
  # subset by season
  sub_data = mod2[mod2$season==season,]
  if(season=="breeding"){col="darkorange2"; sub_data = sub_data[sub_data$mahal<.7,]}else{col="dodgerblue"}
  for (spvar in c("lat","dist","elev")){
      if(spvar=="elev"){sub_data = sub_data[sub_data$IW=="yes",]}
      sub_data$spat = sub_data[,paste0(spvar,"_coef")][[1]] *20
      if(season=="breeding"){xlab=""}else{
        if(spvar=="lat"){xlab="Δ Latitude (°)"}
        if(spvar=="elev"){xlab="Δ Elevation (m)"}
        if(spvar=="dist"){xlab="Δ Distance (km)"}
      }
      low="white"; mid=col; high="white"
      # plot
      gg = ggplot(sub_data, aes(spat,mahal)) +
        geom_point(aes(fill=mahal), size=4, pch=21, col="black") +
        geom_smooth(color="black", method="lm") +
        theme(text = element_text(size=15),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(color = "black"),
              axis.text.y = element_text(color = "black", size=15),
              axis.text.x = element_text(color = "black", size=15),
              legend.position="none") +
        geom_hline(yintercept=0, lty=2) +
        geom_vline(xintercept=0, lty=2) +
        scale_fill_gradient2(low=low, mid=mid, high=high, midpoint=0) +
        xlab(xlab) +
        ylab("Mahalanobis distance")
      gg = ggMarginal(gg, type="density", size=10, fill=col, alpha=.3)
      assign(paste0("sp.mahal.",spvar,".",season), gg)
    }
}
# arrange, save plots
{gga = ggarrange(sp.mahal.lat.breeding, sp.mahal.elev.breeding, 
                 sp.mahal.dist.breeding, sp.mahal.lat.nonbreeding, 
                 sp.mahal.elev.nonbreeding, sp.mahal.dist.nonbreeding, 
                 nrow=2, ncol=3, labels=c("a)","b)","c)","d)","e)","f)"))
ggsave("figs/trends/cross_sp_trends_mahal.jpeg", 
       gga, "jpeg", height=8, width=15, units="in")
}


### baseline vs species comparison - fig 2
# loop seasons
for (season in c("breeding","nonbreeding")){
  # subset seasonally
  sub_data = mod[mod$season==season,]
  # colors and labels
  if(season=="breeding"){col="darkorange2"}else{col="dodgerblue"}
  for (var in c("tmax","modisevi","prec")){ 
    if(var=="tmax"){ylab="Realized Δ Max. temperature (°C)"
    xlab="Δ Max. temperature exposure (°C)"
    low=col; mid="white"; high="gray30"; xlim=c(0,8); ylim=c(-3,9)}
    if(var=="modisevi"){ylab="Realized Δ EVI"
    xlab="Δ EVI exposure"
    low="white"; mid="white"; high="white"; xlim=c(-.02,.05); ylim=c(-.11,.08)}
    if(var=="prec"){ylab="Realized Δ Daily precipitation (mm)"
    xlab="Δ Daily precipitation exposure (mm)"
    low="white"; mid="white"; high="white"; xlim=c(-50,500); ylim=c(-200,650)}
    if(season=="breeding"){xlab=""}
    # adjustments
      sub_data$bl = sub_data[,paste0(var,"_coef_bl")][[1]] *20
      sub_data$sp = sub_data[,paste0(var,"_sp_corr")][[1]] *20
      sub_data$dif = sub_data$sp - sub_data$bl
      # plots
      gg = ggplot(sub_data, aes(bl, sp)) +
        # color by distance from 1:1 line
        geom_point(aes(fill=dif), col="black", size=4, pch=21,
                   show.legend=F) +
        scale_fill_gradient2(low=low, mid=mid, high=high, 
                             midpoint=0) +
        theme(text = element_text(size=16),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(color = "black"),
              axis.text.y = element_text(color = "black", size=16),
              axis.text.x = element_text(color = "black", size=16)) +
        geom_abline(slope=1, lty=3) +
        geom_smooth(color="black", method="lm") +
        geom_hline(yintercept=0, lty=2) +
        geom_vline(xintercept=0, lty=2) +
        xlim(xlim) +
        ylim(ylim) +
        ylab(ylab) +
        xlab(xlab)
      gg = ggMarginal(gg, type="density", size=10, fill=col, alpha=.3)
      assign(paste0("bl.",var,".",season), gg)
      # model
      lm = lm(sp~bl, data=sub_data)
      r2 = summary(lm)$r.squared
      print(paste(season, var, "realized"))
      print(lm[1])
      print(r2)
  }
  
  # baseline vs mitigated comparison - fig 2
  for (var in c("tmax","modisevi","prec")){ 
    if(var=="tmax"){ylab="Mitigated Δ Max. temperature (°C)"
    xlab="Δ Max. temperature exposure (°C)"
    low=col; mid="white"; high="gray30"; xlim=c(0,8); ylim=c(-3,9)}
    if(var=="modisevi"){ylab="Mitigated Δ EVI"
    xlab="Δ EVI exposure"
    low="white"; mid="white"; high="white"; xlim=c(-.02,.05); ylim=c(-.11,.08)}
    if(var=="prec"){ylab="Mitigated Δ Daily precipitation (mm)"
    xlab="Δ Daily precipitation exposure (mm)"
    low="white"; mid="white"; high="white"; xlim=c(-50,500); ylim=c(-200,650)}
    if(season=="breeding"){xlab=""}
    sub_data$bl = sub_data[,paste0(var,"_coef_bl")][[1]] *20
    sub_data$sp = sub_data[,paste0(var,"_sp_mit")][[1]] *20
    sub_data$dif = sub_data$sp - sub_data$bl
    gg = ggplot(sub_data, aes(bl, sp)) +
      # color by distance from 1:1 line
      geom_point(aes(col=sp), fill="black", size=4, pch=21, stroke=1.5, 
                 show.legend=F) +
      scale_colour_gradient2(low=low, mid=mid, high=high, 
                             midpoint=0) +
      theme(text = element_text(size=16),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(color = "black"),
            axis.text.y = element_text(color = "black", size=16),
            axis.text.x = element_text(color = "black", size=16)) +
      geom_hline(yintercept=0, lty=2) +
      geom_vline(xintercept=0, lty=2) +
      xlim(xlim) +
      ylim(ylim) +
      ylab(ylab) +
      xlab(xlab)
    gg = ggMarginal(gg, type="density", size=10, fill=col, alpha=.3)
    assign(paste0("blmit.",var,".",season), gg)
    # model
    lm = lm(sp~bl, data=sub_data)
    r2 = summary(lm)$r.squared
    print(paste(season, var, "mitigated"))
    print(lm[1])
    print(r2)
  }
}
gga = ggarrange(bl.tmax.breeding, blmit.tmax.breeding, 
                 bl.tmax.nonbreeding, blmit.tmax.nonbreeding,
                 ncol=2, nrow=2)
ggsave("figs/trends/cross_sp_trends_bl_v_sp7.jpeg", 
       gga, "jpeg", height=8, width=10, units="in")
gga = ggarrange(bl.modisevi.breeding, bl.prec.breeding,
                 bl.modisevi.nonbreeding, bl.prec.nonbreeding,
                 ncol=2, nrow=2,
                 labels=c("a)","c)","b)","d)"))
ggsave("figs/trends/cross_sp_trends_bl_v_sp_figs1v3.jpeg", 
       gga, "jpeg", height=10, width=12, units="in")


### model relationships between niche and space shift combinations (table s2)
nmetrics = c("tmax","modisevi","prec")
smetrics = c("lat","elev","dist")
output = matrix(nrow=0, ncol=6)
# loop combinations of niche and spatial metrics
for (j in 1:length(smetrics)){
  for (i in 1:length(nmetrics)){
    # assign spatial variable
    svar = mod[,paste0(smetrics[j],"_coef")][[1]]
    if(smetrics[j]=="dist"){svar=abs(svar)} # abs value for distance
    # assign realized and mitigated niche variable
    nvar = mod[,paste0(nmetrics[i],"_sp_corr")][[1]]
    mitvar = mod[,paste0(nmetrics[i],"_sp_mit")][[1]]
    # lm and r squared of realized niche shift/spatial relationship
    lm = lm(svar ~ nvar + mod$season)
    summ1 = c(nmetrics[i], smetrics[j], 
              summary(lm)$coef[2,c(1,2,4)],summary(lm)$r.squared)
    # lm and r squared of mitigated niche/spatial relationship
    lm = lm(svar ~ mitvar + mod$season)
    summ2 = c(paste("Mitagated",nmetrics[i]), smetrics[j], 
              summary(lm)$coef[2,c(1,2,4)],summary(lm)$r.squared)
    # compile
    output = rbind(output, summ1, summ2)
  }
  # mahalanobis distance
  lm = lm(svar ~ mod$mahal + mod$season)
  summ = c("mahal", smetrics[j], 
           summary(lm)$coef[2,c(1,2,4)],summary(lm)$r.squared)
  output = rbind(output, summ)
  
  # seasonal version
  ssvar = mod[mod$season=="breeding",][,paste0(smetrics[j],"_coef")][[1]]
  slm = lm(ssvar ~ mod[mod$season=="breeding",]$mahal)
  summary(slm)
  wsvar = mod[mod$season=="nonbreeding",][,paste0(smetrics[j],"_coef")][[1]]
  wlm = lm(wsvar ~ mod[mod$season=="nonbreeding",]$mahal)
  summary(wlm)
}


# compile coefs and r squared
output = cbind(output[,1:2],round(as.numeric(output[,3]),3),
               round(as.numeric(output[,4]),3),
               round(as.numeric(output[,5]),3),
               round(as.numeric(output[,6]),3))
colnames(output) = c("Niche variable","Spatial variable","Coefficient",
                     "SD","p-value","R2-value")
write_csv(data.frame(output), "outputs/niche_v_sp_models.csv")





##### Shift maps
### Arrow figure (figure 3a-b)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)

# prep data for use in arrow maps (do once)
for (season in c("breeding","nonbreeding")){ 
  # subset data to relevant variables
  arrow = mod[mod$season==season,c("common_name","scientific_name","season",
                  "lat_coef","lon_coef","tmax_sp_corr","tmax_sp_mit")]
  arrow$sci_name = gsub("_"," ",arrow$scientific_name)
  alldata = read_csv(paste0("data/",season,"_processed_20yr.csv"))
  # push december events to next year so years are contiguous
  alldata$month = month(alldata$eventdate)
  alldata$year = ifelse(alldata$month==12 & alldata$year>2009,alldata$year+1,alldata$year)
  # 2000-2004 data
early = alldata[alldata$year<2005,c("sci_name","season","longitude","latitude")] %>%
  group_by(sci_name, season) %>%
  dplyr::summarize(meanlon = mean(longitude), meanlat = mean(latitude)) %>%
  ungroup()
# join data, make adjustments
arrow = left_join(arrow, early, by=c("sci_name","season"))
arrow$latelat = arrow$meanlat + (arrow$lat_coef*20)
arrow$latelon = arrow$meanlon + (arrow$lon_coef*20)
arrow$tmax_sp_corr_adj = arrow$tmax_sp_corr *20
# write data
write_csv(arrow, paste0("outputs/arrow_data_",season,".csv"))
rm(alldata); gc()
}

# map preliminary - get polygons, set theme, palette
theme_set(theme_bw())
world = ne_countries(country=c("united states of america","canada","mexico"),
                      scale = "medium", returnclass = "sf") %>%
  st_crop(extent(-130,-70,20,55))
myPalette = rev(brewer.pal(11, "Spectral"))

# loop seasons
for (season in c("breeding","nonbreeding")){ # season = "nonbreeding"
  # load data, subset to most relevant area of map
 arrow = read_csv(paste0("outputs/arrow_data_",season,".csv"))
 arrow = arrow[arrow$latelat<(55) &
                  arrow$latelat>(23) &
                  arrow$latelon<(-70) &
                  arrow$latelon>(-127),]
 # calculate change between mean and 2016-2020 latitude
 arrow$changelat = arrow$latelat-arrow$meanlat
 if(season == "breeding"){low="#3288bd"; mid="#fee08b"; high="#d53e4f"}else{
   low="#ea9999"; mid="#eeeeee"; high="#0076bf"}
 
 # map for all species, arrows give distance covered, colors give max. temp shift
 map1 = ggplot(data = world) +
   geom_sf(alpha=.5) +
   geom_segment(data = arrow, 
                aes(x = meanlon, xend = latelon, y = meanlat, yend = latelat,
                    col=tmax_sp_corr_adj), size=2,  
                arrow = arrow(length = unit(0.02, "npc"))) + 
   geom_point(data = arrow, aes(x = meanlon, y = meanlat),
              col="black", size=2) +
   scale_colour_gradient2(low = low, mid = mid, high = high,
                          midpoint = 0, limits=c(-3,9.2), name="Realized Δ\nMax. temp.\n(°C)") +
   theme(axis.text=element_blank(),
         axis.ticks=element_blank(),
         legend.position="none") +
   xlab("") +
   ylab("")
 assign(paste0("map1.",season), map1)
 
 # # legend
 # ggsave("figs/trends/shift_map_legend.jpeg", legend, height=6, width=8)
}

# save maps
gga = ggarrange(map1.breeding, 
                 map1.nonbreeding, 
                 nrow=2, ncol=1,
                 font.label = list(size = 20), hjust=-1.5)
ggsave("figs/trends/shift_map7.jpeg", gga, height=10, width=6.5)


