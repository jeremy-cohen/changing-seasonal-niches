# code sheet 3
# functional and phylogenetic trends across species

library(ggplot2)
library(ggpubr)
library(tidyverse)
library(ggExtra)
library(MuMIn)
library(ape)
library(geiger)
library(nlme)
library(phytools)
library(jtools)
#library(treeman)
library(picante)
library(glmnet)
library(visreg)
library(car)
library(ggeffects)
library(pdp)
library(phylolm)

# Import, name species
list = read_csv("~/ranges/species_list.csv")
list$Scientific = gsub(" ","_",list$Scientific)
# Directories
setwd("~/stoat/bg_points/")
dir.create("outputs/cross_sp/")

# Outputs from niche estimation
sp.names = read_csv("~/stoat/bg_points/outputs/mod_slopes.csv")

# Frank's migration data
migration = read.delim("~/stoat/traits/migration.txt")
colnames(migration)[1] = "species"
# cut species outside US
migration = migration[migration$MAX.LAT.b > 27 &
                        migration$LON.b < (-20),]
# resolve names to Clements
sp.names$scientific_name[which(!(sp.names$scientific_name %in% migration$species))]
migration$species = plyr::revalue(
  migration$species,
  c("Calonectris_borealis" = "Calonectris_diomedea",
    "Catharus_swainsoni" = "Catharus_ustulatus",
    "Cistothorus_stellaris" = "Cistothorus_platensis",
    "Cyanecula_svecica" = "Luscinia svecica",
    "Glaucidium_californicum" = "Glaucidium_gnoma",
    "Hesperiphona_vespertina" = "Coccothraustes_vespertinus",
    "Hydrobates_castro" = "Oceanodroma_castro",
    "Hydrobates_furcatus" = "Oceanodroma_furcata",
    "Hydrobates_homochroa" = "Oceanodroma_homochroa",
    "Hydrobates_leucorhous" = "Oceanodroma_leucorhoa",
    "Hydrobates_melania" = "Oceanodroma_melania",
    "Hylatomus_pileatus" = "Dryocopus_pileatus",
    "Oreothlypis_ruficapilla" = "Leiothlypis_ruficapilla",
    "Leuconotopicus_albolarvatus" = "Dryobates_albolarvatus",
    "Leuconotopicus_arizonae" = "Dryobates_arizonae",
    "Leuconotopicus_borealis" = "Dryobates_borealis",
    "Leuconotopicus_albolarvatus" = "Dryobates_albolarvatus",
    "Pachyramphus_uropygialis" = "Melanerpes_uropygialis",
    "Passerculus_bairdii" = "Centronyx_bairdii",
    "Passerculus_henslowii" = "Centronyx_henslowii",
    "Passerella_arborea" = "Spizelloides_arborea",
    "Pica_nutalli" = "Pica_nuttalli",
    "Piranga_hepatica" = "Piranga_flava",
    "Porphyrio_martinicus" = "Porphyrio_martinica",
    "Steganopus_tricolor" = "Phalaropus_tricolor",
    "Sterna_hirundo" = "Sterna_hirundo"))
sp.names = left_join(sp.names, migration, by=c("scientific_name"="species"))

# AVONET
avonet = read_csv("~/avonet/TraitData/AVONET2_eBird.csv") %>%
  dplyr::select(Species2, Mass, Trophic.Niche, Habitat, `Hand-Wing.Index`)
avonet$Species2 = gsub(" ","_",avonet$Species2)

which = which(!(sp.names$scientific_name %in% avonet$Species2))
sp.names$scientific_name[which]

avonet$Species2 = plyr::revalue(
  avonet$Species2,
  c("Caracara_plancus" = "Caracara_cheriway",
    "Canachites_canadensis" = "Falcipennis_canadensis",
    "Nannopterum_auritum" = "Phalacrocorax_auritus",
    "Nannopterum_brasilianum" = "Phalacrocorax_brasilianus" ,
    "Urile_pelagicus" =  "Phalacrocorax_pelagicus",
    "Corthylio_calendula" = "Regulus_calendula"))
sp.names = left_join(sp.names, avonet, by=c("scientific_name"="Species2"))

# Landscape Diversity Index (have to do it seasonally)
ldi_table = read_csv("~/ranges/birds_vs1_sp_level_corrected_with_ldi.csv") %>%
  dplyr::filter(scale==3) %>%
  dplyr::select(common_name, season, ldi)
sp_summer = sp.names[sp.names$season=="breeding",]
ldi_summer = ldi_table[ldi_table$season %in% c("breeding", "resident"),] %>%
  dplyr::select(common_name, ldi)
sp_summer = left_join(sp_summer, ldi_summer, by="common_name")
sp_winter = sp.names[sp.names$season=="nonbreeding",]
ldi_winter = ldi_table[ldi_table$season %in% c("nonbreeding", "resident"),] %>%
  dplyr::select(common_name, ldi)
sp_winter = left_join(sp_winter, ldi_winter, by="common_name")
sp.names = rbind(sp_summer,sp_winter)

# Phylo prep
# get data and bird phylo
data(birds)
birds = as(birds, 'phylo')
# check which are missing
unique=uniques=unique(as.character(sp.names$scientific_name))
sp.names$phylo_names = sp.names$scientific_name
sp.names$phylo_names[which(!(sp.names$phylo_names %in% birds$tip.label))]

{# resolve taxonomy, drop all species other than ours
  sp.names[sp.names$phylo_names=="Spinus_tristis","phylo_names"] = "Carduelis_tristis"
  sp.names[sp.names$phylo_names=="Poecile_atricapillus","phylo_names"] = "Parus_atricapillus"
  sp.names[sp.names$phylo_names=="Setophaga_virens","phylo_names"] = "Dendroica_virens"
  sp.names[sp.names$phylo_names=="Poecile_carolinensis","phylo_names"] = "Parus_carolinensis"
  sp.names[sp.names$phylo_names=="Setophaga_pensylvanica","phylo_names"] = "Dendroica_pensylvanica"
  sp.names[sp.names$phylo_names=="Dryobates_pubescens","phylo_names"] = "Picoides_pubescens"
  sp.names[sp.names$phylo_names=="Ardea_alba","phylo_names"] = "Casmerodius_albus"
  sp.names[sp.names$phylo_names=="Dryobates_villosus","phylo_names"] = "Picoides_villosus"
  sp.names[sp.names$phylo_names=="Setophaga_citrina","phylo_names"] = "Wilsonia_citrina"
  sp.names[sp.names$phylo_names=="Haemorhous_mexicanus","phylo_names"] = "Carpodacus_mexicanus"
  sp.names[sp.names$phylo_names=="Leucophaeus_atricilla","phylo_names"] = "Larus_atricilla"
  sp.names[sp.names$phylo_names=="Setophaga_magnolia","phylo_names"] = "Dendroica_magnolia"
  sp.names[sp.names$phylo_names=="Oreothlypis_ruficapilla","phylo_names"] = "Vermivora_ruficapilla"
  sp.names[sp.names$phylo_names=="Setophaga_americana","phylo_names"] = "Parula_americana"
  sp.names[sp.names$phylo_names=="Setophaga_pinus","phylo_names"] = "Dendroica_pinus"
  sp.names[sp.names$phylo_names=="Setophaga_discolor","phylo_names"] = "Dendroica_discolor"
  sp.names[sp.names$phylo_names=="Haemorhous_purpureus","phylo_names"] = "Carpodacus_purpureus"
  sp.names[sp.names$phylo_names=="Antigone_canadensis","phylo_names"] = "Grus_canadensis"
  sp.names[sp.names$phylo_names=="Tringa_semipalmata","phylo_names"] = "Catoptrophorus_semipalmatus"
  sp.names[sp.names$phylo_names=="Setophaga_coronata","phylo_names"] = "Dendroica_coronata"
  sp.names[sp.names$phylo_names=="Setophaga_petechia","phylo_names"] = "Dendroica_aestiva"
  sp.names[sp.names$phylo_names=="Acanthis_flammea","phylo_names"] = "Carduelis_flammea"
  sp.names[sp.names$phylo_names=="Ammospiza_caudacuta","phylo_names"] = "Ammodramus_caudacutus"
  sp.names[sp.names$phylo_names=="Ammospiza_leconteii","phylo_names"] = "Ammodramus_leconteii"
  sp.names[sp.names$phylo_names=="Ammospiza_maritima","phylo_names"] = "Ammodramus_maritimus"
  sp.names[sp.names$phylo_names=="Ammospiza_nelsoni","phylo_names"] = "Ammodramus_nelsoni"
  sp.names[sp.names$phylo_names=="Anas_diazi","phylo_names"] = "Anas_platyrhynchos"
  sp.names[sp.names$phylo_names=="Anser_caerulescens","phylo_names"] = "Chen_caerulescens"
  sp.names[sp.names$phylo_names=="Anser_rossii","phylo_names"] = "Chen_rossii"
  sp.names[sp.names$phylo_names=="Antrostomus_arizonae","phylo_names"] = "Caprimulgus_carolinensis"
  sp.names[sp.names$phylo_names=="Antrostomus_carolinensis","phylo_names"] = "Caprimulgus_carolinensis"
  sp.names[sp.names$phylo_names=="Antrostomus_vociferus","phylo_names"] = "Caprimulgus_vociferus"
  sp.names[sp.names$phylo_names=="Aphelocoma_wollweberi","phylo_names"] = "Caprimulgus_carolinensis"
  sp.names[sp.names$phylo_names=="Aphelocoma_woodhouseii","phylo_names"] = "Caprimulgus_carolinensis"
  sp.names[sp.names$phylo_names=="Ardenna_creatopus","phylo_names"] = "Puffinus_creatopus"
  sp.names[sp.names$phylo_names=="Ardenna_gravis","phylo_names"] = "Puffinus_gravis"
  sp.names[sp.names$phylo_names=="Ardenna_grisea","phylo_names"] = "Puffinus_gravis"
  sp.names[sp.names$phylo_names=="Ardenna_pacifica","phylo_names"] = "Puffinus_gravis"
  sp.names[sp.names$phylo_names=="Ardenna_tenuirostris","phylo_names"] = "Puffinus_tenuirostris"
  sp.names[sp.names$phylo_names=="Artemisiospiza_belli","phylo_names"] = "Amphispiza_belli"
  sp.names[sp.names$phylo_names=="Artemisiospiza_nevadensis","phylo_names"] = "Amphispiza_belli"
  sp.names[sp.names$phylo_names=="Bubo_scandiacus","phylo_names"] = "Bubo_scandiaca"
  sp.names[sp.names$phylo_names=="Buteo_plagiatus","phylo_names"] = "Buteo_nitidus"
  sp.names[sp.names$phylo_names=="Calidris_subruficollis","phylo_names"] = "Tryngites_subruficollis"
  sp.names[sp.names$phylo_names=="Calidris_virgata","phylo_names"] = "Aphriza_virgata"
  sp.names[sp.names$phylo_names=="Cardellina_canadensis","phylo_names"] = "Wilsonia_canadensis"
  sp.names[sp.names$phylo_names=="Cardellina_pusilla","phylo_names"] = "Wilsonia_pusilla"
  sp.names[sp.names$phylo_names=="Centronyx_bairdii","phylo_names"] = "Ammodramus_bairdii"
  sp.names[sp.names$phylo_names=="Centronyx_henslowii","phylo_names"] = "Ammodramus_henslowii"
  sp.names[sp.names$phylo_names=="Charadrius_nivosus","phylo_names"] = "Charadrius_alexandrinus"
  sp.names[sp.names$phylo_names=="Chasiempis_sclateri","phylo_names"] = "Chasiempis_sandwichensis"
  sp.names[sp.names$phylo_names=="Chlorodrepanis_flava","phylo_names"] = "Hemignathus_flavus"
  sp.names[sp.names$phylo_names=="Chlorodrepanis_stejnegeri","phylo_names"] = "Hemignathus_kauaiensis"
  sp.names[sp.names$phylo_names=="Chlorodrepanis_virens","phylo_names"] = "Hemignathus_virens"
  sp.names[sp.names$phylo_names=="Chroicocephalus_philadelphia","phylo_names"] = "Larus_philadelphia"
  sp.names[sp.names$phylo_names=="Chroicocephalus_ridibundus","phylo_names"] = "Larus_ridibundus"
  sp.names[sp.names$phylo_names=="Circus_hudsonius","phylo_names"] = "Circus_cyaneus"
  sp.names[sp.names$phylo_names=="Drepanis_coccinea","phylo_names"] = "Vestiaria_coccinea"
  sp.names[sp.names$phylo_names=="Dryobates_albolarvatus","phylo_names"] = "Picoides_albolarvatus"
  sp.names[sp.names$phylo_names=="Dryobates_arizonae","phylo_names"] = "Picoides_arizonae"
  sp.names[sp.names$phylo_names=="Dryobates_borealis","phylo_names"] = "Picoides_borealis"
  sp.names[sp.names$phylo_names=="Dryobates_nuttallii","phylo_names"] = "Picoides_nuttallii"
  sp.names[sp.names$phylo_names=="Dryobates_scalaris","phylo_names"] = "Picoides_scalaris"
  sp.names[sp.names$phylo_names=="Falcipennis_canadensis","phylo_names"] = "Dendragapus_canadensis"
  sp.names[sp.names$phylo_names=="Gallinago_delicata","phylo_names"] = "Gallinago_gallinago"
  sp.names[sp.names$phylo_names=="Gallinula_galeata","phylo_names"] = "Gallinula_chloropus"
  sp.names[sp.names$phylo_names=="Gelochelidon_nilotica","phylo_names"] = "Sterna_nilotica"
  sp.names[sp.names$phylo_names=="Geothlypis_formosa","phylo_names"] = "Oporornis_formosus"
  sp.names[sp.names$phylo_names=="Geothlypis_philadelphia","phylo_names"] = "Oporornis_philadelphia"
  sp.names[sp.names$phylo_names=="Geothlypis_tolmiei","phylo_names"] = "Oporornis_tolmiei"
  sp.names[sp.names$phylo_names=="Geranoaetus_albicaudatus","phylo_names"] = "Buteo_albicaudatus"
  sp.names[sp.names$phylo_names=="Haemorhous_cassinii","phylo_names"] = "Carpodacus_cassinii"
  sp.names[sp.names$phylo_names=="Hydrocoloeus_minutus","phylo_names"] = "Larus_minutus"
  sp.names[sp.names$phylo_names=="Hydroprogne_caspia","phylo_names"] = "Sterna_caspia"
  sp.names[sp.names$phylo_names=="Ixoreus_naevius","phylo_names"] = "Zoothera_naevia"
  sp.names[sp.names$phylo_names=="Lanius_borealis","phylo_names"] = "Lanius_excubitor"
  sp.names[sp.names$phylo_names=="Leiothlypis_crissalis","phylo_names"] = "Vermivora_crissalis"
  sp.names[sp.names$phylo_names=="Leiothlypis_luciae","phylo_names"] = "Vermivora_luciae"
  sp.names[sp.names$phylo_names=="Leucophaeus_pipixcan","phylo_names"] = "Larus_pipixcan"
  sp.names[sp.names$phylo_names=="Magumma_parva","phylo_names"] = "Hemignathus_parvus"
  sp.names[sp.names$phylo_names=="Mareca_americana","phylo_names"] = "Anas_americana"
  sp.names[sp.names$phylo_names=="Mareca_penelope","phylo_names"] = "Anas_penelope"
  sp.names[sp.names$phylo_names=="Mareca_strepera","phylo_names"] = "Anas_strepera"
  sp.names[sp.names$phylo_names=="Melanitta_americana","phylo_names"] = "Melanitta_nigra"
  sp.names[sp.names$phylo_names=="Melanitta_deglandi","phylo_names"] = "Melanitta_fusca"
  sp.names[sp.names$phylo_names=="Melozone_aberti","phylo_names"] = "Pipilo_aberti"
  sp.names[sp.names$phylo_names=="Melozone_crissalis","phylo_names"] = "Pipilo_crissalis"
  sp.names[sp.names$phylo_names=="Melozone_fusca","phylo_names"] = "Pipilo_fuscus"
  sp.names[sp.names$phylo_names=="Onychoprion_anaethetus","phylo_names"] = "Sterna_anaethetus"
  sp.names[sp.names$phylo_names=="Onychoprion_fuscatus","phylo_names"] = "Sterna_fuscata"
  sp.names[sp.names$phylo_names=="Parkesia_motacilla","phylo_names"] = "Seiurus_motacilla"
  sp.names[sp.names$phylo_names=="Parkesia_noveboracensis","phylo_names"] = "Seiurus_noveboracensis"
  sp.names[sp.names$phylo_names=="Peucaea_aestivalis","phylo_names"] = "Aimophila_aestivalis"
  sp.names[sp.names$phylo_names=="Peucaea_botterii","phylo_names"] = "Aimophila_botterii"
  sp.names[sp.names$phylo_names=="Peucaea_carpalis","phylo_names"] = "Aimophila_carpalis"
  sp.names[sp.names$phylo_names=="Peucaea_cassinii","phylo_names"] = "Aimophila_cassinii"
  sp.names[sp.names$phylo_names=="Phalaropus_tricolor","phylo_names"] = "Steganopus_tricolor"
  sp.names[sp.names$phylo_names=="Poecile_gambeli","phylo_names"] = "Parus_gambeli"
  sp.names[sp.names$phylo_names=="Poecile_hudsonicus","phylo_names"] = "Parus_hudsonicus"
  sp.names[sp.names$phylo_names=="Poecile_rufescens","phylo_names"] = "Parus_rufescens"
  sp.names[sp.names$phylo_names=="Poecile_sclateri","phylo_names"] = "Parus_sclateri"
  sp.names[sp.names$phylo_names=="Psiloscops_flammeolus","phylo_names"] = "Otus_flammeolus"
  sp.names[sp.names$phylo_names=="Poecile_sclateri","phylo_names"] = "Parus_sclateri"
  sp.names[sp.names$phylo_names=="Rallus_crepitans","phylo_names"] = "Rallus_longirostris"
  sp.names[sp.names$phylo_names=="Rallus_obsoletus","phylo_names"] = "Rallus_longirostris"
  sp.names[sp.names$phylo_names=="Rhynchophanes_mccownii","phylo_names"] = "Calcarius_mccownii"
  sp.names[sp.names$phylo_names=="Selasphorus_calliope","phylo_names"] = "Stellula_calliope"
  sp.names[sp.names$phylo_names=="Setophaga_caerulescens","phylo_names"] = "Dendroica_caerulescens"
  sp.names[sp.names$phylo_names=="Setophaga_castanea","phylo_names"] = "Dendroica_castanea"
  sp.names[sp.names$phylo_names=="Setophaga_cerulea","phylo_names"] = "Dendroica_cerulea"
  sp.names[sp.names$phylo_names=="Setophaga_chrysoparia","phylo_names"] = "Dendroica_chrysoparia"
  sp.names[sp.names$phylo_names=="Setophaga_dominica","phylo_names"] = "Dendroica_dominica"
  sp.names[sp.names$phylo_names=="Setophaga_fusca","phylo_names"] = "Dendroica_fusca"
  sp.names[sp.names$phylo_names=="Setophaga_graciae","phylo_names"] = "Dendroica_graciae"
  sp.names[sp.names$phylo_names=="Setophaga_kirtlandii","phylo_names"] = "Dendroica_kirtlandii"
  sp.names[sp.names$phylo_names=="Setophaga_nigrescens","phylo_names"] = "Dendroica_nigrescens"
  sp.names[sp.names$phylo_names=="Setophaga_occidentalis","phylo_names"] = "Dendroica_occidentalis"
  sp.names[sp.names$phylo_names=="Setophaga_palmarum","phylo_names"] = "Dendroica_palmarum"
  sp.names[sp.names$phylo_names=="Setophaga_striata","phylo_names"] = "Dendroica_striata"
  sp.names[sp.names$phylo_names=="Setophaga_tigrina","phylo_names"] = "Dendroica_tigrina"
  sp.names[sp.names$phylo_names=="Setophaga_townsendi","phylo_names"] = "Dendroica_townsendi"
  sp.names[sp.names$phylo_names=="Spatula_clypeata","phylo_names"] = "Anas_clypeata"
  sp.names[sp.names$phylo_names=="Spatula_cyanoptera","phylo_names"] = "Anas_cyanoptera"
  sp.names[sp.names$phylo_names=="Spatula_discors","phylo_names"] = "Anas_discors"
  sp.names[sp.names$phylo_names=="Spinus_lawrencei","phylo_names"] = "Carduelis_lawrencei"
  sp.names[sp.names$phylo_names=="Spinus_pinus","phylo_names"] = "Carduelis_pinus"
  sp.names[sp.names$phylo_names=="Spinus_psaltria","phylo_names"] = "Carduelis_psaltria"
  sp.names[sp.names$phylo_names=="Sternula_antillarum","phylo_names"] = "Sterna_antillarum"
  sp.names[sp.names$phylo_names=="Synthliboramphus_scrippsi","phylo_names"] = "Synthliboramphus_hypoleucus"
  sp.names[sp.names$phylo_names=="Thalasseus_elegans","phylo_names"] = "Sterna_elegans"
  sp.names[sp.names$phylo_names=="Thalasseus_maximus","phylo_names"] = "Sterna_maxima"
  sp.names[sp.names$phylo_names=="Thalasseus_sandvicensis","phylo_names"] = "Sterna_sandvicensis"
  sp.names[sp.names$phylo_names=="Tringa_incana","phylo_names"] = "Heteroscelus_incanus"
  sp.names[sp.names$phylo_names=="Troglodytes_hiemalis","phylo_names"] = "Troglodytes_troglodytes"
  sp.names[sp.names$phylo_names=="Troglodytes_pacificus","phylo_names"] = "Troglodytes_troglodytes"
  sp.names[sp.names$phylo_names=="Vermivora_cyanoptera","phylo_names"] = "Vermivora_pinus"
  sp.names[sp.names$phylo_names=="Leiothlypis_celata","phylo_names"] = "Vermivora_celata"
  sp.names[sp.names$phylo_names=="Leiothlypis_peregrina","phylo_names"] = "Vermivora_peregrina"
  sp.names[sp.names$phylo_names=="Leiothlypis_ruficapilla","phylo_names"] = "Vermivora_ruficapilla"
  sp.names[sp.names$phylo_names=="Leiothlypis_virginiae","phylo_names"] = "Vermivora_virginiae"
  sp.names[sp.names$phylo_names=="Spizelloides_arborea","phylo_names"] = "Spizella_arborea"
  
  birds.sub=drop.tip(birds, setdiff(birds$tip.label, sp.names$phylo_names))
  save(birds.sub, file="outputs/birds_sub.RData")
}

# Export
write_csv(sp.names, "outputs/mod_slopes_traits.csv")



# IF ALREADY EXTRACTED- BRING OVER TRAITS FOR NEW VERSION
traits = read_csv("species_traits.csv")
spdata = read_csv("outputs/mod_slopes.csv") %>%
  left_join(traits, by=c("common_name", "season")) %>%
  write_csv("outputs/mod_slopes_traits.csv")


##### Cross-species
setwd("~/stoat/bg_points/")

### Functional trait data

# Species names, niche estimation, traits
spdata = read_csv("outputs/mod_slopes_traits.csv") %>%
  # cut duplicates
  distinct()
# phylo object
load("outputs/birds_sub.RData")


### Traits/phylo associations with range, niche shifts 

{# preliminary
  labels = c("Realized Δ\nMax. temperature (°C)",
             "Mitigated Δ\nMax. temperature (°C)",
             "Realized Δ EVI", "Mitigated Δ EVI",
             "Realized Δ Daily precip.", "Mitigated Δ Daily precip.",
             "Δ Latitude", "Δ Elevation (m)", "Δ Distance (km)")
  responses = c("tmax_sp_corr","tmax_sp_mit",
                "modisevi_sp_corr","modisevi_sp_mit",
                "prec_sp_corr","prec_sp_mit",
                "lat_coef","elev_coef","dist_coef")
  
  # set up for models
  predictors = c("Migration.Distance","Body.Mass","LDI","Hand.Wing.Index")
  # log transforming
  spdata$Migration.Distance = log(spdata$distance+1)
  spdata$Body.Mass = log(spdata$Mass)
  spdata$Hand.Wing.Index = log(spdata$`Hand-Wing.Index`)
  spdata$LDI = spdata$ldi
  # cut NAs
  spdata = spdata[complete.cases(spdata),]
  # set up output matrices
  output = matrix(nrow=0, ncol=5)
  phylo = matrix(nrow=0, ncol=6)
}

for (season in c("breeding","nonbreeding")){
  # seasonal data
  sdata = spdata[spdata$season==season,]
  for (i in 1:length(responses)){
    #### phylogenetic least squares models
    sdata$resp = sdata[,responses[i]][[1]] *20
    if(i==9){sdata$resp=abs(sdata$resp)}
    formula1 = as.formula(paste0("resp ~",paste(predictors, collapse="+")))
    pgls = gls(formula1, 
               correlation=corBrownian(phy=birds.sub, form=~phylo_names),
               data=sdata)
    table = round(as.data.frame(summary(pgls)$tTable),3)
    #write.csv(table, paste0("outputs/cross_sp/pgls_",responses[i],"_",season,"_rev1b.csv"))
    # with phylolm
    sdata = distinct(sdata, phylo_names, .keep_all=T)
    rownames(sdata) = sdata$phylo_names
    plm = phylolm(formula1, data = sdata, phy = birds.sub, model = "lambda")
    table = round(as.data.frame(summary(plm)$coefficients),3)
    #write.csv(table, paste0("outputs/cross_sp/pgls_",responses[i],"_",season,"_rev1b.csv"))
    
    # glm
    glm = glm(formula1, data=sdata)
    table = round(as.data.frame(summary(glm)$coefficients),3)
    #write.csv(table, paste0("outputs/cross_sp/pgls_",responses[i],"_",season,"_rev1b.csv"))
    
    # partial residual plots
    prednames = c("Migration distance (log(km))","Body mass (log)",
                  "Habitat Diversity Index", "Hand-wing Index (log)")
    test = data.frame(sdata[,predictors])
    for (pred in 1:4){
      #if(i %in% c(6,9,2)){xlab=prednames[pred]}else{xlab=NULL}
      xlab = prednames[pred]
      #if(pred == 1 & season=="breeding"){ylab=labels[i]}else{ylab=""}
      ylab = labels[i]
      
      par = pdp::partial(plm, pred.var=predictors[pred],
                         trim.outliers = F, chull = T, parallel = T,
                         grid.resolution = 30,  paropts = list(.packages = "phylolm"),
                         type="regression", train=sdata)
      names(par)[1] = "var"
      
      # pr <- predict_response(pgls, "x [all]")
      # pr = ggpredict(pgls)
      # test$resp = as.vector(LMMstar::residuals(
      #   pgls, type="partial", variable=predictors[pred])) + mean(test[,predictors[pred]])
      # test$resp = as.vector(residuals(pgls)) + mean(test[,predictors[pred]])
      
      test$resp = as.vector(predict(plm, sdata, se.fit=T))
      test$var = test[,predictors[pred]]
      
      
      ggpd = ggplot(data=par, aes(x=var, y=yhat)) +
        #geom_smooth(method="lm", level=.99) +
        geom_point(data=test, aes(x=var, y=resp)) + # ADD points manually
        geom_smooth(data=test, aes(x=var, y=resp), method="lm", level=.99) + 
        ylab(ylab) +
        xlab(xlab) +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = 15),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.ticks = element_line(color="gray70"),
              axis.line = element_line(color = "black"),
              axis.text.x = element_text(color = "black", size=15),
              axis.text.y = element_text(color = "black", size=15))
      assign(paste0("gg.pr.",i,".",pred), ggpd)
    }
    
    # univariate models
    for (pred in 1:4){
      formula = as.formula(paste0("resp ~",predictors[pred]))
      # pgls = gls(formula, 
      #            correlation=corBrownian(phy=birds.sub, form=~phylo_names),
      #             data=sdata)
      # table = round(as.data.frame(summary(pgls)$tTable),3)
      plm = phylolm(formula, data = sdata, phy = birds.sub, model = "lambda")
      table = round(as.data.frame(summary(plm)$coefficients),3)
      v = as.numeric(table[2,c(1,4)])
      #write.csv(table, paste0("outputs/cross_sp/uni_pgls_",responses[i],"_",season,"_rev1b.csv"))
      # stats
      name = paste0(responses[i],"_",predictors[pred])
      lm = lm(formula, data=sdata)
      summ = c(season, name, v, summary(lm)$r.squared)
      output = rbind(output, summ)
      
      # plots   
      rownames(sdata) = sdata$phylo_names
      xlab = prednames[pred]
      ylab = labels[i]
      plotdata = sdata
      plotdata$var = sdata[,predictors[pred]][[1]]
      ggpd = ggplot(data=plotdata, aes(x=var, y=resp)) +
        #geom_smooth(method="lm", level=.95) +
        geom_point() + 
        geom_smooth(method="lm", level=.95) + 
        ylab(ylab) +
        xlab(xlab) +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = 15),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.ticks = element_line(color="gray70"),
              axis.line = element_line(color = "black"),
              axis.text.x = element_text(color = "black", size=15),
              axis.text.y = element_text(color = "black", size=15))
      assign(paste0("gg.u.",i,".",pred), ggpd)
      
      
      
    }
    
    #### phylogenetic signal
    response = sdata[,responses[i]][[1]]
    names(response) = sdata$phylo_names
    # lambda
    l = phylosig(birds.sub, response, method="lambda", test=T,
                 nsim=1000, se=NULL, start=NULL, control=list())
    # blomberg's k
    k = phylosig(birds.sub, response, method="K", test=T,
                 nsim=1000, se=NULL, start=NULL, control=list())
    # outputs
    summ = c(season, responses[i], l$lambda, l$P, k$K, k$P)
    phylo = rbind(phylo, summ)
    
  }
  
  # migration, HWI, body mass
  gga = ggarrange(gg.pr.7.1, gg.pr.7.2, gg.pr.7.3, gg.pr.7.4, 
                  gg.pr.1.1,  gg.pr.1.2,  gg.pr.1.3,  gg.pr.1.4,
                  gg.pr.2.1, gg.pr.2.2,  gg.pr.2.3, gg.pr.2.4, 
                  nrow=3, ncol=4
                  #,labels=c("a)","b)","c)","d)","e)","f))"
  )
  # ggsave(paste0("outputs/cross_sp/functional_figure_rev1b_",season,"2.jpeg"), 
  #        gga, "jpeg", height=9, width=14, units="in")
  # # Migration only figure
  # if(season=="breeding"){labs = c("a)","c)","e)")}else{labs = c("b)","d)","f)")}
  # gga = ggarrange(gg.pr.7.1,gg.pr.1.1, gg.pr.2.1,
  #                 nrow=3, ncol=1, labels=labs)
  # ggsave(paste0("outputs/cross_sp/functional_figure_mig_",season,"2.jpeg"), gga, "jpeg",
  #        height=9, width=4, units="in")
  
  # Migration only figure - just realized niche shifts
  ggsave(paste0("outputs/cross_sp/functional_figure_mig_",season,"_realized.jpeg"), 
         gg.pr.1.1, "jpeg", height=4, width=4, units="in")
}

{# univariate models table
colnames(output) = c("Season","Model","Coefficient","p-value","R2-value")
write_csv(data.frame(output), "outputs/cross_sp/univariate_rev1b.csv")

# phylogenetic correlation table
phylo = data.frame(phylo)
phylo[,2] = c("Max. Temperature","Mitigated Max. Temperature",
              "EVI","Mitigated EVI","Precipitation","Mitigated Precipitation",
              "Latitude","Elevation","Distance")
phylo = cbind(phylo[,1:2],round(as.numeric(phylo[,3]),3),
              round(as.numeric(phylo[,4]),3),
              round(as.numeric(phylo[,5]),4), phylo[,6])
colnames(phylo) = c("Season","Metric","Lambda","p-value","K", "P-value")
write_csv(data.frame(phylo), "outputs/cross_sp/phylo_rev1b.csv")
}




