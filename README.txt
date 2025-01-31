## DATA FILES

### nonbreeding_processed_20yr.csv 
### breeding_processed_20yr.csv 
Processed, filtered GBIF data from 2000-2020 for December-February or June-August (respectively, in order of files listed) with spatially linked environmental covariates. Data is structured with rows for each observation. Cells containing NA have no associated environmental information available at a given location.

#### Column IDs
sci_name: scientific name
eventdate: date of observation
latitude: latitude of observation coordinate
longitude: longitude of observation coordinate
tmax: maximum daily temperature over the 30 days prior to the observation
modisevi: Enhanced Vegetation Index over the 30 days prior to the observation
prec: Total precipitation over the 30 days prior to the observation
season: breeding (summer) or nonbreeding (winter)
year: year of observation
lon: longitude rounded to two decimals (for spatial filtering)
lat: latitude rounded to two decimals (for spatial filtering)
week: week of year of observation (1-52; for temporal filtering)
elev: elevation (mm)


### baseline_points_nonbreeding_bl.csv
### baseline_points_breeding_bl.csv
"Baseline" observations, which are pseudo-observations of species based on GBIF data from 2000-2004. Pseudo-observations are created for all years 2000-2020 at the locations of all observations 2000-2004 to estimate climate change exposure in the absence of movement. Baseline observations are for December-February or June-August (respectively, in order of files listed) with spatially linked environmental covariates. Data is structured with rows for each observation. Cells containing NA have no associated environmental information available at a given location.

#### Column IDs
species: scientific name
season: breeding (summer) or nonbreeding (winter)
latitude: latitude of observation coordinate
longitude: longitude of observation coordinate
eventdate: date of observation
ID: ID number
elev: elevation (mm)
tmax: maximum daily temperature over the 30 days prior to the observation
modisevi: Enhanced Vegetation Index over the 30 days prior to the observation
prec: Total precipitation over the 30 days prior to the observation


### birds_vs1_sp_level_corrected_with_ldi.csv 
Species-level estimation of landcover diversity index (a measure of habitat specialization) estimated in a different study (Cohen and Jetz 2025 Global Ecology and Biogeography).

#### Column IDs
common_name: common name of species
sciname: scientific name of species
season: breeding (summer) or nonbreeding (winter)
scale: grain size at which model was fit (not relevant to current study)
ldi: landscape diversity index


### species_list.csv 
Index of all candidate species. Cells containing NA have no information about extinction for a given species.

#### Column IDs
species name: species name
common name: common name
Abbr: four letter species abbreviation
Code: species rarity code- 1 (common) or 2 (uncommon)
ebird_code: species six letter code
TAXON_ORDER_CODE: order code
category: taxonomic category
range: description of range
order: taxonomic order
family: taxonomic family
eBird species group: species group description
waterbird: (not used in manuscript)
extinction: (not used in manuscript)
extinct year: (not used in manuscript)
IW: whether the species is an interior western species (for use in elevation models)

### migration.txt
Results of models fit to estimate trends in species-level workflow, for use in modeling and plotting cross-species trends.

#### Column IDs
scientific: species name
LON.b: median breeding longitude
LAT.b: median breeding latitude
LON.nb" median nonbreeding longitude
LAT.nb: median nonbreeding latitude
MIN.LAT.b: minimum breeding latitude
MAX.LAT.b: maximum breeding latitude
MIN.LAT.nb: minimum nonbreeding latitude
MAX.LAT.nb: maximum nonbreeding latitude
distance: distance between breeding and nonbreeding median latitude
intersect: 1 if breeding and nonbreeding ranges intersect, 0 if not


### mod_slopes.csv 
Results of models fit to estimate trends in species-level workflow, for use in modeling and plotting cross-species trends.

#### Column IDs
common_name: species common name
scientific_name: species name
season: breeding (summer) or nonbreeding (winter)
mahal_slope_niche: Slope of Mahalanobis distances between 3-dimensional realized and baseline niches over time
mahal_se_niche: SE of above
mahal_p_niche: p value of above
dr_slope_niche: Slope of Determinant ratios (measure of niche breadth difference) between 3-dimensional realized and baseline niches over time
dr_se_niche: SE of above
dr_p_niche: p value of above
tmax_coef: slope of max. temperature shift over time across focal species, baseline, and observer datasets
tmax_se: SE of above
tmax_p: p value of above
tmax_coef_bg: slope of max. temperature shift over time for observer dataset (from emmeans)
tmax_se_bg: SE of above
tmax_lcl_bg: lower confidence limit of above at 95%
tmax_ucl_bg: upper confidence limit of above at 95%
tmax_coef_bl: slope of max. temperature shift over time for baseline dataset (from emmeans)
tmax_se_bl: SE of above
tmax_lcl_bl: lower confidence limit of above at 95%
tmax_ucl_bl: upper confidence limit of above at 95%
tmax_coef_sp: slope of max. temperature shift over time for focal species dataset (from emmeans)
tmax_se_sp: SE of above
tmax_lcl_sp: lower confidence limit of above at 95%
tmax_ucl_sp: upper confidence limit of above at 95%
modisevi_coef: slope of Enhanced Vegetation Index over time across focal species, baseline, and observer datasets
modisevi_se: SE of above
modisevi_p: p value of above
modisevi_coef_bg: slope of Enhanced Vegetation Index over time for observer dataset (from emmeans)
modisevi_se_bg: SE of above
modisevi_lcl_bg: lower confidence limit of above at 95%
modisevi_ucl_bg: upper confidence limit of above at 95%
modisevi_coef_bl: slope of Enhanced Vegetation Index over time for baseline dataset (from emmeans)
modisevi_se_bl: SE of above
modisevi_lcl_bl: lower confidence limit of above at 95%
modisevi_ucl_bl: upper confidence limit of above at 95%
modisevi_coef_sp: slope of Enhanced Vegetation Index over time for focal species dataset (from emmeans)
modisevi_se_sp: SE of above
modisevi_lcl_sp: lower confidence limit of above at 95%
modisevi_ucl_sp: upper confidence limit of above at 95%
prec_coef: slope of Total Precipitation over time across focal species, baseline, and observer datasets
prec_se: SE of above
prec_p: p value of above
prec_coef_bg: slope of Total Precipitation over time for observer datasets (from emmeans)
prec_se_bg: SE of above
prec_lcl_bg: lower confidence limit of above at 95%
prec_ucl_bg: upper confidence limit of above at 95%
prec_coef_bl: slope of Total Precipitation over time for baseline datasets (from emmeans)
prec_se_bl: SE of above
prec_lcl_bl: lower confidence limit of above at 95%
prec_ucl_bl: upper confidence limit of above at 95%
prec_coef_sp: slope of Total Precipitation over time for focal species datasets (from emmeans)
prec_se_sp: SE of above
prec_lcl_sp: lower confidence limit of above at 95%
prec_ucl_sp: upper confidence limit of above at 95%
lat_coef: slope of latitude over time across focal species and observer datasets
lat_se: SE of above
lat_p: p value of above
lat_coef_bg: slope of latitude over time for observer dataset (from emmeans)
lat_se_bg: SE of above
lat_lcl_bg: lower confidence limit of above at 95%
lat_ucl_bg: upper confidence limit of above at 95%
lat_coef_sp: slope of latitude over time for focal species dataset (from emmeans)
lat_se_sp: SE of above
lat_lcl_sp: lower confidence limit of above at 95%
lat_ucl_sp: upper confidence limit of above at 95%
elev_coef: slope of elevation over time across focal species and observer datasets
elev_se: SE of above
elev_p: p value of above
elev_coef_bg: slope of elevation over time for observer dataset (from emmeans)
elev_se_bg: SE of above
elev_lcl_bg: lower confidence limit of above at 95%
elev_ucl_bg: upper confidence limit of above at 95%
elev_coef_sp: slope of elevation over time for focal species dataset (from emmeans)
elev_se_sp: SE of above
elev_lcl_sp: lower confidence limit of above at 95%
elev_ucl_sp: upper confidence limit of above at 95%
dist_coef: slope of total distance over time across focal species and observer datasets
dist_se: SE of above
dist_p: p value of above
dist_coef_bg: slope of total distance over time for observer dataset (from emmeans)
dist_se_bg: SE of above
dist_lcl_bg: lower confidence limit of above at 95%
dist_ucl_bg: upper confidence limit of above at 95%
dist_coef_sp: slope of total distance over time for focal species dataset (from emmeans)
dist_se_sp: SE of above
dist_lcl_sp: lower confidence limit of above at 95%
dist_ucl_sp
lon_coef: slope of longitude over time across focal species and observer datasets
lon_se: SE of above
lon_p: p value of above
lon_coef_bg: slope of longitude over time for observer dataset (from emmeans)
lon_se_bg: SE of above
lon_lcl_bg: lower confidence limit of above at 95%
lon_ucl_bg: upper confidence limit of above at 95%
lon_coef_sp: slope of longitude over time for focal species dataset (from emmeans)
lon_se_sp: SE of above
lon_lcl_sp: lower confidence limit of above at 95%
lon_ucl_sp: upper confidence limit of above at 95%
IW: 1 if the species is interior western, 0 if not (for inclusion in elevation comparisons)
tmax_obs_corr: corrective element of max. temperature slope over time (background slope - baseline slope)
tmax_sp_corr: realized max. temperature slope over time, adjusted for observer (focal species slope - corrective element)
tmax_sp_mit: mitigated max. temperature slope over time, adjusted for observer (realized slope - baseline slope)
modisevi_obs_corr: corrective element of enhanced vegetation index over time (background slope - baseline slope)
modisevi_sp_corr: realized enhanced vegetation index slope over time, adjusted for observer (focal species slope - corrective element)
modisevi_sp_mit: mitigated enhanced vegetation index slope over time, adjusted for observer (realized slope - baseline slope)
prec_obs_corr: corrective element of total precipitation slope over time (background slope - baseline slope)
prec_sp_corr: realized total precipitation slope over time, adjusted for observer (focal species slope - corrective element)
prec_sp_mit: mitigated total precipitation slope over time, adjusted for observer (realized slope - baseline slope)





## Sharing/access info
GBIF data can be downloaded via www.gbif.org
Ebird data can be downloaded via www.ebird.org
Breeding bird survey data can be downloaded via https://www.usgs.gov/data/2023-release-north-american-breeding-bird-survey-dataset-1966-2022
Environmental data can be downloaded from CHELSA (https://chelsa-climate.org/) or MODIS (https://lpdaac.usgs.gov/products/mod11a1v006/) 
AVONET trait data is available at https://onlinelibrary.wiley.com/doi/full/10.1111/ele.13898


## CODE FILES

1. baseline_points_github.R - generation of baseline pseudo-observations for GBIF data
2. data_build_ebird_github.R - eBird data compilation and baseline generation
3. data_build_bbs_github.R - Breeding bird survey data compilation and baseline generation
4. bg_points_v12a_github_rev1.R - workflow to estimate trends in niches and distributions from observations
5. trends_github_rev2.R - estimation of cross-species patterns in niche/distribution trends, unrelated to functional traits or phylogeny
6. traits_phylo_github_rev1c.R - estimation of cross-species patterns in niche/distribution trends, related to functional traits or phylogeny
7. fix_stat_ellipse_github.R - code to patch a glitch in ellipse formation for ggplot niche plots.