################################################################################
# TITLE: AIS x Seabirds Analysis
#
# PURPOSE: This script uses the NPPSD v3.0 and satellite AIS vessel tracking data 
# to calculate the overlap between vessel activity and seabird distributions. 
#
# AUTHORS: Ben Sullender & Kelly Kapsar, with guidance from Kathy Kuletz
#
# CREATED: 13 July 2022
# LAST UPDATED ON: 23 February 2023
################################################################################

########################
#### Load Libraries #### 
########################

library(tidyverse)
library(sf)
library(rnaturalearth)
library(biscale)
library(cowplot)
library(ggpubr)
library(grid)
library(gridExtra)

################################
#### Read in Raw Data Files #### 
################################

# Read in vessel traffic data 
hexdir <- "D:/AIS_V2_DayNight_60km6hrgap/Hex_DayNight_Hours/"
hexList <-list.files(hexdir, pattern=".shp")

# Read in blank hexagon template 
# hex <- st_read("../Data_Raw/BlankHexes.shp")

# Read in a custom-created hexagon mask, which masks out all parts of hexagons that are land and calculates
# total marine area. This will be useful in calculating survey effort.
hexMask <- st_read("../Data_Raw/hex_x_ocean/hex_x_ocean.shp") %>%
  select(hexID) %>%
  mutate(AreaKM = c(st_area(.)/1000000)) 

# Load in seabird data from NPPSD v3.0
loc <- read.csv("../Data_Raw/NPPSD_v4.0/Locations_NPPSDv4.csv") 
datobs <- read.csv("../Data_Raw/NPPSD_v4.0/Observations_NPPSDv4.csv")

# Load in basemap for plots
basemap <- st_read("../Data_Raw/AK_CAN_RUS/AK_CAN_RUS.shp") 

# study area
goa <- st_read("../Data_Raw/KodiakPWS.shp")
uni <- st_read("../Data_Raw/Unimak.shp")
berstr <- st_read("../Data_Raw/BeringStrait.shp") %>% filter(Shape_Leng > 3000000) 
berchuk <- st_read("../Data_Raw/BerChukBeau.shp")
npac <- st_as_sf(data.frame(name="All Alaska", geometry= st_as_sfc(st_bbox(hexMask))))

###################################
#### Specify Input Information #### 
###################################

# Specify metric used to calculate vessel activity 
metric <- "Hours" # MUST be "OperatingDays", "Ships", or "Hours"
metricName <-"Hours"
## OperatingDays = number of ship days per month (i.e., the same ship in the same hex each day for a month = 30)
## Ships = number of unique ships per month (i.e., the same ship int h same hex each day for a month = 1)

timeofday <- "Night" #If true, will only calculate nighttime vessel traffic (ignoring daytime)
# If false, will calculate all vessel traffic, including both day and night

# months <- c(9:11) # Numeric value(s) of months to be included in the analysis
months <- c(6:8)

# monthsname <- "Fall" # Text label describing the numeric months in the analysis (e.g., "Summer", "Annual")
monthsname <- "Summer"

startyear <- 2006 # Earliest year for which bird observations will be included in the analysis 
# NOTE: start year is for seabird observations only. Vessel traffic includes all data from 2015-2020

effortThreshold <- 0.01 # Percentage area of each hex that has to be observed in order to include the hex in the analysis

# Save folder
# Specify file path for folder where results will be saved
savefolder <- "../Data_Processed/"

#####################################
#### Seabird Species of Interest #### 
#####################################

allSpp <- read.csv("../Data_Raw/NPPSD_Bird_Codes_Only_Revised.csv")
totalBirds <- allSpp$Code
seaducks <- c("WWSC","SUSC","UNSC","BLSC","KIEI","STEI","SPEI","UNEI","LTDU")
shear <- c("STSH","SOSH","UNSH","UNDS")

## Lower priority bird taxa groups

murre <- c("COMU","TBMU","UNMU")
auklet <- c("CRAU","LEAU","PAAU","UNAU","USDA")
phal <- c("REPH","RNPH","UNPH")
eider <- c("KIEI","STEI","SPEI","UNEI")
murrelet <- c("KIMU","MAMU","BRMU","ANMU","UNML")
kitti <- c("BLKI","RLKI","UNKI")
gull <- c("GWGU","MEGU","GLGU","UNGU","ROGU","SAGU","IVGU","HEGU",
         "SBIG","BOGU","HERG","THGU","ICGU","SBGU","SBAG")
puffin <- c("TUPU","HOPU","UNPU")
guill <- c("PIGU","BLGU","UNGI")
scoter <- c("WWSC","SUSC","UNSC","BLSC")
corm <- c("PECO","DCCO","RFCO","UNCO")
loon <- c("YBLO","PALO","ARLO","RTLO","COLO","UNLO")
alba <- c("UALB","BFAL","LAAL","STAL")
duckswangoose <- c("LTDU","HADU","UNGO","BAGO","COGO","UNDU","CAGO","GWFG","SNGO",
"ROGO","BRAN","BLBR","CANG","UNGO","TRUS","TUNS")
stormpet <- c("FTSP","LESP","UNSP")
nofu <- c("NOFU")
kiei <- c("KIEI")

larid <- c(kitti, gull)
alcid <- c(murre, auklet, guill, puffin)


AllAKList <- list(totalBirds, seaducks, shear,
                  murre, auklet, phal, eider,
                  murrelet, kitti, gull, puffin,
                  guill, scoter, corm, loon,
                  alba, duckswangoose, stormpet, nofu,
                  kiei, larid, alcid)


names(AllAKList) <- c("Seabirds","Seaducks", "Shearwaters",
                      "Murres", "Auklets", "Phalaropes", "Eiders",
                      "Murrelets", "Kittiwakes", "Gulls", "Puffins",
                      "Guillemots", "Scoters","Cormorants","Loons",
                      "Albatross", "DuckSwanGoose","Storm Petrels", "Northern Fulmars",
                      "King Eiders", "Larids", "Alcids")

# AlaskaList <- list(totalBirds)
# names(AlaskaList <- "Seabirds")

AleutList <- list(alba)
names(AleutList) <- c("Albatross")

GOAList <- list(shear, stormpet)
names(GOAList) <- c("Shearwaters", "Storm Petrels")

BerChukList <- list(seaducks, shear, auklet)
names(BerChukList) <- c("Seaducks", "Shearwaters", "Auklets")

###########################
#### Load in Functions #### 
###########################

source("./Seabird_Vessel_Analysis_Functions.R")


##########################
#### Run Analysis #### 
##########################

# For one species 
# birdHexesByEffort(taxaNames= aethia, 
#                  taxaLabel= "Aethia",
#                  hexMask=hexMask, 
#                  effortThreshold=effortThreshold, 
#                  loc=loc,
#                  mnths=months,
#                  mnthsnam=monthsname,
#                  startyr=startyear,
#                  savefolder=savefolder,
#                  filedir=hexList, 
#                  metric=metric, 
#                  timeofday=timeofday)

# For a list of species
## ALEUTIANS
lapply(1:length(AllAKList), function(x){birdHexesByEffort(taxaNames= AllAKList[[x]],
                                                         taxaLabel= names(AllAKList[x]),
                                                         hexMask=hexMask,
                                                         effortThreshold=effortThreshold,
                                                         loc=loc,
                                                         mnths=months,
                                                         mnthsnam=monthsname,
                                                         studyarea=uni, 
                                                         studyareaname = "Eastern Aleutians",
                                                         figfolder="../Figures/Aleutians/",
                                                         startyr=startyear,
                                                         savefolder=savefolder,
                                                         filedir=hexList,
                                                         metric=metric,
                                                         timeofday=timeofday)})
## GULF OF ALASKA
lapply(1:length(AllAKList), function(x){birdHexesByEffort(taxaNames= AllAKList[[x]],
                                                         taxaLabel= names(AllAKList[x]),
                                                         hexMask=hexMask,
                                                         effortThreshold=effortThreshold,
                                                         loc=loc,
                                                         mnths=months,
                                                         mnthsnam=monthsname,
                                                         studyarea=goa, 
                                                         studyareaname = "Gulf of Alaska",
                                                         figfolder="../Figures/GOA/",
                                                         startyr=startyear,
                                                         savefolder=savefolder,
                                                         filedir=hexList,
                                                         metric=metric,
                                                         timeofday=timeofday)})
# NORTHERN BERING & CHUKCHI
lapply(1:length(AllAKList), function(x){birdHexesByEffort(taxaNames= AllAKList[[x]],
                                                         taxaLabel= names(AllAKList[x]),
                                                         hexMask=hexMask,
                                                         effortThreshold=effortThreshold,
                                                         loc=loc,
                                                         mnths=months,
                                                         mnthsnam=monthsname,
                                                         studyarea=berchuk, 
                                                         studyareaname = "Northern Bering & Chukchi Seas",
                                                         figfolder="../Figures/NorthBeringAndChukchi/",
                                                         startyr=startyear,
                                                         savefolder=savefolder,
                                                         filedir=hexList,
                                                         metric=metric,
                                                         timeofday=timeofday)})
## ALL ALASKA 
lapply(1:length(AllAKList), function(x){birdHexesByEffort(taxaNames= AllAKList[[x]],
                                                         taxaLabel= names(AllAKList[x]),
                                                         hexMask=hexMask,
                                                         effortThreshold=effortThreshold,
                                                         loc=loc,
                                                         mnths=months,
                                                         mnthsnam=monthsname,
                                                         startyr=startyear,
                                                         studyarea=npac, 
                                                         studyareaname = "All Alaska",
                                                         figfolder="../Figures/AllAlaska/",
                                                         savefolder=savefolder,
                                                         filedir=hexList,
                                                         metric=metric,
                                                         timeofday=timeofday)})

################################# RESULTS PLOTS #################################
# Plots for All Alaska  
lapply(1:length(AleutList), function(x){tryCatch(plotResults(studyarea=uni,
                                                            studyareaname = "Eastern Aleutians",
                                                            basemap=basemap,
                                                            figfolder="../Figures/Aleutians/",
                                                            savefolder=savefolder,
                                                            monthsname=monthsname,
                                                            taxaLabel=names(AleutList[x]),
                                                            timeofday=timeofday,
                                                            metricName=metricName),error=function(e) NULL)})


lapply(1:length(GOAList), function(x){tryCatch(plotResults(studyarea=goa, 
                                                             studyareaname = "Gulf of Alaska",
                                                             basemap=basemap,
                                                             figfolder="../Figures/GOA/",
                                                             savefolder=savefolder,
                                                             monthsname=monthsname,
                                                             taxaLabel=names(GOAList[x]),
                                                             timeofday=timeofday,
                                                             metricName=metricName),error=function(e) NULL)})

lapply(1:length(BerChukList), function(x){tryCatch(plotResults(studyarea=berchuk, 
                                                           studyareaname = "Northern Bering & Chukchi Seas",
                                                           basemap=basemap,
                                                           figfolder="../Figures/NorthBeringAndChukchi/",
                                                           savefolder=savefolder,
                                                           monthsname=monthsname,
                                                           taxaLabel=names(BerChukList[x]),
                                                           timeofday=timeofday,
                                                           metricName=metricName),error=function(e) NULL)})

lapply(1:length(AllAKList), function(x){tryCatch(plotResults(studyarea=npac, 
                                                             studyareaname = "All Alaska",
                                                             basemap=basemap,
                                                             figfolder="../Figures/AllAlaska/",
                                                             savefolder=savefolder,
                                                             monthsname=monthsname,
                                                             taxaLabel=names(AllAKList[x]),
                                                             timeofday=timeofday,
                                                             metricName=metricName),error=function(e) NULL)})

# browseURL("https://www.youtube.com/watch?v=K1b8AhIsSYQ")

# bird <- st_read("../Data_Processed/ObsInHexes_Fall.shp")
# ship <- read.csv("../Data_Processed/TraffInHexes_Fall_NightOnlyTRUE.csv")
# finaldf <- st_read("../Data_Processed/FinalDF_Albatross_Fall_NightOnlyTRUE.shp")

