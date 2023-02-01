################################################################################
# TITLE: AIS x Seabirds Analysis
#
# PURPOSE: This script uses the NPPSD v3.0 and satellite AIS vessel tracking data 
# to calculate the overlap between vessel activity and seabird distributions. 
#
# AUTHORS: Ben Sullender & Kelly Kapsar, with guidance from Kathy Kuletz
#
# CREATED: 13 July 2022
# LAST UPDATED ON: 1 December 2022
################################################################################

########################
#### Load Libraries #### 
########################

library(tidyverse)
library(sf)
library(rnaturalearth)

################################
#### Read in Raw Data Files #### 
################################

# Save folder
# Specify file path for folder where results will be saved
savefolder <- "../Data_Processed/"
figfolder <- "../Figures/"

# Read in vessel traffic data 
hexdir <- "D:/AIS_V2_DayNight_60km6hrgap/Hex/"
hexList <-list.files(hexdir, pattern=".shp")

# Read in blank hexagon template 
# hex <- st_read("../Data_Raw/BlankHexes.shp")

# Read in a custom-created hexagon mask, which masks out all parts of hexagons that are land and calculates
# total marine area. This will be useful in calculating survey effort.
hexMask <- st_read("../Data_Raw/hex_x_ocean/hex_x_ocean.shp") %>%
  select(hexID) %>%
  mutate(AreaKM = c(st_area(.)/1000000)) 

# Load in seabird data from NPPSD v3.0
loc <- read.csv("../Data_Raw/NPPSD_v3.0/tbl_LOCATION.csv") 
datobs <- read.csv("../Data_Raw/NPPSD_v3.0/tbl_DATA_OBS.csv")

# Load in basemap for plots
studyarea <- st_read("../Data_Raw/AK_CAN_RUS/AK_CAN_RUS.shp") 

###################################
#### Specify Input Information #### 
###################################

# Specify metric used to calculate vessel activity 
metric <- "OperatingDays" # MUST be "OperatingDays" or "Ships"
## OperatingDays = number of ship days per month (i.e., the same ship in the same hex each day for a month = 30)
## Ships = number of unique ships per month (i.e., the same ship int h same hex each day for a month = 1)

nightonly <- TRUE #If true, will only calculate nighttime vessel traffic (ignoring daytime)
# If false, will calculate all vessel traffic, including both day and night

months <- c(9:11) # Numeric value(s) of months to be included in the analysis
# months <- c(6:8)

monthsname <- "Fall" # Text label describing the numeric months in the analysis (e.g., "Summer", "Annual")
# monthsname <- "Summer" 

startyear <- 2006 # Earliest year for which bird observations will be included in the analysis 
# NOTE: start year is for seabird observations only. Vessel traffic includes all data from 2015-2020

effortThreshold <- 0.01 # Percentage area of each hex that has to be observed in order to include the hex in the analysis


#### Seabird Species of Interest #### 
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

taxaList <- list(totalBirds, seaducks, shear, 
                 murre, auklet, phal, eider,
                 murrelet, kitti, gull, puffin, 
                 guill, scoter, corm, loon, 
                 alba, duckswangoose, stormpet, nofu, 
                 kiei, larid, alcid)
names(taxaList) <- c("AllBirds","Seaducks", "Shearwaters", 
                     "Murre", "Auklet", "Phalarope", "Eider",
                     "Murrelet", "Kittiwake", "Gull", "Puffin",
                     "Guillemot", "Scoter","Cormorant","Loon", 
                     "Albatross", "DuckSwanGoose","StormPetrel", "NorthernFulmar",
                     "KingEider", "Larids", "Alcids")

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
#                  night=nightonly)

# For a list of species 
lapply(1:length(taxaList), function(x){birdHexesByEffort(taxaNames= taxaList[[x]], 
                              taxaLabel= names(taxaList[x]),
                              hexMask=hexMask, 
                              effortThreshold=effortThreshold, 
                              loc=loc,
                              mnths=months,
                              mnthsnam=monthsname,
                              startyr=startyear,
                              savefolder=savefolder,
                              figfolder=figfolder,
                              filedir=hexList, 
                              studyarea=studyarea,
                              metric=metric, 
                              night=nightonly)})

# night2 <- st_read("../Data_Processed/FinalDF_AllBirds_Fall_NightOnlyTRUE.shp")
# allday2 <- st_read("../Data_Processed/FinalDF_AllBirds_Fall_NightOnlyFALSE.shp")
# all.equal(night2, allday2)
# all.equal(night, allday)

browseURL("https://www.youtube.com/watch?v=K1b8AhIsSYQ")

# bird <- st_read("../Data_Processed/ObsInHexes_Fall.shp")
# ship <- read.csv("../Data_Processed/TraffInHexes_Fall_NightOnlyTRUE.csv")
# finaldf <- st_read("../Data_Processed/FinalDF_Albatross_Fall_NightOnlyTRUE.shp")

