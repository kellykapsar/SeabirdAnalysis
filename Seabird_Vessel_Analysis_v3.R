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

################################
#### Read in Raw Data Files #### 
################################

# Save folder
# Specify file path for folder where results will be saved
savefolder <- "../Data_Processed/"

# Read in vessel traffic data 
hexdir <- "E:/AIS_V2_DayNight_60km6hrgap/Hex/"
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

monthsname <- "Fall" # Text label describing the numeric months in the analysis (e.g., "Summer", "Annual")

startyear <- 2006 # Earliest year for which bird observations will be included in the analysis 
# NOTE: start year is for seabird observations only. Vessel traffic includes all data from 2015-2020

effortThreshold <- 0.01 # Percentage area of each hex that has to be observed in order to include the hex in the analysis


#### Seabird Species of Interest #### 
allSpp <- read.csv("../Data_Raw/NPPSD_Bird_Codes_Only_Revised.csv")
totalBirds <- allSpp$Code
seaducks <- c("WWSC","SUSC","UNSC","BLSC","KIEI","STEI","SPEI","UNEI","LTDU")
aethia <- c("CRAU","LEAU","PAAU","UNAU","USDA")
shear <- c("STSH","SOSH","UNSH","UNDS")

taxaList <- list(totalBirds, seaducks, aethia, shear)
names(taxaList) <- c("AllBirds","Seaducks","Aethia","Shearwaters")


## Lower priority bird taxa groups

#murre <- c("COMU","TBMU","UNMU")
#auklet <- c("CRAU","LEAU","PAAU","UNAU","USDA")
#phal <- c("REPH","RNPH","UNPH")
#shear <- c("STSH","SOSH","UNSH")
#eider <- c("KIEI","STEI","SPEI","UNEI")
#murrelet <- c("KIMU","MAMU","BRMU","ANMU","UNML")
#kitti <- c("BLKI","RLKI","UNKI")
#nofu <- c("NOFU")
#gull <- c("GWGU","MEGU","GLGU","UNGU","ROGU","SAGU","IVGU","HEGU",
#          "SBIG","BOGU","HERG","THGU","ICGU","SBGU","SBAG")
#puffin <- c("TUPU","HOPU","UNPU")
#guill <- c("PIGU","BLGU","UNGI")
#scoter <- c("WWSC","SUSC","UNSC","BLSC")
#eidscot <- c(eider,scoter)
#corm <- c("PECO","DCCO","RFCO","UNCO")
#loon <- c("YBLO","PALO","ARLO","RTLO","COLO","UNLO")
#alba <- c("UALB","BFAL","LAAL","STAL")
#duckswangoose <- c("LTDU","HADU","UNGO","BAGO","COGO","UNDU","CAGO","GWFG","SNGO",
#            "ROGO","BRAN","BLBR","CANG","UNGO","TRUS","TUNS")
#stormpet <- c("FTSP","LESP","UNSP")

##########################
#### Define Functions #### 
##########################

#### Observations and Survey Effort #### 

# Function to spatially intersect at sea bird observations with vessel traffic hex grid
# Also calculates survey effort within each hex during the time period 

surveyEffort <- function(loc, datobs, hex, startyr, mnths, survfilename){
  # Drop observations that are too old or off-transect
  # Based on Kathy Kuletz's feedback, prior to 2007 is too old.
  loc <- loc %>%
    filter(Year > startyr) %>%
    # Removing off-transect from location data
    filter(Modified.Survey.Type != "Off Transect Observation") %>%
    droplevels()
  
  datobs <- datobs %>%
    # Removing off-transect from observations
    filter(OT.OBS != "Off") %>%
    left_join(x=datobs,y=loc,by="Master.Key") %>%
    filter(Year > startyr) %>%
    droplevels()
  
  # Identify unique 4-digit codes for each bird spp 
  birdNames <- unique(datobs$NPPSD.4.Letter.Code)
  
  # convert to SF and transform to Alaska Albers (epsg 3338), otherwise hexagons over -180:180 get warped
  datSF <- st_as_sf(datobs,coords = c("Lon","Lat"),crs=4326) %>%
    st_transform(crs=3338)
  locSF <- st_as_sf(datobs,coords = c("Lon","Lat"),crs=4326) %>%
    st_transform(crs=3338)
  
  # Set up one df for bird obs and another df for survey effort
  locIn <- locSF %>%
    # Remove observations from months outside of window of interest
    filter(Month %in% mnths) %>%
    # Spatially join with hexagon data 
    st_join(hex) %>%
    st_drop_geometry()
  
  obsIn <- datSF %>%
    filter(Month %in% mnths)  %>%
    st_join(hex) %>%
    st_drop_geometry() %>%
    filter(Number!="")
  
  res <- hex
  sppNout <- as.list(NA)
  
  # Iterate through each spp and calculate number of observations per hex
  for (i in 1:length(birdNames)){
    obsIn2 <- obsIn[obsIn$NPPSD.4.Letter.Code==birdNames[i],]
    sppNout[[i]] <- obsIn2 %>%
      group_by(hexID) %>%
      summarize(birdN = sum(as.numeric(Number))) %>%
      as.data.frame()
    names(sppNout[[i]]) <- c("hexID",birdNames[i])
  }
  
  for (j in 1:length(birdNames)){
    res <- left_join(res,sppNout[[j]],by="hexID")
  }
  
  # Calculate survey effort within each hex 
  surveyEff <- locIn %>%
    group_by(hexID) %>%
    summarize(survEff = sum(as.numeric(Sample.Area.y))) %>%
    as.data.frame()
  
  res <- left_join(res,surveyEff,by="hexID")
  
  # switch NAs to 0s
  res[is.na(res)] <- 0
  
  # Save output 
  st_write(res,survfilename)
}

#### Vessel Activity #### 

hexStack <- function(filedir, mnths, metric, night=TRUE, trafffilename){
  
  # Isolate month values of interest from file names 
  hexes <- filedir[as.numeric(substr(filedir,15,16)) %in% mnths]
  
  # Read in data 
  temp <- lapply(hexes, function(x){st_read(paste0(hexdir,x)) %>% st_drop_geometry()})
  hexAll <- do.call(rbind, temp)
  
  # Select metrics based on inputs 
  lab <- ifelse(metric=="OperatingDays", "OpD", 
                ifelse(metric=="Ships","Shp", NA))
  lab <- ifelse(night==TRUE, paste0("N_",lab),lab)
  
  test <- hexAll[,grep(lab, colnames(hexAll))]
  test$hexID <- hexAll$hexID
  
  # Replace NA values with zero 
  test[is.na(test)] <- 0
  
  # Calculate total vessel traffic across study period for months of interest 
  hexRes <- test %>%
    group_by(hexID) %>%
    summarise(across(everything(), sum))
  
  # Only use data for "all" vessels instead of subsetting by vessel type 
  allcol <- colnames(hexRes[grepl("_A",colnames(hexRes))])
  
  calcdis <- hexRes %>% dplyr::select(all_of(allcol)) %>% pull()
  
  hexRes$AllShip <- calcdis
  
  # Calculate SD categories for each hex's vessel activity  
  hexRes$QuantShip <- ecdf(calcdis)(calcdis)
  
  hexRes$ClassShip <- ifelse(hexRes$AllShip<mean(hexRes$AllShip),1,
         ifelse(hexRes$AllShip>=c(mean(hexRes$AllShip)+sd(hexRes$AllShip)),3,2))
  
  # Isolate only columns of interest
  hexFinal <- hexRes %>% dplyr::select(hexID, AllShip, QuantShip, ClassShip)
  
  # Save results 
  write.csv(hexFinal,trafffilename)
}

#### Seabird Metrics #### 
# NOTE: This function is dependent on the surveyEffort and hexStack functions 

birdHexesByEffort <- function(dataobs,
                              loc,
                              taxaNames, 
                              taxaLabel,
                              hexMask, 
                              effortThreshold, 
                              mnths,
                              mnthsnam,
                              startyr,
                              savefolder,
                              filedir, 
                              metric, 
                              night=TRUE){
  
  # Make sure appropriate vessel data exist and if not, generate it
  trafffilename <- paste0(savefolder,"TraffInHexes_",mnthsnam,".csv")
  
  if(!file.exists(trafffilename)){
    hexStack(filedir = hexList, mnths = mnths, metric = metric, night = night, trafffilename)
  }
  
  hexFinal <- read.csv(trafffilename)
  
  # Make sure survey effort has been calculated and saved for appropriate months and if not, generate it
  survfilename <- paste0(savefolder,"ObsInHexes_",monthsname,".shp")
  
  if(!file.exists(survfilename)){
    surveyEffort(loc=loc, datobs=datobs, hex=hexMask, startyr=startyr, mnths=mnths, survfilename=survfilename)
  }
  res <- st_read(survfilename)
  
  # Make sure this function has not already been run with these exact parameter specifications 
  finaldfname <- paste0(savefolder,"FinalDF_",taxaLabel,"_",monthsname,"_night",night[1],".shp")
  
  if(!file.exists(finaldfname)){
    
    # Select only hexes with sufficient survey effort
    resGuild <- res %>%
      filter(as.numeric(survEff)>c(effortThreshold*as.numeric(AreaKM)))
    
    # Isolate taxa of interest and calculate total number of observations 
    resGuildAll <- resGuild[,c(colnames(resGuild) %in% taxaNames)] %>%
      st_drop_geometry()
    resGuild$AllBird <- rowSums(resGuildAll)
    
    # Stack just the relevant guilds, remove extraneous columns
    resStacked <- resGuild %>%
      select(hexID,AllBird,survEff,AreaKM)
    
    # Calculate SD categories for each hex's seabird observations 
    resFinal <- resStacked %>%
      mutate(QuantBird = ecdf(AllBird/survEff)(AllBird/survEff))
    
    # Calculate effort-weighted classes for the number of observations 
    resFinal$ClassBird <- ifelse(c(resFinal$AllBird/resFinal$survEff)<mean(c(resFinal$AllBird/resFinal$survEff)),1,
                             ifelse(c(resFinal$AllBird/resFinal$survEff)>=c(mean(c(resFinal$AllBird/resFinal$survEff))+sd(c(resFinal$AllBird/resFinal$survEff))),3,2))
    
    # Create column to specify taxa of interest
    resFinal$taxa <- taxaLabel
    
    # Join to vessel traffic data 
    finalCombined <- resFinal %>%
      left_join(y=hexFinal,by="hexID")
    
    # Save results 
    st_write(finalCombined, paste0(savefolder,"FinalDF_",taxaLabel,"_",monthsname,"_night",night[1],".shp"))
    }
}

##########################
#### Run Analysis #### 
##########################

# For one species 
birdHexesByEffort(taxaNames= totalBirds, 
                 taxaLabel= "AllBirds",
                 hexMask=hexMask, 
                 effortThreshold=effortThreshold, 
                 loc=loc,
                 mnths=months,
                 mnthsnam=monthsname,
                 startyr=startyear,
                 savefolder=savefolder,
                 filedir=hexList, 
                 metric=metric, 
                 night=nightonly)

# For a list of species 
lapply(1:length(taxaList), function(x){birdHexesByEffort(taxaNames= taxaList[x], 
                              taxaLabel= names(taxaList[x]),
                              hexMask=hexMask, 
                              effortThreshold=effortThreshold, 
                              loc=loc,
                              mnths=months,
                              mnthsnam=monthsname,
                              startyr=startyear,
                              savefolder=savefolder,
                              filedir=hexList, 
                              metric=metric, 
                              night=nightonly)})