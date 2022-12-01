# AIS x Seabirds Pipeline
# Based on NPPSD v3.0
#
# Ben Sullender, with guidance from Kathy Kuletz
# July 13, 2022

library(tidyverse)
library(sf)

# Save folder
# Specify file path for folder where results will be saved
savefolder <- "../Data_Processed/"


# read in just hexagons, drop all attributes except hex ID
# setwd("/Users/bensullender/Documents/KickstepApproaches_Projects/FWS_AIS/Seabird_Analysis")

hexdir <- "E:/AIS_V2_DayNight_60km6hrgap/Hex/"
hexList <-list.files(hexdir, pattern=".shp")

hex <- st_read("../Data_Raw/BlankHexes.shp")
metric <- "OperatingDays" # MUST be "OperatingDays" or "Ships"
## OperatingDays = number of ship days per month (i.e., the same ship in the same hex each day for a month = 30)
## Ships = number of unique ships per month (i.e., the same ship int h same hex each day for a month = 1)

nightonly <- TRUE #If true, will only calculate nighttime vessel traffic (ignoring daytime)
# If false, will calculate all vessel traffic, including both day and night

# From downloaded NPPSD:
loc <- read.csv("../Data_Raw/NPPSD_v3.0/tbl_LOCATION.csv")
datobs <- read.csv("../Data_Raw/NPPSD_v3.0/tbl_DATA_OBS.csv")

# loc has 460,285 rows; 169,283 are >= 2007; 169,273 are >= 2007 and on-transect
# datobs has 823,326 rows; 305,556 are >= 2007; 305,478 are >= 2007 and on-transect

# And we'll also read in a custom-created hexagon mask, which masks out all parts of hexagons that are land and calculates
#         total marine area. This will be useful in calculating survey effort.
hexMask <- st_read("../Data_Raw/hex_x_ocean/hex_x_ocean.shp") %>%
  select(hexID) %>%
  mutate(AreaKM = c(st_area(.)/1000000)) 

months <- c(9:11) # Numeric value(s) of months to be included in the analysis 
monthsname <- "Fall" # Text label describing the numeric months in the analysis (e.g., "Summer", "Annual")
startyear <- 2006 # Earliest year for which bird observations will be included in the analysis 
# NOTE: start year is for seabird observations only. Vessel traffic includes all data from 2015-2020
effortThreshold <- 0.01 # Percentage area of each hex that has to be observed in order to include the hex in the analysis

##############################################################################################################
# Function to spatially intersect at sea bird observations with vessel traffic hex grid
# Also calculates survey effort within each hex during the time period 
# startyr = oldest year of data collection to be included in the analysis
# mnths = months of observations to be included in the anlaysis
# mnthsname = name used to save file and identify month grouping (e.g., "Summer")
surveyEffort <- function(loc, datobs, hex, startyr, mnths, survfilename){
  # Drop ones that are too old - based on Kathy Kuletz's feedback, prior to 2007 is too old.
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
  
  birdNames <- unique(datobs$NPPSD.4.Letter.Code)
  
  
  # convert to SF and transform to Alaska Albers (epsg 3338), otherwise hexagons over -180:180 get warped
  datSF <- st_as_sf(datobs,coords = c("Lon","Lat"),crs=4326) %>%
    st_transform(crs=3338)
  locSF <- st_as_sf(datobs,coords = c("Lon","Lat"),crs=4326) %>%
    st_transform(crs=3338)
  
  # Set up one df for bird obs and another df for survey effort
  locIn <- locSF %>%
    filter(Month %in% mnths) %>%
    st_join(hex) %>%
    st_drop_geometry()
  
  obsIn <- datSF %>%
    filter(Month %in% mnths)  %>%
    st_join(hex) %>%
    st_drop_geometry() %>%
    filter(Number!="")
  
  res <- hex
  sppNout <- as.list(NA)
  
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
  
  surveyEff <- locIn %>%
    group_by(hexID) %>%
    summarize(survEff = sum(as.numeric(Sample.Area.y))) %>%
    as.data.frame()
  
  res <- left_join(res,surveyEff,by="hexID")
  
  # switch NAs to 0s
  res[is.na(res)] <- 0
  
  # SET UP AUTO NAME GENERATION
  st_write(res,survfilename)
  
}

##############################################################################################################
#
# Stack Vessel Traffic
#
# hexFall <- hexList %>%
#   filter(c(substr(name,15,16) %in% c("09","10","11")))

# Metric Options
## OperatingDays = number of ship days per month (i.e., the same ship in the same hex each day for a month = 30)
## Ships = number of unique ships per month (i.e., the same ship int h same hex each day for a month = 1)

hexStack <- function(filedir, mnths, metric, night=TRUE, trafffilename){
  
  hexes <- filedir[as.numeric(substr(filedir,15,16)) %in% mnths]
  
  # h = hex; S = summer; 1 = first in list
  temp <- lapply(hexes, function(x){st_read(paste0(hexdir,x)) %>% st_drop_geometry()})
  hexAll <- do.call(rbind, temp)
  
  lab <- ifelse(metric=="OperatingDays", "OpD", 
                ifelse(metric=="Ships","Shp", NA))
  lab <- ifelse(night==TRUE, paste0("N_",lab),lab)
  
  test <- hexAll[,grep(lab, colnames(hexAll))]
  test$hexID <- hexAll$hexID
  
  test[is.na(test)] <- 0
  
  hexRes <- test %>%
    group_by(hexID) %>%
    summarise(across(everything(), sum))
  
      
  allcol <- colnames(hexRes[grepl("_A",colnames(hexRes))])
  
  calcdis <- hexRes %>% dplyr::select(all_of(allcol)) %>% pull()
  
  hexRes$AllShip <- calcdis
  hexRes$QuantShip <- ecdf(calcdis)(calcdis)
  
  hexRes$ClassShip <- ifelse(hexRes$AllShip<mean(hexRes$AllShip),1,
         ifelse(hexRes$AllShip>=c(mean(hexRes$AllShip)+sd(hexRes$AllShip)),3,2))
  
  hexFinal <- hexRes %>% dplyr::select(hexID, AllShip, QuantShip, ClassShip)
  
  write.csv(hexFinal,trafffilename)
}

##############################################################################################################
##
##          Guilding bird spp
##
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





##
##      Sample runthrough
##

# Inputs
# effortThreshold = ???? 
# mnths = months of observations to be included in the anlaysis
# taxa = list of 4 digit codes for the groups of seabirds to be included in each analysis 

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
  
  finaldfname <- paste0(savefolder,"FinalDF_",taxaLabel,"_",monthsname,"_night",night[1],".shp")
  
  if(!file.exists(finaldfname)){
    # First, let's select only hexes with sufficient survey effort
    resGuild <- res %>%
      # and filter by 1% of hex area must be surveyed
      filter(as.numeric(survEff)>c(effortThreshold*as.numeric(AreaKM)))
    
    # now for each taxonomic group, create one column with sum of all obs
    resGuildAll <- resGuild[,c(colnames(resGuild) %in% taxaNames)] %>%
      st_drop_geometry()
    resGuild$AllBird <- rowSums(resGuildAll)
    
    # Stack just the relevant guilds, remove extraneous columns
    resStacked <- resGuild %>%
      select(hexID,AllBird,survEff,AreaKM)
    
    resFinal <- resStacked %>%
      mutate(QuantBird = ecdf(AllBird/survEff)(AllBird/survEff))
    
    resFinal$ClassBird <- ifelse(c(resFinal$AllBird/resFinal$survEff)<mean(c(resFinal$AllBird/resFinal$survEff)),1,
                             ifelse(c(resFinal$AllBird/resFinal$survEff)>=c(mean(c(resFinal$AllBird/resFinal$survEff))+sd(c(resFinal$AllBird/resFinal$survEff))),3,2))
    
    resFinal$taxa <- taxaLabel
    
    finalCombined <- resFinal %>%
      # join in the actual hex marine area sf
      left_join(y=hexFinal,by="hexID")
    
    st_write(finalCombined, paste0(savefolder,"FinalDF_",taxaLabel,"_",monthsname,"_night",night[1],".shp"))
    }
}

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

temp <- list.files("../Data_Processed/", pattern="FinalDF",full.names=TRUE)
temp2 <- temp[grep(".shp", temp)]
temp3 <- lapply(temp2, st_read)

lapply(temp3,dim)
