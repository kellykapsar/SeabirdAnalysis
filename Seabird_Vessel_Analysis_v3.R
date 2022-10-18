# AIS x Seabirds Pipeline
# Based on NPPSD v3.0
#
# Ben Sullender, with guidance from Kathy Kuletz
# July 13, 2022

library(tidyverse)
library(sf)

# read in just hexagons, drop all attributes except hex ID
# setwd("/Users/bensullender/Documents/KickstepApproaches_Projects/FWS_AIS/Seabird_Analysis")

hexdir <- "D:/AlaskaConservation_AIS_20210225/Data_Processed_Hex"
hex <- st_read(paste0(hexdir,"/SpeedHex_2015-01.shp")) %>%
  select(hexID)

# From downloaded NPPSD:
loc <- read.csv("../Data_Raw/NPPSD_v3.0/tbl_LOCATION.csv")
datobs <- read.csv("../Data_Raw/NPPSD_v3.0/tbl_DATA_OBS.csv")

# loc has 460,285 rows; 169,283 are >= 2007; 169,273 are >= 2007 and on-transect
# datobs has 823,326 rows; 305,556 are >= 2007; 305,478 are >= 2007 and on-transect

# And we'll also read in a custom-created hexagon mask, which masks out all parts of hexagons that are land and calculates
#         total marine area. This will be useful in calculating survey effort.
hexMask <- st_read("../Data_Raw/hex_x_ocean/hex_x_ocean.shp") %>%
  select(hexID) %>%
  mutate(AreaKM = c(st_area(.)/1000000)) %>%
  st_drop_geometry()

# Function to spatially intersect at sea bird observations with vessel traffic hex grid
# Also calculates survey effort within each hex during the time period 
# startyr = oldest year of data collection to be included in the analysis
# mnths = months of observations to be included in the anlaysis
# mnthsname = name used to save file and identify month grouping (e.g., "Summer")
surveyeffort <- function(loc, dataobs, hex, startyr = 2006, mnths = c(1:12),mnthsname = "Annual"){
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
  
  # Now split into 2 time periods: summer = June, July, Aug;
  #                                 fall = Sept, Oct, Nov
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
  
  # SET UP AUTO NAME GENERATION
  st_write(res,paste0("../Data_Processed/ObsInHexes_",mnthsname,".shp"))
  
}

surveyeffort(loc, dataobs, hex, startyr=2006, mnths=c(6,7,8), mnthsname = "Summer")
surveyeffort(loc, dataobs, hex, startyr=2006, mnths=c(9,10,11), mnthsname = "Fall")

##############################################################################################################
##
##          Guilding bird spp
##
allSpp <- read.csv("NPPSD_Bird_Codes_Only_Revised.csv")
totalBirds <- allSpp$Code
seaducks <- c("WWSC","SUSC","UNSC","BLSC","KIEI","STEI","SPEI","UNEI","LTDU")
aethia <- c("CRAU","LEAU","PAAU","UNAU","USDA")
shear <- c("STSH","SOSH","UNSH","UNDS")


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

effortThreshold <- c(0.01)

# First, let's select only hexes with sufficient survey effort and turn observation NAs into 0s

resSummGuild <- resSumm %>%
  # join in the actual hex marine area sf
  left_join(y=hexMask,by="hexID") %>%
  # and filter by 1% of hex area must be surveyed
  filter(as.numeric(survEff)>c(effortThreshold*as.numeric(AreaKM)))
# switch NAs to 0s
resSummGuild[is.na(resSummGuild)] <- 0

resFallGuild <- resFall %>%
  # join in the actual hex marine area sf
  left_join(y=hexMask,by="hexID") %>%
  # and filter by 1% of hex area must be surveyed
  filter(as.numeric(survEff)>c(effortThreshold*as.numeric(AreaKM)))
# switch NAs to 0s
resFallGuild[is.na(resFallGuild)] <- 0


# now for each taxonomic group, create one column with sum of all obs
resSummGuildAll <- resSummGuild[,c(colnames(resSummGuild) %in% totalBirds)] %>%
  st_drop_geometry()
resSummGuild$SummAll <- rowSums(resSummGuildAll)

resSummGuildDucks <- resSummGuild[,c(colnames(resSummGuild) %in% seaducks)] %>%
  st_drop_geometry()
resSummGuild$SummSeaducks <- rowSums(resSummGuildDucks)

resSummGuildAethia <- resSummGuild[,c(colnames(resSummGuild) %in% aethia)] %>%
  st_drop_geometry()
resSummGuild$SummAethia <- rowSums(resSummGuildAethia)

resSummGuildShear <- resSummGuild[,c(colnames(resSummGuild) %in% shear)] %>%
  st_drop_geometry()
resSummGuild$SummShear <- rowSums(resSummGuildShear)


# Stack just the relevant guilds, remove extraneous columns
resSummStacked <- resSummGuild %>%
  select(hexID,SummAll,SummSeaducks,SummAethia,SummShear,survEff,AreaKM)



# and the same for fall
resFallGuildAll <- resFallGuild[,c(colnames(resFallGuild) %in% totalBirds)] %>%
  st_drop_geometry()
resFallGuild$FallAll <- rowSums(resFallGuildAll)

resFallGuildDucks <- resFallGuild[,c(colnames(resFallGuild) %in% seaducks)] %>%
  st_drop_geometry()
resFallGuild$FallSeaducks <- rowSums(resFallGuildDucks)

resFallGuildAethia <- resFallGuild[,c(colnames(resFallGuild) %in% aethia)] %>%
  st_drop_geometry()
resFallGuild$FallAethia <- rowSums(resFallGuildAethia)

resFallGuildShear <- resFallGuild[,c(colnames(resFallGuild) %in% shear)] %>%
  st_drop_geometry()
resFallGuild$FallShear <- rowSums(resFallGuildShear)


# Stack just the relevant guilds, remove extraneous columns
resFallStacked <- resFallGuild %>%
  select(hexID,FallAll,FallSeaducks,FallAethia,FallShear,survEff,AreaKM)


# Now switch from absolute counts into quantiles of effort-weighted densities
resFallFinal <- resFallStacked %>%
  mutate(FallAllQ = ecdf(FallAll/survEff)(FallAll/survEff),
         FallSeaducksQ = ecdf(FallSeaducks/survEff)(FallSeaducks/survEff),
         FallAethiaQ = ecdf(FallAethia/survEff)(FallAethia/survEff),
         FallShearQ = ecdf(FallShear/survEff)(FallShear/survEff))

resSummFinal <- resSummStacked %>%
  mutate(SummAllQ = ecdf(SummAll/survEff)(SummAll/survEff),
         SummSeaducksQ = ecdf(SummSeaducks/survEff)(SummSeaducks/survEff),
         SummAethiaQ = ecdf(SummAethia/survEff)(SummAethia/survEff),
         SummShearQ = ecdf(SummShear/survEff)(SummShear/survEff))


#
# Stack Vessel Traffic
#




hexList <- as.data.frame(list.files(hexdir, pattern=".shp"))
names(hexList) <- c("name")

hexSumm <- hexList %>%
  filter(c(substr(name,15,16) %in% c("06","07","08")))

hexFall <- hexList %>%
  filter(c(substr(name,15,16) %in% c("09","10","11")))

# h = hex; S = summer; 1 = first in list
hS1 <- st_read(getwd(),substr(hexSumm[1,1],1,16)) %>%
  st_drop_geometry()
hexAll <- hS1
for (h in 2:nrow(hexSumm)){
  hSn <- st_read(getwd(),substr(hexSumm[h,1],1,16)) %>%
    st_drop_geometry()
  hexAll <- rbind(hexAll,hSn)
}

hexSummRes <- hexAll %>%
  group_by(hexID) %>%
  summarize(DaysAll = sum(nOprD_A,na.rm=T),DaysCargo = sum(nOprD_C,na.rm=T),DaysTanker = sum(nOprD_T,na.rm=T),
            DaysFishing = sum(nOprD_F,na.rm=T),DaysOther = sum(nOprD_O,na.rm=T))

# h = hex; F = fall; 1 = first in list
hF1 <- st_read(getwd(),substr(hexFall[1,1],1,16)) %>%
  st_drop_geometry()
hexAll <- hF1
for (h in 2:nrow(hexFall)){
  hFn <- st_read(getwd(),substr(hexFall[h,1],1,16)) %>%
    st_drop_geometry()
  hexAll <- rbind(hexAll,hFn)
}

hexFallRes <- hexAll %>%
  group_by(hexID) %>%
  summarize(DaysAll = sum(nOprD_A,na.rm=T),DaysCargo = sum(nOprD_C,na.rm=T),DaysTanker = sum(nOprD_T,na.rm=T),
            DaysFishing = sum(nOprD_F,na.rm=T),DaysOther = sum(nOprD_O,na.rm=T))

# Set up quantiles

hexFallFinal <- hexFallRes %>%
  mutate(DaysAllQuant = ecdf(DaysAll)(DaysAll)) %>%
  select(hexID,DaysAll,DaysAllQuant)

hexSummFinal <- hexSummRes %>%
  mutate(DaysAllQuant = ecdf(DaysAll)(DaysAll)) %>%
  select(hexID,DaysAll,DaysAllQuant)




# Use Renner & Kuletz (2015) classification: < avg, > avg, and highest (> (avg + 1 SD))

hexSummFinal$Class <- ifelse(hexSummFinal$DaysAll<mean(hexSummFinal$DaysAll),1,
                             ifelse(hexSummFinal$DaysAll>=c(mean(hexSummFinal$DaysAll)+sd(hexSummFinal$DaysAll)),3,2))

hexFallFinal$Class <- ifelse(hexFallFinal$DaysAll<mean(hexFallFinal$DaysAll),1,
                             ifelse(hexFallFinal$DaysAll>=c(mean(hexFallFinal$DaysAll)+sd(hexFallFinal$DaysAll)),3,2))

hexFallFinal %>% count(Class)






resSummFinal$AllClass <- ifelse(c(resSummFinal$SummAll/resSummFinal$survEff)<mean(c(resSummFinal$SummAll/resSummFinal$survEff)),1,
                                ifelse(c(resSummFinal$SummAll/resSummFinal$survEff)>=c(mean(c(resSummFinal$SummAll/resSummFinal$survEff))+sd(c(resSummFinal$SummAll/resSummFinal$survEff))),3,2))

resSummFinal$DuckClass <- ifelse(c(resSummFinal$SummSeaducks/resSummFinal$survEff)<mean(c(resSummFinal$SummSeaducks/resSummFinal$survEff)),1,
                                 ifelse(c(resSummFinal$SummSeaducks/resSummFinal$survEff)>=c(mean(c(resSummFinal$SummSeaducks/resSummFinal$survEff))+sd(c(resSummFinal$SummSeaducks/resSummFinal$survEff))),3,2))

resSummFinal$AethClass <- ifelse(c(resSummFinal$SummAethia/resSummFinal$survEff)<mean(c(resSummFinal$SummAethia/resSummFinal$survEff)),1,
                                 ifelse(c(resSummFinal$SummAethia/resSummFinal$survEff)>=c(mean(c(resSummFinal$SummAethia/resSummFinal$survEff))+sd(c(resSummFinal$SummAethia/resSummFinal$survEff))),3,2))

resSummFinal$ShearClass <- ifelse(c(resSummFinal$SummShear/resSummFinal$survEff)<mean(c(resSummFinal$SummShear/resSummFinal$survEff)),1,
                                  ifelse(c(resSummFinal$SummShear/resSummFinal$survEff)>=c(mean(c(resSummFinal$SummShear/resSummFinal$survEff))+sd(c(resSummFinal$SummShear/resSummFinal$survEff))),3,2))




resFallFinal$AllClass <- ifelse(c(resFallFinal$FallAll/resFallFinal$survEff)<mean(c(resFallFinal$FallAll/resFallFinal$survEff)),1,
                                ifelse(c(resFallFinal$FallAll/resFallFinal$survEff)>=c(mean(c(resFallFinal$FallAll/resFallFinal$survEff))+sd(c(resFallFinal$FallAll/resFallFinal$survEff))),3,2))

resFallFinal$DuckClass <- ifelse(c(resFallFinal$FallSeaducks/resFallFinal$survEff)<mean(c(resFallFinal$FallSeaducks/resFallFinal$survEff)),1,
                                 ifelse(c(resFallFinal$FallSeaducks/resFallFinal$survEff)>=c(mean(c(resFallFinal$FallSeaducks/resFallFinal$survEff))+sd(c(resFallFinal$FallSeaducks/resFallFinal$survEff))),3,2))

resFallFinal$AethClass <- ifelse(c(resFallFinal$FallAethia/resFallFinal$survEff)<mean(c(resFallFinal$FallAethia/resFallFinal$survEff)),1,
                                 ifelse(c(resFallFinal$FallAethia/resFallFinal$survEff)>=c(mean(c(resFallFinal$FallAethia/resFallFinal$survEff))+sd(c(resFallFinal$FallAethia/resFallFinal$survEff))),3,2))

resFallFinal$ShearClass <- ifelse(c(resFallFinal$FallShear/resFallFinal$survEff)<mean(c(resFallFinal$FallShear/resFallFinal$survEff)),1,
                                  ifelse(c(resFallFinal$FallShear/resFallFinal$survEff)>=c(mean(c(resFallFinal$FallShear/resFallFinal$survEff))+sd(c(resFallFinal$FallShear/resFallFinal$survEff))),3,2))





fallCombined <- resFallFinal %>%
  # join in the actual hex marine area sf
  left_join(y=hexFallFinal,by="hexID")

summCombined <- resSummFinal %>%
  # join in the actual hex marine area sf
  left_join(y=hexSummFinal,by="hexID")


st_write(summCombined,paste0(getwd(),"/summCombined.shp"))
st_write(fallCombined,paste0(getwd(),"/fallCombined.shp"))






#
# For 60 N for greater: manual select in GIS, then export selected
#
hexSummFinal2 <- st_read(paste0(getwd(),"/summCombined_60N.shp"))
hexFallFinal2 <- st_read(paste0(getwd(),"/fallCombined_60N.shp"))
hex2 <- st_read(paste0(getwd(),"/summCombined.shp"))

setwd("/Users/bensullender/Documents/KickstepApproaches_Projects/FWS_AIS/Seabird_Analysis/Output_v2")
library(sf)
# Set this here for 60N vs. all
hexSummFinal <- hexSummFinal2 %>% st_drop_geometry()
hexFallFinal <- hexFallFinal2 %>% st_drop_geometry()
resSummFinal <- hexSummFinal2
resFallFinal <- hexFallFinal2


# Use Renner & Kuletz (2015) classification: < avg, > avg, and highest (> (avg + 1 SD))

hexSummFinal$Class <- ifelse(hexSummFinal$DaysAll<mean(hexSummFinal$DaysAll),1,
                             ifelse(hexSummFinal$DaysAll>=c(mean(hexSummFinal$DaysAll)+sd(hexSummFinal$DaysAll)),3,2))

hexFallFinal$Class <- ifelse(hexFallFinal$DaysAll<mean(hexFallFinal$DaysAll),1,
                             ifelse(hexFallFinal$DaysAll>=c(mean(hexFallFinal$DaysAll)+sd(hexFallFinal$DaysAll)),3,2))

hexFallFinal %>% count(Class)






resSummFinal$AllClass <- ifelse(c(resSummFinal$SummAll/resSummFinal$survEff)<mean(c(resSummFinal$SummAll/resSummFinal$survEff)),1,
                                ifelse(c(resSummFinal$SummAll/resSummFinal$survEff)>=c(mean(c(resSummFinal$SummAll/resSummFinal$survEff))+sd(c(resSummFinal$SummAll/resSummFinal$survEff))),3,2))

resSummFinal$DuckClass <- ifelse(c(resSummFinal$SmmSdck/resSummFinal$survEff)<mean(c(resSummFinal$SmmSdck/resSummFinal$survEff)),1,
                                 ifelse(c(resSummFinal$SmmSdck/resSummFinal$survEff)>=c(mean(c(resSummFinal$SmmSdck/resSummFinal$survEff))+sd(c(resSummFinal$SmmSdck/resSummFinal$survEff))),3,2))

resSummFinal$AethClass <- ifelse(c(resSummFinal$SummAth/resSummFinal$survEff)<mean(c(resSummFinal$SummAth/resSummFinal$survEff)),1,
                                 ifelse(c(resSummFinal$SummAth/resSummFinal$survEff)>=c(mean(c(resSummFinal$SummAth/resSummFinal$survEff))+sd(c(resSummFinal$SummAth/resSummFinal$survEff))),3,2))

resSummFinal$ShearClass <- ifelse(c(resSummFinal$SummShr/resSummFinal$survEff)<mean(c(resSummFinal$SummShr/resSummFinal$survEff)),1,
                                  ifelse(c(resSummFinal$SummShr/resSummFinal$survEff)>=c(mean(c(resSummFinal$SummShr/resSummFinal$survEff))+sd(c(resSummFinal$SummShr/resSummFinal$survEff))),3,2))




resFallFinal$AllClass <- ifelse(c(resFallFinal$FallAll/resFallFinal$survEff)<mean(c(resFallFinal$FallAll/resFallFinal$survEff)),1,
                                ifelse(c(resFallFinal$FallAll/resFallFinal$survEff)>=c(mean(c(resFallFinal$FallAll/resFallFinal$survEff))+sd(c(resFallFinal$FallAll/resFallFinal$survEff))),3,2))

resFallFinal$DuckClass <- ifelse(c(resFallFinal$FllSdck/resFallFinal$survEff)<mean(c(resFallFinal$FllSdck/resFallFinal$survEff)),1,
                                 ifelse(c(resFallFinal$FllSdck/resFallFinal$survEff)>=c(mean(c(resFallFinal$FllSdck/resFallFinal$survEff))+sd(c(resFallFinal$FllSdck/resFallFinal$survEff))),3,2))

resFallFinal$AethClass <- ifelse(c(resFallFinal$FallAth/resFallFinal$survEff)<mean(c(resFallFinal$FallAth/resFallFinal$survEff)),1,
                                 ifelse(c(resFallFinal$FallAth/resFallFinal$survEff)>=c(mean(c(resFallFinal$FallAth/resFallFinal$survEff))+sd(c(resFallFinal$FallAth/resFallFinal$survEff))),3,2))

resFallFinal$ShearClass <- ifelse(c(resFallFinal$FallShr/resFallFinal$survEff)<mean(c(resFallFinal$FallShr/resFallFinal$survEff)),1,
                                  ifelse(c(resFallFinal$FallShr/resFallFinal$survEff)>=c(mean(c(resFallFinal$FallShr/resFallFinal$survEff))+sd(c(resFallFinal$FallShr/resFallFinal$survEff))),3,2))





fallCombined <- resFallFinal %>%
  # join in the actual hex marine area sf
  left_join(y=hexFallFinal,by="hexID")

summCombined <- resSummFinal %>%
  # join in the actual hex marine area sf
  left_join(y=hexSummFinal,by="hexID")


st_write(summCombined,paste0(getwd(),"/summCombined_60N_res.shp"))
st_write(fallCombined,paste0(getwd(),"/fallCombined_60N_res.shp"))





#
#       Bonus code
#



##
## Calculate survey effort
##

# Calculating transect length (km) 

surveySumm <- locSumm %>%
  group_by(hexID) %>%
  # /1000 since width is in m and area is sq km
  summarize(survEff2 = sum(as.numeric(Sample.Area.y)/as.numeric(Transect.Width/1000))) %>%
  as.data.frame()

surveyFall <- locFall %>%
  group_by(hexID) %>%
  # /1000 since width is in m and area is sq km
  summarize(survEff2 = sum(as.numeric(Sample.Area.y)/as.numeric(Transect.Width/1000))) %>%
  as.data.frame()

# only looking at survey length >= 20km from Kuletz et al. 2015
nrow(surveySumm[surveySumm$survEff2>=20,])
nrow(surveyFall[surveyFall$survEff2>=20,])


effortThreshold <- c(0.01)
surveySumm2 <- resSumm %>%
  # join in the actual hex marine area sf
  left_join(y=hexMask,by="hexID") %>%
  # and filter by 1% of hex area must be surveyed
  filter(as.numeric(survEff)>c(effortThreshold*as.numeric(AreaKM)))
nrow(surveySumm2)


surveyFall2 <- resFall %>%
  # join in the actual hex marine area sf
  left_join(y=hexMask,by="hexID") %>%
  # and filter by 1% of hex area must be surveyed
  filter(as.numeric(survEff)>c(effortThreshold*as.numeric(AreaKM)))
nrow(surveyFall2)
# 1712 with >20km; 1687 with 1%
# 1333 with >20km; 1296 with 1%





library(rnaturalearth)
US <- st_as_sf(ne_countries(country = 'united states of america'))
RUS <- st_as_sf(ne_countries(country = 'russia'))
CAN <- st_as_sf(ne_countries(country = 'canada'))


ggplot() +
  geom_sf(data=resMguild,aes(fill = EiderQuant)) +
  geom_sf(data=US,colour='gray') +
  geom_sf(data=RUS,colour='gray') +
  coord_sf(xlim = st_bbox(resMguild)[c(1, 3)], st_bbox(resMguild)[c(2, 4)],
           expand = FALSE)


# select just top 90% quantiles for guild and plot by vessel traffic


resMmurr <- resMguild %>%
  filter(MurreQuant>=0.90)

resMeid <- resMguild %>%
  filter(EiderQuant>=0.95)

ggplot() +
  geom_sf(data=resMmurr,aes(fill = DaysQuant)) +
  geom_sf(data=US,colour='gray') +
  geom_sf(data=RUS,colour='gray') +
  coord_sf(xlim = st_bbox(resMmurr)[c(1, 3)], st_bbox(resMmurr)[c(2, 4)],
           expand = FALSE) +
  labs(title="Vessel Traffic Overlap with Top 90% of Murre Areas")

ggplot() +
  geom_sf(data=resMguild,aes(fill = DaysQuant)) +
  geom_sf(data=US,colour='gray') +
  geom_sf(data=RUS,colour='gray') +
  coord_sf(xlim = st_bbox(resMeid)[c(1, 3)], st_bbox(resMeid)[c(2, 4)],
           expand = FALSE) +
  labs(title="Vessel Traffic Overlap with Top 95% of Eider Areas")


ggplot()+
  geom_point(data=resMguild,aes(x=EiderQuant,y=DaysAll))




ggplot() +
  geom_sf(data=resEguild,aes(fill = MurreQuant)) +
  geom_sf(data=US,colour='gray') +
  geom_sf(data=RUS,colour='gray') +
  coord_sf(xlim = st_bbox(resMguild)[c(1, 3)], st_bbox(resMguild)[c(2, 4)],
           expand = FALSE)



##
##          Stack Vessel Traffic
##

library(zip)
setwd("/Users/bensullender/Documents/KickstepApproaches_Projects/FWS_AIS/Seabird_Analysis/Hex")
hexList <- as.list(list.files())
lapply(hexList,unzip)


setwd("/Users/bensullender/Documents/KickstepApproaches_Projects/FWS_AIS/Seabird_Analysis/Hex")
hexList <- as.data.frame(list.files(pattern=".shp"))
names(hexList) <- c("name")

hexE <- hexList %>%
  filter(c(substr(name,15,16) %in% c("05","06")))

hexM <- hexList %>%
  filter(c(substr(name,15,16) %in% c("07","08")))

hexL <- hexList %>%
  filter(c(substr(name,15,16) %in% c("09","10","11")))

hM1 <- st_read(getwd(),substr(hexM[1,1],1,16)) %>%
  st_drop_geometry()
hexAll <- hM1
for (h in 2:nrow(hexM)){
  hMn <- st_read(getwd(),substr(hexM[h,1],1,16)) %>%
    st_drop_geometry()
  hexAll <- rbind(hexAll,hMn)
}

hexMres <- hexAll %>%
  group_by(hexID) %>%
  summarize(DaysAll = sum(nOprD_A,na.rm=T),DaysCargo = sum(nOprD_C,na.rm=T),DaysTanker = sum(nOprD_T,na.rm=T),
            DaysFishing = sum(nOprD_F,na.rm=T),DaysOther = sum(nOprD_O,na.rm=T))


##
##      Survey Effort
##


## Omit hexes with too low survey effort
hexMask <- st_read("/Users/bensullender/Documents/KickstepApproaches_Projects/FWS_AIS/Seabird_Analysis","hex_x_ocean") %>%
  select(hexID) %>%
  mutate(AreaKM = c(st_area(.)/1000000)) %>%
  st_drop_geometry()

effortThreshold <- c(0.0)

resMeff <- resM %>%
  left_join(y=hexMask,by="hexID") %>%
  filter(as.numeric(survEff)>c(effortThreshold*as.numeric(AreaKM)))

print(paste0(c(effortThreshold*100),"% gives ",nrow(resMeff)," hexes."))


##
##            Some charts
##


# Calculating sum of birds per species per season
obsEna <- obsE %>%
  filter(!is.na(hexID)) %>%
  group_by(NPPSD.4.Letter.Code) %>%
  summarize(totalObs = sum(as.numeric(Number)))

obsMna <- obsM %>%
  filter(!is.na(hexID)) %>%
  group_by(NPPSD.4.Letter.Code) %>%
  summarize(totalObs = sum(as.numeric(Number)))

obsLna <- obsL %>%
  filter(!is.na(hexID)) %>%
  group_by(NPPSD.4.Letter.Code) %>%
  summarize(totalObs = sum(as.numeric(Number)))

#write.csv(obsEna,"obsEna.csv")
#write.csv(obsMna,"obsMna.csv")
#write.csv(obsLna,"obsLna.csv")

