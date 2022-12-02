################################################################################
# TITLE: AIS x Seabirds Functions 
#
# PURPOSE: This script contains the functions used to calculate the overlap 
# between vessel activity and seabird distributions. 
#
# AUTHORS: Ben Sullender & Kelly Kapsar, with guidance from Kathy Kuletz
#
# CREATED: 13 July 2022
# LAST UPDATED ON: 2 December 2022
################################################################################

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

