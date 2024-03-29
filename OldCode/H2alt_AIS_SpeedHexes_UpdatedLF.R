################################################################################
# TITLE: Speed Hex Data Processing
# PURPOSE: This script takes raw AIS data and processes it into summary statistics
  # aggregated into hexagonal polygons (generated by H1_HexGridCreation.R script).
  # Summary statistics include number of ships, number of operating days, average
  # speed, and sd of speed on a monthly scale for cargo, fishing, tanker, other,
  # long-fast, and all vessels combined.
# AUTHOR: Ben Sullender & Kelly Kapsar
# CREATED: 2021
# LAST UPDATED ON:
# 
# NOTE: This script contains an alternate version of the long, fast ship calculation.
  # This version calculates long, fast ships as those with over 50% of signals
  # within a given hex on a given day at a speed of >10 knots (according to SOG
  # field in AIS transmission). THIS CODE HAS NOT YET BEEN RUN/IMPLEMENTED.
# NOTE 2: CODE DESIGNED TO RUN ON HPCC. 
################################################################################

# Start timer
start <- proc.time()

# Load libraries 
library(maptools)
library(rgdal)
library(dplyr)
library(tidyr)
library(tibble)
library(sf)
library(foreach)
library(doParallel)

####################################################################
##################### AIS PROCESSING FUNCTION ######################
####################################################################

# INPUTS: A list of lists containing all daily csv file names for one year of AIS data.  
## Inner list = file paths/names for daily AIS csvs
## Outer list = a list of months 

# OUTPUTS: 
## Monthly vector files (hex grid) containing ship traffic summary statistics 
## by ship type and ship length. 

FWS.AIS.SpeedHex <- function(csvList, hexgrid){
  # start timer 
  starttime <- proc.time()
  
  # clear variables
  AIScsv <- NA
  AIScsvDF <- NA
  AISlookup <- NA
  temp <- NA
  
  start <- proc.time()
  # read in csv (specifying only columns that we want based on position in dataframe)
  temp <-  lapply(csvList, read.csv, header=TRUE, na.strings=c("","NA"),
                  colClasses = c(rep("character", 2), "NULL", "character", "NULL", "NULL", 
                                 "character", rep("NULL", 6), "character", "NULL", rep("character",8),
                                 "NULL", "character", "NULL", "character", "NULL", rep("character", 2), rep("NULL", 109)
                  ))
  
  
  AIScsv <- do.call(rbind , temp)
  
  # Count number of rows
  all_pts <- length(AIScsv$Latitude)
  
  importtime <- (proc.time() - start)[[3]]/60
  start <- proc.time()
  
  # Create AIS_ID field
  AIScsv <- AIScsv %>% add_column(AIS_ID = paste0(AIScsv$MMSI,"-",substr(AIScsv$Time,1,8)))

  # Convert character columns to numeric as needed
  numcols <- c(1:2, 6:12, 14:17)
  AIScsv[,numcols] <- lapply(AIScsv[,numcols], as.numeric)
  
  # this will come in handy later. chars 28 to 34 = "yyyy-mm"
  MoName <- substr(csvList[[1]][1],45, 51)
  yr <- substr(MoName, 1, 4) 
  mnth <- substr(MoName, 6, 7)
  print(paste0("Processing",yr, mnth))
  
  # create df
  # we only care about lookup col, time, lat + long, and non-static messages
  AIScsvDF_step1 <- AIScsv %>%
    dplyr::select(MMSI,Latitude,Longitude,Time,Message_ID,AIS_ID,SOG) %>%
    filter(Message_ID!=c(5,24)) 
  
  # Calculate original number of ships with position information
  # Excludes ships for which we just received static information 
  posit_aisids <- length(unique(AIScsvDF_step1$AIS_ID))
  posit_mmsis <- length(unique(AIScsvDF_step1$MMSI))
  posit_pts <- length(AIScsvDF_step1$AIS_ID)
  
  # Remove positions with incorrect lat/long/mmsi
  AIScsvDF_step2 <- AIScsvDF_step1 %>% 
    filter(!is.na(Latitude)) %>%
    filter(!is.na(Longitude)) %>%
    filter(nchar(trunc(abs(MMSI))) > 8)
  
  
  # Convert to spatial data, get coords in projected format, and intersect with hex grid 
  AIScsvDF <- AIScsvDF_step2 %>%
    st_as_sf(coords=c("Longitude","Latitude"),crs=4326) %>%
    # project into Alaska Albers (or other CRS that doesn't create huge gap in mid-Bering with -180W and 180E)
    st_transform(crs=3338) %>%
    st_join(hexgrid["hexID"]) %>% 
    mutate(long =  st_coordinates(.)[,"X"], lat =  st_coordinates(.)[,"Y"]) %>% 
    st_drop_geometry()
  
  # Fix time stamp 
  AIScsvDF <- AIScsvDF %>% 
    mutate(Time = as.POSIXct(Time, format="%Y%m%d_%H%M%OS")) %>% 
    arrange(Time) 
  
  dftime <- (proc.time() - start)[[3]]/60
  start <- proc.time()
  
  # Remove remaining duplicate rows of data 
  AIScsvDF <- AIScsvDF %>% distinct(.keep_all=TRUE)
  
  # Calculate original number of ships with accurate position information (i.e. valid lat/lon/mmsu)
  # Excludes ships for which we just received static information 
  filt_aisids <- length(unique(AIScsvDF$AIS_ID))
  filt_mmsis <- length(unique(AIScsvDF$MMSI))
  filt_pts <- length(AIScsvDF$MMSI)
  
  distincttime <- (proc.time() - start)[[3]]/60
  start <- proc.time()
  
  # Remove points outside of hex grid 
  # system.time(joinTest <- joinTest_step1[!is.na(joinTest_step1$hexID),]) # Slow
  # system.time(na.omit(joinTest_step1[,c("hexID")])) # Slowest
  # system.time(joinTest_step1[-which(is.na(joinTest_step1$hexID)),]) # Slower 
  AISspeed <- subset(AIScsvDF, !is.na(AIScsvDF$hexID)) # Fastest
  
  
  inhex_aisids <- length(unique(AISspeed$AIS_ID))
  inhex_mmsis <- length(unique(AISspeed$MMSI))
  inhex_pts <- length(AISspeed$MMSI)
  
  nothextime <- (proc.time() - start)[[3]]/60
  start <- proc.time()
  # Function from: https://www.reddit.com/r/rstats/comments/8czqni/i_have_spatiotemporal_movement_data_how_would_i/
  # Yes. I found it on reddit. Please don't judge. I'm pretty sure it works. (See SpeedFilterTesting_20210433.R for verification
  # and comparison with previous method)  
  euclidean_speed <- function(lat2, lat1, long2, long1, time2, time1) {
    latdiff <- lat2 - lat1
    longdiff <- long2 - long1
    distance <- sqrt(latdiff^2 + longdiff^2)/1000
    timediff <- as.numeric(difftime(time2,time1,units=c("hours")))
    return(distance / timediff)
  }
  
  # Calculate euclidean speed############################################################
  AISspeed <- AISspeed %>% 
    group_by(AIS_ID) %>%
    arrange(AIS_ID, Time) %>% 
    mutate(speed = euclidean_speed(lat, lag(lat), long, lag(long), Time, lag(Time)))
  
  # Filter out points with euclidean speed > 100 km/hr 
  AISspeed <- AISspeed %>% filter(speed < 100 | is.na(speed))
  
  # Metadata collection
  speed_aisids <- length(unique(AISspeed$AIS_ID))
  speed_mmsis <- length(unique(AISspeed$MMSI))
  speed_pts <- length(AISspeed$MMSI)
  
  speedtime <- (proc.time() - start)[[3]]/60
  start <- proc.time()
  
  # create lookup table
  # we only care about 7 columns in total: lookup col (MMSI + date), name + IMO (in case we have duplicates / want to do an IMO-based lookup in the future),
  #     ship type, and size (in 3 cols: width, length, Draught)
  #     I'm including Destination and Country because that would be dope! We could do stuff with innocent passage if that's well populated.
  AISlookup <- AIScsv %>%
    # add_column(DimLength = AIScsv$Dimension_to_Bow+AIScsv$Dimension_to_stern, DimWidth = AIScsv$Dimension_to_port+AIScsv$Dimension_to_starboard) %>%
    dplyr::select(MMSI, Message_ID, Country, Vessel_Name, IMO, Ship_Type, Draught, Destination, Navigational_status,
                  SOG, AIS_ID, Dimension_to_Bow, Dimension_to_stern, Dimension_to_port, Dimension_to_starboard) %>%
    filter(Message_ID==c(5,24)) %>%
    filter(nchar(trunc(abs(MMSI))) > 8) %>% 
    distinct(AIS_ID, .keep_all=TRUE)
  
  # Calculate number of ships with 0 values in dimensions and convert to NA
  nolength <- length(which(is.na(AISlookup$Dimension_to_Bow | AISlookup$Dimension_to_stern)))
  zerolength <- which(AISlookup$Dimension_to_Bow == 0 | AISlookup$Dimension_to_stern == 0)
  pctzerolength <- round((nolength + length(zerolength))/length(AISlookup$Dimension_to_Bow)*100,2)
  
  nowidth <- length(which(is.na(AISlookup$Dimension_to_port | AISlookup$Dimension_to_starboard)))
  zerowidth <- which(AISlookup$Dimension_to_port == 0 | AISlookup$Dimension_to_starboard == 0)
  pctzerowidth <- round((nowidth + length(zerowidth))/length(AISlookup$Dimension_to_port)*100,2)

  # Remove zero value rows for consideration of respective measurement
  # (i.e. if either bow or stern is zero, then both get NA 
  # and if either port or starboard is zero then both get NA)
  AISlookup$Dimension_to_Bow[zerolength] <- NA
  AISlookup$Dimension_to_stern[zerolength] <- NA
  AISlookup$Dimension_to_port[zerowidth] <- NA
  AISlookup$Dimension_to_starboard[zerowidth] <- NA
  
  # Create length and width values from dimensions 
  AISlookup <- AISlookup %>% 
    add_column(DimLength = AISlookup$Dimension_to_Bow+AISlookup$Dimension_to_stern, 
             DimWidth = AISlookup$Dimension_to_port+AISlookup$Dimension_to_starboard) %>% 
    dplyr::select(MMSI, Message_ID, Country, Vessel_Name, IMO, Ship_Type, Draught, Destination, Navigational_status,
                  SOG, AIS_ID, DimLength, DimWidth)
  
  # Calculate how many ships missing attributes and vice versa
  InSpeedNotLookup_mmsis <- length(unique(AISspeed$MMSI[!(AISspeed$MMSI %in% AISlookup$MMSI)]))
  InSpeedNotLookup_npts <- length(AISspeed$MMSI[!(AISspeed$MMSI %in% AISlookup$MMSI)])
  InLookupNotSpeed_mmsis <- length(unique(AISlookup$MMSI[!(AISlookup$MMSI %in% AISspeed$MMSI)]))
  # InLookupNotSpeed_npts <- length(unique(AISlookup$MMSI[!(AISlookup$MMSI %in% AISspeed$MMSI)]))
  
  # step 2: join lookup table to the points 
  AISjoined <- AISspeed %>%
    left_join(AISlookup,by="AIS_ID")
  
  # step 3: split lines by ship type
  # link to ship type/numbers table: 
  # https://help.marinetraffic.com/hc/en-us/articles/205579997-What-is-the-significance-of-the-AIS-Shiptype-number-

  # Remove NA ship type - mostly signals from static platforms (e.g., terrestrial AIS receivers)
  AISjoined <- AISjoined[!is.na(AISjoined$Ship_Type),]

  #  ID other ship types
  AISjoined$AIS_Type <- ifelse(substr(AISjoined$Ship_Type,1,1)==7, "Cargo",
                               ifelse(substr(AISjoined$Ship_Type,1,1)==8, "Tanker", 
                                      ifelse(substr(AISjoined$Ship_Type,1,2)==30, "Fishing", "Other")))
  
  jointime <- (proc.time() - start)[[3]]/60
  start <- proc.time()
  
  # Number of points by ship type
  ntank_pts <- length(which(AISjoined$AIS_Type == "Tanker"))
  ntank_mmsis <- length(unique(AISjoined$MMSI.x[which(AISjoined$AIS_Type == "Tanker")]))
  ntank_aisids <- length(unique(AISjoined$AIS_ID[which(AISjoined$AIS_Type == "Tanker")]))
  
  nfish_pts <- length(which(AISjoined$AIS_Type == "Fishing"))
  nfish_mmsis <- length(unique(AISjoined$MMSI.x[which(AISjoined$AIS_Type == "Fishing")]))
  nfish_aisids <- length(unique(AISjoined$AIS_ID[which(AISjoined$AIS_Type == "Fishing")]))
  
  ncargo_pts <- length(which(AISjoined$AIS_Type == "Cargo"))
  ncargo_mmsis <- length(unique(AISjoined$MMSI.x[which(AISjoined$AIS_Type == "Cargo")]))
  ncargo_aisids <- length(unique(AISjoined$AIS_ID[which(AISjoined$AIS_Type == "Cargo")]))
  
  nother_pts <- length(which(AISjoined$AIS_Type == "Other"))
  nother_mmsis <- length(unique(AISjoined$MMSI.x[which(AISjoined$AIS_Type == "Other")]))
  nother_aisids <- length(unique(AISjoined$AIS_ID[which(AISjoined$AIS_Type == "Other")]))
  
  ntotal_pts <- length(AISjoined$AIS_Type)
  ntotal_mmsis <- length(unique(AISjoined$MMSI.x))
  ntotal_aisids <- length(unique(AISjoined$AIS_ID))
  
  # Calculate total % position data removed from initial to final data set
  pctremoved <- round((length(unique(AIScsvDF_step1$MMSI)) - length(unique(AISjoined$MMSI.x)))/length(unique(AIScsvDF_step1$MMSI))*100, 2)
  
  # Thrown out for missing/incorrect lat/lon/MMSI, duplicate points, speed > 100 km/hr, outisde hex grid
  pctmissingwidth <- round(sum(is.na(AISjoined$DimWidth))/length(AISjoined$DimWidth)*100, 2)
  pctmissinglength <- round(sum(is.na(AISjoined$DimLength))/length(AISjoined$DimLength)*100, 2)
  pctmissingSOG <-  round(sum(is.na(AISjoined$SOG.x))/length(AISjoined$DimLength)*100, 2)
  
  # # Loop through each ship type and calculate summary statistics
  allTypes <- unique(AISjoined$AIS_Type)

  # Calculate summary stats for each ship type
  for (k in 1:length(allTypes)){
    # Select ship type
    AISfilteredType <- AISjoined %>%
      filter(AIS_Type==allTypes[k])
    
    # Calculate average speed within hex grid 
    joinOut <- AISfilteredType %>%
      group_by(hexID) %>%
      summarize(SOG_kts=mean(SOG.x, na.rm=T),
                SOGsd = sd(SOG.x, na.rm=T),
                EucSpeed_kts = mean(speed, na.rm=T)/1.852, # converted to knots 
                EucSpeedsd = sd(speed, na.rm=T),
                ShipLength = mean(DimLength, na.rm=T), 
                nPts=n(), 
                nMMSI=length(unique(MMSI.x)),
                nOperDays=length(unique(AIS_ID)))
                
    colnames(joinOut)[2:ncol(joinOut)] <- paste0(colnames(joinOut)[2:ncol(joinOut)],"_",substring(allTypes[k],1,1))
    
    hexgrid <- left_join(hexgrid, joinOut, by="hexID")
  }

  # Calculate summary stats for ships > 65 feet length with at least 50% of points per hex at > 10 knots
  # Create id for each unique ship in each unique hex in each day 
  AISjoined$aishex_ID <- paste0(AISjoined$AIS_ID, "-", AISjoined$hexID)
  
  longfast_aishex_IDs <- AISjoined %>% 
    # Group all points by unique ship/hex/day combo 
    group_by(aishex_ID) %>% 
    # Remove NA SOG values
    filter(!is.na(SOG.x)) %>% 
    # Calculate percentage of transmissions from each ship in each hex in each day that are > 10 knots
    summarize(pctfast = round(sum(SOG.x > 10)/length(SOG.x), 2)) %>% 
    # Isolate out ids from only those with >50% of transmission at > 10 knots speed
    filter(pctfast > 0.5)
  
  # Isolate all points from the identified long fast ids from the main data set 
  fastlong <- AISjoined[which(AISjoined$aishex_ID %in% longfast_aishex_IDs$aishex_ID),]
  
  nfastlong_pts <- length(fastlong$AIS_Type)
  nfastlong_mmsis <- length(unique(fastlong$MMSI.x))
  nfastlong_aisids <- length(unique(fastlong$AIS_ID))
  
  # Calculate total number of points, ships, and operating days from long, fast ships for each 
  fastlong <- fastlong %>% 
    group_by(hexID) %>%
    summarize(SOG_kts=mean(SOG.x, na.rm=T),
              SOGsd = sd(SOG.x, na.rm=T),
              ShipLength = mean(DimLength, na.rm=T), 
              nPts=n(), 
              nMMSI=length(unique(MMSI.x)),
              nOperDays=length(unique(AIS_ID)))
  
  colnames(fastlong)[2:ncol(fastlong)] <- paste0(colnames(fastlong)[2:ncol(fastlong)],"_LongFast")
  
  hexgrid <- left_join(hexgrid, fastlong, by="hexID")
    
  # Calculate summary stats for all ship types in aggregate
  # Calculate average speed within hex grid 
  allships <- AISjoined %>%
    group_by(hexID) %>%
    summarize(SOG_kts=mean(SOG.x, na.rm=T),
              SOGsd = sd(SOG.x, na.rm=T),
              Speed_kmh = mean(speed, na.rm=T),
              Speedsd = sd(speed, na.rm=T),
              ShipLength = mean(DimLength, na.rm=T),  
              nPts=n(), 
              nMMSI=length(unique(MMSI.x)),
              nOperDays=length(unique(AIS_ID)))
  
  colnames(allships)[2:ncol(allships)] <- paste0(colnames(allships)[2:ncol(allships)],"_All")
  
  hexgrid <- left_join(hexgrid, allships, by="hexID")
  
  hexpts <- st_as_sf(AISjoined, coords = c("long", "lat"), crs = 3338)
  
  # Save data in vector format
  # write_sf(hexgrid, paste0("../Data_Processed_TEST/Hex/SpeedHex_",MoName,"_",ndays,".shp"))
  write_sf(hexgrid, paste0("../Data_Processed_TEST/Hex/SpeedHex_",MoName,".shp"))
  # write_sf(AISjoined, paste0("../Data_Processed_TEST/Hex/SpeedPts_",MoName,"_",ndays,".shp"))
  # write_sf(hexpts, paste0("../Data_Processed_TEST/Hex/SpeedPts_",MoName,".shp"))
  
  
  # Save processing info to text file 
  hextime <- (proc.time() - start)[[3]]/60
  
  runtime <- proc.time() - starttime 
  runtime_min <- runtime[[3]]/60 
  summarystats <- data.frame(cbind(yr, mnth, runtime_min,all_pts, 
                                   posit_aisids, posit_mmsis, posit_pts,
                                   filt_aisids, filt_mmsis, filt_pts,
                                   inhex_aisids, inhex_mmsis, inhex_pts,
                                   speed_aisids, speed_mmsis, speed_pts,
                                   InSpeedNotLookup_mmsis, InSpeedNotLookup_npts, InLookupNotSpeed_mmsis,
                                   ntank_pts, ntank_aisids, ntank_mmsis,
                                   nfish_pts, nfish_aisids, nfish_mmsis,
                                   ncargo_pts, ncargo_aisids, ncargo_mmsis,
                                   nother_pts, nother_aisids, nother_mmsis,
                                   nfastlong_pts, nfastlong_aisids, nfastlong_mmsis,
                                   ntotal_pts,ntotal_mmsis, ntotal_aisids,
                                   pctremoved, pctmissingwidth, pctmissinglength, pctmissingSOG))
  write.csv(summarystats, paste0("../Data_Processed_TEST/Hex/Metadata_SpeedHex_",MoName,".csv"))
  
  
  runtimes <- data.frame(cbind(yr, mnth, runtime_min, importtime, dftime, distincttime, nothextime, speedtime, jointime, hextime)) 
  write.csv(runtimes, paste0("../Data_Processed_TEST/Hex/Runtimes_SpeedHex_",MoName,".csv"))
  # write.csv(runtimes, paste0("../Data_Processed_TEST/Hex/Runtimes_SpeedHex_",MoName,"_",ndays,".csv"))
  # print(runtimes)
  # return(runtimes)
}

####################################################################
####################### RUNNING SPEED SCRIPT ####################### 
####################################################################


# # Import hex grid
hexgrid <- st_read("../Data_Raw/HexGrid/RevisedHexes.shp")

# Pull up list of AIS files
files <- paste0("../Data_Raw/2020/", list.files("../Data_Raw/2020", pattern='.csv'))

# Separate file names into monthly lists
jan <- files[grepl("-01-", files)]

csvList <- jan[1]

# Run the speed hex creation script
starthex1 <- proc.time()
FWS.AIS.SpeedHex(csvList, hexgrid)
endhex1 <- proc.time() - starthex1


browseURL("https://www.youtube.com/watch?v=K1b8AhIsSYQ&ab_channel=RHINO")


temp <- st_read("../Data_Processed_TEST/Hex/SpeedHex_2020-01.shp")

####################################################################
####################### PARALLELIZATION CODE ####################### 
####################################################################

# Pull up list of AIS files
files <- paste0("../Data_Raw/2015/", list.files("../Data_Raw/2015", pattern='.csv'))

# Separate file names into monthly lists
jan <- files[grepl("-01-", files)]
feb <- files[grepl("-02-", files)]
mar <- files[grepl("-03-", files)]
apr <- files[grepl("-04-", files)]
may <- files[grepl("-05-", files)]
jun <- files[grepl("-06-", files)]
jul <- files[grepl("-07-", files)]
aug <- files[grepl("-08-", files)]
sep <- files[grepl("-09-", files)]
oct <- files[grepl("-10-", files)]
nov <- files[grepl("-11-", files)]
dec <- files[grepl("-12-", files)]

# Create a list of lists of all csv file names grouped by month
csvsByMonth <- list(jan, feb, mar, apr, may, jun, jul, aug, sep, oct, nov, dec)

# Pull up hex grid
hexgrid <- st_read("../Data_Raw/HexGrid/RevisedHexes.shp")

## MSU HPCC: https://wiki.hpcc.msu.edu/display/ITH/R+workshop+tutorial#Rworkshoptutorial-Submittingparalleljobstotheclusterusing{doParallel}:singlenode,multiplecores
# Request a single node (this uses the "multicore" functionality)
registerDoParallel(cores=as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE")[1]))

# create a blank list to store the results (I truncated the code before the ship-type coding, and just returned the sf of all that day's tracks so I didn't 
#       have to debug the raster part. If we're writing all results within the function - as written here and as I think we should do - the format of the blank list won't really matter.)
res=list()

# foreach and %dopar% work together to implement the parallelization
# note that you have to tell each core what packages you need (another reason to minimize library use), so it can pull those over
# I'm using tidyverse since it combines dplyr and tidyr into one library (I think)
res=foreach(i=1:12,.packages=c("maptools", "rgdal", "dplyr", "tidyr", "tibble", "sf", "foreach", "doParallel"),
            .errorhandling='pass',.verbose=T,.multicombine=TRUE) %dopar% FWS.AIS.SpeedHex(csvList=csvsByMonth[[i]], hexgrid=hexgrid)
# lapply(csvsByMonth, FWS.AIS)

# Elapsed time and running information
tottime <- proc.time() - start
tottime_min <- tottime[[3]]/60

cat("Time elapsed:", tottime_min, "\n")
cat("Currently registered backend:", getDoParName(), "\n")
cat("Number of workers used:", getDoParWorkers(), "\n")


########################################################################## 
############################# Post-Processing ############################
####################### Simplification of Hex Grid ####################### 
########################################################################## 

# Read in all hex files
filenames <- list.files("../Data_Processed_TEST/Hex/", pattern='.shp')
filepaths <-paste0("../Data_Processed_TEST/Hex/", 
                   list.files("../Data_Processed_TEST/Hex/", pattern='.shp'))

# Identify year and month for each file
yr <- substr(filenames, 10, 13)
mo <- substr(filenames, 15, 16)

# Read in shape fiels 
temphpcc <- lapply(filepaths, st_read)

# Add year and month as attributes 
for(i in 1:length(temphpcc)){
  print(i)
  t <- temphpcc[[i]]
  t <- t %>% mutate(year = yr[i], month = mo[i])
  temphpcc[[i]] <- t
}

# Select columns of interest and rename to standard names
newshp <- lapply(temphpcc, function(x){x %>% 
    select(hexID, year, month,
          SOG_k_C, SOGsd_C, nMMSI_C, nOprD_C,
          SOG_k_T, SOGsd_T, nMMSI_T, nOprD_T,
          SOG_k_O, SOGsd_O, nMMSI_O, nOprD_O,
          SOG_k_F, SOGsd_F, nMMSI_F, nOprD_F,
          SOG__LF, SOGs_LF, nMMSI_L, nOpD_LF,
          SOG_k_A, SOGsd_A, nMMSI_A, nOprD_A) %>% 
    rename(SOG_C = SOG_k_C,
           SOG_T = SOG_k_T,
           SOG_O = SOG_k_O, 
           SOG_F = SOG_k_F, 
           SOG_L = SOG__LF,
           SOG_A = SOG_k_A, 
           SOGsd_L = SOGs_LF, 
           nOPrD_L = nOpD_LF)})
# Save new files 
for(i in 1:length(newshp)){
  t <- newshp[[i]]
  st_write(t, paste0("../Data_Processed_TEST/Hex/Simplified/SpeedHex_",yr[i],"-",mo[i],".shp"))
}

