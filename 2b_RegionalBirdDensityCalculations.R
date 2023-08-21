RegionalDensity <- function(loc, datobs, hex, startyr, mnths, survfilename, taxaNames, studyareahexIDs, effortThreshold){
  # Drop observations that are too old or off-transect
  # Based on Kathy Kuletz's feedback, prior to 2007 is too old.
  newloc <- loc %>%
    mutate(date = lubridate::as_datetime(local_date_time), 
           year = lubridate::year(date), 
           month = lubridate::month(date)) %>% 
    filter(year > startyr) %>%
    # Removing off-transect from location data
    filter(modified_survey_type != "Off Transect Observation") %>%
    droplevels()
  
  # newdatobs <- datobs %>%
  #   # Removing off-transect from observations
  #   filter(on_off_tx == "ON") %>%
  #   left_join(.,y=newloc,by="master_key") %>%
  #   filter(year > startyr) %>%
  #   droplevels()
  
  # Identify unique 4-digit codes for each bird spp 
  birdNames <- unique(newdatobs$species_code)
  
  # convert to SF and transform to Alaska Albers (epsg 3338), otherwise hexagons over -180:180 get warped
  # datSF <- st_as_sf(newdatobs,coords = c("longitude","latitude"),crs=4326) %>%
  #   st_transform(crs=3338)
  locSF <- st_as_sf(newloc,coords = c("longitude","latitude"),crs=4326) %>%
    st_transform(crs=3338)
  
  # Set up one df for bird obs and another df for survey effort
  locIn <- locSF %>%
    # Spatially join with hexagon data 
    st_join(hex) %>%
    filter(!is.na(hexID)) %>% 
  
  summloc <- locIn %>% 
    filter(month %in% c(6,7,8)) %>% 
    group_by(hexID, AreaKM) %>% 
    summarize(survEff = sum(sample_area)) %>% 
    filter(survEff > (as.numeric(AreaKM)*effortThreshold))
  fallloc <- locIn %>% filter(month %in% c(9,10,11)) %>% 
    group_by(hexID, AreaKM) %>% 
    summarize(survEff = sum(sample_area)) %>% 
    filter(survEff > (as.numeric(AreaKM)*effortThreshold))
  
  sharedIDs <- intersect(unique(summloc$hexID), unique(fallloc$hexID))
  
  plot(summloc$geometry[summloc$hexID %in% sharedIDs])
  plot(hex$geometry[hex$hexID %in% sharedIDs])
  
  # obsIn <- datSF %>%
  #   filter(month %in% mnths)  %>%
  #   st_join(hex) %>%
  #   filter(number!="") %>% 
  #   filter(species_code %in% taxaNames) 
  
  test <- obsIn[st_intersects(obsIn, studyarea, sparse=FALSE),]
  
  sum(st_intersects(summloc[summloc$hexID %in% sharedIDs,], studyarea, sparse=FALSE))
  
  
}