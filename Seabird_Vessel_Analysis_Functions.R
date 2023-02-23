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
    mutate(date = lubridate::as_datetime(local_date_time), 
           year = lubridate::year(date), 
           month = lubridate::month(date)) %>% 
    filter(year > startyr) %>%
    # Removing off-transect from location data
    filter(modified_survey_type != "Off Transect Observation") %>%
    droplevels()
  
  datobs <- datobs %>%
    # Removing off-transect from observations
    filter(on_of_tx != "Off") %>%
    left_join(x=datobs,y=loc,by="master_key") %>%
    filter(year > startyr) %>%
    droplevels()
  
  # Identify unique 4-digit codes for each bird spp 
  birdNames <- unique(datobs$species_code)
  
  # convert to SF and transform to Alaska Albers (epsg 3338), otherwise hexagons over -180:180 get warped
  datSF <- st_as_sf(datobs,coords = c("longitude","latitude"),crs=4326) %>%
    st_transform(crs=3338)
  locSF <- st_as_sf(loc,coords = c("longitude","latitude"),crs=4326) %>%
    st_transform(crs=3338)
  
  # Set up one df for bird obs and another df for survey effort
  locIn <- locSF %>%
    # Remove observations from months outside of window of interest
    filter(month %in% mnths) %>%
    # Spatially join with hexagon data 
    st_join(hex) %>%
    st_drop_geometry()
  
  obsIn <- datSF %>%
    filter(month %in% mnths)  %>%
    st_join(hex) %>%
    st_drop_geometry() %>%
    filter(number!="")
  
  res <- hex
  sppNout <- as.list(NA)
  
  # Iterate through each spp and calculate number of observations per hex
  for (i in 1:length(birdNames)){
    obsIn2 <- obsIn[obsIn$species_code==birdNames[i],]
    sppNout[[i]] <- obsIn2 %>%
      group_by(hexID) %>%
      summarize(birdN = sum(as.numeric(number))) %>%
      as.data.frame()
    names(sppNout[[i]]) <- c("hexID",birdNames[i])
  }
  
  for (j in 1:length(birdNames)){
    res <- left_join(res,sppNout[[j]],by="hexID")
  }
  
  # Calculate survey effort within each hex 
  surveyEff <- locIn %>%
    group_by(hexID) %>%
    summarize(survEff = sum(as.numeric(sample_area))) %>%
    as.data.frame()
  
  res <- left_join(res,surveyEff,by="hexID")
  
  # switch NAs to 0s
  res[is.na(res)] <- 0
  
  # Save output 
  st_write(res,survfilename)
}

#### Vessel Activity #### 

hexStack <- function(filedir, mnths, metric, night, trafffilename){
  
  # Isolate month values of interest from file names 
  hexes <- filedir[as.numeric(substr(filedir,10,11)) %in% mnths]
  
  # Read in data 
  temp <- lapply(hexes, function(x){st_read(paste0(hexdir,x)) %>% st_drop_geometry()})
  hexAll <- do.call(rbind, temp)
  
  # Select metrics based on inputs 
  lab <- ifelse(metric=="OperatingDays", "OpD", 
                ifelse(metric=="Ships","Shp", NA))
  lab <- ifelse(night==TRUE, paste0("N_",lab, "_Al"),paste0(lab, "_Al"))
  
  test <- hexAll %>% dplyr::select(as.character(lab))
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
                              figfolder,
                              filedir,
                              studyarea,
                              metric, 
                              night){
  
  # Make sure appropriate vessel data exist and if not, generate it
  trafffilename <- paste0(savefolder,"TraffInHexes_",mnthsnam,"_NightOnly", night,".csv")
  
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
  finaldfname <- paste0(savefolder,"FinalDF_",taxaLabel,"_",monthsname,"_NightOnly",night[1],".shp")
  
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
    
    # Calculate SD categories for each hex's effort-weighted seabird observations 
    resFinal <- resStacked %>%
      mutate(QuantBird = ecdf(AllBird/survEff)(AllBird/survEff))
    
    # Calculate effort-weighted densities of seabirds (# birds per square km of surveyed area)
    resFinal$DensBird <- c(resFinal$AllBird/resFinal$survEff)
    
    # Calculate effort-weighted classes for the number of observations 
    resFinal$ClassBird <- ifelse(c(resFinal$AllBird/resFinal$survEff)<mean(c(resFinal$AllBird/resFinal$survEff)),1,
                             ifelse(c(resFinal$AllBird/resFinal$survEff)>=c(mean(c(resFinal$AllBird/resFinal$survEff))+sd(c(resFinal$AllBird/resFinal$survEff))),3,2))
    
    # Create column to specify taxa of interest
    resFinal$taxa <- taxaLabel
    
    # Join to vessel traffic data 
    finalCombined <- resFinal %>%
      left_join(y=hexFinal,by="hexID")
    
    # Evaluate risk levels 
    finalCombined$risk <- ifelse(finalCombined$ClassBird == 1 & finalCombined$ClassShip == 1, "low",
                ifelse(finalCombined$ClassBird == 1 & finalCombined$ClassShip == 2 | finalCombined$ClassBird == 2 & finalCombined$ClassShip == 1, "medium",
                ifelse(finalCombined$ClassBird == 1 & finalCombined$ClassShip == 3 | finalCombined$ClassBird == 3 & finalCombined$ClassShip == 1, "medium",
                ifelse(finalCombined$ClassBird == 3 & finalCombined$ClassShip == 2 | finalCombined$ClassBird == 2 & finalCombined$ClassShip == 3, "high",
                ifelse(finalCombined$ClassBird == 2 & finalCombined$ClassShip == 2, "high",
                ifelse(finalCombined$ClassBird == 3 & finalCombined$ClassShip == 3, "veryhigh", NA))))))
    
    finalCombined$risk <- factor(finalCombined$risk,  c("low","medium","high","veryhigh"))
    
    # Save data 
    st_write(finalCombined, paste0(savefolder,"FinalDF_",taxaLabel,"_",monthsname,"_NightOnly",night[1],".shp"))
    
    ######## Plot results
    studyareanew <- studyarea %>% st_crop(st_buffer(finalCombined, 100000))
    
    p1 <- ggplot() +
      geom_sf(data=studyareanew, fill="#fefeff", lwd=0) +
      geom_sf(data=finalCombined,aes(fill = risk)) +
      scale_fill_manual(values = c("low" = "#73b2ff",
                                   "medium" = "#55fe01",
                                   "high" = "#ffff01",
                                   "veryhigh" = "#e31a1c"), 
                        name="Risk", 
                        labels = c("Low", "Medium", "High", "Very High")) +
      xlab("") +
      ylab("") +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      ggtitle(ifelse(night == FALSE, paste0(taxaLabel, " - ", monthsname, " - All Traffic"), 
                     paste0(taxaLabel, " - ", monthsname, " - Night Traffic"))) +
      theme_bw() +
      theme(text = element_text(size = 25),
            axis.text=element_blank(),
            panel.background = element_rect(fill = "#bcc7dd"),
            panel.border =  element_rect(colour = "black"),
            panel.grid.major = element_line(colour = "transparent"))
    
    # Save figure 
    ggsave(filename = paste0(figfolder,"Map_",taxaLabel,"_",monthsname,"_NightOnly",night[1],".png"), 
           plot = p1, width=12, height=8, units="in")
    
    
    ######### Bivariate color palette plot
    finalNoBird <- finalCombined %>% dplyr::filter(DensBird == 0)
    finalBird <- finalCombined %>% dplyr::filter(DensBird > 0)
    
    finalBird <- bi_class(finalBird, x = DensBird, y = AllShip, style = "quantile", dim = 3)
    
    p2 <- ggplot() +
      geom_sf(data=studyareanew, fill="#fefeff", lwd=0) +
      geom_sf(data=finalNoBird, fill="black") + 
      geom_sf(data=finalBird,aes(fill = bi_class), show.legend = FALSE) +
      bi_scale_fill(pal = "DkViolet2", dim = 3, flip_axes = TRUE) +
      # xlab("") +
      # ylab("") +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      ggtitle(ifelse(night == FALSE, paste0(taxaLabel, " - ", monthsname, " - All Traffic"), 
                     paste0(taxaLabel, " - ", monthsname, " - Night Traffic"))) +
      bi_theme() +
      theme_bw() +
      theme(text = element_text(size = 25),
            axis.text=element_blank(),
            panel.background = element_rect(fill = "black"),
            panel.border =  element_rect(colour = "black"),
            panel.grid.major = element_line(colour = "transparent")) 
    
    legend <- bi_legend(pal = "DkViolet2",
                        dim = 3,
                        xlab = "Bird Density",
                        ylab = "Vessel Traffic",
                        size = 14, flip_axes = TRUE)
    
    # Create white rectangles to obscure error boxes where arrows should be on legend
    rect <- rectGrob(
      x=0.57, 
      y=.61,
      width=0.015, height=0.015,
      gp = gpar(fill = "white", alpha = 1, col="white")
    )
    rect2 <- rectGrob(
      x=0.675, 
      y=.425,
      width=0.015, height=0.02,
      gp = gpar(fill = "white", alpha = 1, col="white")
    )
    
    p2Final <- ggdraw() +
      draw_plot(p2, 0, 0, 1, 1) +
      draw_plot(legend, x=0.50, y=.4, width=0.25, height=0.25) +
      draw_grob(rect) +
      draw_grob(rect2)
    
    # Save figure 
    ggsave(filename = paste0(figfolder,"BiVariateMap_",taxaLabel,"_",monthsname,"_NightOnly",night,".png"), 
           plot = p2Final, width=10, height=8, units="in")
    
    }
}
