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

birdDensity <- function(loc, datobs, hex, startyr, mnths, survfilename){
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
  
  newdatobs <- datobs %>%
    # Removing off-transect from observations
    filter(on_off_tx == "ON") %>%
    left_join(.,y=newloc,by="master_key") %>%
    filter(year > startyr) %>%
    droplevels()
  
  # Identify unique 4-digit codes for each bird spp 
  birdNames <- unique(newdatobs$species_code)
  
  # convert to SF and transform to Alaska Albers (epsg 3338), otherwise hexagons over -180:180 get warped
  datSF <- st_as_sf(newdatobs,coords = c("longitude","latitude"),crs=4326) %>%
    st_transform(crs=3338)
  locSF <- st_as_sf(newloc,coords = c("longitude","latitude"),crs=4326) %>%
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
  # for (i in 1:length(birdNames)){
  #   obsIn2 <- obsIn[obsIn$species_code==birdNames[i],]
  #    temp <- obsIn2 %>%
  #     group_by(hexID, master_key) %>%
  #     summarize(nBirds = sum(number), sample_area = first(sample_area)) %>%
  #     as.data.frame()
  #   sppNout[[i]] <- temp %>% 
  #     group_by(hexID) %>% 
  #     summarize(birdDen = sum(nBirds)/sum(sample_area))
  #   names(sppNout[[i]]) <- c("hexID",birdNames[i])
  # }
  for (i in 1:length(birdNames)){
    obsIn2 <- obsIn[obsIn$species_code==birdNames[i],]
    sppNout[[i]] <- obsIn2 %>% 
      group_by(hexID) %>% 
      summarize(AllBird = sum(number))
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

traffDensity <- function(filedir, mnths, metric, timeofday, trafffilename){
  
  # Isolate month values of interest from file names 
  hexes <- filedir[as.numeric(substr(filedir,10,11)) %in% mnths]
  
  # Read in data 
  temp <- lapply(hexes, function(x){st_read(paste0(hexdir,x)) %>% st_drop_geometry()})
  hexAll <- do.call(rbind, temp)
  
  # Select metrics based on inputs 
  lab <- ifelse(metric=="OperatingDays", "OpD", 
         ifelse(metric=="Ships","Shp", 
         ifelse(metric=="Hours", "Hrs", NA)))
  
  if(timeofday=="Night"){
    lab <- paste0("N_",lab, "_Al")
  }
  if(timeofday=="Day"){
    lab <- paste0("D_",lab, "_Al")
  }
  if(timeofday=="All"){
    lab <- paste0(lab, "_Al")
  }
  
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
  
  # Save results 
  write.csv(hexRes,trafffilename)
}

#### Seabird Metrics #### 
# NOTE: This function is dependent on the birdDensity and traffDensity functions 

birdHexesByEffort <- function(dataobs,
                              loc,
                              taxaNames, 
                              taxaLabel,
                              hexMask, 
                              effortThreshold, 
                              mnths,
                              mnthsnam,
                              studyarea, 
                              studyareaname,
                              startyr,
                              savefolder,
                              filedir,
                              metric, 
                              timeofday){
  
  # Make sure appropriate vessel data exist and if not, generate it
  trafffilename <- paste0(savefolder,"Traffic_Cleaned/TraffInHexes_",metric,"_",mnthsnam,"_", timeofday,".csv")
  
  if(!file.exists(trafffilename)){
    traffDensity(filedir = hexList, mnths = mnths, metric = metric, timeofday = timeofday, trafffilename)
  }
  
  traffdf <- read.csv(trafffilename)
  
  # Make sure survey effort has been calculated and saved for appropriate months and if not, generate it
  survfilename <- paste0(savefolder,"Seabird_Cleaned/ObsInHexes_",monthsname,".shp")
  
  if(!file.exists(survfilename)){
    birdDensity(loc=loc, datobs=datobs, hex=hexMask, startyr=startyr, mnths=mnths, survfilename=survfilename)
  }
  birddf <- st_read(survfilename)
  
  # Make sure this function has not already been run with these exact parameter specifications 
  finaldfname <-paste0(savefolder,"Region_Taxa_Season_TimeOfDay_DFs/FinalDF_",studyareaname,"_",taxaLabel,"_",monthsname,"_", timeofday,".shp")
  
  if(!file.exists(finaldfname)){

    # Isolate hexes within study area
    birddf <- birddf[st_contains(studyarea, birddf, sparse = FALSE),]

    ## Bird Risk Calculations
    # Select only hexes with sufficient survey effort
    birdGuild <- birddf %>%
      filter(as.numeric(survEff)>c(effortThreshold*as.numeric(AreaKM)))

    # Isolate taxa of interest and calculate total number of observations
    birdGuildAll <- birdGuild[,c(colnames(birdGuild) %in% taxaNames)] %>%
      st_drop_geometry()
    birdGuild$AllBird <- rowSums(birdGuildAll)

    # Stack just the relevant guilds, remove extraneous columns
    birdStacked <- birdGuild %>%
      select(hexID,AllBird,survEff,AreaKM)

    # Calculate effort-weighted densities of seabirds (# birds per square km of surveyed area)
    birdFinal <- birdStacked %>% mutate(DensBird = AllBird/survEff)
    
    # Calculate SD categories for each hex's effort-weighted seabird observations
    birdFinal <- birdFinal %>%
      mutate(QuantBird = ecdf(DensBird)(DensBird))

    # Calculate effort-weighted classes for the number of observations
    birdFinal$ClassBird <- ifelse(birdFinal$AllBird == 0, 0,
                            ifelse(birdFinal$DensBird<mean(birdFinal$DensBird),1,
                             ifelse(birdFinal$DensBird>=c(mean(birdFinal$DensBird)+sd(birdFinal$DensBird)),3,2)))
    # Standardize values
    # birdFinal$DensBird_rescale <- birdFinal$DensBird/sd(birdFinal$DensBird)

    # Create column to specify taxa of interest
    birdFinal$taxa <- taxaLabel

    # Select traffic data columns of interest
    traffFinal <- traffdf %>% dplyr::select(hexID, AllShip)

    # Join to vessel traffic data
    df <- birdFinal %>%
      left_join(y=traffFinal,by="hexID")

    ## Vessel Risk Calculations
    # Calculate SD categories for each hex's vessel activity
    df$QuantShip <- ecdf(df$AllShip)(df$AllShip)

    df$ClassShip <- ifelse(df$AllShip<mean(df$AllShip),1,
                                ifelse(df$AllShip>=c(mean(df$AllShip)+sd(df$AllShip)),3,2))

    # Standardize values
    # df$AllShip_rescale <- df$AllShip/sd(df$AllShip)
    
    # Evaluate risk levels
    df$risk <- ifelse(df$ClassBird == 1 & df$ClassShip == 1, "low",
                ifelse(df$ClassBird == 1 & df$ClassShip == 2 | df$ClassBird == 1 & df$ClassShip == 3, "mediumSHIP",
                ifelse(df$ClassBird == 2 & df$ClassShip == 1 | df$ClassBird == 3 & df$ClassShip == 1, "mediumBIRD",
                ifelse(df$ClassBird == 3 & df$ClassShip == 2 | df$ClassBird == 2 & df$ClassShip == 3, "high",
                ifelse(df$ClassBird == 2 & df$ClassShip == 2, "high",
                ifelse(df$ClassBird == 3 & df$ClassShip == 3, "veryhigh", NA))))))


    df$risk <- factor(df$risk,  c("low","mediumSHIP","mediumBIRD","high","veryhigh"))
    
    # df$risk_cont <- df$AllShip_rescale/df$DensBird_rescale

    df$season <- mnthsnam
    df$timeofday <- timeofday

    # Save data
    st_write(df, finaldfname)
  }
}


#### Region/Taxa Summaries 
summstats <- function(taxaNames, 
                      taxaLabel,
                      studyarea, 
                      studyareaname,
                      savefolder,
                      figfolder, 
                      basemap){

  
  dfs <- list.files(path=paste0(savefolder,"Region_Taxa_Season_TimeOfDay_DFs/"), 
                    pattern=paste0("FinalDF_", studyareaname,"_", taxaLabel), full.names=TRUE)
  dfs <- dfs[grep(".shp", dfs)]
  
  dfs <- lapply(dfs, st_read)
  df <- do.call(rbind, dfs)
  
  df <- df %>% dplyr::select(-ClassBird, -QuantBird, -ClassShip, -QuantShip, -risk)
  
  bothseasons <- df %>% st_drop_geometry() %>% group_by(hexID) %>% summarize(n=n()) %>% filter(n == 6)
  
  subdf <- df %>% filter(hexID %in% bothseasons$hexID) %>% mutate(subset = paste0(season, "_", timeofday))
  
  timeofday <- c("All", "Night")
  filtdf <- data.frame()
  
  for(i in 1:2){
    temp <- subdf[subdf$timeofday == timeofday[i],]

    temp$ClassBird <- ifelse(temp$AllBird == 0, 0, 
                      ifelse(c(temp$AllBird/temp$survEff)<mean(c(temp$AllBird/temp$survEff)),1,
                      ifelse(c(temp$AllBird/temp$survEff)>=c(mean(c(temp$AllBird/temp$survEff))+sd(c(temp$AllBird/temp$survEff))),3,2)))

    temp$ClassShip <- ifelse(temp$AllShip<mean(temp$AllShip),1,
                               ifelse(temp$AllShip>=c(mean(temp$AllShip)+sd(temp$AllShip)),3,2))
    
    # Evaluate risk levels 
    temp$riskcat <- ifelse(temp$ClassBird == 1 & temp$ClassShip == 1, "low",
                    ifelse(temp$ClassBird == 1 & temp$ClassShip == 2 | temp$ClassBird == 1 & temp$ClassShip == 3, "mediumSHIP",
                    ifelse(temp$ClassBird == 2 & temp$ClassShip == 1 | temp$ClassBird == 3 & temp$ClassShip == 1, "mediumBIRD",
                    ifelse(temp$ClassBird == 3 & temp$ClassShip == 2 | temp$ClassBird == 2 & temp$ClassShip == 3, "high",
                    ifelse(temp$ClassBird == 2 & temp$ClassShip == 2, "high",
                    ifelse(temp$ClassBird == 3 & temp$ClassShip == 3, "veryhigh", NA))))))
    
    # Evaluate continuous risk 
    temp$QuantBird <- ecdf(temp$DensBird)(temp$DensBird)
    # BirdVal99 <- min(temp$DensBird[which(temp$QuantBird > 0.99)])
    # temp$DensBird_rescale1 <- ifelse(temp$QuantBird > 0.99, BirdVal99, temp$DensBird)
    
    temp$QuantShip <- ecdf(temp$AllShip)(temp$AllShip)
    # ShipVal99 <- min(temp$AllShip[which(temp$QuantShip > 0.99)])
    # temp$AllShip_rescale1 <-  ifelse(temp$QuantShip > 0.99, ShipVal99, temp$AllShip)
    
    temp$DensBird_rescale <- log(temp$DensBird+1)
    temp$AllShip_rescale <- log(temp$AllShip+1)

    # temp$DensBird_rescale3 <- rescale(log(temp$DensBird+1), to = c(0,10))
    # temp$AllShip_rescale3 <- rescale(log(temp$AllShip+1), to = c(0, 10))
    
    # temp$riskcont1 <- temp$DensBird_rescale1*temp$AllShip_rescale1
    temp$riskcont <- temp$DensBird_rescale*temp$AllShip_rescale
    # temp$riskcont3 <- temp$DensBird_rescale3*temp$AllShip_rescale3
    
    filtdf <- rbind(filtdf, temp)
  }
  
  filtdf$riskcat <- ordered(filtdf$riskcat,  levels=c("low","mediumSHIP","mediumBIRD","high","veryhigh"))
  
  # ggplot(subset(filtdf, !is.na(filtdf$riskcat)), aes(fill=riskcat, x=subset))+
  #   geom_bar(position="dodge")
  
  st_write(filtdf, paste0(savefolder, "FinalShapefiles/AllSeasonsAllTimeOfDay_", studyareaname, "_", taxaLabel, ".shp"))
  # filtdf$relDensBird <- filtdf$DensBird/max(filtdf$DensBird)
  # filtdf$relAllShip <- filtdf$AllShip/max(filtdf$AllShip)
  # 
  # filtdf$relriskcatprenorm <- filtdf$relDensBird*filtdf$relAllShip
  # filtdf <- filtdf %>% 
  #   mutate(relriskcatpostnorm = filtdf$DensBird*filtdf$AllShip) %>% 
  #   mutate(relriskcatpostnorm = relriskpostnorm/max(relriskpostnorm))
  
  widedf <- filtdf %>% dplyr::select(hexID, riskcat, subset, taxa) %>% spread(key=subset, value=riskcat)
  widecont <- filtdf %>% dplyr::select(hexID, riskcont, subset, taxa) %>% spread(key=subset, value=riskcont)
  
  allvsnightrisk <- data.frame(taxa = taxaLabel, 
                               studyarea = studyareaname, 
                               numhexes = length(widedf$Fall_Night),
                               nighthighrisk_Fall = sum(widedf$Fall_Night %in% c("high", "veryhigh")),
                               nighthighrisk_Summer = sum(widedf$Summer_Night %in% c("high", "veryhigh")),
                               # nightcontrisk_Summer_mean = mean(widecont$Summer_Night),
                               # nightcontrisk_Summer_median = median(widecont$Summer_Night),
                               # nightcontrisk_Summer_sd = sd(widecont$Summer_Night),
                               # nightcontrisk_Fall_mean = mean(widecont$Fall_Night),
                               # nightcontrisk_Fall_median = median(widecont$Fall_Night),
                               # nightcontrisk_Fall_sd = sd(widecont$Fall_Night),
                               
                               
                               allhighrisk_Fall = sum(widedf$Fall_All %in% c("high", "veryhigh")),
                               allhighrisk_Summer = sum(widedf$Summer_All %in% c("high", "veryhigh")),
                               # allcontrisk_Summer = mean(widecont$Summer_All),
                               # allcontrisk_Summer_median = median(widecont$Summer_All),
                               # allcontrisk_Summer_sd = sd(widecont$Summer_All),
                               # allcontrisk_Fall = mean(widecont$Fall_All),
                               # allcontrisk_Fall_median = median(widecont$Fall_All),
                               # allcontrisk_Fall_sd = sd(widecont$Fall_All),
                               
                               
                               nightmorethanall_Fall = sum(widedf$Fall_All < widedf$Fall_Night, na.rm=T), 
                               nightequaltoall_Fall= sum(widedf$Fall_All == widedf$Fall_Night, na.rm=T),
                               nightlessthanall_Fall= sum(widedf$Fall_All > widedf$Fall_Night, na.rm=T),
                               nightmorethanall_Summer = sum(widedf$Summer_All < widedf$Summer_Night, na.rm=T), 
                               nightequaltoall_Summer= sum(widedf$Summer_All == widedf$Summer_Night, na.rm=T),
                               nightlessthanall_Summer= sum(widedf$Summer_All > widedf$Summer_Night, na.rm=T))
  
  write.csv(allvsnightrisk, paste0("../Data_Processed/AllVsNightRisk/AllVsNightRisk_",studyareaname, "_", taxaLabel, ".csv" ))

  
  summfiltdf <- filter(st_drop_geometry(filtdf), filtdf$season == "Summer")
  fallfiltdf <- filter(st_drop_geometry(filtdf), filtdf$season == "Fall")
  
  summ <- data.frame(region = studyareaname, 
                     taxa = taxaLabel,
                     nhexes = length(unique(filtdf$hexID)),
                     surveyarea_summ = sum(filter(st_drop_geometry(filtdf), filtdf$season == "Summer" & filtdf$timeofday == "All") %>% dplyr::select(survEff)), 
                     surveyarea_fall = sum(filter(st_drop_geometry(filtdf), filtdf$season == "Fall" & filtdf$timeofday == "All") %>% dplyr::select(survEff)), 
                     totalstudyarea = sum(filter(st_drop_geometry(filtdf), filtdf$season == "Summer" & filtdf$timeofday == "Night") %>% dplyr::select(AreaKM)),
                     
                     densbird_summ = sum(summfiltdf$AllBird[summfiltdf$timeofday == "All"])/sum(summfiltdf$survEff[summfiltdf$timeofday == "All"]), 
                     alltraff_summ = sum(summfiltdf$AllShip[summfiltdf$timeofday == "All"]),
                     nighttraff_summ =  sum(summfiltdf$AllShip[summfiltdf$timeofday == "Night"]),
                     
                     densbird_mean_fall = sum(fallfiltdf$AllBird[fallfiltdf$timeofday == "All"])/sum(fallfiltdf$survEff[fallfiltdf$timeofday == "All"]),
                     alltraff_fall = sum(fallfiltdf$AllShip[fallfiltdf$timeofday == "All"]), 
                     nighttraff_fall = sum(fallfiltdf$AllShip[fallfiltdf$timeofday == "Night"])
                     )
  
  st_write(summ, paste0(savefolder, "SummaryStatistics/SummaryStatistics_", studyareaname, "_", taxaLabel, ".csv"))
  
  
  ### Plot prep 
  plottheme <- list(
    scale_x_longitude(ticks = 5, expand = c(0, 0)),
    scale_y_latitude(ticks = 5, expand = c(0, 0)),
    theme_bw(),
    theme(text = element_text(size = 18),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          plot.caption = element_text(size = 8, hjust=0),
          plot.margin = margin(t=0.5, r=0.5, b=0.5, l=0.5, unit="cm"),
          # panel.background = element_rect(fill = "#73b2ff"),
          panel.border =  element_rect(colour = "black"),
          axis.text = element_text(colour = "darkgray", size=8)))
  
  combos <- data.frame(id = 1:4, 
                       season = c("Summer", "Fall", "Summer", "Fall"), 
                       tod = c("All", "All", "Night", "Night"))
  tod <- c("All", "Night")
  seasons <- c("Summer","Fall")

  basemapnew <- basemap %>% st_crop(st_buffer(st_as_sfc(st_bbox(filtdf)), 10000))
  
  box <- st_as_sf(name="boundary", st_buffer(st_as_sfc(st_bbox(filtdf)), 10000))
  
  
  # Get color scale that works well for skewed data 
  cols <- c(colorRampPalette(c("#e7f0fa", "#c9e2f6", "#95cbee", "#0099dc", "#4ab04a", "#ffd73e"))(25),
            colorRampPalette(c("#eec73a", "#e29421", "#e29421", "#f05336","#ce472e"), bias=2)(25))
  
  ######################################################### 
  ##################### Risk Plot: Continuous #################
  ######################################################### 
  
  # for(i in 1:2){
  #   dfnew <- filtdf %>% filter(timeofday == tod[i])
  #   riskplotcontname <- paste0(figfolder, "Risk_Continuous_", tod[i],"_",taxaLabel, "_", studyareaname, ".png")
  #   
  #   for(j in 1:length(unique(combos$season))){
  #     dfsub <- filtdf %>% filter(season == unique(combos$season)[j])
  #     plt <- ggplot() +
  #       geom_sf(data=box, fill=NA, color=NA,lwd=0) +   
  #       geom_sf(data=basemapnew, fill="lightgray", lwd=0) +
  #       geom_sf(data=filter(dfsub, riskcont == 0), color="darkgray", fill=NA) +
  #       geom_sf(data=filter(dfsub, riskcont != 0),aes(fill = riskcont), color="darkgray") +
  #       scale_fill_gradientn(colours=cols, labels=scales::label_number(),name="Risk Index", na.value="white") +
  #       # scale_fill_steps(trans="log",low = "yellow", high = "red",nice.breaks=TRUE, labels=scales::label_number(), 
  #       #                  name="Density \n(Ind'ls/km\u00b2)", guide = guide_coloursteps(show.limits = TRUE)) + 
  #       # labs(caption = paste0("*Empty hexes were surveyed, but no ", taxaLabel, " were sighted during study period.")) + 
  #       guides(fill = guide_colourbar(barwidth = 25, 
  #                                     barheight = 1, 
  #                                     title.hjust = 0.5,
  #                                     ticks.colour="black", 
  #                                     frame.colour = "black", 
  #                                     ticks.linewidth = 0.75)) +
  #       plottheme +
  #       ggtitle(unique(combos$season)[j])
  #     assign(paste0("b", j), plt)
  #   }
  #   comboplot <- ggarrange(b1, b2, ncol=2, nrow=1, common.legend=TRUE, legend = "bottom")
  #   
  #   comboplot <- annotate_figure(comboplot, top = text_grob(dfsub$taxa[1], face = "bold", size = 30)) + 
  #     theme(panel.background = element_rect(fill = "white"))
  #   
  #   
  #   # Save figure
  #   if(!file.exists(riskplotcontname)){
  #     ifelse(studyareaname == "Eastern Aleutians",
  #            ggsave(filename=riskplotcontname,
  #                   plot= comboplot,
  #                   width=12, height=4, units="in"),
  #            ggsave(filename=riskplotcontname,
  #                   plot= comboplot,
  #                   width=10, height=6, units="in"))
  #   }
  # }
  
  ################################################# 
  ##################### Risk plot: Categorical #################
  ################################################# 

  for(i in 1:2){
    dfnew <- filtdf %>% filter(timeofday == tod[i])
    riskplotallname <- paste0(figfolder, "Risk_Categorical_", tod[i],"_",taxaLabel, "_", studyareaname, ".png")
    for(j in 1:2){
      dfsub <- dfnew %>% filter(season == seasons[j])
      dfsub$riskcat <- factor(dfsub$riskcat, 
                           levels = c("low", "mediumBIRD", "mediumSHIP", "high", "veryhigh"))
      plt <- ggplot() +
        geom_sf(data=box, fill=NA, color=NA,lwd=0) +      
        geom_sf(data=basemapnew, fill="lightgray",lwd=0) +
        geom_sf(data=dfsub,aes(fill = riskcat), color="darkgray") +
        scale_fill_manual(values = c("low" = "#73b2ff",
                                     "mediumBIRD" = "#90ee90",
                                     "mediumSHIP" = "#14ff96",
                                     "high" = "#f7f570",
                                     "veryhigh" = "#e31a1c"),
                          na.value = "white",
                          name="Risk", 
                          drop=F,
                          labels = c("Low", "Moderate (bird)", "Moderate (ship)", "High", "Very High")) +
        plottheme +
        ggtitle(paste0(dfsub$season[1]))
      assign(paste0("p", j), plt)
    }
    comboplot <- ggarrange(p1, p2, ncol=2, nrow=1, common.legend=TRUE, legend = "bottom")
    comboplot <- annotate_figure(comboplot, top = text_grob(dfsub$taxa[1], face = "bold", size = 30)) + 
      theme(panel.background = element_rect(fill = "white"))
    
      # Save figure
    if(!file.exists(riskplotallname)){
      ifelse(studyareaname == "Eastern Aleutians",
        ggsave(filename=riskplotallname,
               plot= comboplot,
               width=12, height=4, units="in"),
        ggsave(filename=riskplotallname,
               plot= comboplot,
               width=10, height=6, units="in"))
    }
  }

  #################################################### 
  ##################### Traffic plot #################
  #################################################### 
  
  if(studyareaname %in% c("All Alaska", "Gulf of Alaska", "Eastern Aleutians")){
    br <- c(0,10,100, 1000, 10000, 100000, 100000)
    lims <- c(0,150000)
  }
  if(studyareaname == "Northern Bering & Chukchi Seas"){
    br <- c(0,10,100, 1000, 10000)
    lims <- c(0,50000)
  }
  
  for(i in 1:2){
    dfnew <- filtdf %>% filter(timeofday == unique(combos$tod)[i])
    traffplotname <- paste0(figfolder,"TrafficCombos_", unique(combos$tod)[i],"_", studyareaname,".png")
    
    for(j in 1:2){
      dfsub <- dfnew %>% filter(season == unique(combos$season)[j])
      plt <- ggplot() +
        geom_sf(data=box, fill=NA, color=NA,lwd=0) +   
        geom_sf(data=basemapnew, fill="lightgray", lwd=0) +
        geom_sf(data=dfsub,aes(fill = AllShip), color="darkgray") +
        # scale_fill_steps(trans="log",n.breaks=4, low = "yellow", high = "red",nice.breaks=TRUE, labels=scales::comma,
        #                  name="Total Hours\nof Vessel Traffic", guide = guide_coloursteps(show.limits = TRUE)) +
        scale_fill_gradientn(colours=cols, trans="pseudo_log", breaks=br, labels=scales::label_comma(),name="Vessel Activity\n(Hours)", na.value="white", limits=lims) +
        # labs(caption = paste0("*One operating day is equal to one vessel present in a hex on a given day.")) + 
        guides(fill = guide_colourbar(barwidth = 25, 
                                      barheight = 1, 
                                      title.hjust = 0.5,
                                      ticks.colour="black", 
                                      frame.colour = "black", 
                                      ticks.linewidth = 0.75)) +
        plottheme +
        ggtitle(dfsub$season[1])
      assign(paste0("t", j), plt)
  }
    traffcombo <- ggarrange(t1, t2, ncol=2, nrow=1, common.legend=TRUE, legend = "bottom")
    
    traffcombo <- traffcombo + theme(panel.background = element_rect(fill = "white"))
    
  
  # Save figure 
  if(!file.exists(traffplotname)){
    ifelse(studyareaname == "Eastern Aleutians", 
           ggsave(filename=traffplotname, 
                  plot= traffcombo, 
                  width=12, height=4, units="in"), 
           ggsave(filename=traffplotname, 
                  plot= traffcombo, 
                  width=10, height=6, units="in"))
    }
  }
  
  ######################################################### 
  ##################### Bird density plot #################
  ######################################################### 
  
  birdplotname <- paste0(figfolder,"BirdDensityCombos_",studyareaname,"_",taxaLabel, ".png")
  
  for(i in 1:length(unique(combos$season))){
    dfsub <- filtdf %>% filter(season == unique(combos$season)[i])
    plt <- ggplot() +
      geom_sf(data=box, fill=NA, color=NA,lwd=0) +   
      geom_sf(data=basemapnew, fill="lightgray", lwd=0) +
      geom_sf(data=filter(dfsub, AllBird == 0), color="darkgray", fill=NA) +
      geom_sf(data=filter(dfsub, AllBird != 0),aes(fill = DensBird), color="darkgray") +
      scale_fill_gradientn(colours=cols, trans="log10", labels=scales::label_number(),name="Density \n(Ind'ls/km\u00b2)", na.value="white") +
      # scale_fill_steps(trans="log",low = "yellow", high = "red",nice.breaks=TRUE, labels=scales::label_number(), 
      #                  name="Density \n(Ind'ls/km\u00b2)", guide = guide_coloursteps(show.limits = TRUE)) + 
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      # labs(caption = paste0("*Empty hexes were surveyed, but no ", taxaLabel, " were sighted during study period.")) + 
      guides(fill = guide_colourbar(barwidth = 25, 
                                    barheight = 1, 
                                    title.hjust = 0.5,
                                    ticks.colour="black", 
                                    frame.colour = "black", 
                                    ticks.linewidth = 0.75)) +
      plottheme +
      ggtitle(unique(combos$season)[i])
    assign(paste0("b", i), plt)
  }
  birdcombo <- ggarrange(b1, b2, ncol=2, nrow=1, common.legend=TRUE, legend = "bottom")
  
  birdcombo <- annotate_figure(birdcombo, top = text_grob(dfsub$taxa[1], face = "bold", size = 30)) + 
    theme(panel.background = element_rect(fill = "white"))
  
  # Save figure
  if(!file.exists(birdplotname)){
    ifelse(studyareaname == "Eastern Aleutians", 
           ggsave(filename=birdplotname, 
                  plot= birdcombo, 
                  width=12, height=4, units="in"), 
           ggsave(filename=birdplotname, 
                  plot= birdcombo,
                  width=10, height=6, units="in"))
  }
  
}
