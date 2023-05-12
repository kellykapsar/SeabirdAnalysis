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
                              figfolder,
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
    
    # Calculate SD categories for each hex's effort-weighted seabird observations 
    birdFinal <- birdStacked %>%
      mutate(QuantBird = ecdf(AllBird/survEff)(AllBird/survEff))
    
    # Calculate effort-weighted densities of seabirds (# birds per square km of surveyed area)
    birdFinal$DensBird <- c(birdFinal$AllBird/birdFinal$survEff)
    
    # Calculate effort-weighted classes for the number of observations
    birdFinal$ClassBird <- ifelse(birdFinal$AllBird == 0, 0, 
                            ifelse(c(birdFinal$AllBird/birdFinal$survEff)<mean(c(birdFinal$AllBird/birdFinal$survEff)),1,
                             ifelse(c(birdFinal$AllBird/birdFinal$survEff)>=c(mean(c(birdFinal$AllBird/birdFinal$survEff))+sd(c(birdFinal$AllBird/birdFinal$survEff))),3,2)))
    
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

    # Evaluate risk levels 
    df$risk <- ifelse(df$ClassBird == 1 & df$ClassShip == 1, "low",
                ifelse(df$ClassBird == 1 & df$ClassShip == 2 | df$ClassBird == 2 & df$ClassShip == 1, "medium",
                ifelse(df$ClassBird == 1 & df$ClassShip == 3 | df$ClassBird == 3 & df$ClassShip == 1, "medium",
                ifelse(df$ClassBird == 3 & df$ClassShip == 2 | df$ClassBird == 2 & df$ClassShip == 3, "high",
                ifelse(df$ClassBird == 2 & df$ClassShip == 2, "high",
                ifelse(df$ClassBird == 3 & df$ClassShip == 3, "veryhigh", NA))))))


    df$risk <- factor(df$risk,  c("low","medium","high","veryhigh"))
    
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
  
  bothseasons <- df %>% st_drop_geometry() %>% group_by(hexID) %>% summarize(n=n()) %>% filter(n == 4)
  
  filtdf <- df %>% filter(hexID %in% bothseasons$hexID) %>% mutate(subset = paste0(season, "_", timeofday))

  filtdf$ClassBird <- ifelse(filtdf$AllBird == 0, 0, 
                      ifelse(c(filtdf$AllBird/filtdf$survEff)<mean(c(filtdf$AllBird/filtdf$survEff)),1,
                      ifelse(c(filtdf$AllBird/filtdf$survEff)>=c(mean(c(filtdf$AllBird/filtdf$survEff))+sd(c(filtdf$AllBird/filtdf$survEff))),3,2)))
  
  filtdf$ClassShip <- ifelse(filtdf$AllShip<mean(filtdf$AllShip),1,
                             ifelse(filtdf$AllShip>=c(mean(filtdf$AllShip)+sd(filtdf$AllShip)),3,2))
  
  # Evaluate risk levels 
  filtdf$risk <- ifelse(filtdf$ClassBird == 1 & filtdf$ClassShip == 1, "low",
                    ifelse(filtdf$ClassBird == 1 & filtdf$ClassShip == 2 | filtdf$ClassBird == 2 & filtdf$ClassShip == 1, "medium",
                    ifelse(filtdf$ClassBird == 1 & filtdf$ClassShip == 3 | filtdf$ClassBird == 3 & filtdf$ClassShip == 1, "medium",
                    ifelse(filtdf$ClassBird == 3 & filtdf$ClassShip == 2 | filtdf$ClassBird == 2 & filtdf$ClassShip == 3, "high",
                    ifelse(filtdf$ClassBird == 2 & filtdf$ClassShip == 2, "high",
                    ifelse(filtdf$ClassBird == 3 & filtdf$ClassShip == 3, "veryhigh", NA))))))
  
  filtdf$risk <- factor(filtdf$risk,  c("low","medium","high","veryhigh"))
  # filtdf$relDensBird <- filtdf$DensBird/max(filtdf$DensBird)
  # filtdf$relAllShip <- filtdf$AllShip/max(filtdf$AllShip)
  # 
  # filtdf$relriskprenorm <- filtdf$relDensBird*filtdf$relAllShip
  # filtdf <- filtdf %>% 
  #   mutate(relriskpostnorm = filtdf$DensBird*filtdf$AllShip) %>% 
  #   mutate(relriskpostnorm = relriskpostnorm/max(relriskpostnorm))
  
  widedf <- filtdf %>% dplyr::select(hexID, risk, subset, taxa) %>% spread(key=subset, value=risk)

  summ <- data.frame(region = studyareaname, 
                     taxa = taxaLabel,
                     surveyarea = sum(df$survEff), 
                     
                     totbird_summ = sum(filter(st_drop_geometry(df), df$season == "Summer") %>% dplyr::select(AllBird)), 
                     traff_summ = sum(filter(st_drop_geometry(df), df$season == "Summer") %>% dplyr::select(AllShip)), 
                     daytraff_summ = sum(filter(st_drop_geometry(df), df$timeofday == "Day" & df$season == "Summer") %>% dplyr::select(AllShip)), 
                     nighttraff_summ = sum(filter(st_drop_geometry(df), df$timeofday == "Night" & df$season == "Summer") %>% dplyr::select(AllShip)), 
                     
                     totbird_fall = sum(filter(st_drop_geometry(df), df$season == "Fall") %>% dplyr::select(AllBird)), 
                     traff_fall = sum(filter(st_drop_geometry(df), df$season == "Fall") %>% dplyr::select(AllShip)), 
                     daytraff_fall = sum(filter(st_drop_geometry(df), df$timeofday == "Day" & df$season == "Fall") %>% dplyr::select(AllShip)), 
                     nighttraff_fall = sum(filter(st_drop_geometry(df), df$timeofday == "Night" & df$season == "Fall") %>% dplyr::select(AllShip)) 
                     )
  
  write.csv(summ, paste0(savefolder, "SummaryStatistics/SummaryStatistics_", studyareaname, "_", taxaLabel, ".csv"))
  
  combos <- data.frame(id = 1:4, 
                       season = c("Summer", "Summer", "Fall", "Fall"), 
                       tod = c("Day", "Night", "Day", "Night"))

  basemapnew <- basemap %>% st_crop(st_buffer(st_as_sfc(st_bbox(filtdf)), 10000))
  
  box <- st_as_sf(name="boundary", st_buffer(st_as_sfc(st_bbox(filtdf)), 10000))
  
  ################################################# 
  ##################### Risk plot #################
  ################################################# 
  riskplotname <- paste0(figfolder, "RiskCombos_", taxaLabel, "_", studyareaname, ".png")
  
  for(i in 1:length(combos$id)){
    dfsub <- filtdf %>% filter(season == combos$season[i] & timeofday == combos$tod[i])
    plt <- ggplot() +
      geom_sf(data=box, fill=NA, color=NA,lwd=0) +      
      geom_sf(data=basemapnew, fill="lightgray",lwd=0) +
      geom_sf(data=dfsub,aes(fill = risk), color="darkgray") +
      scale_fill_manual(values = c("low" = "#73b2ff",
                                   "medium" = "#55fe01",
                                   "high" = "#ffff01",
                                   "veryhigh" = "#e31a1c"),
                        na.value = "white",
                        name="Risk",
                        labels = c("Low", "Medium", "High", "Very High")) +
      xlab("") +
      ylab("") +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      theme_bw() +
      theme(text = element_text(size = 18),
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            plot.caption = element_text(size = 8, hjust=0),
            axis.ticks = element_blank(),
            axis.text=element_blank(),
            panel.border =  element_rect(colour = "black"),
            panel.grid.major = element_line(colour = "transparent"), 
            panel.background = element_rect(fill = "white")) +
      ggtitle(paste0(dfsub$season[1], ": ", dfsub$timeofday[1]))
    assign(paste0("p", i), plt)
  }
  comboplot <- ggarrange(p1, p3, p2, p4, ncol=2, nrow=2, common.legend=TRUE, legend = "bottom")
  comboplot <- annotate_figure(comboplot, top = text_grob(dfsub$taxa[1], face = "bold", size = 30)) + 
    theme(panel.background = element_rect(fill = "white"))
  
    # Save figure
  if(!file.exists(riskplotname)){
    ifelse(studyareaname == "Eastern Aleutians", 
      ggsave(filename=riskplotname, 
             plot= comboplot, 
             width=12, height=6, units="in"), 
      ggsave(filename=riskplotname, 
             plot= comboplot, 
             width=12, height=12, units="in"))
  }

  #################################################### 
  ##################### Traffic plot #################
  #################################################### 
  traffplotname <- paste0(figfolder,"TrafficCombos_", studyareaname,".png")
  
  # Get color scale that works well for skewed data 
  cols <- c(colorRampPalette(c("#e7f0fa", "#c9e2f6", "#95cbee", "#0099dc", "#4ab04a", "#ffd73e"))(25),
            colorRampPalette(c("#eec73a", "#e29421", "#e29421", "#f05336","#ce472e"), bias=2)(25))
  
  if(studyareaname %in% c("All Alaska", "Gulf of Alaska", "Eastern Aleutians")){
    br <- c(0,10,100, 1000, 10000, 100000, 100000)
    lims <- c(0,150000)
  }
  if(studyareaname == "Northern Bering & Chukchi Seas"){
    br <- c(0,10,100, 1000, 10000)
    lims <- c(0,50000)
  }
  
  for(i in 1:length(combos$id)){
    dfsub <- filtdf %>% filter(season == combos$season[i] & timeofday == combos$tod[i])
    plt <- ggplot() +
      geom_sf(data=box, fill=NA, color=NA,lwd=0) +   
      geom_sf(data=basemapnew, fill="lightgray", lwd=0) +
      geom_sf(data=dfsub,aes(fill = AllShip), color="darkgray") +
      # scale_fill_steps(trans="log",n.breaks=4, low = "yellow", high = "red",nice.breaks=TRUE, labels=scales::comma,
      #                  name="Total Hours\nof Vessel Traffic", guide = guide_coloursteps(show.limits = TRUE)) +
      scale_fill_gradientn(colours=cols, trans="pseudo_log", breaks=br, labels=scales::label_comma(),name="Vessel Activity\n(Hours)", na.value="white", limits=lims) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      # labs(caption = paste0("*One operating day is equal to one vessel present in a hex on a given day.")) + 
      guides(fill = guide_colourbar(barwidth = 25, 
                                    barheight = 1, 
                                    title.hjust = 0.5,
                                    ticks.colour="black", 
                                    frame.colour = "black", 
                                    ticks.linewidth = 0.75)) +
      theme_bw() +
      theme(text = element_text(size = 18),
            axis.ticks=element_blank(), 
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5), 
            plot.caption = element_text(size = 8, hjust=0),
            axis.text=element_blank(),
            panel.border =  element_rect(colour = "black"),
            panel.grid.major = element_line(colour = "transparent")) +
    ggtitle(paste0(dfsub$season[1], ": ", dfsub$timeofday[1]))
    assign(paste0("t", i), plt)
  }
  traffcombo <- ggarrange(t1, t3, t2, t4, ncol=2, nrow=2, common.legend=TRUE, legend = "bottom")
  traffcombo <- traffcombo + theme(panel.background = element_rect(fill = "white"))
  
  # Save figure 
  if(!file.exists(traffplotname)){
    ifelse(studyareaname == "Eastern Aleutians", 
           ggsave(filename=traffplotname, 
                  plot= traffcombo, 
                  width=12, height=6, units="in"), 
           ggsave(filename=traffplotname, 
                  plot= traffcombo, 
                  width=12, height=12, units="in"))
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
      theme_bw() +
      theme(text = element_text(size = 18),
            axis.ticks=element_blank(), 
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5), 
            plot.caption = element_text(size = 8, hjust=0),
            plot.margin = margin(t=0.5, r=0.5, b=0.5, l=0.5, unit="cm"), 
            axis.text=element_blank(),
            # panel.background = element_rect(fill = "#73b2ff"),
            panel.border =  element_rect(colour = "black"),
            panel.grid.major = element_line(colour = "transparent")) +
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