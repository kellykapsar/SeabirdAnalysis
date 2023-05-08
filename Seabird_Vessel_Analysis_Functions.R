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

traffDensity <- function(filedir, mnths, metric, night, trafffilename){
  
  # Isolate month values of interest from file names 
  hexes <- filedir[as.numeric(substr(filedir,10,11)) %in% mnths]
  
  # Read in data 
  temp <- lapply(hexes, function(x){st_read(paste0(hexdir,x)) %>% st_drop_geometry()})
  hexAll <- do.call(rbind, temp)
  
  # Select metrics based on inputs 
  lab <- ifelse(metric=="OperatingDays", "OpD", 
         ifelse(metric=="Ships","Shp", 
         ifelse(metric=="Hours", "Hrs", NA)))
  
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
  
  # Save results 
  write.csv(hexRes,trafffilename)
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
                              studyarea, 
                              studyareaname,
                              startyr,
                              savefolder,
                              figfolder,
                              filedir,
                              metric, 
                              night){
  
  # Make sure appropriate vessel data exist and if not, generate it
  trafffilename <- paste0(savefolder,"TraffInHexes_",metric,"_",mnthsnam,"_NightOnly", night,".csv")
  
  if(!file.exists(trafffilename)){
    traffDensity(filedir = hexList, mnths = mnths, metric = metric, night = night, trafffilename)
  }
  
  traffdf <- read.csv(trafffilename)
  
  # Make sure survey effort has been calculated and saved for appropriate months and if not, generate it
  survfilename <- paste0(savefolder,"ObsInHexes_",monthsname,".shp")
  
  if(!file.exists(survfilename)){
    birdDensity(loc=loc, datobs=datobs, hex=hexMask, startyr=startyr, mnths=mnths, survfilename=survfilename)
  }
  birddf <- st_read(survfilename)
  
  # Make sure this function has not already been run with these exact parameter specifications 
  finaldfname <- paste0(savefolder,"FinalDF_",studyareaname,"_",taxaLabel,"_",monthsname,"_NightOnly",night[1],".shp")
  
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
    birdFinal$ClassBird <- ifelse(c(birdFinal$AllBird/birdFinal$survEff)<mean(c(birdFinal$AllBird/birdFinal$survEff)),1,
                             ifelse(c(birdFinal$AllBird/birdFinal$survEff)>=c(mean(c(birdFinal$AllBird/birdFinal$survEff))+sd(c(birdFinal$AllBird/birdFinal$survEff))),3,2))
    
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
    
    # Save data 
    st_write(df, paste0(savefolder,"FinalDF_",studyareaname,"_",taxaLabel,"_",monthsname,"_NightOnly",night[1],".shp"))
  }
}

#### Results Plots #### 
plotResults <- function(basemap,
                        studyarea,
                        studyareaname, 
                        figfolder,
                        savefolder,
                        monthsname,
                        taxaLabel,
                        night,
                        metricName){

  # Load in data file
  dfname <- paste0(savefolder,"FinalDF_",studyareaname,"_",taxaLabel,"_",monthsname,"_NightOnly",night[1],".shp")
  
  if(!file.exists(dfname)){
    stop("Combined dataframe does not exist. Try running birdHexesByEffort function first.")
  }
  
  df <- st_read(dfname)  %>% st_crop(studyarea)
  
  # Prep basemap  
  basemapnew <- basemap %>% st_crop(studyarea, 100000)
  
  finalNoBird <- df %>% dplyr::filter(DensBird == 0)
  finalBird <- df %>% dplyr::filter(DensBird > 0)
  
  # Get color scale that works well for skewed data 
  cols <- c(colorRampPalette(c("#e7f0fa", "#c9e2f6", "#95cbee", "#0099dc", "#4ab04a", "#ffd73e"))(25),
            colorRampPalette(c("#eec73a", "#e29421", "#e29421", "#f05336","#ce472e"), bias=2)(25))
  
  
  if(length(finalBird$hexID) == 0){stop(paste0("No ", taxaLabel, " were found in this study area."))}
    
  #### Traffic plot ####
  traffplotname <- paste0(figfolder,"TrafficPlots/TrafficDensity_", studyareaname,"_",monthsname,"_NightOnly",night,".png")
  
  p0 <- ggplot() +
    geom_sf(data=basemapnew, fill="lightgray", lwd=0) +
    geom_sf(data=df,aes(fill = AllShip), color="lightgray") +
    # scale_fill_steps(trans="log",low = "yellow", high = "red",nice.breaks=TRUE, labels=scales::comma,
    #                  name="Total Hours\nof Vessel Traffic", guide = guide_coloursteps(show.limits = TRUE)) +
    scale_fill_gradientn(colours=cols, trans="log10", labels=scales::label_number(),name="Vessel Activity\n(Hours)") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    # labs(caption = paste0("*One operating day is equal to one vessel present in a hex on a given day.")) + 
    theme_bw() +
    theme(text = element_text(size = 18),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5), 
          plot.caption = element_text(size = 8, hjust=0),
          axis.text=element_blank(),
          panel.border =  element_rect(colour = "black"),
          panel.grid.major = element_line(colour = "transparent")) 
    
    ifelse(night == TRUE, 
      p0alone <- p0 + ggtitle(label = paste0("Nighttime Vessel Traffic Intensity\n", studyareaname), 
                   subtitle = paste0("(",monthsname," 2015-2022)")), 
      p0alone <- p0 + ggtitle(label = paste0("Vessel Traffic Intensity\n", studyareaname), 
                         subtitle = paste0("(",monthsname," 2015-2022)")) 
    )
    
    # Save figure 
    if(!file.exists(traffplotname)){
      ggsave(filename = traffplotname, plot = p0alone, width=10, height=8, units="in")
    }

  #### Risk plot ####
  riskplotname <- paste0(figfolder,"RiskPlots/RiskMap_", studyareaname,"_",taxaLabel,"_",monthsname,"_NightOnly",night,".png")
  
  df$risk <- factor(df$risk, levels=c("low", "medium", "high", "veryhigh"))
  
  p1 <- ggplot() +
    geom_sf(data=basemapnew, fill="#fefeff", lwd=0) +
    geom_sf(data=df,aes(fill = risk)) +
    scale_fill_manual(values = c("low" = "#73b2ff",
                                 "medium" = "#55fe01",
                                 "high" = "#ffff01",
                                 "veryhigh" = "#e31a1c", 
                                 na.value = "#bcc7dd"), 
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
          axis.text=element_blank(),
          panel.background = element_rect(fill = "#bcc7dd"),
          panel.border =  element_rect(colour = "black"),
          panel.grid.major = element_line(colour = "transparent"))
  
  ifelse(night == TRUE, 
         p1alone <- p1 + ggtitle(label = paste0(taxaLabel, " and Nighttime Vessel Traffic\n", studyareaname), 
                            subtitle = paste0("(",monthsname," 2015-2022)")), 
         p1alone <- p1 + ggtitle(label = paste0(taxaLabel, " and Vessel Traffic\n", studyareaname), 
                            subtitle = paste0("(",monthsname," 2015-2022)")))
  
  # Save figure 
  if(!file.exists(riskplotname)){
    ggsave(filename = riskplotname,
           plot = p1alone, width=12, height=8, units="in")
  }
  
  #### Bivariate plot ####
  bivarplotname <- paste0(figfolder,"BivariatePlots/BiVariateMap_",studyareaname,"_",taxaLabel,"_",monthsname,"_NightOnly",night,".png")
  bivarlegendname <- paste0(figfolder,"BivariatePlots/BiVariateMap_",studyareaname,"_",taxaLabel,"_",monthsname,"_NightOnly",night,"_Legend.png")
  
  # Possibly create custom color palette
  # 
  # custom_pal3 <- c(
  #   "1-1" = "#d3d3d3", # low x, low y
  #   "2-1" = "#A1D08A",
  #   "3-1" = "#66CC33", # high x, low y
  #   "1-2" = "#E7E773",
  #   "2-2" = "#FF9966", # medium x, medium y
  #   "3-2" = "#ACB54A",
  #   "1-3" = "#FFFF00", # low x, high y
  #   "2-3" = "#FFC738",
  #   "3-3" = "#FF1C13" # high x, high y
  # )
    # And custom breaks using https://stackoverflow.com/questions/46750635/cut-and-quantile-in-r-in-not-including-zero
  # But also we've already developed our own breaks for the risk, so this just makes it confusing... 

  df <- bi_class(df, x = DensBird , y = AllShip, style = "jenks", dim = 3)
  bs <- bi_class_breaks(df, x = DensBird, y = AllShip, style = "jenks", dim = 3, split=TRUE)
  
  # round(mean(c(min(df$DensBird[df$ClassBird == 2]), max(df$DensBird[df$ClassBird == 1]))),1)
  # round(mean(c(min(df$DensBird[df$ClassBird == 3]), max(df$DensBird[df$ClassBird == 2]))),1)
  # 
  # plyr::round_any(mean(c(min(df$AllShip[df$ClassShip == 2]), max(df$AllShip[df$ClassShip == 1]))),100, f=ceiling)
  # plyr::round_any(mean(c(min(df$AllShip[df$ClassShip == 3]), max(df$AllShip[df$ClassShip == 2]))),100, f=ceiling)

  p2 <- ggplot() +
    geom_sf(data=basemapnew, fill="#8ba761", lwd=0) +
    # geom_sf(data=finalNoBird, fill="#73b2ff", color="lightgray") + 
    geom_sf(data=df,aes(fill = bi_class), color="#73b2ff", show.legend = FALSE) +
    bi_scale_fill(pal = "PurpleOr", dim = 3) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(caption = paste0("*Empty hexes were surveyed, but no ", taxaLabel, " were sighted during study period.")) + 
    bi_theme() +
    theme_bw() +
    theme(text = element_text(size = 18),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5), 
          plot.caption = element_text(size = 8, hjust=0),
          axis.text=element_blank(),
          panel.background = element_rect(fill = "#73b2ff"),
          panel.border =  element_rect(colour = "black"),
          panel.grid.major = element_line(colour = "transparent")) 
  
  ifelse(night == TRUE, 
         p2alone <- p2 + ggtitle(label = paste0(taxaLabel, " and Nighttime Vessel Traffic\n", studyareaname), 
                            subtitle = paste0("(",monthsname," 2015-2022)")), 
         p2alone <- p2 + ggtitle(label = paste0(taxaLabel, " and Vessel Traffic\n", studyareaname), 
                            subtitle = paste0("(",monthsname," 2015-2022)")))

  legend <- bi_legend(pal = "PurpleOr",
                      dim = 3,
                      xlab = "Bird Density",
                      ylab = "Vessel Traffic",
                      size = 10, 
                      breaks = bs, 
                      arrows=TRUE)
  leggrob <- ggplotGrob(legend)
  
  # Can't figure out how to draw legend onto plot for different study areas so going to disable this and 
  # save legend separately for now. 
  
  # p2Final <- ggdraw() +
  #   draw_plot(p2, 0, 0, 1, 1) 
  # 
  # p2Final <- p2Final + draw_plot(legend + theme(plot.background = element_rect(fill = "#8ba761", color = NA), 
  #                                                 text = element_text(color = "black")), x=0.65, y=.38, width=0.25, height=0.265) # +


  # Save figure 
  if(!file.exists(bivarplotname)){
    ggsave(filename = bivarplotname,
           plot = p2alone, width=10, height=8, units="in")
  
    ggsave(filename = bivarlegendname,
           plot = legend, width=2, height=2, units="in")
  }
  
  #### Bird density plot ####
  birdplotname <- paste0(figfolder,"BirdDensityPlots/DensityMap_",studyareaname,"_",taxaLabel,"_",monthsname,"_NightOnly",night,".png")
  
  p3 <- ggplot() +
    geom_sf(data=basemapnew, fill="lightgray", lwd=0) +
    # geom_sf(data=finalNoBird, fill="#73b2ff", color="lightgray") + 
    geom_sf(data=finalBird,aes(fill = DensBird), color="lightgray") +
    scale_fill_gradientn(colours=cols, trans="log10", labels=scales::label_number(),name="Density \n(Ind'ls/km\u00b2)") +
    # scale_fill_steps(trans="log",low = "yellow", high = "red",nice.breaks=TRUE, labels=scales::label_number(), 
    #                  name="Density \n(Ind'ls/km\u00b2)", guide = guide_coloursteps(show.limits = TRUE)) + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(caption = paste0("*Empty hexes were surveyed, but no ", taxaLabel, " were sighted during study period.")) + 
    theme_bw() +
    theme(text = element_text(size = 18),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5), 
          plot.caption = element_text(size = 8, hjust=0),
          axis.text=element_blank(),
          # panel.background = element_rect(fill = "#73b2ff"),
          panel.border =  element_rect(colour = "black"),
          panel.grid.major = element_line(colour = "transparent")) 
  
  p3alone <- p3 + ggtitle(label = paste0("Effort-weighted density of ", taxaLabel), 
                              subtitle = paste0("(",monthsname,")"))
  
  # Save figure
  if(!file.exists(birdplotname)){
    ggsave(filename = birdplotname,
           plot = p3alone, width=10, height=8, units="in")
  }

  #### Combo plot ####
  
  comboname <- paste0(figfolder,"ComboPlot_",studyareaname,"_",taxaLabel,"_",monthsname,"_NightOnly",night,".png")
  
  # Set matrix for layout of plots 1-5
  lay <- rbind(c(1,1,1,2,3,3,3,3),
               c(4,4,4,4,5,5,5,5))
  
  # Turn plots into "graphical objects" (aka grobs)
  p2g <- as_grob(p2)
  p1g <- as_grob(p1)
  p3g <- as_grob(p3)
  p0g <- as_grob(p0)
  
  # Arrange grobs by layout 
  combo <- grid.arrange(grobs=gList(p2g, leggrob, p1g, p3g, p0g), layout_matrix=lay)
  
  # Add custom title 
  titletext <-  ifelse(night == TRUE, 
        paste0(taxaLabel, " and Nighttime Vessel Traffic\n", studyareaname, " (",monthsname," 2015-2022)"), 
        paste0(taxaLabel, " and Vessel Traffic\n", studyareaname, " (",monthsname," 2015-2022)"))
  
  combofinal <- annotate_figure(combo, top=text_grob(titletext, color = "black", face = "bold", size = 20))
  
  # Make panel background white
  combofinal <- combofinal + theme(panel.background = element_rect(fill = "white"))
  
  # Save output 
  ggsave(filename = comboname,
         plot = combofinal, width=18, height=10, units="in")
}

  
