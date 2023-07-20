
library(tidyverse)
library(sf)
library(scales)
library(ggpubr)
library(metR)


taxaNames <- c("Total Seabirds", "Auklets", "Northern Fulmars", "Seaducks", "Shearwaters", "Storm Petrels")

ContinuousRiskMaps <- function(taxaName){
  
  filelist <- list.files("../Data_Processed/FinalShapefiles/", pattern="_All Alaska_", full.names=T)
  filelist <- filelist[grep(".shp", filelist)]
  file <- filelist[grep(taxaName, filelist)]
  df <- st_read(file)
  
  daydf <- df %>% 
    filter(timefdy == "All") %>% 
    mutate(DensBird_rescale = rescale(log(DensBrd+1), to = c(0,10)),
           AllShip_rescale = rescale(log(AllShip+1), to = c(0, 10)))
  nightdf <-  df %>% 
    filter(timefdy == "Night") %>%
    mutate(DensBird_rescale = rescale(log(DensBrd+1), to = c(0,10)),
           AllShip_rescale = rescale(log(AllShip+1), to = c(0, 10)))
  
  dfnew <- rbind(daydf, nightdf)
  
  dfnew$riskscaled <- dfnew$DensBird_rescale*dfnew$AllShip_rescale
  
  tod <- c("All", "Night")
  seasons <- c("Summer","Fall")
  
  basemapnew <- st_read("../Data_Raw/AK_CAN_RUS/AK_CAN_RUS.shp") %>% st_crop(st_buffer(st_as_sfc(st_bbox(dfnew)), 10000))
  
  box <- st_as_sf(name="boundary", st_buffer(st_as_sfc(st_bbox(dfnew)), 10000))
  
  
  # Get color scale that works well for skewed data 
  cols <- c(colorRampPalette(c("#e7f0fa", "#c9e2f6", "#95cbee", "#0099dc", "#4ab04a", "#ffd73e"))(25),
            colorRampPalette(c("#eec73a", "#e29421", "#e29421", "#f05336","#ce472e"), bias=2)(25))
  
  
  for(i in 1:2){
    dffilt <- dfnew %>% filter(timefdy == tod[i])
    riskplotcontname <- paste0("../Figures/Risk_Continuous_", tod[i],"_",taxaName, "_All Alaska.png")
  
    for(j in 1:length(seasons)){
      dfsub <- dffilt %>% filter(season == seasons[j])
      plt <- ggplot() +
        geom_sf(data=box, fill=NA, color=NA,lwd=0) +
        geom_sf(data=basemapnew, fill="lightgray", lwd=0) +
        geom_sf(data=filter(dfsub, riskscaled == 0), color="darkgray", fill=NA) +
        geom_sf(data=filter(dfsub, riskscaled != 0),aes(fill = riskscaled), color="darkgray") +
        scale_fill_gradientn(colours=cols, labels=scales::label_number(),name="Risk Index", na.value="white") +
        # scale_fill_steps(trans="log",low = "yellow", high = "red",nice.breaks=TRUE, labels=scales::label_number(),
        #                  name="Density \n(Ind'ls/km\u00b2)", guide = guide_coloursteps(show.limits = TRUE)) +
        # labs(caption = paste0("*Empty hexes were surveyed, but no ", taxaLabel, " were sighted during study period.")) +
        guides(fill = guide_colourbar(barwidth = 25,
                                      barheight = 1,
                                      title.hjust = 0.5,
                                      ticks.colour="black",
                                      frame.colour = "black",
                                      ticks.linewidth = 0.75)) +
        scale_x_longitude(ticks = 5, expand = c(0, 0)) +
        scale_y_latitude(ticks = 5, expand = c(0, 0)) +
        theme_bw() +
        theme(text = element_text(size = 18),
              plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5),
              plot.caption = element_text(size = 8, hjust=0),
              plot.margin = margin(t=0.5, r=0.5, b=0.5, l=0.5, unit="cm"),
              # panel.background = element_rect(fill = "#73b2ff"),
              panel.border =  element_rect(colour = "black"),
              axis.text = element_text(colour = "darkgray", size=8)) +
        ggtitle(seasons[j])
      assign(paste0("b", j), plt)
    }
    comboplot <- ggarrange(b1, b2, ncol=2, nrow=1, common.legend=TRUE, legend = "bottom")
  
    comboplot <- annotate_figure(comboplot, top = text_grob(dfsub$taxa[1], face = "bold", size = 30)) +
      theme(panel.background = element_rect(fill = "white"))
  
  
    # Save figure
    if(!file.exists(riskplotcontname)){
      ggsave(filename=riskplotcontname,
             plot= comboplot,
             width=12, height=6, units="in")
    }
  }
  
  goa <- st_read("../Data_Raw/KodiakPWS.shp") 
  uni <- st_read("../Data_Raw/Unimak.shp")
  berchuk <- st_read("../Data_Raw/BerChukBeau.shp")
  hexMask <- st_read("../Data_Raw/hex_x_ocean/hex_x_ocean.shp") %>%
    select(hexID) %>%
    mutate(AreaKM = c(st_area(.)/1000000)) 
  npac <- st_as_sf(data.frame(name="All Alaska", geometry= st_as_sfc(st_bbox(hexMask))))
  
  areas <- list(npac, goa, uni, berchuk)
  names(areas) <- c( "All Alaska", "Gulf of Alaska", "Eastern Aleutians", "Northern Bering & Chukchi Seas")
  
  summall <- dfnew %>% filter(season == "Summer", timefdy == "All")
  summnight <- dfnew %>% filter(season == "Summer", timefdy == "Night")
  fallall <- dfnew %>% filter(season == "Fall", timefdy == "All")
  fallnight <- dfnew %>% filter(season == "Fall", timefdy == "Night")
  
  riskdf <- data.frame()
  
  for(i in 1:length(areas)){
    loc <- areas[[i]]
    
    temp <- data.frame(
    taxa =taxaName,
    loc =names(areas[i]),
    
    risk_summ_all_mean =mean(summall$riskscaled[st_contains(loc, summall, sparse=F)]),
    risk_summ_all_median =median(summall$riskscaled[st_contains(loc, summall, sparse=F)]),
    risk_summ_all_sd =sd(summall$riskscaled[st_contains(loc, summall, sparse=F)]),
    risk_summ_all_max =max(summall$riskscaled[st_contains(loc, summall, sparse=F)]),
    
    risk_fall_all_mean =mean(fallall$riskscaled[st_contains(loc, fallall, sparse=F)]),
    risk_fall_all_median =median(fallall$riskscaled[st_contains(loc, fallall, sparse=F)]),
    risk_fall_all_sd=sd(fallall$riskscaled[st_contains(loc, fallall, sparse=F)]),
    risk_fall_all_max=max(fallall$riskscaled[st_contains(loc, fallall, sparse=F)]),
    
    risk_summ_night_mean =mean(summnight$riskscaled[st_contains(loc, summnight, sparse=F)]),
    risk_summ_night_median =median(summnight$riskscaled[st_contains(loc, summnight, sparse=F)]),
    risk_summ_night_sd =sd(summnight$riskscaled[st_contains(loc, summnight, sparse=F)]),
    risk_summ_night_max =max(summnight$riskscaled[st_contains(loc, summnight, sparse=F)]),
    
    risk_fall_night_mean =mean(fallnight$riskscaled[st_contains(loc, fallnight, sparse=F)]),
    risk_fall_night_median =median(fallnight$riskscaled[st_contains(loc, fallnight, sparse=F)]),
    risk_fall_night_sd =sd(fallnight$riskscaled[st_contains(loc, fallnight, sparse=F)]),
    risk_fall_night_max =max(fallnight$riskscaled[st_contains(loc, fallnight, sparse=F)])
    )
    riskdf <- rbind(riskdf, temp)
    riskdf[,3:ncol(riskdf)] <- round(riskdf[,3:ncol(riskdf)], 2)
  }
  write.csv(riskdf, paste0("../Data_Processed/Risk_Scaled/Risk_Scaled_", taxaName, ".csv"))
}
lapply(taxaNames, ContinuousRiskMaps)

# Master summary statistics csv file 
filelist <- list.files("../Data_Processed/Risk_Scaled/", pattern=".csv", full.names=T)

files <- lapply(filelist, read.csv)
summstats <- do.call(rbind, files)
write.csv(summstats, "../Data_Processed/Master_Risk_Scaled.csv")
