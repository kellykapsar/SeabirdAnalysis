

library(sf)
library(tidyverse)
library(ggpubr)

# Master summary statistics csv file 
filelist <- list.files("../Data_Processed/SummaryStatistics/", pattern=".csv", full.names=T)

files <- lapply(filelist, read.csv)
summstats <- do.call(rbind, files)
write.csv(summstats, "../Data_Processed/Master_SummaryStatistics.csv")
rm(filelist)
rm(files)

# Load in basemap for plots
basemap <- st_read("../Data_Raw/AK_CAN_RUS/AK_CAN_RUS.shp") 

hexMask <- st_read("../Data_Raw/hex_x_ocean/hex_x_ocean.shp") %>%
  select(hexID) %>%
  mutate(AreaKM = c(st_area(.)/1000000)) 



filelist <- list.files("../Data_Processed/AllVsNightRisk/", pattern="AllVsNight", full.names=T)

files <- lapply(filelist, read.csv)
hexes <- do.call(rbind, files)
write.csv(hexes, "../Data_Processed/Master_AllVsNightRisk.csv")

hexes$pct_hi_summ_all <- round(hexes$allhighrisk_Summer/hexes$numhexes*100,2)
hexes$pct_hi_summ_night <- round(hexes$nighthighrisk_Summer/hexes$numhexes*100,2)

hexes$pct_hi_fall_all <- round(hexes$allhighrisk_Fall/hexes$numhexes*100,2)
hexes$pct_hi_fall_night <- round(hexes$nighthighrisk_Fall/hexes$numhexes*100,2)

h <- hexes %>% select(pct_hi_summ_all, pct_hi_summ_night, pct_hi_fall_all, pct_hi_fall_night, taxa, studyarea) %>% gather(key=subset,value=riskpct, -taxa, -studyarea)

test <- strsplit(h$subset,split = "_")
h$season <- sapply(test, "[[", 3)
h$timefdy <-  sapply(test, "[[", 4)

h$taxa <- factor(h$taxa, 
                 levels = c("Total Seabirds", "Auklets", "Northern Fulmars", "Seaducks", "Shearwaters", "Storm Petrels"),
                 labels = c("Total Seabirds", "Auklets", "Northern Fulmars", "Seaducks", "Shearwaters", "Storm Petrels"))
h$season <- factor(h$season, levels = c("summ", "fall"), labels = c("Summer", "Fall"))
##################################################

############################################
hnight <- h %>% filter(timefdy == "night", taxa == "Total Seabirds")

p <- ggplot(hnight, aes(x=season, y=riskpct)) +
  geom_bar(position="dodge", stat="identity") +
  # scale_fill_manual(values=c("all" = "#fee227", 
  #                            "night" = "#191970"), 
  #                   labels = c("All", "Night Only"),
  #                   name="Vessel Traffic") +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  facet_wrap(~studyarea) + 
  ylab("Percent of Study Area at High or Very High Risk") + 
  xlab("") +
  theme_bw() + 
  theme(text = element_text(size=20))
p
ggsave(filename = paste0("../Figures/RiskBarGraph_TotalSeabirds_Night_ByRegion.png"),
       plot = p, width=10, height=8, units="in")


hall <- h %>% filter(taxa == "Total Seabirds")

p <- ggplot(hall, aes(x=season, y=riskpct, fill=timefdy)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("all" = "#fee227",
                             "night" = "#191970"),
                    labels = c("All", "Night Only"),
                    name="Vessel Traffic") +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  facet_wrap(~studyarea) + 
  ylab("Percent of Study Area at High or Very High Risk") + 
  xlab("") +
  theme_bw() + 
  theme(text = element_text(size=20))
p
ggsave(filename = paste0("../Figures/RiskBarGraph_TotalSeabirds_ByRegion.png"),
       plot = p, width=12, height=8, units="in")

##################################################
# List of regions 
loclist <- c("Northern Bering & Chukchi Seas", "Gulf of Alaska", "All Alaska", "Eastern Aleutians")

for(i in 1:length(loclist)){
  hnew <- h %>% filter(studyarea == loclist[[i]])

  
  p <- ggplot(hnew, aes(x=season, y=riskpct, fill=timefdy)) +
    geom_bar(position="dodge", stat="identity") +
    scale_fill_manual(values=c("all" = "#fee227", 
                               "night" = "#191970"), 
                      labels = c("All", "Night Only"),
                      name="Vessel Traffic") +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    facet_wrap(~taxa) + 
    ylab("Percent of Study Area at High or Very High Risk") + 
    xlab("") +
    theme_bw() + 
    theme(text = element_text(size=20))
  p
  ggsave(filename = paste0("../Figures/RiskBarGraph_", loclist[i],".png"),
         plot = p, width=16, height=8, units="in")

  ############################################
  
  filelist <- list.files("../Data_Processed/FinalShapefiles/", pattern=loclist[i], full.names=T)
  filelist <- filelist[grep(".shp", filelist)]
  filelist <- filelist[-grep("Seabird", filelist)]
  
  
  files <- lapply(filelist, st_read)
  hexes <- do.call(rbind, files)
  
  test <- hexes %>% st_drop_geometry() %>% filter(riskcat %in% c("high", "veryhigh")) %>% group_by(hexID, season, timefdy) %>% summarize(n=length(unique(taxa)))
  
  hexMaskCrop <- hexMask[hexMask$hexID %in% unique(hexes$hexID),]
  
  combos <- data.frame(id = 1:4, 
                       season = c("Summer", "Fall", "Summer", "Fall"), 
                       tod = c("All", "All", "Night", "Night"))
  
  basemapnew <- basemap %>% st_crop(st_buffer(st_as_sfc(st_bbox(hexMaskCrop)), 10000))

  box <- st_as_sf(name="boundary", st_buffer(st_as_sfc(st_bbox(hexMaskCrop)), 10000))
  
  ################################################# 
  ##################### Risk plot: All Traffic #################
  ################################################# 
  taxariskcounts <- paste0("../Figures/RiskCounts_All_", loclist[i], ".png")
  
  for(j in 1:2){
    t <- test %>% filter(season == combos$season[j] & timefdy == combos$tod[j])
    t$n <- factor(t$n)
    tnew <- left_join(hexMaskCrop, t, by=c("hexID"))
    
    plt <- ggplot() +
      geom_sf(data=box, fill=NA, color=NA,lwd=0) +      
      geom_sf(data=basemapnew, fill="lightgray",lwd=0) +
      geom_sf(data=tnew,aes(fill = n), color="darkgray") +
      scale_fill_manual(values = c("1" = "#ffae52",
                                   "2" = "#a83b00",
                                   "3" = "#281863",
                                   "4" = "black"),
                        na.value = "white",
                        name="Number of Taxa", 
                        drop=F) +
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
      ggtitle(paste0(t$season[1]))
    assign(paste0("p", j), plt)
  }
  comboplot <- ggarrange(p1, p2, ncol=2, nrow=1, common.legend=TRUE, legend = "bottom")
  comboplot <- annotate_figure(comboplot, top = text_grob(tnew$taxa[1], face = "bold", size = 30)) + 
    theme(panel.background = element_rect(fill = "white"))
  
  # Save figure
  if(!file.exists(taxariskcounts)){
    ifelse(loclist[i] == "Eastern Aleutians", 
           ggsave(filename=taxariskcounts, 
                  plot= comboplot, 
                  width=12, height=4, units="in"),
           ggsave(filename=taxariskcounts, 
                  plot= comboplot, 
                  width=12, height=6, units="in"))
  }
  
  ################################################# 
  ##################### Risk plot: Night Traffic #################
  ################################################# 
  taxanightriskcounts <- paste0("../Figures/RiskCounts_Night_", loclist[i], ".png")
  
  for(k in 3:4){
    t <- test %>% filter(season == combos$season[k] & timefdy == combos$tod[k])
    t$n <- factor(t$n)
    tnew <- left_join(hexMaskCrop, t, by=c("hexID"))
    
    plt <- ggplot() +
      geom_sf(data=box, fill=NA, color=NA,lwd=0) +      
      geom_sf(data=basemapnew, fill="lightgray",lwd=0) +
      geom_sf(data=tnew,aes(fill = n), color="darkgray") +
      scale_fill_manual(values = c("1" = "#ffae52",
                                   "2" = "#a83b00",
                                   "3" = "#281863",
                                   "4" = "black"),
                        na.value = "white",
                        name="Number of Taxa", 
                        drop=F) +
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
      ggtitle(paste0(t$season[1]))
    assign(paste0("p", k), plt)
  }
  comboplot <- ggarrange(p1, p2, ncol=2, nrow=1, common.legend=TRUE, legend = "bottom")
  comboplot <- annotate_figure(comboplot, top = text_grob(tnew$taxa[1], face = "bold", size = 30)) + 
    theme(panel.background = element_rect(fill = "white"))
  
  # Save figure
  if(!file.exists(taxariskcounts)){
    ifelse(loclist[i] == "Eastern Aleutians", 
           ggsave(filename=taxariskcounts, 
                  plot= comboplot, 
                  width=12, height=4, units="in"),
           ggsave(filename=taxariskcounts, 
                  plot= comboplot, 
                  width=12, height=6, units="in"))
  }
  
}


# hall <- h %>% filter(timefdy == "all")
# 
# p <- ggplot(hall, aes(x=taxa, y=riskpct, fill=season)) +
#   geom_bar(position="dodge", stat="identity") +
#   scale_fill_manual(values=c("#fee227", "#191970")) +
#   theme(axis.text.x = element_text(angle=45, hjust=1)) 
# p
# ggsave(filename = "../Figures/RiskBarGraph.png",
#        plot = p, width=16, height=8, units="in")



##################
# Log bird density vs. log vessel traffic - Risk Categories
# akall <- st_read("../Data_Processed/FinalShapefiles/AllSeasonsAlltimefdy_All Alaska_Seabirds.shp") %>% st_drop_geometry()
# akallsumm <- akall %>% filter(subset == "Summer_All")
# akallfall <- akall %>% filter(subset == "Fall_All")
# 
# akallsumm$risk <- factor(akallsumm$risk, levels = c("low", "medium", "high", "veryhigh"), labels = c("Low", "Medium", "High", "Very High"))
# akallfall$risk <- factor(akallfall$risk, levels = c("low", "medium", "high", "veryhigh"), labels = c("Low", "Medium", "High", "Very High"))
# 
# risk_pal <- c("#73b2ff", "#55fe01", "#ffff01", "#e31a1c")
# names(risk_pal) <- levels(akallsumm$risk)
# 
# ggplot(akallsumm, aes(x = DensBird, y = AllShip)) +
#   geom_point(aes(color=risk)) +
#   scale_color_manual(values = risk_pal) +
#   scale_x_continuous(trans=scales::pseudo_log_trans(base = 10)) +
#   scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
#   theme(text = element_text(size=30))
# 
# ggplot(akallfall, aes(x = DensBird, y = AllShip)) +
#   geom_point(aes(color=risk)) +
#   scale_color_manual(values = risk_pal) +
#   scale_x_continuous(trans=scales::pseudo_log_trans(base = 10)) +
#   scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
#   theme(text = element_text(size=30))
# 
