################################################################################
# TITLE: Summary Plots from Hex Metadata
# PURPOSE: This code uses the metadata files developed in the H2_AIS_SpeedHexes.R
#   script to create summary plots of vessel traffic over the study period. 
# AUTHOR: Kelly Kapsar
# CREATED: ??
# LAST UPDATED ON: 2022-06-09
# 
# NOTE:
################################################################################

# Load Libraries 
library(tidyr)
library(dplyr)
library(ggplot2)
library(sf)
library(paletteer)

# Combine all metadata files
temphpcc <-paste0("../Data_Processed_Hex/Metadata/", 
              list.files("../Data_Processed_Hex/Metadata", pattern='Metadata'))
temphpcc <- lapply(temphpcc, read.csv)
meta <- do.call(rbind, temphpcc)
# st_write(meta, "../Data_Processed_Hex/Metadata/Metadata_SpeedHex_All.csv")


# Combine all runtimes files
temphpcc2 <-paste0("../Data_Processed_Hex/Metadata/", 
                  list.files("../Data_Processed_Hex/Metadata", pattern='Runtimes'))
temphpcc2 <- lapply(temphpcc2, read.csv)
runs <- do.call(rbind, temphpcc2)
# st_write(runs, "../Data_Processed_Hex/Metadata/Runtimes_SpeedHex_All.csv")

###############################################################################################################

meta <- read.csv("../Data_Processed_Hex/Metadata/Metadata_SpeedHex_ALL.csv")

# Turn year and month columns into a date format
meta$yrmnth <- as.Date(paste(meta$yr,meta$mnth,"01",sep="-"), format="%Y-%m-%d")

# Number of unique vessels over time
meta_long <- meta %>% 
             dplyr::select(ntotal_pts, ntotal_mmsis, ntotal_aisids, yrmnth) %>% 
             gather(key="type", value="count", -yrmnth)

type.labs <- c("Ships", "Operating Days", "Points")
names(type.labs) <- c("ntotal_mmsis", "ntotal_aisids", "ntotal_pts")


meta_long$count[which(meta_long$type == "ntotal_pts")] <- meta_long$count[which(meta_long$type == "ntotal_pts")]/1000000


p1 <- ggplot(meta_long, aes(x=yrmnth, y=count)) +
  geom_line(lwd=1) +
  xlab("Year") +
  scale_x_date(date_labels = "%Y", date_breaks="1 year", expand=c(0,0)) +
  theme_bw() +
  theme(text = element_text(size=30)) +
  facet_grid(rows=vars(type), scales="free_y",
             labeller = as_labeller(c(ntotal_aisids = "Operating\nDays", 
                                      ntotal_mmsis = "Unique\nShips", 
                                      ntotal_pts = "AIS\nSignals\n(Millions)") ) )  +
  ylab(NULL) +
  theme(strip.background = element_blank(),
        strip.placement = "outside", text = element_text(size=35))

p1

ggsave("D:/AlaskaConservation_AIS_20210225/Figures/NumberShipsPerMonth.png",width=20,height=10,units='in',dpi=300)


# Shipping days by type over time
longmeta_naisids <- meta %>% 
  dplyr::select(yrmnth, ntank_aisids, nfish_aisids, ncargo_aisids, nother_aisids) %>% 
  gather(key=type,value=naisids, ntank_aisids:nother_aisids)

ggplot(longmeta_naisids, aes(x=yrmnth, y=naisids)) +
  geom_line(aes(color=type), lwd=1)+
  scale_color_brewer(palette = "Dark2", name="Ship Type", 
                     breaks=c("nfish_aisids", "nother_aisids", "ncargo_aisids", "ntank_aisids"), 
                     labels=c("Fishing", "Other", "Cargo", "Tanker")) +
  xlab("Year") +
  scale_x_date(date_labels = "%Y", date_breaks="1 year", expand=c(0,0)) +
  scale_y_continuous(expand=c(0,1000)) +
  ylab("Operating days") + 
  theme_bw() +
  theme(text = element_text(size=30))
# ggsave("D:/AlaskaConservation_AIS_20210225/Figures/OperatingDaysPerMonth_ByType.png",width=25,height=10,units='in',dpi=300)


longmeta_npts <- meta %>% 
  dplyr::select(yrmnth, ntank_pts, nfish_pts, ncargo_pts, nother_pts) %>% 
  gather(key=type,value=npts, ntank_pts:nother_pts)

ggplot(longmeta_npts, aes(x=yrmnth, y=npts)) +
  geom_line(aes(color=type), lwd=1)+
  scale_color_brewer(palette = "Dark2", name="Ship Type", 
                     breaks=c("nfish_pts", "nother_pts", "ncargo_pts", "ntank_pts"), 
                     labels=c("Fishing", "Other", "Cargo", "Tanker")) +
  xlab("Year") +
  scale_x_date(date_labels = "%Y", date_breaks="1 year", expand=c(0,0)) +
  scale_y_continuous(expand=c(0,1000)) +
  ylab("Points") + 
  theme_bw() + 
  theme(text = element_text(size=30))



