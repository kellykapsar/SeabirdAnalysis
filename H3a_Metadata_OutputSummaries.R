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
hexdir <- "D:/AIS_V2_DayNight_60km6hrgap/Hex_DayNight/"
hexList <-list.files(hexdir, pattern="Metadata", full.names=TRUE)

# hexList <- lapply(hexList, read.csv)
# meta <- do.call(bind_rows, hexList)
# write.csv(meta, "D:/AIS_V2_DayNight_60km6hrgap/Hex_DayNight/Metadata_SpeedHex_All.csv")
meta <- read.csv("D:/AIS_V2_DayNight_60km6hrgap/Hex_DayNight/Metadata_SpeedHex_All.csv")

# Combine all runtimes files
# hexList2 <-list.files(hexdir, pattern="Runtimes", full.names=TRUE)
# hexList2 <- lapply(hexList2, read.csv)
# runs <- do.call(bind_rows, hexList2)
# write.csv(runs, "D:/AIS_V2_DayNight_60km6hrgap/Hex_DayNight/Runtimes_SpeedHex_All.csv")

# Read in data from individual hexes 
hexdat <-readRDS("./hexdat.rds") # Calculated in NightVsDayTraff.R script. 


###############################################################################################################

# Turn year and month columns into a date format
meta$yrmnth <- as.Date(paste(meta$yr,meta$mnth,"01",sep="-"), format="%Y-%m-%d")

# Number of unique vessels over time
meta_long <- meta %>% 
             dplyr::select(ntotal_mmsis, yrmnth) %>% 
             gather(key="type", value="count", -yrmnth)

type.labs <- c("Ships")
names(type.labs) <- c("ntotal_mmsis")


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
  dplyr::select(yrmnth, ntank_mmsis, nfish_mmsis, ncargo_mmsis, ntug_mmsis, npass_mmsis, nsail_mmsis, npleas_mmsis, nother_mmsis) %>% 
  gather(key=type,value=nmmsis, ntank_mmsis:nother_mmsis)

ggplot(longmeta_naisids, aes(x=yrmnth, y=nmmsis)) +
  geom_line(aes(color=type), lwd=1)+
  scale_color_brewer(palette = "Dark2", name="Ship Type", 
                     breaks=c("ntank_mmsis", "nfish_mmsis", "ncargo_mmsis", "ntug_mmsis", "npass_mmsis", "nsail_mmsis", "npleas_mmsis", "nother_mmsis"), 
                     labels=c("Tanker", "Fishing", "Cargo", "Other", "Tug", "Passenger", "Sailing", "Pleasure")) +
  xlab("Year") +
  scale_x_date(date_labels = "%Y", date_breaks="1 year", expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
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

# Shipping days by hex over time 
ggplot(hexdat, aes(x=date, y=OpD_Al, group=hexID)) +
  geom_line(alpha=0.5)

ggplot(hexdat, aes(x=date, y=nShp_Al, group=hexID)) +
  geom_line(alpha=0.5)



