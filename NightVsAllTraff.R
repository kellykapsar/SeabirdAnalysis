################################################################################
# TITLE: AIS Day vs. Night 
#
# PURPOSE: This script compares vessel traffic occurring at night to all vessel 
# traffic in the North Pacific with the goal of determiniing whether nighttime 
# vessel traffic is just a subset of all vessel traffic or if there are areas 
# where vessels are more or less likely to travel at night. 
#
# AUTHORS: Kelly Kapsar 
#
# CREATED: 17 January 2023
# LAST UPDATED ON: 
################################################################################


# Import libraries
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)

################################
#### Read in Raw Data Files #### 
################################

# Save folder
# Specify file path for folder where results will be saved
savefolder <- "../Data_Processed/"
figfolder <- "../Figures/"

# Read in vessel traffic data 
hexdir <- "D:/AIS_V2_DayNight_60km6hrgap/Hex/"
hexList <-list.files(hexdir, pattern=".shp", full.names=TRUE)

# Read in blank hexagon template 
# hex <- st_read("../Data_Raw/BlankHexes.shp")

# Read in a custom-created hexagon mask, which masks out all parts of hexagons that are land and calculates
# total marine area. This will be useful in calculating survey effort.
hexMask <- st_read("../Data_Raw/hex_x_ocean/hex_x_ocean.shp") %>%
  select(hexID) %>%
  mutate(AreaKM = c(st_area(.)/1000000)) 

# Load in basemap for plots
studyarea <- st_read("../Data_Raw/AK_CAN_RUS/AK_CAN_RUS.shp")


###################################
#### Specify Input Information #### 
###################################

# Specify metric used to calculate vessel activity 
metric <- "OperatingDays" # MUST be "OperatingDays" or "Ships"
## OperatingDays = number of ship days per month (i.e., the same ship in the same hex each day for a month = 30)
## Ships = number of unique ships per month (i.e., the same ship int h same hex each day for a month = 1)

nightonly <- TRUE #If true, will only calculate nighttime vessel traffic (ignoring daytime)
# If false, will calculate all vessel traffic, including both day and night

months <- c(9:11) # Numeric value(s) of months to be included in the analysis
# months <- c(6:8)

# monthsname <- "Fall" # Text label describing the numeric months in the analysis (e.g., "Summer", "Annual")
monthsname <- "Summer"

###########################
#### Load in Functions #### 
###########################

source("./Seabird_Vessel_Analysis_Functions.R")


jan <- st_read(hexList[[1]]) %>% dplyr::select(hexID, year, month, N_OpD_Al, OpD_Al)

# Read in data 
temp <- lapply(hexList, function(x){st_read(x) %>% st_drop_geometry()})
hexAll <- do.call(rbind, temp)

# Select metrics 

test <- hexAll %>% dplyr::select(hexID, year, month, N_OpD_Al, OpD_Al)

# Replace NA values with zero 
test[is.na(test)] <- 0

# Calculate proportion of vessel traffic at night for each hex in each month
test$propnight <- test$N_OpD_Al/test$OpD_Al

test$date <- as.POSIXct(strptime(paste0(test$year, "-", test$month, "-01"),format="%Y-%m-%d", tz="UTC"))


testwide <- test %>% dplyr::select(date, hexID, propnight) %>% spread(., key = date, value=propnight)

hexpropnight <- left_join(hexMask, testwide, by=c("hexID" = "hexID"))

hexpropnight$lat <- st_coordinates(st_centroid(hexpropnight))[,"Y"]

testnona <- test[test$OpD_Al != 0,]

testnona <- left_join(testnona, st_drop_geometry(hexpropnight[, c("hexID", "lat")]), by =c("hexID"="hexID"))


ggplot(testnona, aes(x = month, y = propnight)) +
  geom_point(aes(color = lat, size=OpD_Al)) +
  scale_size(trans="log10")


ggplot(testnona, aes(x = month, y = lat)) +
  geom_jitter(aes(color = propnight, size=OpD_Al)) +
  scale_color_gradient2(low="yellow", mid="gray", high="black", 
                        limits = c(0,1)) +
  scale_size(trans="log10") 


ggplot() +
  geom_sf(data=hexpropnight, aes(fill=2015-01-01), lwd=0) +
  # geom_sf(data=finalCombined,aes(fill = risk)) +
  # scale_fill_manual(values = c("low" = "#73b2ff",
  #                              "medium" = "#55fe01",
  #                              "high" = "#ffff01",
  #                              "veryhigh" = "#e31a1c"), 
  #                   name="Risk", 
  #                   labels = c("Low", "Medium", "High", "Very High")) +
  xlab("") +
  ylab("") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  # ggtitle(ifelse(night == FALSE, paste0(taxaLabel, " - ", monthsname, " - All Traffic"), 
  #                paste0(taxaLabel, " - ", monthsname, " - Night Traffic"))) +
  theme_bw() +
  theme(text = element_text(size = 25),
        axis.text=element_blank(),
        panel.background = element_rect(fill = "#bcc7dd"),
        panel.border =  element_rect(colour = "black"),
        panel.grid.major = element_line(colour = "transparent"))


st_write(hexpropnight, "../Sandbox/PropTraffAtNight_20230207.shp")
# Save figure 
ggsave(filename = paste0(figfolder,"Map_",taxaLabel,"_",monthsname,"_NightOnly",night[1],".png"), 
       plot = p1, width=12, height=8, units="in")
