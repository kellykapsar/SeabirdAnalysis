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
library(chillR)

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

# Calculate the latitude for the centroid of each hex (to be used in daylight calculations)
hexcent <- st_read(hexList[[1]]) %>% dplyr::select(hexID) %>% st_transform(4326)
hexcent$lat <- st_coordinates(st_centroid(hexcent))[,"Y"] 
hexcent <- st_drop_geometry(hexcent)

# Select metrics and calculate julian dates of month for sunlight calculations
nightprop <- hexAll %>% 
  dplyr::select(hexID, year, month, N_OpD_Al, OpD_Al) %>% 
  left_join(., hexcent, by=c("hexID" = "hexID")) %>% 
  mutate(date =  as.POSIXct(strptime(paste0(year, "-", month, "-01"),format="%Y-%m-%d", tz="UTC")), 
         jdatestart = as.numeric(as.character(format(date, "%j"))), 
         ndayspermnth = lubridate::days_in_month(date), 
         jdateend = (jdatestart + ndayspermnth-1)) 

# Calculate the average proportion of daytime in each hex in each month 
nightprop$propnight <- NA
for(i in 1:length(nightprop$hexID)){
  print(i)
  nightprop$propnight[i] <- mean(chillR::daylength(latitude = nightprop$lat[i], 
                                                   JDay = nightprop$jdatestart[i]:nightprop$jdateend[i])$Daylength)
}


# Replace NA values with zero 
nightprop[is.na(nightprop)] <- 0

# Convert to proportion night 
nightprop$propnight <- 1-nightprop$propnight/24

# Calculate proportion of vessel traffic at night 
nightprop$propnighttraff <- nightprop$N_OpD_Al/nightprop$OpD_Al

# Save results (for loop above takes a while to run)
saveRDS(nightprop, "./nightprop.rds")
nightprop <-readRDS("./nightprop.rds")

# Remove hexes with no vessel traffic and combine years to look at monthly data 
nightpropall <- nightprop %>%  
  filter(OpD_Al > 0) %>% 
  group_by(hexID, month) %>% 
  summarize(N_OpD_Al = sum(N_OpD_Al), 
            OpD_Al = sum(OpD_Al), 
            propnighttraff = N_OpD_Al/OpD_Al,
            propnight=mean(propnight), 
            lat = mean(lat)) %>% 
  ungroup()

# Check to make sure date calculations worked 
unique(nightprop[c("month", "ndayspermnth", "jdatestart", "jdateend")])

#### Remove hexes with insufficient vessel traffic that would skew proportions 
nightpropnew <- nightpropall %>% filter(OpD_Al > 100)

#### Plots 

# Histograms 
ggplot(nightpropall, aes(x = propnight, color=month)) + geom_histogram()
ggplot(nightpropnew, aes(x = propnight, color=month)) + geom_histogram()

ggplot(nightpropall, aes(x = propnighttraff, color=month)) + geom_histogram()
ggplot(nightpropnew, aes(x = propnighttraff, color=month)) + geom_histogram()

# Proportion night vs. proportion nighttime vessel traffic 
# Expect 1:1 correlation if no selection 
ggplot(nightpropnew, aes(x = propnight, y = propnighttraff)) +
  geom_point(aes(color=month, size=OpD_Al)) +
  stat_smooth(method="lm", se=TRUE) + 
  geom_abline(color="red")

# Proportion night vs. proportion nighttime vessel traffic colored by month
ggplot(nightpropnew, aes(x = propnight, y = propnighttraff)) +
  geom_point(aes(color = month)) 

# Proportion night vs. latitude colored by month
ggplot(nightpropnew, aes(x = propnight, y = lat)) +
  geom_point(aes(color = month)) 

# Proportion night vs. proportion nighttime vessel traffic colored by hexID
ggplot(nightpropnew, aes(x = propnight, y = propnighttraff)) +
  geom_point(aes(color = hexID)) 

# Models 
m1 <- lm(propnighttraff~propnight, data=nightpropnew)
summary(m1)

hist(residuals(m1))
ggplot(m1, aes(x = resdiuals, color=month)) + geom_histogram()

qqnorm(residuals(m1))
qqline(residuals(m1))

m2 <- t.test(x=nightpropnew$propnighttraff, y=nightpropnew$propnight)
m2


