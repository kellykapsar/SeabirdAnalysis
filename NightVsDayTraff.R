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
hexdir <- "D:/AIS_V2_DayNight_60km6hrgap/Hex_DayNight/"
hexList <-list.files(hexdir, pattern=".shp", full.names=TRUE)

# Read in data 
temp <- lapply(hexList, function(x){st_read(x) %>% st_drop_geometry()})
hexAll <- do.call(rbind, temp)

# Read in a custom-created hexagon mask, which masks out all parts of hexagons that are land and calculates
# total marine area. This will be useful in calculating survey effort.
hexMask <- st_read("../Data_Raw/hex_x_ocean/hex_x_ocean.shp") %>%
  select(hexID) %>%
  mutate(AreaKM = c(st_area(.)/1000000)) 

# Load in basemap for plots
studyarea <- st_read("../Data_Raw/AK_CAN_RUS/AK_CAN_RUS.shp") %>% st_crop(st_buffer(hexMask, 100000))

##########################################################
# Calculate proportion of night per hex per month  
##########################################################

# # Calculate the latitude for the centroid of each hex (to be used in daylight calculations)
# hexcent <- st_read(hexList[[1]]) %>% dplyr::select(hexID) %>% st_transform(4326)
# hexcent$lat <- st_coordinates(st_centroid(hexcent))[,"Y"]
# hexcent <- st_drop_geometry(hexcent)
# 
# # Select metrics and calculate julian dates of month for sunlight calculations
# hexdat <- hexAll %>%
#   # dplyr::select(hexID, year, month, N_OpD_Al, D_OpD_Al, OpD_Al) %>%
#   left_join(., hexcent, by=c("hexID" = "hexID")) %>%
#   mutate(date =  as.POSIXct(strptime(paste0(year, "-", month, "-01"),format="%Y-%m-%d", tz="UTC")),
#          jdatestart = as.numeric(as.character(format(date, "%j"))),
#          ndayspermnth = lubridate::days_in_month(date),
#          jdateend = (jdatestart + ndayspermnth-1))
# 
# # Calculate the average proportion of daytime in each hex in each month
# hexdat$propnight <- NA
# for(i in 1:length(hexdat$hexID)){
#   print(i)
#   hexdat$propnight[i] <- mean(chillR::daylength(latitude = hexdat$lat[i],
#                                                    JDay = hexdat$jdatestart[i]:hexdat$jdateend[i])$Daylength)
# }
# 
# # Replace NA values with zero
# hexdat[is.na(hexdat)] <- 0
# 
# # Convert to proportion night
# hexdat$propnight <- 1-hexdat$propnight/24
# 
# # Calculate proportion of vessel traffic at night
# hexdat$propnighttraff <- hexdat$N_OpD_Al/hexdat$OpD_Al
# 
# # Save results (for loop above takes a while to run)
# saveRDS(hexdat, "./hexdat.rds")

#########################################################################
# Evaluate differences between daytime and nighttime vessel traffic  
#########################################################################
hexdat <-readRDS("./hexdat.rds")

# Remove hexes with no vessel traffic and combine years to look at monthly data
hexdatall <- hexdat %>%
  filter(OpD_Al > 0) %>%
  group_by(hexID) %>%
  summarize(N_OpD_Al = sum(N_OpD_Al),
            D_OpD_Al = sum(D_OpD_Al),
            OpD_Al = sum(OpD_Al),
            propnighttraff = N_OpD_Al/OpD_Al,
            propnight=mean(propnight),
            lat = mean(lat)) %>%
  ungroup()

# Check to make sure date calculations worked
unique(hexdat[c("month", "ndayspermnth", "jdatestart", "jdateend")])

#### Remove hexes with insufficient vessel traffic that would skew proportions 
hexdatnew <- hexdatall %>% filter(OpD_Al > 100, !propnight %in% c(0,1))

# Calculate ratio of day to night traffic
hexdatnew$daynight <- hexdatnew$D_OpD_Al/hexdatnew$N_OpD_Al
hexdatnew$subtract <- hexdatnew$D_OpD_Al - hexdatnew$N_OpD_Al
summary(hexdatnew$daynight)
summary(hexdatnew$subtract)
hist(hexdatnew$subtract)

# Get hexIDs for top 10% of values for day-night traff
topday <- hexdatnew %>% 
  arrange(desc(subtract)) %>% 
  slice_head(prop=0.05) 

topnight <- hexdatnew %>% 
  arrange(subtract) %>% 
  slice_head(prop=0.05) 

dayonly <- topday[!topday$hexID %in% topnight$hexID,]
nightonly <- topnight[!topnight$hexID %in% topday$hexID,]

length(unique(nightonly$hexID))
length(unique(dayonly$hexID))

length(which(topday$hexID %in% topnight$hexID))
length(which(topnight$hexID %in% topday$hexID))

topdaysf <- topday  %>% left_join(hexMask) %>% st_as_sf()
topnightsf <- topnight  %>% left_join(hexMask) %>% st_as_sf()

hexdatsf <- hexdatall  %>% left_join(hexMask) %>% st_as_sf()


# Cells with the greatest difference in daytime and nighttime traffic 
ggplot() +
  geom_sf(data=studyarea, fill="#8ba761", lwd=0) +
  geom_sf(data=hexdatsf,aes(fill = propnight), color="lightgray") +
  # scale_fill_continuous(trans="log",low = "yellow", high = "red", labels=scales::comma, name="Total \nOperating Days") + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(caption = paste0("*One operating day is equal to one vessel present in a hex on a given day.")) + 
  theme_bw() +
  theme(text = element_text(size = 18),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5), 
        plot.caption = element_text(size = 8, hjust=0),
        axis.text=element_blank(),
        panel.background = element_rect(fill = "#73b2ff"),
        panel.border =  element_rect(colour = "black"),
        panel.grid.major = element_line(colour = "transparent")) 

# INDIVIDUAL CELL MAPS 
ggplot() +
  geom_sf(data=studyarea, fill="#8ba761", lwd=0) +
  geom_sf(data=hexMask[hexMask$hexID == 1984,], fill="red", color="red") + # PUT HEXID HERE... 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(text = element_text(size = 18),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5), 
        plot.caption = element_text(size = 8, hjust=0),
        axis.text=element_blank(),
        panel.background = element_rect(fill = "#73b2ff"),
        panel.border =  element_rect(colour = "black"),
        panel.grid.major = element_line(colour = "transparent")) 


#### Plots 

# Histograms 
ggplot(hexdatall, aes(x = propnight)) + geom_histogram()
ggplot(hexdatnew, aes(x = propnight)) + geom_histogram()

ggplot(hexdatall, aes(x = propnighttraff)) + geom_histogram()
ggplot(hexdatnew, aes(x = propnighttraff)) + geom_histogram()

ggplot(hexdatnew, aes(x = subtract)) + geom_histogram()

# Proportion night vs. proportion nighttime vessel traffic 
# Expect 1:1 correlation if no selection 
ggplot(hexdatnew, aes(x = propnight, y = daynight)) +
  geom_point(aes(size=OpD_Al))

# Proportion night vs. proportion nighttime vessel traffic colored by month
ggplot(hexdatnew, aes(x = D_OpD_Al, y = N_OpD_Al)) +
  geom_point(aes(color=lat)) +
  stat_smooth(method="lm", se=TRUE) + 
  geom_abline(color="black", lwd=1)

ggplot(hexdatnew, aes(x = D_OpD_Al, y = N_OpD_Al)) +
  geom_point(aes(color=propnight)) +
  stat_smooth(method="lm", se=TRUE) + 
  geom_abline(color="black", lwd=1)

ggplot(hexdatnew %>% filter(OpD_Al < 2000), aes(x = D_OpD_Al, y = N_OpD_Al)) +
  geom_point(aes(color=lat)) +
  stat_smooth(method="lm", se=TRUE) + 
  geom_abline(color="black", lwd=1)

ggplot(hexdatnew %>% filter(OpD_Al < 500), aes(x = D_OpD_Al, y = N_OpD_Al)) +
  geom_point(aes(color=propnight)) +
  stat_smooth(method="lm", se=TRUE) + 
  geom_abline(color="black", lwd=1)

# Proportion night vs. latitude colored by month
ggplot(hexdatnew, aes(x = propnight, y = lat)) +
  geom_point() 

# Proportion night vs. proportion nighttime vessel traffic colored by hexID
ggplot(hexdatnew, aes(x = propnight, y = propnighttraff)) +
  geom_point(aes(color = hexID)) 

# Models 
m1 <- lm(N_OpD_Al~D_OpD_Al, data=hexdatnew)
summary(m1)

cor.test(~N_OpD_Al + D_OpD_Al, data=hexdatnew, method="spearman")
cor(hexdatnew)

hist(residuals(m1),  nbin=30)
ggplot(m1, aes(x = residuals, color=month)) + geom_histogram()

qqnorm(residuals(m1))
qqline(residuals(m1))

m2 <- t.test(x=hexdatnew$propnighttraff, y=hexdatnew$propnight)
m2


