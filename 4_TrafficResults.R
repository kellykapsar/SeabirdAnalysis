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
library(ggpmisc)

################################
#### Read in Raw Data Files #### 
################################

# Save folder
# Specify file path for folder where results will be saved
savefolder <- "../Data_Processed/"
figfolder <- "../Figures/"

# Read in vessel traffic data 
# hexdir <- "D:/AIS_V2_DayNight_60km6hrgap/Hex_DayNight_Hours/"
# hexList <-list.files(hexdir, pattern=".shp", full.names=TRUE)

# Read in data 
# temp <- lapply(hexList, function(x){st_read(x) %>% st_drop_geometry()})
# hexAll <- do.call(rbind, temp)

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
#   # dplyr::select(hexID, year, month, N_Hrs_Al, D_Hrs_Al, Hrs_Al) %>%
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
# hexdat$propnighttraff <- hexdat$N_Hrs_Al/hexdat$Hrs_Al
# 
# # Save results (for loop above takes a while to run)
# saveRDS(hexdat, "./hexdat.rds")

#########################################################################
# Summary statistics for hexes in the All Alaska analysis 
#########################################################################
hexdat <-readRDS("./hexdat.rds")


# Load in bird data to get hex ids included in the final data set 
birdhex <- st_read("../Data_Processed/FinalShapefiles/AllSeasonsAllTimeOfDay_All Alaska_Seabirds.shp")
ids <- data.frame(hexID = unique(birdhex$hexID))

ids <- left_join(ids, dplyr::select(birdhex, hexID, AreaKM), multiple = "first")

#########################################################
#################### TRAFFIC RESULTS ####################
#########################################################

# Filter out hexes without enough survey effort in both seasons
# Select columns and months of interest 
hexdatnew <- hexdat %>% 
  filter(hexID %in% ids$hexID) %>% 
  select(hexID, year, month, lat, date, Hrs_Al, D_Hrs_Al, N_Hrs_Al, propnight, propnighttraff, ndayspermnth) %>% 
  mutate(month = as.numeric(month)) %>% 
  mutate(N_Hrs_Expected = Hrs_Al*propnight) %>% 
  mutate(nightratio = N_Hrs_Al/N_Hrs_Expected)

hexsumm <- hexdatnew %>% filter(month %in% c(6,7,8))
hexfall <- hexdatnew %>% filter(month %in% c(9,10,11))

#################### Results statements ####################

# Number of hexes included and total area 
length(ids$hexID)
sum(ids$AreaKM)

# Hours of daylight at minimum and maximum latitudes for All Alaska and by season 
max(hexfall$lat, na.rm=TRUE)
mean(24- (hexsumm$propnight[hexsumm$lat == max(hexsumm$lat, na.rm=TRUE)] * 24))
mean(24 - (hexfall$propnight[hexfall$lat == max(hexfall$lat, na.rm=TRUE)] * 24))

min(hexfall$lat, na.rm=TRUE)
mean(24 - (hexsumm$propnight[hexsumm$lat == min(hexsumm$lat, na.rm=TRUE)] * 24))
mean(24 - (hexfall$propnight[hexfall$lat == min(hexfall$lat, na.rm=TRUE)] * 24))


# Total amount of vessel traffic 
sum(hexfall$Hrs_Al)
sum(hexfall$ndayspermnth[hexfall$hexID == first(hexfall$hexID)])
sum(hexfall$Hrs_Al)/sum(hexfall$ndayspermnth[hexfall$hexID == first(hexfall$hexID)])

sum(hexsumm$Hrs_Al)
sum(hexsumm$ndayspermnth[hexsumm$hexID == first(hexsumm$hexID)])
sum(hexsumm$Hrs_Al)/sum(hexsumm$ndayspermnth[hexsumm$hexID == first(hexsumm$hexID)])


mod <- glm(propnighttraff~propnight, data=hexdatnew)
mod <- glm(N_Hrs_Al~N_Hrs_Expected, data=hexdatnew)

summary(mod)


ggplot(hexdatnew, aes(x = N_Hrs_Expected, y = N_Hrs_Al)) +
  geom_point(aes(color=lat)) +
  scale_y_continuous(trans='pseudo_log') +
  scale_x_continuous(trans='pseudo_log') +
  geom_smooth(method="lm", se=TRUE) + 
  geom_abline(color="black", lwd=1) +
  stat_poly_eq(parse=T, aes(label = ..eq.label..), formula=y~x)

ggplot(hexdatnew, aes(x = Hrs_Al, y = N_Hrs_Al)) +
  geom_point(aes(color=lat)) +
  scale_y_continuous(trans='pseudo_log') +
  scale_x_continuous(trans='pseudo_log') +
  geom_smooth(method="lm", se=TRUE) + 
  geom_abline(color="black", lwd=1) +
  stat_poly_eq(parse=T, aes(label = ..eq.label..), formula=y~x)

ggplot(hexdatnew, aes(x = propnight, y = propnighttraff)) +
  geom_point(aes(color=lat)) +
  scale_y_continuous(trans='pseudo_log') +
  scale_x_continuous(trans='pseudo_log') +
  geom_smooth(method="lm", se=TRUE) + 
  geom_abline(color="black", lwd=1) +
  stat_poly_eq(parse=T, aes(label = ..eq.label..), formula=y~x)


#########################################################
#################### SEABIRD RESULTS ####################
#########################################################
# Bird data summaries...
birdseason <-  birdhex %>% 
    st_drop_geometry() %>% 
    filter(timeofday == "All") %>% 
    group_by(season) %>% 
    summarize(DensBird = mean(DensBird), 
              nBird = sum(AllBird), 
              survEff=sum(survEff),
              AreaKM = sum(AreaKM))

# total survey effort
sum(birdseason$nBird)

sum(birdseason$survEff)
sum(birdseason$AreaKM[birdseason$season == "Fall"])
sum(birdseason$survEff)/sum(birdseason$AreaKM[birdseason$season == "Fall"])


######################################################
#################### RISK RESULTS ####################
######################################################



############################################################################

temp <- hexdatnew %>% group_by(hexID) %>% summarize(nightratio = mean(nightratio, na.rm=TRUE)) %>% mutate(nightratio_new = 1-nightratio)
hexdatsf <- left_join(hexMask, temp) 
hexdatsf <- hexdatsf[!is.na(hexdatsf$nightratio_new),]
# Ratio of expected vs. observed nighttime vessel traffic
ggplot() +
  geom_sf(data=studyarea, fill="#8ba761", lwd=0) +
  geom_sf(data=hexdatsf, aes(fill = nightratio_new), color="lightgray") +
  scale_fill_gradient2() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  # labs(caption = paste0("*One operating day is equal to one vessel present in a hex on a given day.")) + 
  theme_bw() +
  theme(text = element_text(size = 18),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5), 
        plot.caption = element_text(size = 8, hjust=0),
        axis.text=element_blank(),
        panel.background = element_rect(fill = "#73b2ff"),
        panel.border =  element_rect(colour = "black"),
        panel.grid.major = element_line(colour = "transparent")) 


# Cells with the greatest difference in daytime and nighttime traffic 
ggplot() +
  geom_sf(data=studyarea, fill="#8ba761", lwd=0) +
  geom_sf(data=hexdatsf,aes(fill = propnight), color="lightgray") +
  # scale_fill_continuous(trans="log",low = "yellow", high = "red", labels=scales::comma, name="Total \nOperating Days") + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  # labs(caption = paste0("*One operating day is equal to one vessel present in a hex on a given day.")) + 
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
ggplot(hexdatnew, aes(x = propnight, y = propnighttraff)) +
  geom_point(aes(size=Hrs_Al))

# Proportion night vs. proportion nighttime vessel traffic colored by month
ggplot(hexdatnew, aes(x = propnight, y = propnighttraff)) +
  geom_point(aes(color=lat)) +
  stat_smooth(method="lm", se=TRUE) + 
  geom_abline(color="black", lwd=1)

ggplot(hexdatnew, aes(x = N_Hrs_Expected, y = N_Hrs_Al)) +
  geom_point(aes(color=lat)) +
  stat_smooth(method="lm", se=TRUE) + 
  geom_abline(color="black", lwd=1)




ggplot(hexdatnew, aes(x = lat, y = N_Hrs_Al)) +
  geom_point(aes(color=propnight)) +
  stat_smooth(method="lm", se=TRUE) + 
  geom_abline(color="black", lwd=1)+
  stat_poly_eq(parse=T, aes(label = ..eq.label..), formula=y~x)

ggplot(hexdatnew %>% filter(Hrs_Al < 2000), aes(x = N_Hrs_Al_Expected, y = N_Hrs_Al)) +
  geom_point(aes(color=lat)) +
  stat_smooth(method="lm", se=TRUE) + 
  geom_abline(color="black", lwd=1)+
  stat_poly_eq(parse=T, aes(label = ..eq.label..), formula=y~x)


ggplot(hexdatnew %>% filter(Hrs_Al < 500), aes(x = N_Hrs_Al_Expected, y = N_Hrs_Al)) +
  geom_point(aes(color=propnight)) +
  geom_smooth(method="lm", se=TRUE, ) + 
  geom_abline(color="black", lwd=1) +
  stat_poly_eq(parse=T, aes(label = ..eq.label..), formula=y~x)

# Proportion night vs. latitude colored by month
ggplot(hexdatnew, aes(x = propnight, y = lat)) +
  geom_point(aes(color=month)) 

# Proportion night vs. proportion nighttime vessel traffic colored by hexID
ggplot(hexdatnew, aes(x = propnight, y = propnighttraff)) +
  geom_point() 

# Models 
m1 <- lm(N_Hrs_Al~propnight*Hrs_Al, data=hexdatnew)
summary(m1)

cor.test(~N_Hrs_Al + N_Hrs_Al_Expected, data=hexdatnew, method="spearman")
cor.test(~N_Hrs_Al + Hrs_Al, data=hexdatnew, method="spearman")
cor.test(~N_Hrs_Al + propnight, data=hexdatnew, method="spearman")

cor(hexdatnew)

hist(residuals(m1),  nbin=30)
ggplot(m1, aes(x = residuals, color=month)) + geom_histogram()

qqnorm(residuals(m1))
qqline(residuals(m1))

m2 <- t.test(x=hexdatnew$propnighttraff, y=hexdatnew$propnight)
m2



# Shipping days by hex over time 
ggplot(hexdat, aes(x=date, y=OpD_Al, group=hexID)) +
  geom_line(alpha=0.5)

ggplot(hexdat, aes(x=date, y=nShp_Al, group=hexID)) +
  geom_line(alpha=0.5)
