

library(sf)
library(tidyverse)

loc <- "Northern Bering"
loc <- "Gulf of Alaska"
loc <- "All Alaska"
loc <- "Aleutian"


filelist <- list.files("../Data_Processed/AllVsNightRisk/", pattern=loc, full.names=T)

files <- lapply(filelist, read.csv)
hexes <- do.call(rbind, files)

hexes$pct_hi_summ_all <- round(hexes$allhighrisk_Summer/hexes$numhexes*100,2)
hexes$pct_hi_summ_night <- round(hexes$nighthighrisk_Summer/hexes$numhexes*100,2)

hexes$pct_hi_fall_all <- round(hexes$allhighrisk_Fall/hexes$numhexes*100,2)
hexes$pct_hi_fall_night <- round(hexes$nighthighrisk_Fall/hexes$numhexes*100,2)

h <- hexes %>% select(pct_hi_summ_all, pct_hi_summ_night, pct_hi_fall_all, pct_hi_fall_night, taxa) %>% gather(key=subset,value=riskpct, -taxa)

test <- strsplit(h$subset,split = "_")
h$season <- sapply(test, "[[", 3)
h$timeofday <-  sapply(test, "[[", 4)

p <- ggplot(h, aes(x=season, y=riskpct, fill=timeofday)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("#fee227", "#191970")) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  facet_wrap(~taxa)
  # annotate("text", x = 1:4, y = - 400,
  #          label = rep(c("Variety 1", "Variety 2"), 2)) +
  # annotate("text", c(1.5, 3.5), y = -0.24, label = c("Irrigated", "Dry")) +
  # theme_classic() +
  # theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
  #       axis.title.x = element_blank(),
  #       axis.text.x = element_blank())
p
ggsave(filename = "../Figures/RiskBarGraph.png",
       plot = p, width=16, height=8, units="in")


##################
# Log bird density vs. log vessel traffic - Risk Categories
akall <- st_read("../Data_Processed/FinalShapefiles/AllSeasonsAllTimeOfDay_All Alaska_Seabirds.shp") %>% st_drop_geometry()
akallsumm <- akall %>% filter(subset == "Summer_All")
akallfall <- akall %>% filter(subset == "Fall_All")

akallsumm$risk <- factor(akallsumm$risk, levels = c("low", "medium", "high", "veryhigh"), labels = c("Low", "Medium", "High", "Very High"))
akallfall$risk <- factor(akallfall$risk, levels = c("low", "medium", "high", "veryhigh"), labels = c("Low", "Medium", "High", "Very High"))

risk_pal <- c("#73b2ff", "#55fe01", "#ffff01", "#e31a1c")
names(risk_pal) <- levels(akallsumm$risk)

ggplot(akallsumm, aes(x = DensBird, y = AllShip)) +
  geom_point(aes(color=risk)) +
  scale_color_manual(values = risk_pal) +
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  theme(text = element_text(size=30))

ggplot(akallfall, aes(x = DensBird, y = AllShip)) +
  geom_point(aes(color=risk)) +
  scale_color_manual(values = risk_pal) +
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  theme(text = element_text(size=30))

