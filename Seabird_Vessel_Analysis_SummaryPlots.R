

library(sf)
library(tidyverse)

loc <- "Gulf of Alaska"


files <- intersect(list.files("../Data_Processed/", pattern=loc, full.names=T), 
                   list.files("../Data_Processed/", pattern=".shp", full.names = T))




df <- data.frame(do.call(rbind, strsplit(files, "_"))) %>% select(-X1, -X2, -X6)
colnames(df) <- c("location", "taxa", "season")
df$season <- factor(df$season, levels=c("Summer", "Fall"))
df$nhex <- NA
df$riskpct <- NA
df$risk <- NA
df$den <- NA
df$denpct <- NA

for(i in 1:length(files)){
  t <- st_read(files[i], quiet=T) %>% st_drop_geometry()
  df$nhex[i] <- length(t$hexID)
  df$risk[i] <- sum(t$risk %in% c("high", "veryhigh"))
  df$riskpct[i] <- round(df$risk[i]/df$nhex[i]*100, 2)
  df$den[i] <- length(which(t$ClassBird == 3))
  df$denpct[i] <- round(df$den[i]/df$nhex[i]*100, 2)
}

akallsumm <- st_read("../Data_Processed/FinalDF_All Alaska_Seabirds_Summer_NightOnlyFALSE.shp") %>% st_drop_geometry()
akallfall <- st_read("../Data_Processed/FinalDF_All Alaska_Seabirds_Fall_NightOnlyFALSE.shp") %>% st_drop_geometry()

# Log bird density vs. log vessel traffic
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

######

ggplot(df, aes(x = denpct, y = riskpct)) +
  geom_point(aes(color=season)) +
  theme(text = element_text(size=30)) 

ggplot(df, aes(x = denpct, y = riskpct)) +
  geom_point(aes(color=season)) +
  theme(text = element_text(size=30)) +
  stat_smooth(method="lm", se=TRUE)  +
  geom_abline(color="black", lwd=1)

dfnew <- df %>% filter(taxa %in% c("Albatross", "Shearwaters", "Auklets", "Murres", "Kittiwakes", "Gulls", "Storm Petrels", "Northern Fulmars"))

p <- ggplot(df, aes(x=reorder(taxa, -riskpct), y=riskpct, fill=season)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("#7ab845", "#44828f")) +
  theme(axis.text.x = element_text(angle=45, hjust=1))

ggsave(filename = "../Figures/RiskBarGraph.png",
       plot = p, width=16, height=8, units="in")

