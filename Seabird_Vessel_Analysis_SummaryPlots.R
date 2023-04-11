

library(sf)
library(tidyverse)

loc <- "All Alaska"


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
  t <- st_read(files[37], quiet=T) %>% st_drop_geometry()
  df$nhex[i] <- length(t$hexID)
  df$risk[i] <- sum(t$risk %in% c("high", "veryhigh"))
  df$riskpct[i] <- round(df$risk[i]/df$nhex[i]*100, 2)
  df$den[i] <- length(which(t$ClassBird == 3))
  df$denpct[i] <- round(df$den[i]/df$nhex[i]*100, 2)
}

tnext <- 

ggplot(t, aes(x = log(DensBird)+1, y = log(AllShip)+1)) +
  geom_point(aes(color=risk)) +
  theme(text = element_text(size=30)) 

ggplot(df, aes(x = denpct, y = riskpct)) +
  geom_point(aes(color=season)) +
  theme(text = element_text(size=30)) 

ggplot(df, aes(x = denpct, y = riskpct)) +
  geom_point(aes(color=season)) +
  theme(text = element_text(size=30)) +
  stat_smooth(method="lm", se=TRUE)  +
  geom_abline(color="black", lwd=1)

p <- ggplot(df, aes(x=taxa, y=riskpct, fill=season)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("#44828f", "#7ab845")) +
  theme(axis.text.x = element_text(angle=45, hjust=1))

ggsave(filename = "../Figures/RiskBarGraph.png",
       plot = p, width=16, height=8, units="in")


17/20
