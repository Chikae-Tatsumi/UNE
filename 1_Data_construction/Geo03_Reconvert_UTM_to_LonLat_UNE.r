# https://rdrr.io/cran/PBSmapping/man/convUL.html
library(oce)

# Import files
setwd("~/R/Analysis/2_UNE/Autocorrelation")
data <- read.csv("UTM_UNE_rep.separated.csv",header=T)

# for UTM_zone = 19
subset_19 <- subset(data, data$UTM_zone==19)
lonlat <- utm2lonlat(subset_19$Corrected_UTM_X, subset_19$Corrected_UTM_Y, zone = 19, hemisphere = "N", km = TRUE)
lonlat_19 <- cbind(subset_19, lonlat)

# for UTM_zone = 18
subset_18 <- subset(data, data$UTM_zone==18)
lonlat <- utm2lonlat(subset_18$Corrected_UTM_X, subset_18$Corrected_UTM_Y, zone = 18, hemisphere = "N", km = TRUE)
lonlat_18 <- cbind(subset_18, lonlat)

# Save
bind <- rbind (lonlat_19,lonlat_18)
bind <- bind[order(bind$Site,decreasing = FALSE),]
bind <- bind[,-1:-2]

colnames(bind)[(ncol(bind)-1):ncol(bind)] <- c("Corrected_Longitude","Corrected_Latitude")
bind$UTM_zone <- data$UTM_zone
write.csv (bind, "lon_lat_UNE_rep.separated.csv")
