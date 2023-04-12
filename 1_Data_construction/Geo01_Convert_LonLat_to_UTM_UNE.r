# https://rdrr.io/cran/PBSmapping/man/convUL.html
library(PBSmapping)

# Import files
setwd("~/R/Analysis/2_UNE/Autocorrelation")
data <- read.csv("lon_lat_UNE.csv",header=T)

# for UTM_zone = 19
subset_19 <- subset(data, data$UTM_zone==19)
xydata <- cbind(subset_19$Longitude, subset_19$Latitude)
colnames(xydata) <- c("X","Y")
xydata <- data.frame(xydata)
attr(xydata,"projection")="LL"

UTM <- convUL (xydata, km=TRUE, southern=NULL)
UTM_19 <- cbind(subset_19, UTM)

# for UTM_zone = 18
subset_18 <- subset(data, data$UTM_zone==18)
xydata <- cbind(subset_18$Longitude, subset_18$Latitude)
colnames(xydata) <- c("X","Y")
xydata <- data.frame(xydata)
attr(xydata,"projection")="LL"

UTM <- convUL (xydata, km=TRUE, southern=NULL)
UTM_18 <- cbind(subset_18, UTM)

# Save
bind <- rbind (UTM_19,UTM_18)
bind <- bind[order(bind$Site,decreasing = FALSE),]
colnames(bind)[(ncol(bind)-1):ncol(bind)] <- c("UTM_X","UTM_Y")
bind$UTM_zone <- data$UTM_zone
write.csv (bind, "UTM_UNE.csv")
