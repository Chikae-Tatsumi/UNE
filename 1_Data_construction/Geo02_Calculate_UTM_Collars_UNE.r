# Import files
setwd("~/R/Analysis/2_UNE/Autocorrelation")
UTM <- read.csv(file = "UTM_UNE.csv",header=T)

AA01 <- subset(UTM, UTM$Site=="AA01")
BH02 <- subset(UTM, UTM$Site=="BH02")
BM03 <- subset(UTM, UTM$Site=="BM03")
HF04 <- subset(UTM, UTM$Site=="HF04")
HF05 <- subset(UTM, UTM$Site=="HF05")
HF06 <- subset(UTM, UTM$Site=="HF06")
HW07 <- subset(UTM, UTM$Site=="HW07")
SW08 <- subset(UTM, UTM$Site=="SW08")

sample.list <- list()
sample.list[[1]] <- AA01
sample.list[[2]] <- BH02
sample.list[[3]] <- BM03
sample.list[[4]] <- HF04
sample.list[[5]] <- HF05
sample.list[[6]] <- HF06
sample.list[[7]] <- HW07
sample.list[[8]] <- SW08

NAME <- c("AA01","BH02","BM03","HF04","HF05","HF06","HW07","SW08")

Result <- NULL
for (i in 1:8){
data <- sample.list[[i]]
lm <- lm(UTM_Y~UTM_X, data=data)
slope <- abs(lm$coefficients[2])
if (data$UTM_X[1]>data$UTM_X[9] & data$UTM_Y[1]>data$UTM_Y[9]){
    data$Corrected_UTM_X <- data$UTM_X + 0.005*slope / sqrt(slope^2+1)
    data$Corrected_UTM_Y <- data$UTM_Y - 0.005 / sqrt(slope^2+1)
    for(k in 1:10){
    if(k%%2 == 0){
    data$Corrected_UTM_X[k] <- data$UTM_X[k] - 0.005*slope / sqrt(slope^2+1)
    data$Corrected_UTM_Y[k] <- data$UTM_Y[k] + 0.005 / sqrt(slope^2+1)
}}} else if (data$UTM_X[1]<data$UTM_X[9] & data$UTM_Y[1]>data$UTM_Y[9]){
    data$Corrected_UTM_X <- data$UTM_X + 0.005*slope / sqrt(slope^2+1)
    data$Corrected_UTM_Y <- data$UTM_Y + 0.005 / sqrt(slope^2+1)
    for(k in 1:10){
    if(k%%2 == 0){
    data$Corrected_UTM_X[k] <- data$UTM_X[k] - 0.005*slope / sqrt(slope^2+1)
    data$Corrected_UTM_Y[k] <- data$UTM_Y[k] - 0.005 / sqrt(slope^2+1)
}}} else if (data$UTM_X[1]<data$UTM_X[9] & data$UTM_Y[1]<data$UTM_Y[9]){
    data$Corrected_UTM_X <- data$UTM_X - 0.005*slope / sqrt(slope^2+1)
    data$Corrected_UTM_Y <- data$UTM_Y + 0.005 / sqrt(slope^2+1)
    for(k in 1:10){
    if(k%%2 == 0){
    data$Corrected_UTM_X[k] <- data$UTM_X[k] + 0.005*slope / sqrt(slope^2+1)
    data$Corrected_UTM_Y[k] <- data$UTM_Y[k] - 0.005 / sqrt(slope^2+1)
}}} else if (data$UTM_X[1]>data$UTM_X[9] & data$UTM_Y[1]<data$UTM_Y[9]){
    data$Corrected_UTM_X <- data$UTM_X - 0.005*slope / sqrt(slope^2+1)
    data$Corrected_UTM_Y <- data$UTM_Y - 0.005 / sqrt(slope^2+1)
    for(k in 1:10){
    if(k%%2 == 0){
    data$Corrected_UTM_X[k] <- data$UTM_X[k] + 0.005*slope / sqrt(slope^2+1)
    data$Corrected_UTM_Y[k] <- data$UTM_Y[k] + 0.005 / sqrt(slope^2+1)
}}} else {}
Result <- rbind(Result, data)}

Result$UTM_zone <- UTM$UTM_zone
write.csv(Result, "UTM_UNE_rep.separated.csv")
