library(fields)

# Import files
setwd("~/R/Analysis/2_UNE")
fieldData<-read.csv("experimental_design.csv",header=TRUE)

dec.degrees.mat<-as.matrix(cbind(fieldData$Corrected_Longitude, fieldData$Corrected_Latitude)); rownames(dec.degrees.mat)<-fieldData$Sample.ID;
distance.matrix<-rdist.earth(dec.degrees.mat[,1:2],miles=FALSE)

write.csv(distance.matrix,"Autocorrelation/dd.distances.csv",row.names=TRUE)
