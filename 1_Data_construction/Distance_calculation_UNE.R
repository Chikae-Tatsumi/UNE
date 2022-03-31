#read in data
setwd("~/R/Analysis/2_UNE")
fieldData<-read.csv("experimental_design.csv",header=TRUE)

require(fields)
dec.degrees.mat<-as.matrix(cbind(fieldData$Longitude, fieldData$Latitude)); rownames(dec.degrees.mat)<-fieldData$Sample.ID;
distance.matrix<-rdist.earth(dec.degrees.mat[,1:2],miles=FALSE)

dir.create("Autocorrelation")
write.csv(distance.matrix,"Autocorrelation/dd.distances.csv",row.names=TRUE)
