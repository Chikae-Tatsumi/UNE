# Import files
setwd("~/R/Analysis/2_UNE/qPCR")
result <- read.csv(file = "copies.ul.csv",header=T)
soil <- read.csv(file = "soil data.csv",header=T)

# Add info
sample.dilt <- 20 # dilution factor
final.nfw <- 50 # ul, the final extraction volume

# Calculation
conc <- result*sample.dilt  # copies/ul
conc.extraction <- conc*final.nfw  # copies in extraction
soil.dry.weight <- soil$weight_gram*(1-(soil$water_content_percent/100)) #g soil

data <- NULL
for (i in 1:6){
copies.gsoil <- conc.extraction[,i]/soil.dry.weight # copies/g soil
log.copies.gsoil <- log10(copies.gsoil) # log

copies.ngDNA <- conc[,i]/(soil$DNA_ngperml/1000) # copies/ng DNA
log.copies.ngDNA <- log10(copies.ngDNA) # log

bind <- cbind(copies.gsoil,copies.ngDNA)
name <- colnames(result)[i]
colnames(bind) <- c(paste(name,colnames(bind)[1], sep="."),paste(name,colnames(bind)[2], sep="."))
data <- cbind(data,bind)}

# Save
write.csv(data, "qPCR.calc.result.UNE.csv")