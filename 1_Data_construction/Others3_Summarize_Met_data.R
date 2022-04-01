setwd("~/R/Analysis/2_UNE/Met_data")

# aggregate
data <- cbind(logger_daily_mean$below_0.5C,logger_daily_mean$freeze_to_thaw_transition_0.5C)
Result <- aggregate(data,by=list(logger_daily_mean$dfe,logger_daily_mean$site_code,logger_daily_mean$Year),FUN = sum,na.rm=T)
colnames(Result)<-c("DFE","site_code","Year","below_0.5C","freeze_to_thaw_transition_0.5C")

data_moist <- cbind(logger_daily_mean_moist$below_10perc,logger_daily_mean_moist$wet_to_dry_transition_10perc)
Result_moist <- aggregate(data_moist,by=list(logger_daily_mean_moist$dfe,logger_daily_mean_moist$site_code,logger_daily_mean_moist$Year),FUN = sum,na.rm=T)
colnames(Result_moist)<-c("DFE","site_code","Year","below_10perc","wet_to_dry_transition_10perc")

# Save
bind <- cbind(Result,Result_moist)
bind<-bind[,-6:-8]
write.csv(bind,"Met.result.csv")
