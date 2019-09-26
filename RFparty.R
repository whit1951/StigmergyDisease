# Random Forest w/ party package ------------------------------------------

#Load libraries
library(party)

#Load data
merged_data<-read.csv("~/StigmergyDisease/stig_summary.csv")
merged_data<-merged_data[,-1]
merged_data$outbreak<-ifelse(merged_data$max_I>1, 1,0)

set.seed(123)

tic=Sys.time()
fit.cf <- cforest(outbreak ~ n.initial+rec_rate+scent_load+pathogen_load+ scent_decay+ inf_decay +dir_move, data=merged_data, controls=cforest_unbiased(ntree=500))
v <- varimp(fit.cf, conditional= TRUE)
v<-v[order(v)]
write.csv(v, "partyRF_logit500.csv")
dotchart(v[order(v)])
print(difftime(Sys.time(),tic,units="mins"))

merged_data<-subset(merged_data, merged_data$outbreak==1) #analyze only those outbreaks that spread beyond initial individual

tic=Sys.time()
fit.cf1 <- cforest(max_prevalence ~ n.initial+rec_rate+scent_load+pathogen_load+ scent_decay+ inf_decay +dir_move, data=merged_data, controls=cforest_unbiased(ntree=500))
v1 <- varimp(fit.cf1, conditional= TRUE)
v1<-v1[order(v1)]
write.csv(v1, "partyRF_logitprev500.csv")
dotchart(v1[order(v1)])

print(difftime(Sys.time(),tic,units="mins"))

tic=Sys.time()
fit.cf2 <- cforest(duration ~ n.initial+rec_rate+scent_load+pathogen_load+ scent_decay+ inf_decay +dir_move, data=merged_data, controls=cforest_unbiased(ntree=500))
v2 <- varimp(fit.cf2, conditional= TRUE)
v2<-v2[order(v2)]
write.csv(v2, "partyRF_logitdur500.csv")
dotchart(v2[order(v2)])

print(difftime(Sys.time(),tic,units="mins"))

