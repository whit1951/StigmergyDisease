#' Random forest analysis of stigmergy simulations
#' @date May 29, 2019
#' @author Lauren White

#Load libraries
library(parallel)
library(doParallel)
library(foreach)
library(randomForest)

#Testing
#Detect cores
ncores <- parallel::detectCores()
doParallel::registerDoParallel(ncores)


#Load data
merged_data<-read.csv("~/StigmergyDisease/stig_summary.csv")
merged_data<-merged_data[,-1]
# merged_data<-subset(merged_data, rec_rate==0.01)
# merged_data<-merged_data[which(merged_data$n.initial==125),]
merged_data$outbreak<-ifelse(merged_data$max_I>1, 1,0)


x<-merged_data[,c(2:7)] #covariates
y<-merged_data$outbreak #response variable, outbreak (0 or 1)


tic=Sys.time()
set.seed(123)
rf_logit <- foreach(ntree=rep(100, 100), .combine=combine, .multicombine=TRUE,
              .packages='randomForest') %dopar% {
                randomForest(x, y, ntree=ntree, importance=TRUE)
              }
print(difftime(Sys.time(),tic,units="mins"))
# save(rf_logit, file = "RF_logit.RData")

#varImpPlot(rf, type=2, main= "Duration")

setwd("~/StigmergyDisease")
imp_logit<-rf_logit$importance
write.csv(imp_logit, "imp_logit.csv")


merged_data<-subset(merged_data, merged_data$outbreak==1) #analyze only those outbreaks that spread beyond initial individual
x<-merged_data[,c(2:7)] #covariates

tic=Sys.time()
set.seed(123)
# y<-merged_data$max_I/mean(merged_data$maxI) #response variable, maxI
y<-merged_data$max_prevalence
rf_prev <- foreach(ntree=rep(100, 100), .combine=combine, .multicombine=TRUE,
                    .packages='randomForest') %dopar% {
                      randomForest(x, y, ntree=ntree, importance=TRUE)
                    }
print(difftime(Sys.time(),tic,units="mins"))
# save(rf_prev, file = "RF_prev.RData")

imp_prev<-rf_prev$importance
write.csv(imp_prev, "imp_prev.csv")


tic=Sys.time()
set.seed(123)
y<-merged_data$duration/mean(merged_data$duration) #response variable, duration
rf_dur <- foreach(ntree=rep(100, 100), .combine=combine, .multicombine=TRUE,
                   .packages='randomForest') %dopar% {
                     randomForest(x, y, ntree=ntree, importance=TRUE)
                   }
print(difftime(Sys.time(),tic,units="mins"))
# save(rf_dur, file = "RF_dur.RData")

imp_dur<-rf_dur$importance
write.csv(imp_dur, "imp_dur.csv")


