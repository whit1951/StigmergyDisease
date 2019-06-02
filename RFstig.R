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
merged_data<-read.csv("stig_summary.csv")
merged_data<-merged_data[,-1]
merged_data<-merged_data[which(merged_data$n.initial==125),]
merged_data$outbreak<-ifelse(merged_data$max_I>1, 1,0)


x<-merged_data[,c(2:8)] #covariates
y<-merged_data[,15] #response variable, outbreak (0 or 1)


tic=Sys.time()
set.seed(123)
rf_logit <- foreach(ntree=rep(100, 100), .combine=combine, .multicombine=TRUE,
              .packages='randomForest') %dopar% {
                randomForest(x, y, ntree=ntree, importance=TRUE)
              }
print(difftime(Sys.time(),tic,units="mins"))
# save(rf_logit, file = "RF_logit.RData")

#varImpPlot(rf, type=2, main= "Duration")

imp_logit<-rf_logit$importance
write.csv(imp_logit, "imp_logit.csv")


tic=Sys.time()
set.seed(123)
y<-merged_data[,11]/mean(merged_data[,11]) #response variable, maxI
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
y<-merged_data[,9]/mean(merged_data[,9]) #response variable, duration
rf_dur <- foreach(ntree=rep(100, 100), .combine=combine, .multicombine=TRUE,
                   .packages='randomForest') %dopar% {
                     randomForest(x, y, ntree=ntree, importance=TRUE)
                   }
print(difftime(Sys.time(),tic,units="mins"))
# save(rf_dur, file = "RF_dur.RData")

imp_dur<-rf_dur$importance
write.csv(imp_dur, "imp_dur.csv")

# rm(list=ls(all=T)) #clear workspace