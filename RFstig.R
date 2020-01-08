#' Random forest analysis of stigmergy simulations
#' @date May 29, 2019
#' @author Lauren White

#Load libraries and set seed
library(parallel)
library(doParallel)
library(foreach)
library(randomForest)
library(caret)

set.seed(123)
setwd("~/StigmergyDisease")

#Testing
#Detect cores
ncores <- parallel::detectCores()
doParallel::registerDoParallel(ncores)

#Load data
merged_data<-read.csv("~/StigmergyDisease/stig_summary.csv")
merged_data<-merged_data[,-1]

merged_data$outbreak<-ifelse(merged_data$max_I>1, 1, 0) #create an outbreak success column, if successful=1, if unsuccessful=0
merged_data$outbreak<-as.factor(merged_data$outbreak)

samp <- sample(nrow(merged_data), 0.8 * nrow(merged_data))
train <- merged_data[samp, ]
test <- merged_data[-samp, ]


# Outbreak success --------------------------------------------------------

x_train<-train[,c(2:8)] #covariates
y_train<-train$outbreak #response variable, categorical, outbreak (0 or 1)

tic=Sys.time()
# (rf_logit <- foreach(ntree=rep(10, 10), .combine=combine, .multicombine=TRUE,
#               .packages='randomForest') %dopar% {
#                 randomForest(x_train, y_train, ntree=ntree, importance=TRUE)
#               })
(rf_logit<-randomForest(x_train, y_train, ntree=1000, importance=TRUE))
print(difftime(Sys.time(),tic,units="mins"))
plot(rf_logit)

OOB_mat<-table(predict(rf_logit),train$outbreak)

(OOB_error<-(OOB_mat[1,2]+OOB_mat[2,1])/(OOB_mat[1,2]+OOB_mat[2,1]+ OOB_mat[1,1]+OOB_mat[2,2])) #OOB error estimate

pred <- predict(rf_logit, newdata = test)
conf_mat<-table(pred, test$outbreak)
(acc<-(conf_mat[1,1]+conf_mat[2,2])/(conf_mat[1,2]+conf_mat[2,1]+ conf_mat[1,1]+conf_mat[2,2]))
# save(rf_logit, file = "RF_logit.RData")

# rf_fit <- train(x_train, y_train, method = "rf")

varImpPlot(rf_logit, type=1, main= "Outbreak Success")


imp_logit<-rf_logit$importance
write.csv(imp_logit, "imp_logit.csv")


# Maximum prevalence ------------------------------------------------------

merged_data<-subset(merged_data, merged_data$outbreak==1) #analyze only those outbreaks that spread beyond initial individual
merged_data$duration<-merged_data$duration/mean(merged_data$duration) #normalize duration, since not contained between 0 and 1

samp <- sample(nrow(merged_data), 0.8 * nrow(merged_data))
train <- merged_data[samp, ]
test <- merged_data[-samp, ]

x_train<-train[,c(2:8)] #covariates
y_train<-train$max_prevalence

tic=Sys.time()
# (rf_prev <- foreach(ntree=rep(10, 10), .combine=combine, .multicombine=TRUE,
#                     .packages='randomForest') %dopar% {
#                       randomForest(x_train, y_train, ntree=ntree, importance=TRUE)
#                     })
(rf_prev<-randomForest(x_train, y_train, ntree=1000, importance=TRUE))
print(difftime(Sys.time(),tic,units="mins"))
plot(rf_prev)

predicted<-unname(predict(rf_prev, test))
actual <- test$max_prevalence

R2 <- 1 - (sum((actual-predicted)^2, na.rm=TRUE)/sum((actual-mean(actual))^2, na.rm=TRUE))
# save(rf_prev, file = "RF_prev.RData")

varImpPlot(rf_prev, type=1, main= "Prevalence")

imp_prev<-rf_prev$importance
write.csv(imp_prev, "imp_prev.csv")


# Duration ----------------------------------------------------------------

y_train<-train$duration

tic=Sys.time()
# rf_dur <- foreach(ntree=rep(10, 10), .combine=combine, .multicombine=TRUE,
#                    .packages='randomForest') %dopar% {
#                      randomForest(x, y, ntree=ntree, importance=TRUE)
#                    }
(rf_dur<-randomForest(x_train, y_train, ntree=1000, importance=TRUE))
print(difftime(Sys.time(),tic,units="mins"))
# save(rf_dur, file = "RF_dur.RData")

plot(rf_dur)
predicted<-unname(predict(rf_dur, test))
actual <- test$duration

(R2 <- 1 - (sum((actual-predicted)^2, na.rm=TRUE)/sum((actual-mean(actual))^2, na.rm=TRUE)))
# save(rf_prev, file = "RF_prev.RData")

varImpPlot(rf_dur, type=1, main= "Duration")
imp_dur<-rf_dur$importance
write.csv(imp_dur, "imp_dur.csv")


# Test calculation of R2 for combined randomforest ------------------------
# 
# tic=Sys.time()
# (rf_prev<-randomForest(x_train, y_train, ntree=1000, importance=TRUE))
# print(difftime(Sys.time(),tic,units="mins"))
# 
# rf_prev$mse
# rf_prev$rsq*100 #percent variance explained
# 
# (predicted<-unname(predict(rf_prev, train)))
# predicted<-unname(rf_prev$predicted)
# actual <- train$max_prevalence
# # predicted <- unname(predict(rfPar, dat))
# 
# R2 <- 1 - (sum((actual-predicted)^2, na.rm=TRUE)/sum((actual-mean(actual))^2, na.rm=TRUE))
