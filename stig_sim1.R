#Run Stigmergy_loop_function in parallel via clusters
#May 15, 2019
#Lauren White

library(parallel)
library(foreach)
library(doParallel)
library(SDMTools)

#load functions
source('~/StigmergyDisease/StigmergyFunctions.R')
source('~/StigmergyDisease/StigmergyLoopFunction.R')

#specify parameters to run
maxT<-1000
nsim<-100
lsize<-100 #dimensions of landscape= 2^k+1
n.initial<-c(25, 50)
inf_prob<-0.1
rec_rate<-0.01
dur_scent<- c(1,10) 
initial_load<-c(1,10) 
scent_decay<-c(0.1, 0.5)
inf_decay<-c(0.1, 0.5)


params<-expand.grid(maxT=maxT, nsim=nsim, lsize=lsize, n.initial=n.initial, inf_prob=inf_prob, rec_rate=rec_rate, dur_scent=dur_scent, initial_load=initial_load, scent_decay=scent_decay, inf_decay=inf_decay)
#params<-params[1:4,]
#params<-params[-which(params$beta2==5 & params$beta3==0),]
#summary_data<-list(NA)

## Apply the declared function in parallel
#StigLoop(maxT=1000, nsim=1, lsize=100, n.initial=50, inf_prob=0.5, rec_rate=0.01, dur_scent=10, initial_load=10, scent_decay=0.1, inf_decay=0.1)

ncores <- parallel::detectCores()
print(ncores)
doParallel::registerDoParallel(ncores)
tic=Sys.time()
foreach(n = 1:nrow(params)) %dopar% StigLoop(params[n,]$maxT, params[n,]$nsim, params[n,]$lsize, params[n,]$n.initial, params[n,]$inf_prob, params[n,]$rec_rate, params[n,]$dur_scent, params[n,]$initial_load, params[n,]$scent_decay, params[n,]$inf_decay)
print(difftime(Sys.time(),tic,units="mins"))

#save(summary_data, file="move_sim1.RData")
#save(params, file="parameters_sim1.RData")

rm(list=ls(all=T)) #clear workspace