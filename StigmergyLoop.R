#' Building up ToyModelStigmergy.R to run in a loop, completing `N` simulations per parameter set
#' @author Lauren White
#' @date May 8, 2019

#' Run toy model with stigmergy functions
#' May 6, 2019
#' Lauren White

#load required functions to run simulation
source('~/StigmergyDisease/StigmergyFunctions.R')

#Set up initial conditions/parameters to loop for simulations
lsize<-50 #landscape size
n.initial <- 150 # inital population size
n.offset<-1 #neighborhood of nine cells, including current cell
rowcol.delta <- expand.grid(-n.offset:n.offset,-n.offset:n.offset) #possible moves for given neighborhood size
dist_weights<- get.dist.weight(rowcol.delta)
scent_load<-10 #strength of scent mark upon initial deposition
pathogen_load<-1 #initial pathogen load deposited into environment upon visiting a cell
maxT<-200 #duration of simulation
lxy<-longxy(lsize) #convenience data frame with x, y coordinates for number system of matrices in R
inf_prob<-0.2 #probability of infection per interaction per time step
rec_rate<-0.02
scent_decay<-0.5 #rate at which scent cues decay from the environment (N0*exp(-scent_decay*t))
inf_decay<-0.01 #rate at which infectious agents decay from the environment (N0*exp(-inf_decay*t))
prob_mat<-create.prob(dist_weights) #probability matrix governing movement choices after scent encounter

nsim<-2 #number of simulations to iterate

##Initiate data structure storage across simulations
mydata<-list(inds=NULL, landscape=NULL, N=NULL, movedat=NULL, infdat=NULL)
filename<-paste("lsize", lsize, "n.initial", n.initial, "inf", inf_prob, "rec", rec_rate,"scent_load",scent_load,"pathogen_load",pathogen_load, "scent_decay", scent_decay, "inf_decay", inf_decay, sep="")
summary<-data.frame(lsize=rep(lsize, times=nsim), n.initial=rep(n.initial, times=nsim), inf_prob=rep(inf_prob, times=nsim), rec_rate=rep(rec_rate, times=nsim), scent_load=rep(scent_load, times=nsim), pathogen_load=rep(pathogen_load, times=nsim), scent_decay=rep(scent_decay, times=nsim), inf_decay= rep(inf_decay, times=nsim), duration=NaN, max_I=NaN, max_prevalence=NaN)
infected<-matrix(NaN, nrow=maxT, ncol=nsim)


for(count in 1:nsim){
  #create dataframe of individuals with one infectious individual
  inds<-make.inds(n.initial, lsize, nI=1) 

  # Create scent landscapes for each individual and for infection
  landscape<-createland(lsize=lsize, inds=inds, scent_load=scent_load, infected= FALSE) #create array of landscapes to track scent mark location and strength for each animal
  Num<-calc.dens(lsize=lsize, inds=inds)
  inf_landscape<-Num[[3]]*pathogen_load

  #ignore which animals have been exposed to scent at t= 1, because no prior direction information to inform movement choices

  #record individual characteristics
  mydata$inds[[count]]<-inds

  #Initiate storage for movement and infection data
  movedat<-inds$vec
  infdat<-inds$status
  
  # Object for storing results for single simulation (disease dynamics through time)
  N <-data.frame(S=NaN, I=NaN, R=NaN)
  N[1,] <- c(sum(inds$status=="S"), sum(inds$status=="I"), sum(inds$status=="R"))

  for(t in 2:maxT){
    
    #generate list of possible cells that each animal *could* select
    possible_loc<-get.neighbors(inds[,2:3], mapdim=c(lsize,lsize), rowcol.delta=rowcol.delta, n.offset= n.offset, torus=TRUE, na.val=0)
    
    #sample among possible new cells and udate `inds` df
    if(dir_move==TRUE){
      inds<-new.cell.directed(possible_loc,inds, prob_mat=prob_mat, lsize, scent_load=scent_load)
    }
    if(dir_move==FALSE){
      inds<-new.cell(possible_loc,inds,lsize, dist_weights=dist_weights) #no aversion to scent marks
    }
    movedat<-cbind(movedat, inds$vec)
    
    
    #Disease processes
    #Which individuals become directly infected by other conspecfics?
    # Num<-calc.dens(lsize=lsize, inds=inds)
    # inds<-infect_direct(inds, nS=as.numeric(Num[[2]]), nI=as.numeric(Num[[3]]), transProb=inf_prob, lxy=lxy)
    
    #Which individuals becom indirectly infected by environmental exposure?
    Num<-calc.dens(lsize=lsize, inds=inds) #recalculate S vs. I after direct transmission
    inds<-infect_indirect(inds, nS=as.numeric(Num[[2]]), nI=as.numeric(Num[[3]]), lxy=lxy, inf_landscape=inf_landscape)
    
    #Which animals recover?
    inds<-recover.inds(inds, gamma=rec_rate)
    
    #Update disease status lists
    infdat<-cbind(infdat,inds$status)
    N[t,] <- c(sum(inds$status=="S"), sum(inds$status=="I"), sum(inds$status=="R"))
    
    #Update scent mark landscapes (if time remaining on landscape is greater than one, subtract one time step)
    landscape<-landscape*exp(-scent_decay*1)
    
    for(i in 1:n.initial) #for each individual in simulation
    {
      #array[row, col, layer]
      landscape[inds$y[i], inds$x[i], i]<-scent_load #at new location, deposit scent mark of initial strength
    }
    
    #Exposed to scent mark in this latest move? if sum of scent layers > individual deposition, exposure = Y
    for(j in 1:n.initial)
    {
      scent_sum<-apply(landscape[,,-j], MARGIN=c(1, 2), sum) #sum scent landscapes minus the current individual's profile
      inds$scent_exp[j]<-scent_sum[inds$vec[j]] #reassign scent exposure value (can be greater than one)
        #ifelse(scent_sum[inds$vec[j]]>0, 1, 0) #if scent lingers in the new cell, assign a value of 1, otherwise 0
    }
    
    # Update infection landscape (assume that deposit of infectious agents is additive/cumulative)
    inf_landscape<-inf_landscape*exp(-inf_decay*1)
    inf_landscape<-Num[[3]]*pathogen_load+inf_landscape
    
    if(sum(inds$status=="I")==0){ #If number of infected in individuals--> 0, stop running
      summary$duration[count]<-t-1
      for (i in t:maxT){
        N[i,]<-N[t,] #automatically fill remaining time slots with current N values
      }
      break
    }
  
  }
  mydata$N[[count]]<-N
  mydata$movedat[[count]]<-movedat
  mydata$infdat[[count]]<-infdat
  summary$max_I[count]<-max(N$I)
  summary$max_prevalence[count]<-max(N$I)/n.initial
  infected[,count]<-N$I
}

#Record output
write.csv(summary,file=paste(filename,"summary.csv", sep="_"))
write.csv(infected,file=paste(filename,"infected.csv", sep="_"))
print(filename)
