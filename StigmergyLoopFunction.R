#' StigmergyLoopFunction.R
#' Setting up StigmergyLoop.R to run as a function as a precursor to running scripts in parallel
#' @author Lauren White
#' @date May 8, 2019

#load dependent functions
source('~/StigmergyDisease/StigmergyFunctions.R')


#'@param maxT - maximum duration of simulation length
#'@param nsim - number of simulations to be repeated per parameter set
#'@param lsize - landscape size (length of one dimension)
#'@param n.initial - number of simulated agents on the landscape
#'@param rec_rate - recovery rate
#'@param scent_load - how long scent marks last in the environment
#'@param pathogen_load - initial pathogen load deposited into environment upon visiting a cell
#'@param scent_decay - rate at which scent cues decay from the environment (N0*exp(-scent_decay*t))
#'@param inf_decay - rate at which scent cues decay from the environment (N0*exp(-scent_decay*t))
#'@param dir_move - boolean (TRUE OR FALSE)- is movement directed in response to scent cues or random walk?
StigLoop<-function(maxT, nsim, lsize, n.initial, rec_rate, scent_load, pathogen_load, scent_decay, inf_decay, dir_move){
  # Reproducible random numbers in parallel 
  # https://bookdown.org/rdpeng/rprogdatascience/parallel-computation.html#generating-random-numbers-1
  RNGkind("L'Ecuyer-CMRG")
  set.seed(2718)
  
  #Set up initial conditions/parameters to loop for simulations
  n.offset<-1 #neighborhood of nine cells, including current cell
  rowcol.delta <- expand.grid(-n.offset:n.offset,-n.offset:n.offset) #possible moves for given neighborhood size
  dist_weights<- get.dist.weight(rowcol.delta)
  prob_mat<-create.prob(dist_weights) #probability matrix governing movement choices after scent encounter
  lxy<-longxy(lsize) #convenience data frame with x, y coordinates for number system of matrices in R
  
  ##Initiate data structure storage across simulations
  mydata<-list(inds=NULL, landscape=NULL, N=NULL, movedat=NULL, infdat=NULL)
  filename<-paste("lsize", lsize, "n.initial", n.initial, "rec", rec_rate,"scent_load",scent_load,"pathogen_load",pathogen_load, "scent_decay", scent_decay, "inf_decay", inf_decay, "dir_move", dir_move, sep="")
  summary<-data.frame(lsize=rep(lsize, times=nsim), n.initial=rep(n.initial, times=nsim), rec_rate=rep(rec_rate, times=nsim), scent_load=rep(scent_load, times=nsim), pathogen_load=rep(pathogen_load, times=nsim), scent_decay=rep(scent_decay, times=nsim), inf_decay= rep(inf_decay, times=nsim), dir_move=rep(dir_move, times=nsim), duration=NaN, max_I=NaN, max_prevalence=NaN)
  infected<-matrix(NaN, nrow=maxT, ncol=nsim)
  
  
  for(count in 1:nsim){
    #create dataframe of individuals with one infectious individual
    inds<-make.inds(n.initial, lsize, nI=1) 
    
    # Create scent landscapes for each individual and for infection
    landscape<-createland(lsize=lsize, inds=inds, scent_load=scent_load, infected= FALSE) #create array of landscapes to track scent mark location and strength for each animal
    Num<-calc.dens(lsize=lsize, inds=inds)
    inf_landscape<-Num[[3]]*pathogen_load
    
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
      #
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
      
      for(i in 1:n.initial)
      {
        #array[row, col, layer]
        landscape[inds$y[i], inds$x[i], i]<-scent_load #at new location, deposit scent mark of initial strength
      }
      
      #Exposed to scent mark in this latest move? if sum of scent layers > individual deposition, exposure = Y
      for(j in 1:n.initial)
      {
        scent_sum<-apply(landscape[,,-j], MARGIN=c(1, 2), sum) #sum scent landscapes minus the current individual's profile
        inds$scent_exp[j]<-scent_sum[inds$vec[j]] #reassign scent exposure value (can be greater than one)
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
}


#Test

# tic=Sys.time()
# StigLoop(maxT=200, nsim=1, lsize=50, n.initial=50, rec_rate=0.01, scent_load=0.1, pathogen_load=0.5, scent_decay=0.01, inf_decay=0.01, dir_move=TRUE)
# print(difftime(Sys.time(),tic,units="mins"))
