#' Run toy model with stigmergy functions
#' May 6, 2019
#' Lauren White

#load required functions to run simulation
source('~/StigmergyDisease/StigmergyFunctions.R')

#Set up initial conditions for all simulations
lsize<-50
n.initial <- 20 # inital population size
n.offset<-1 #neighborhood of nine cells, including current cell
rowcol.delta <- expand.grid(-n.offset:n.offset,-n.offset:n.offset) #possible moves for given neighborhood size
dur_scent<-50 #how long scent marks last in the environment
initial_load<-1 #initial pathogen load deposited into environment upon visiting a cell
T<-100 #duration of simulation
lxy<-longxy(lsize) #convenience data frame with x, y coordinates for number system of matrices in R
inf_prob<-0.2 #probability of infection per interaction per time step
rec_rate<-0.02
scent_decay<-0.5 #rate at which scent cues decay from the environment (N0*exp(-scent_decay*t))
inf_decay<-0.5 #rate at which infectious agents decay from the environment (N0*exp(-inf_decay*t))
prob_mat<-create.prob() #probability matrix governing movement choices after scent encounter

inds<-make.inds(n.initial, lsize, nI=1) #create dataframe of individuals with one infectious individual

# Create scent landscapes for each individual and for infection
landscape<-createland(lsize=lsize, inds=inds, dur_scent=dur_scent, infected= FALSE) #create array of landscapes to track scent mark location and strength for each animal
#inf_landscape<-createland(lsize=lsize, N=n.initial, inds=inds, dur_scent=dur_scent, infected=TRUE) #create a separate landscape to keep track of infected cells
#inf_landscape<- matrix(0, nrow= lsize, ncol= lsize) #y, x, one layer for each individual
Num<-calc.dens(lsize=lsize, inds=inds)
inf_landscape<-Num[[3]]*initial_load

#update which animals have been exposed to scent at t= 1
scented <-which(inds$vec %in% which(Num[[1]]>1))
inds$scent_exp[scented]<-1

#Initiate longterm storage for movement and infection data
movedat<-inds$vec
infdat<-inds$status

# Object for storing results for single simulation (disease dynamics through time)
N <-data.frame(S=NaN, I=NaN, R=NaN)
N[1,] <- c(sum(inds$status=="S"), sum(inds$status=="I"), sum(inds$status=="R"))

for(t in 2:T){
  
  
#generate list of possible cells that each animal *could* select
possible_loc<-get.neighbors(inds[,2:3], mapdim=c(lsize,lsize), rowcol.delta=rowcol.delta, n.offset= n.offset, torus=TRUE, na.val=0)

#sample among possible new cells and udate `inds` df
#newloc.vec<-apply(possible_loc, 1, sampfun) 
inds<-new.cell(possible_loc,inds)
movedat<-cbind(movedat,inds$vec)


#Disease processes
#Which individuals become directly infected by other conspecfics?
Num<-calc.dens(lsize=lsize, inds=inds)
inds<-infect_direct(inds, nS=as.numeric(Num[[2]]), nI=as.numeric(Num[[3]]), transProb=inf_prob, lxy=lxy)

#Which individuals becom indirectly infected by environmental exposure?
Num<-calc.dens(lsize=lsize, inds=inds) #recalculate S vs. I after direct transmission
inds<-infect_indirect(inds, nS=as.numeric(Num[[2]]), nI=as.numeric(Num[[3]]), transProb=inf_prob, lxy=lxy, inf_landscape=inf_landscape)

#Which animals recover?
inds<-recover.inds(inds, gamma=rec_rate)

#Update disease status lists
infdat<-cbind(infdat,inds$status)
N[t,] <- c(sum(inds$status=="S"), sum(inds$status=="I"), sum(inds$status=="R"))

#Update scent mark landscapes (if time remaining on landscape is greater than one, subtract one time step)
#landscape[which(landscape>1)]<-landscape[which(landscape>1)]-1
landscape[which(landscape>1)]<-landscape[which(landscape>1)]*exp(-scent_decay*1)

for(i in 1:n.initial)
{
  #array[row, col, layer]
  landscape[inds$y[i], inds$x[i], i]<-dur_scent #at new location, deposit scent mark of initial strength
}

#Exposed to scent mark in this latest move? if sum of scent layers > individual deposition, exposure = Y
#scent_sum<-apply(landscape, MARGIN=c(1, 2), sum)
for(j in 1:n.initial)
{
  scent_sum<-apply(landscape[,,-j], MARGIN=c(1, 2), sum) #sum scent landscapes minus the current individual's profile
  inds$scent_exp[j]<-ifelse(scent_sum[inds$vec[j]]>0, 1, 0) #if scent lingers in the new cell, assign a value of 1, otherwise 0
}

# Update infection landscape (assume that deposit of infectious agents is additive/cumulative)
#inf_landscape[which(inf_landscape>1)]<-inf_landscape[which(inf_landscape>1)]-1
inf_landscape[which(inf_landscape>1)]<-inf_landscape[which(inf_landscape>1)]*exp(-inf_decay*1)
inf_landscape<-Num[[3]]*initial_load+inf_landscape

if(sum(inds$status=="I")==0){ #If number of infected in individuals--> 0, stop running
  #summary$duration[count]<-t-1
  for (i in t:T){
    N[i,]<-N[t,] #automatically fill remaining time slots with current N values
  }
  break
}

}

