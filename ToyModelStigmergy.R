#' Run toy model with stigmergy functions
#' May 6, 2019
#' Lauren White

#load functions
source('~/StigmergyDisease/StigmergyFunctions.R')

#Set up initial conditions for all simulations
lsize<-5
n.initial <- 5 # inital population size
n.offset<-1 #neighborhood of nine cells, including current cell
rowcol.delta <- expand.grid(-n.offset:n.offset,-n.offset:n.offset) #possible moves given neighborhoodsize
dur_scent<-10 #how long scent marks last in the environment
T<-50 #duration of simulation

inds<-make.inds(lsize, n.initial, nI=1) #create dataframe of individuals with one infectious individual
landscape<-createland(lsize=lsize, N=n.initial, inds=inds, dur_scent=dur_scent) #create array of landscapes to track scent mark location and strength for each animal

##Record initial locations
startloc<-inds$vec
newloc.vec<-startloc

##Convert vector notation to matrix coordinates
newloc<-as.data.frame(Rmatrix(newloc.vec-1,lsize))

#Initiate longterm storage for movement and infection data
movedat<-startloc
infdat<-inds$status

# Object for storing results for single simulation (disease dynamics through time)
N <-data.frame(S=NaN, I=NaN, R=NaN)
N[1,] <- c(sum(inds$status=="S"), sum(inds$status=="I"), sum(inds$status=="R"))

for(t in 1:T){
#generate list of possible cells that each animal *could* select
possible_loc<-get.neighbors(newloc, mapdim=c(lsize,lsize), rowcol.delta=rowcol.delta, n.offset= n.offset, torus=TRUE, na.val=0)

#sample along the rows of the possible locations and update `newloc` vector
newloc.vec<-apply(possible_loc, 1, sampfun) 
newloc<-as.data.frame(Rmatrix(newloc.vec-1,lsize))

##Update location lists
inds$xloc<-newloc[,1] #col
inds$yloc<-newloc[,2] #row
movedat<-cbind(movedat,newloc.vec)
inds$vec<-newloc.vec

#Update scent marks
landscape[which(landscape>1)]<-landscape[which(landscape>1)]-1
for(i in 1:n.initial)
{
  #array[row, col, layer]
  landscape[ newloc$row[i], newloc$col[i], i]<-dur_scent
}


}
