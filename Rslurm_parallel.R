#' Run code in parallel using Rslurm
#' @author Lauren White
#' @date May 17, 2019
#' REMEMBER: Make sure extra libraries are detached
#' Helpful tutorial: https://cran.r-project.org/web/packages/rslurm/vignettes/rslurm.html#adding-auxiliary-data-and-functions

source('~/StigmergyDisease/StigmergyLoopFunction.R')
source('~/StigmergyDisease/MoveLoopFunction.R')
source('~/StigmergyDisease/StigmergyFunctions.R')

# add functions from `Stigmergy Functions` file
add_objects<-c("calc.dens", "Cmatrix", "create.prob", "createland", "get.neighbors", "index_val", "infect_direct", "infect_indirect", 
              "longxy", "make.inds", "new.cell", "new.cell.directed", "recover.inds", "Rmatrix", "samp_direction", "sampfun"
)

#specify parameters to run
maxT<-5000
nsim<-100
lsize<-50 
n.initial<-c(125, 250)
inf_prob<-c(0.10, 0.25, 0.50)
rec_rate<-c(0.10, 0.05, 0.01)
dur_scent<- c(1,10) 
initial_load<-c(1,10) 
scent_decay<-c(0.1, 0.5)
inf_decay<-c(0.1, 0.5)

params<-expand.grid(maxT=maxT, nsim=nsim, lsize=lsize, n.initial=n.initial, inf_prob=inf_prob, rec_rate=rec_rate, dur_scent=dur_scent, initial_load=initial_load, scent_decay=scent_decay, inf_decay=inf_decay)
# params<-params[1:8,] #subset of parameters


library(rslurm)
#options for slurm cluster
sopt <- list(time = '120:00:00', share = TRUE, "mail-user"= "lwhite@sesync.org", "mail-type"="ALL", partition="sesync")

sjob1 <- slurm_apply(StigLoop, params, jobname = 'stig_sim',
                    nodes = 24, cpus_per_node = 8, submit = TRUE, add_objects=add_objects, slurm_options=sopt)

sjob2 <- slurm_apply(MoveLoop, params, jobname = 'move_sim',
                     nodes = 24, cpus_per_node = 8, submit = TRUE, add_objects=add_objects, slurm_options=sopt)

# Check on or cancel jobs:
# print_job_status(sjob)
# cancel_slurm(sjob)