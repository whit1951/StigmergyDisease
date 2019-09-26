#' Run code in parallel using Rslurm
#' @author Lauren White
#' @date May 17, 2019
#' REMEMBER: Make sure extra libraries are detached
#' Helpful tutorial: https://cran.r-project.org/web/packages/rslurm/vignettes/rslurm.html#adding-auxiliary-data-and-functions

source('~/StigmergyDisease/StigmergyLoopFunction.R')
source('~/StigmergyDisease/StigmergyFunctions.R')
# source('~/StigmergyDisease/StigmergyLoopRepFunction.R')


# add functions from `Stigmergy Functions` file
add_objects<-c("calc.dens", "Cmatrix", "create.prob", "createland", "get.dist.weight", "get.neighbors", "infect_direct", "infect_indirect", 
              "longxy", "make.inds", "new.cell", "new.cell.directed", "recover.inds", "Rmatrix"
)

#specify parameters to run
maxT<-5000
nsim<-100
lsize<-50 
n.initial<-c(50, 100)
rec_rate<-c(0.10, 0.05, 0.01)
scent_load<- c(0.5, 1, 10) 
pathogen_load<-c(0.5, 1, 10) 
scent_decay<-c(0.01, 0.1, 1)
inf_decay<-c(0.01, 0.1, 1)
dir_move<-c(TRUE, FALSE)

params<-expand.grid(maxT=maxT, nsim=nsim, lsize=lsize, n.initial=n.initial, rec_rate=rec_rate, scent_load=scent_load, pathogen_load=pathogen_load, scent_decay=scent_decay, inf_decay=inf_decay, dir_move=dir_move)
# params<-params[c(14, 16, 18, 32, 34, 36, 50,  52,  54,  68,  70,  72,  86,  88,  90, 104, 106, 108, 122, 124, 126, 140, 142, 144, 158, 160,
#                   162, 176, 178, 180, 194, 196, 198, 212, 214, 216, 230, 232, 234, 248, 250, 252, 266, 268, 270, 284, 286, 288),]
# 
# nreps<-10
# params_reps<-data.frame()
# params_reps<-params[rep(seq_len(nrow(params)), each=nreps),]
# params_reps$nsim<-rep(10, times=nrow(params_reps))
# params_reps$nrep<-rep(1:10, times=nrow(params_reps)/nreps)

library(rslurm)
#options for slurm cluster
sopt <- list(time = '192:00:00', share = TRUE, "mail-user"= "lwhite@sesync.org", "mail-type"="ALL", partition="sesync")
# 
sjob1 <- slurm_apply(StigLoop, params, jobname = 'stig_sim_2019_07_18',
                    nodes = 24, cpus_per_node = 8, submit = TRUE, add_objects=add_objects, slurm_options=sopt)

# sopt2 <- list(time = '48:00:00', share = TRUE, "mail-user"= "lwhite@sesync.org", "mail-type"="ALL", partition="sesync")
# 
# sjob2 <- slurm_apply(StigLoopRep, params_reps, jobname = 'stig_rep',
#                      nodes = 24, cpus_per_node = 8, submit = TRUE, add_objects=add_objects, slurm_options=sopt2)
# 
# Check on or cancel jobs:
# print_job_status(sjob)
# cancel_slurm(sjob)

#specify parameters to run
maxT<-5000
nsim<-100
lsize<-50 
n.initial<-150
rec_rate<-c(0.10, 0.05, 0.01)
scent_load<- c(0.5, 1, 10) 
pathogen_load<-c(0.5, 1, 10) 
scent_decay<-c(0.01, 0.1, 1)
inf_decay<-c(0.01, 0.1, 1)
dir_move<-c(TRUE, FALSE)

params<-expand.grid(maxT=maxT, nsim=nsim, lsize=lsize, n.initial=n.initial, rec_rate=rec_rate, scent_load=scent_load, pathogen_load=pathogen_load, scent_decay=scent_decay, inf_decay=inf_decay, dir_move=dir_move)

sopt <- list(time = '192:00:00', share = TRUE, "mail-user"= "lwhite@sesync.org", "mail-type"="ALL", partition="sesync")
# 
sjob2 <- slurm_apply(StigLoop, params, jobname = 'stig_sim_2019_09_20',
                     nodes = 24, cpus_per_node = 8, submit = TRUE, add_objects=add_objects, slurm_options=sopt)

