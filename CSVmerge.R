#'May 16, 2019
#'For analysis of stigmergy simulations
#'Merge read ins from data files

library(doBy)
library(ggplot2)
library(viridis)
library(ggthemes)
library(scales)
library(lsr)

##Read in file names from directory--> n.initial=150
filenames <- list.files(path = "~/StigmergyDisease/_rslurm_stig_sim_2019_09_20")
summaries<-filenames[grep("summary", filenames)] #summary data
infecteds<-filenames[grep("infected", filenames)] #I data

##Check complete list
#specify parameters to run
maxT<-5000
nsim<-100
lsize<-50 
n.initial<-c(150)
rec_rate<-c(0.10, 0.05, 0.01)
scent_load<- c(0.5, 1, 10) 
pathogen_load<-c(0.5, 1, 10) 
scent_decay<-c(0.01, 0.1, 1)
inf_decay<-c(0.01, 0.1, 1)
dir_move<-c(TRUE, FALSE)

complete_list<-vector(mode="character")
params<-expand.grid(maxT=maxT, nsim=nsim, lsize=lsize, n.initial=n.initial, rec_rate=rec_rate, scent_load=scent_load, pathogen_load=pathogen_load, scent_decay=scent_decay, inf_decay=inf_decay, dir_move=dir_move)
for(i in 1:nrow(params)){
complete_list[i]<-paste("lsize", params$lsize[i], "n.initial", params$n.initial[i], "rec", params$rec_rate[i],"scent_load", params$scent_load[i],"pathogen_load",params$pathogen_load[i], "scent_decay", params$scent_decay[i], "inf_decay", params$inf_decay[i], "dir_move", params$dir_move[i], "_summary.csv", sep="")
}

miss_sim<-logical(length(complete_list))
for(i in 1:length(complete_list)){
if(complete_list[i] %in% summaries)
  miss_sim[i]<-TRUE
}
which(miss_sim==FALSE)  

##Make a merged dataset from all available .csv files
setwd("~/StigmergyDisease/_rslurm_stig_sim_2019_09_20")
filenames <- list.files(path = "~/StigmergyDisease/_rslurm_stig_sim_2019_09_20")
summaries<-filenames[grep("summary", filenames)] #summary data
infecteds<-filenames[grep("infected", filenames)] #I data
merged_data150<- do.call("rbind", lapply(summaries, read.csv, header = TRUE))
merged_data150$X<-NULL

##Read in file names from directory--> n.initial=50 & 100
setwd("~/StigmergyDisease/_rslurm_stig_sim_2019_07_18")
filenames <- list.files(path = "~/StigmergyDisease/_rslurm_stig_sim_2019_07_18")
summaries<-filenames[grep("summary", filenames)] #summary data
infecteds<-filenames[grep("infected", filenames)] #I data

merged_data<- do.call("rbind", lapply(summaries, read.csv, header = TRUE))
merged_data$X<-NULL
merged_data[which(is.na(merged_data$duration)),]

merged_data<-rbind(merged_data, merged_data150)



#Create character strings for landscape structure (Hurst exponent and proportion available habitat)
merged_data$PL<-NA #pathogen load
merged_data$SL<-NA #scent load
merged_data$PL[which(merged_data$pathogen_load==0.5)]<-"lowPL"
merged_data$PL[which(merged_data$pathogen_load==1)]<-"medPL"
merged_data$PL[which(merged_data$pathogen_load==10)]<-"highPL"
merged_data$SL[which(merged_data$scent_load==0.5)]<-"lowSL"
merged_data$SL[which(merged_data$scent_load==1)]<-"medSL"
merged_data$SL[which(merged_data$scent_load==10)]<-"highSL"
merged_data$decay <- paste(merged_data$scent_decay,merged_data$inf_decay)
setwd("~/StigmergyDisease/")
# write.csv(merged_data, "stig_summary.csv")


merged_data<-read.csv("stig_summary.csv")
#Order factors
# merged_data$betas<- factor(merged_data$betas, levels = c("0 0 0", "0 0 -1", "0 1 -1", "0 1 -0.5", "0 2 -1", "0 2 -0.5", "3 0 0", "3 0 -1", "3 1 -1", "3 1 -0.5", "3 2 -1", "3 2 -0.5", "6 0 0", "6 0 -1", "6 1 -1", "6 1 -0.5", "6 2 -1", "6 2 -0.5"))
# merged_data$lstructure<- factor(merged_data$lstructure, levels = c("L/HP", "L/MP", "L/LP", "M/HP", "M/MP", "M/LP", "H/HP", "H/MP", "H/LP"))


##Make some sub data sets by density and recovery rate
size50<-merged_data[which(merged_data$n.initial==50),]
size100<-merged_data[which(merged_data$n.initial==100),]
size150<-merged_data[which(merged_data$n.initial==150),]

sub_dT_r0.01<-merged_data[which(merged_data$dir_move==TRUE & merged_data$rec_rate==0.01),]
sub_dT_r0.05<-merged_data[which(merged_data$dir_move==TRUE & merged_data$rec_rate==0.05),]
sub_dT_r0.1<-merged_data[which(merged_data$dir_move==TRUE & merged_data$rec_rate==0.1),]

sub_n50_dT_r0.01<-merged_data[which(merged_data$dir_move==TRUE & merged_data$rec_rate==0.01 & merged_data$n.initial==50),]
sub_n50_dT_r0.05<-merged_data[which(merged_data$dir_move==TRUE & merged_data$rec_rate==0.05 & merged_data$n.initial==50),]
sub_n50_dT_r0.1<-merged_data[which(merged_data$dir_move==TRUE & merged_data$rec_rate==0.1 & merged_data$n.initial==50),]

sub_n100_dT_r0.01<-merged_data[which(merged_data$dir_move==TRUE & merged_data$rec_rate==0.01 & merged_data$n.initial==100),]
sub_n100_dT_r0.05<-merged_data[which(merged_data$dir_move==TRUE & merged_data$rec_rate==0.05 & merged_data$n.initial==100),]
sub_n100_dT_r0.1<-merged_data[which(merged_data$dir_move==TRUE & merged_data$rec_rate==0.1 & merged_data$n.initial==100),]

sub_n150_dT_r0.01<-merged_data[which(merged_data$dir_move==TRUE & merged_data$rec_rate==0.01 & merged_data$n.initial==150),]
sub_n150_dT_r0.05<-merged_data[which(merged_data$dir_move==TRUE & merged_data$rec_rate==0.05 & merged_data$n.initial==150),]
sub_n150_dT_r0.1<-merged_data[which(merged_data$dir_move==TRUE & merged_data$rec_rate==0.1 & merged_data$n.initial==150),]

sub_n50_dF_r0.01<-merged_data[which(merged_data$dir_move==FALSE & merged_data$rec_rate==0.01 & merged_data$n.initial==50),]
sub_n50_dF_r0.05<-merged_data[which(merged_data$dir_move==FALSE & merged_data$rec_rate==0.05 & merged_data$n.initial==50),]
sub_n50_dF_r0.1<-merged_data[which(merged_data$dir_move==FALSE & merged_data$rec_rate==0.1 & merged_data$n.initial==50),]

sub_n100_dF_r0.01<-merged_data[which(merged_data$dir_move==FALSE & merged_data$rec_rate==0.01 & merged_data$n.initial==100),]
sub_n100_dF_r0.05<-merged_data[which(merged_data$dir_move==FALSE & merged_data$rec_rate==0.05 & merged_data$n.initial==100),]
sub_n100_dF_r0.1<-merged_data[which(merged_data$dir_move==FALSE & merged_data$rec_rate==0.1 & merged_data$n.initial==100),]

sub_n150_dF_r0.01<-merged_data[which(merged_data$dir_move==FALSE & merged_data$rec_rate==0.01 & merged_data$n.initial==150),]
sub_n150_dF_r0.05<-merged_data[which(merged_data$dir_move==FALSE & merged_data$rec_rate==0.05 & merged_data$n.initial==150),]
sub_n150_dF_r0.1<-merged_data[which(merged_data$dir_move==FALSE & merged_data$rec_rate==0.1 & merged_data$n.initial==150),]



# Delta Prevalence and Duration -------------------------------------------

median(sub_n150_dF_r0.01$max_prevalence)
median(sub_n150_dT_r0.01$max_prevalence)

mean(sub_n150_dF_r0.01$max_prevalence)
mean(sub_n150_dT_r0.01$max_prevalence)

max(sub_n150_dF_r0.01$max_prevalence)
max(sub_n150_dT_r0.01$max_prevalence)

median(sub_n150_dF_r0.01$duration)
median(sub_n150_dT_r0.01$duration)

mean(sub_n150_dF_r0.01$duration)
mean(sub_n150_dT_r0.01$duration)

max(sub_n150_dF_r0.01$duration)
max(sub_n150_dT_r0.01$duration)

median(sub_n150_dF_r0.05$max_prevalence)
median(sub_n150_dT_r0.05$max_prevalence)

mean(sub_n150_dF_r0.05$max_prevalence)
mean(sub_n150_dT_r0.05$max_prevalence)

max(sub_n150_dF_r0.05$max_prevalence)
max(sub_n150_dT_r0.05$max_prevalence)

median(sub_n150_dF_r0.1$max_prevalence)
median(sub_n150_dT_r0.1$max_prevalence)

mean(sub_n150_dF_r0.1$max_prevalence)
mean(sub_n150_dT_r0.1$max_prevalence)

max(sub_n150_dF_r0.1$max_prevalence)
max(sub_n150_dT_r0.1$max_prevalence)
# HEATMAPS ----------------------------------------------------------------

##Heatmap of max_prevelance for n.initial=50 (Decay rates)
sdf <- summaryBy(max_I~rec_rate+dur_scent +initial_load +scent_decay + inf_decay +dir_move, data=size50, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(dir_move), y=as.factor(rec_rate), fill=max_I.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Maximum number of \n infected individuals")
gg <- gg + coord_equal()
gg <- gg + facet_grid(scent_decay~ inf_decay, labeller = labeller(
  scent_decay = c("0.01" = "Slow scent decay", "0.1"="Med. scent decay", "1"= "Fast scent decay"),
  inf_decay = c("0.01" = "Slow pathogen decay", "0.1"="Med. pathogen decay", "1" = "Fast pathogen decay")
))
gg<- gg+ labs(x ="Directed movement", y="Recovery rate")
gg<- gg+theme(legend.position="bottom")
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=8)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg


##Heatmap of max_prevelance for n.initial=100 (Decay rates)
sdf <- summaryBy(max_I~rec_rate+dur_scent +initial_load +scent_decay + inf_decay +dir_move, data=size100, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(dir_move), y=as.factor(rec_rate), fill=max_I.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Maximum number of \n infected individuals")
gg <- gg + coord_equal()
gg <- gg + facet_grid(scent_decay~ inf_decay, labeller = labeller(
  scent_decay = c("0.01" = "Slow scent decay", "0.1"="Med. scent decay", "1"= "Fast scent decay"),
  inf_decay = c("0.01" = "Slow pathogen decay", "0.1"="Med. pathogen decay", "1" = "Fast pathogen decay")
))
gg<- gg+ labs(x ="Directed movement", y="Recovery rate")
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=8)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg

##Heatmap of max_prevelance for n.initial=150 (Decay rates)
sdf <- summaryBy(max_I~rec_rate+dur_scent +initial_load +scent_decay + inf_decay +dir_move, data=size150, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(dir_move), y=as.factor(rec_rate), fill=max_I.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Maximum number of \n infected individuals")
gg <- gg + coord_equal()
gg <- gg + facet_grid(scent_decay~ inf_decay, labeller = labeller(
  scent_decay = c("0.01" = "Slow scent decay", "0.1"="Med. scent decay", "1"= "Fast scent decay"),
  inf_decay = c("0.01" = "Slow pathogen decay", "0.1"="Med. pathogen decay", "1" = "Fast pathogen decay")
))
gg<- gg+ labs(x ="Directed movement", y="Recovery rate")
gg<- gg+theme(legend.position="bottom")
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=8)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg

##Heatmap duration for n.initial=50 (Decay rates)
sdf <- summaryBy(duration~rec_rate+dur_scent +initial_load +scent_decay + inf_decay + dir_move, data=size50, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(dir_move), y=as.factor(rec_rate), fill=duration.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Duration")
gg <- gg + coord_equal()
gg <- gg + facet_grid(scent_decay~ inf_decay, labeller = labeller(
  scent_decay = c("0.01" = "Slow scent decay", "0.1"="Med. scent decay", "1"= "Fast scent decay"),
  inf_decay = c("0.01" = "Slow pathogen decay", "0.1"="Med. pathogen decay", "1" = "Fast pathogen decay")
))
gg<- gg+ labs(x ="Directed movement", y="Recovery rate")
gg<- gg+theme(legend.position="bottom")
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=8)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg


##Heatmap duration for n.initial=100 (Decay rates)
sdf <- summaryBy(duration~rec_rate+dur_scent +initial_load +scent_decay + inf_decay + dir_move, data=size100, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(dir_move), y=as.factor(rec_rate), fill=duration.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Duration")
gg <- gg + coord_equal()
gg <- gg + facet_grid(scent_decay~ inf_decay, labeller = labeller(
  scent_decay = c("0.01" = "Slow scent decay", "0.1"="Med. scent decay", "1"= "Fast scent decay"),
  inf_decay = c("0.01" = "Slow pathogen decay", "0.1"="Med. pathogen decay", "1" = "Fast pathogen decay")
))
gg<- gg+ labs(x ="Directed movement", y="Recovery rate")
gg<- gg+theme(legend.position="bottom")
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg

##Heatmap duration for n.initial=150 (Decay rates)
sdf <- summaryBy(duration~rec_rate+dur_scent +initial_load +scent_decay + inf_decay + dir_move, data=size150, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(dir_move), y=as.factor(rec_rate), fill=duration.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Duration")
gg <- gg + coord_equal()
gg <- gg + facet_grid(scent_decay~ inf_decay, labeller = labeller(
  scent_decay = c("0.01" = "Slow scent decay", "0.1"="Med. scent decay", "1"= "Fast scent decay"),
  inf_decay = c("0.01" = "Slow pathogen decay", "0.1"="Med. pathogen decay", "1" = "Fast pathogen decay")
))
gg<- gg+ labs(x ="Directed movement", y="Recovery rate")
gg<- gg+theme(legend.position="bottom")
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg

##Heatmap of max_prevelance for n.initial=50 (Initial deposits)
sdf <- summaryBy(max_I~rec_rate+scent_load +pathogen_load +scent_decay + inf_decay +dir_move, data=size50, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(dir_move), y=as.factor(rec_rate), fill=max_I.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Maximum number of \n infected individuals")
gg <- gg + coord_equal()
gg <- gg + facet_grid(pathogen_load~scent_load, labeller = labeller(
  initial_load = c("0.5" = "Low pathogen load", "1"= "Med. pathogen load", "10"="High pathogen load"),
  dur_scent = c("1" = "Low scent load", "1" = "Med. scent load", "10" = "High scent load")
))
gg<- gg+ labs(x ="Directed movement", y="Recovery rate")
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg


##Heatmap of max_prevelance for n.initial=100 (Initial deposits)
sdf <- summaryBy(max_I~rec_rate+scent_load +pathogen_load +scent_decay + inf_decay +dir_move, data=size100, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(dir_move), y=as.factor(rec_rate), fill=max_I.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Maximum number of \n infected individuals")
gg <- gg + coord_equal()
gg <- gg + facet_grid(pathogen_load~scent_load, labeller = labeller(
  initial_load = c("0.5" = "Low pathogen load", "1"= "Med. pathogen load", "10"="High pathogen load"),
  dur_scent = c("0.5" = "Low scent load", "1" = "Med. scent load", "10" = "High scent load")
))
gg<- gg+ labs(x ="Directed movement", y="Recovery rate")
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg

##Heatmap of max_prevelance for n.initial=150 (Initial deposits)
sdf <- summaryBy(max_I~rec_rate+scent_load +pathogen_load +scent_decay + inf_decay +dir_move, data=size150, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(dir_move), y=as.factor(rec_rate), fill=max_I.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Maximum number of \n infected individuals")
gg <- gg + coord_equal()
gg <- gg + facet_grid(pathogen_load~scent_load, labeller = labeller(
  initial_load = c("0.5" = "Low pathogen load", "1"= "Med. pathogen load", "10"="High pathogen load"),
  dur_scent = c("0.5" = "Low scent load", "1" = "Med. scent load", "10" = "High scent load")
))
gg<- gg+ labs(x ="Directed movement", y="Recovery rate")
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg


##Heatmap duration for n.initial=50 (Initial deposits)
sdf <- summaryBy(duration~rec_rate+scent_load +pathogen_load +scent_decay + inf_decay +dir_move, data=size50, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(dir_move), y=as.factor(rec_rate), fill=duration.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Duration")
gg <- gg + coord_equal()
gg <- gg + facet_grid(pathogen_load~scent_load, labeller = labeller(
  initial_load = c("0.5" = "Low pathogen load", "1"= "Med. pathogen load", "10"="High pathogen load"),
  dur_scent = c("0.5" = "Low scent load", "1" = "Med. scent load", "10" = "High scent load")
))
gg<- gg+ labs(x ="Directed movement", y="Recovery rate")
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg

##Heatmap duration for n.initial=100 (Initial deposits)
sdf <- summaryBy(duration~rec_rate+scent_load +pathogen_load +scent_decay + inf_decay +dir_move, data=size100, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(dir_move), y=as.factor(rec_rate), fill=duration.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Duration")
gg <- gg + coord_equal()
gg <- gg + facet_grid(pathogen_load~scent_load, labeller = labeller(
  initial_load = c("0.5" = "Low pathogen load", "1"= "Med. pathogen load", "10"="High pathogen load"),
  dur_scent = c("0.5" = "Low scent load", "1" = "Med. scent load", "10" = "High scent load")
))
gg<- gg+ labs(x ="Directed movement", y="Recovery rate")
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg

##Heatmap duration for n.initial=150 (Initial deposits)
sdf <- summaryBy(duration~rec_rate+scent_load +pathogen_load +scent_decay + inf_decay +dir_move, data=size150, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(dir_move), y=as.factor(rec_rate), fill=duration.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Duration")
gg <- gg + coord_equal()
gg <- gg + facet_grid(pathogen_load~scent_load, labeller = labeller(
  initial_load = c("0.5" = "Low pathogen load", "1"= "Med. pathogen load", "10"="High pathogen load"),
  dur_scent = c("0.5" = "Low scent load", "1" = "Med. scent load", "10" = "High scent load")
))
gg<- gg+ labs(x ="Directed movement", y="Recovery rate")
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg


##Heatmap duration for n.initial=100, dir_move=T and recovery rate= 0.01
sdf <- summaryBy(duration~ rec_rate+scent_load +pathogen_load +scent_decay + inf_decay +dir_move, data=sub_n100_dT_r0.01, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(pathogen_load), y=as.factor(scent_load), fill=duration.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Mean\nDuration")
gg <- gg + coord_equal()
gg <- gg + facet_grid(scent_decay~ inf_decay, labeller = labeller(
  scent_decay = c("0.01" = "Slow scent decay", "0.1"="Med. scent decay", "1"= "Fast scent decay"),
  inf_decay = c("0.01" = "Slow pathogen decay", "0.1"="Med. pathogen decay", "1" = "Fast pathogen decay")
))
gg<- gg+ labs(x ="Initial Pathogen Load", y="Initial Scent Load")
gg<- gg+ theme(axis.title = element_text(face="bold", size=20))
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=12)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg

sdf <- summaryBy(max_prevalence~ rec_rate+scent_load +pathogen_load +scent_decay + inf_decay +dir_move, data=sub_n100_dT_r0.01, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(pathogen_load), y=as.factor(scent_load), fill=max_prevalence.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Mean\nMaximum\nPrevalence")
gg <- gg + coord_equal()
gg <- gg + facet_grid(scent_decay~ inf_decay, labeller = labeller(
  scent_decay = c("0.01" = "Slow scent decay", "0.1"="Med. scent decay", "1"= "Fast scent decay"),
  inf_decay = c("0.01" = "Slow pathogen decay", "0.1"="Med. pathogen decay", "1" = "Fast pathogen decay")
))
gg<- gg+ labs(x ="Initial Pathogen Load", y="Initial Scent Load")
gg<- gg+ theme(axis.title = element_text(face="bold", size=20))
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=14)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg

##Heatmap duration for n.initial=150, dir_move=T and recovery rate= 0.01
sdf <- summaryBy(duration~ rec_rate+scent_load +pathogen_load +scent_decay + inf_decay +dir_move, data=sub_n150_dT_r0.01, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(pathogen_load), y=as.factor(scent_load), fill=duration.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Mean\nDuration")
gg <- gg + coord_equal()
gg <- gg + facet_grid(scent_decay~ inf_decay, labeller = labeller(
  scent_decay = c("0.01" = "Slow scent decay", "0.1"="Med. scent decay", "1"= "Fast scent decay"),
  inf_decay = c("0.01" = "Slow pathogen decay", "0.1"="Med. pathogen decay", "1" = "Fast pathogen decay")
))
gg<- gg+ labs(x ="Initial Pathogen Load", y="Initial Scent Load")
gg<- gg+ theme(axis.title = element_text(face="bold", size=20))
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=12)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg

sdf <- summaryBy(max_prevalence~ rec_rate+scent_load +pathogen_load +scent_decay + inf_decay +dir_move, data=sub_n150_dT_r0.01, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(pathogen_load), y=as.factor(scent_load), fill=max_prevalence.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Mean\nMaximum\nPrevalence")
gg <- gg + coord_equal()
gg <- gg + facet_grid(scent_decay~ inf_decay, labeller = labeller(
  scent_decay = c("0.01" = "Slow scent decay", "0.1"="Med. scent decay", "1"= "Fast scent decay"),
  inf_decay = c("0.01" = "Slow pathogen decay", "0.1"="Med. pathogen decay", "1" = "Fast pathogen decay")
))
gg<- gg+ labs(x ="Initial Pathogen Load", y="Initial Scent Load")
gg<- gg+ theme(axis.title = element_text(face="bold", size=20))
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=14)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg


# BOX PLOTS ---------------------------------------------------------------

sub_data<-merged_data[which(merged_data$max_I>1),] #successful outbreaks
# colnames(sub_data)[3]<-"d" #rename density as "d"
# colnames(sub_data)[11]<-"r" #rename perceptual range as "r"
#means <- aggregate(max_prevalence ~ density+ rec_rate +percep, sub_data, mean)
#means_dur <- aggregate(duration ~ density+ rec_rate +percep, sub_data, mean)
#tiff("boxplot_prev.tiff", height = 8, width = 12, units = "in", compression = "lzw", res = 300)
gg<- ggplot(data=size150, aes(y=max_prevalence, x=as.factor(rec_rate)))
#gg<- gg+ geom_bar(stat="identity", position=position_dodge(), colour="black")+ facet_grid(Rec.prob. ~ percep, scales= "free_x", labeller=label_both)
gg<-gg+ geom_boxplot()+ facet_grid(  inf_decay~dir_move, scales= "free_x", labeller= labeller(
  dir_move = c("FALSE" = "Random", "TRUE"="Stigmergy"),
  inf_decay = c("0.01" = "Slow pathogen decay", "0.1"="Med. pathogen decay", "1" = "Fast pathogen decay")))
gg<- gg+ labs(x = expression(Recovery~rate~(gamma)), y = "Maximum prevalence") #,title="(A)")
#gg<- gg+ theme(axis.title.x = element_text(face="bold", size=10),title=element_text(face="bold", size=10))
gg<- gg + theme(axis.text = element_text(size=20))
gg<- gg + theme(axis.text.x =element_text(vjust= 0, angle=-90))
gg<- gg + theme(axis.title = element_text(vjust=0.5, size=18))
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=18))
gg<- gg + theme(plot.title = element_text(size = 18))
gg<- gg + theme(plot.margin=unit(c(0.1,0.1,0,0.1), "cm"))
#gg<-gg + geom_text(data = means , aes(label = duration, y = duration + 0.08))
gg

gg<- ggplot(data=size150, aes(y=duration, x=as.factor(rec_rate)))
#gg<- gg+ geom_bar(stat="identity", position=position_dodge(), colour="black")+ facet_grid(Rec.prob. ~ percep, scales= "free_x", labeller=label_both)
gg<-gg+ geom_boxplot()+ facet_grid(  inf_decay~dir_move, scales= "free_x", labeller= labeller(
  dir_move = c("FALSE" = "Random", "TRUE"="Stigmergy"),
  inf_decay = c("0.01" = "Slow pathogen decay", "0.1"="Med. pathogen decay", "1" = "Fast pathogen decay")))
gg<- gg+ labs(x = expression(Recovery~rate~(gamma)), y = "Duration") #,title="(A)")
#gg<- gg+ theme(axis.title.x = element_text(face="bold", size=10),title=element_text(face="bold", size=10))
gg<- gg + theme(axis.text = element_text(size=20))
gg<- gg + theme(axis.text.x =element_text(vjust= 0, angle=-90))
gg<- gg + theme(axis.title = element_text(vjust=0.5, size=18))
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=18))
gg<- gg + theme(plot.title = element_text(size = 18))
gg<- gg + theme(plot.margin=unit(c(0.1,0.1,0,0.1), "cm"))
#gg<-gg + geom_text(data = means , aes(label = duration, y = duration + 0.08))
gg



ii<- ggplot(data=sub_n150_dT_r0.01, aes(x=as.factor(scent_load), y=max_prevalence, fill=as.factor(pathogen_load)))
ii<- ii+ geom_boxplot() #+ geom_jitter(alpha=0.5)
ii <- ii + facet_grid(scent_decay~ inf_decay, labeller = labeller(
  scent_decay = c("0.01" = "Slow scent decay", "0.1"="Med. scent decay", "1"= "Fast scent decay"),
  inf_decay = c("0.01" = "Slow pathogen decay", "0.1"="Med. pathogen decay", "1" = "Fast pathogen decay")
))
ii<- ii+ labs(x = "Scent Load", y = "Maximum prevalence", fill= "Pathogen Load")
ii<- ii+ theme(plot.title = element_text(size=12))
ii<- ii+ theme(axis.title = element_text(size=11))
ii<- ii + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=0, size=10))
ii<- ii + theme(axis.text.y = element_text(vjust=0.5, size=10)) 
ii<- ii + theme(strip.text = element_text(vjust=0.5, size=11)) 
ii<- ii + theme(legend.text = element_text(vjust=0.5, size=11), legend.title=element_text(vjust=0.5, size=11))
ii<- ii + geom_vline(xintercept=c(6.5,12.5), linetype="dashed")
ii<- ii + theme(plot.margin=unit(c(0.1,0.1,0,0.1), "cm"))
ii

ii<- ggplot(data=sub_n150_dT_r0.01, aes(x=as.factor(scent_load), y=duration, fill=as.factor(pathogen_load)))
ii<- ii+ geom_boxplot() #+ geom_jitter(alpha=0.5)
ii <- ii + facet_grid(scent_decay~ inf_decay, labeller = labeller(
  scent_decay = c("0.01" = "Slow scent decay", "0.1"="Med. scent decay", "1"= "Fast scent decay"),
  inf_decay = c("0.01" = "Slow pathogen decay", "0.1"="Med. pathogen decay", "1" = "Fast pathogen decay")
))
ii<- ii+ labs(x = "Scent Load", y = "Duration", fill= "Pathogen Load")
ii<- ii+ theme(plot.title = element_text(size=12))
ii<- ii+ theme(axis.title = element_text(size=11))
ii<- ii + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=0, size=10))
ii<- ii + theme(axis.text.y = element_text(vjust=0.5, size=10)) 
ii<- ii + theme(strip.text = element_text(vjust=0.5, size=11)) 
ii<- ii + theme(legend.text = element_text(vjust=0.5, size=11), legend.title=element_text(vjust=0.5, size=11))
ii<- ii + geom_vline(xintercept=c(6.5,12.5), linetype="dashed")
ii<- ii + theme(plot.margin=unit(c(0.1,0.1,0,0.1), "cm"))
ii


ii<- ggplot(data=sub_n150_dT_r0.01, aes(x=as.factor(inf_decay), y=max_prevalence, fill=as.factor(scent_decay)))
# ii<- ii+ geom_boxplot() #+ geom_jitter(alpha=0.5)
# ii<- ii+ stat_summary(fun.y=mean, geom="point", col="blue")
ii <- ii+ geom_boxplot() + facet_grid(scent_load~ pathogen_load, labeller = labeller(
  scent_load = c("0.5" = "Low SL", "1"= "Med SL", "10"= "High SL"),
  pathogen_load = c("0.5" = "Low PL","1" = "Med PL", "10" = "High PL")))
ii<- ii+ labs(x = "Infection decay rate", y = "Maximum prevalence", fill= "Scent\ndecay rate")
ii<- ii+ theme(plot.title = element_text(size=12))
ii<- ii+ theme(axis.title = element_text(size=11))
ii<- ii + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=0, size=10))
ii<- ii + theme(axis.text.y = element_text(vjust=0.5, size=10)) 
ii<- ii + theme(strip.text = element_text(vjust=0.5, size=11)) 
ii<- ii + theme(legend.text = element_text(vjust=0.5, size=11), legend.title=element_text(vjust=0.5, size=11))
ii<- ii + geom_vline(xintercept=c(6.5,12.5), linetype="dashed")
ii<- ii + theme(plot.margin=unit(c(0.1,0.1,0,0.1), "cm"))
ii

ii<- ggplot(data=sub_n150_dT_r0.01, aes(x=as.factor(inf_decay), y=duration, fill=as.factor(scent_decay)))
# ii<- ii+ geom_boxplot() #+ geom_jitter(alpha=0.5)
# ii<- ii+ stat_summary(fun.y=mean, geom="point", col="blue")
ii <- ii+ geom_boxplot() + facet_grid(scent_load~ pathogen_load, labeller = labeller(
  scent_load = c("0.5" = "Low SL", "1"= "Med SL", "10"= "High SL"),
  pathogen_load = c("0.5" = "Low PL","1" = "Med PL", "10" = "High PL")))
ii<- ii+ labs(x = "Infection decay rate", y = "Duration", fill= "Scent\ndecay rate")
ii<- ii+ theme(plot.title = element_text(size=12))
ii<- ii+ theme(axis.title = element_text(size=11))
ii<- ii + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=0, size=10))
ii<- ii + theme(axis.text.y = element_text(vjust=0.5, size=10)) 
ii<- ii + theme(strip.text = element_text(vjust=0.5, size=11)) 
ii<- ii + theme(legend.text = element_text(vjust=0.5, size=11), legend.title=element_text(vjust=0.5, size=11))
ii<- ii + geom_vline(xintercept=c(6.5,12.5), linetype="dashed")
ii<- ii + theme(plot.margin=unit(c(0.1,0.1,0,0.1), "cm"))
ii
