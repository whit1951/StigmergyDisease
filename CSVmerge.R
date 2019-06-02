#'May 16, 2019
#'For analysis of stigmergy simulations
#'Merge read ins from data files

library(doBy)
library(ggplot2)
library(viridis)
library(ggthemes)
library(scales)
library(lsr)

##Read in file names from directory
filenames <- list.files(path = "~/StigmergyDisease/_rslurm_stig_sim")
summaries<-filenames[grep("summary", filenames)] #summary data
infecteds<-filenames[grep("infected", filenames)] #I data


##Check complete list
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

complete_list<-vector(mode="character")
params<-expand.grid(maxT=maxT, nsim=nsim, lsize=lsize, n.initial=n.initial, inf_prob=inf_prob, rec_rate=rec_rate, dur_scent=dur_scent, initial_load=initial_load, scent_decay=scent_decay, inf_decay=inf_decay)
for(i in 1:nrow(params)){
complete_list[i]<-paste("lsize", params$lsize[i], "n.initial", params$n.initial[i], "inf", params$inf_prob[i], "rec", params$rec_rate[i],"dur_scent", params$dur_scent[i],"initial_load",params$initial_load[i], "scent_decay", params$scent_decay[i], "inf_decay", params$inf_decay[i], "_summary.csv", sep="")
}

miss_sim<-logical(length(complete_list))
for(i in 1:length(complete_list)){
if(complete_list[i] %in% summaries)
  miss_sim[i]<-TRUE
}
which(miss_sim==FALSE)  

##Make a merged dataset from all available .csv files
setwd("~/StigmergyDisease/_rslurm_stig_sim")
merged_data<- do.call("rbind", lapply(summaries, read.csv, header = TRUE))
merged_data$X<-NULL
merged_data$inf_prob[which(merged_data$inf_prob==0.1)]<-0.10 #Same number of precision points
#merged_data$betas <- paste(merged_data$beta1,merged_data$beta2)
#merged_data$betas <- ordered(merged_data$betas, levels = c("0 -6", "0 -3", "0 0", "3 -6", "3 -3", "3 0", "6 -6", "6 -3", "6 0"))
#merged_data$pbyH <- paste(merged_data$p,merged_data$H) #, sep=',')
merged_data[which(is.na(merged_data$duration)),]
# merged_data$duration[which(is.na(merged_data$duration))]<-1000

#Create character strings for landscape structure (Hurst exponent and proportion available habitat)
merged_data$PL<-NA #pathogen load
merged_data$SL<-NA #scent load
merged_data$PL[which(merged_data$initial_load==1)]<-"lowPL"
merged_data$PL[which(merged_data$initial_load==10)]<-"highPL"
merged_data$SL[which(merged_data$dur_scent==1)]<-"lowSL"
merged_data$SL[which(merged_data$dur_scent==10)]<-"highSL"
merged_data$decay <- paste(merged_data$scent_decay,merged_data$inf_decay)
# write.csv(merged_data, "stig_summary.csv")
#Order factors
# merged_data$betas<- factor(merged_data$betas, levels = c("0 0 0", "0 0 -1", "0 1 -1", "0 1 -0.5", "0 2 -1", "0 2 -0.5", "3 0 0", "3 0 -1", "3 1 -1", "3 1 -0.5", "3 2 -1", "3 2 -0.5", "6 0 0", "6 0 -1", "6 1 -1", "6 1 -0.5", "6 2 -1", "6 2 -0.5"))
# merged_data$lstructure<- factor(merged_data$lstructure, levels = c("L/HP", "L/MP", "L/LP", "M/HP", "M/MP", "M/LP", "H/HP", "H/MP", "H/LP"))


##Make some sub data sets by density and recovery rate
size125<-merged_data[which(merged_data$n.initial==125),]
size250<-merged_data[which(merged_data$n.initial==250),]

sub_i0.1_r0.05<-merged_data[which(merged_data$inf_prob==0.1 & merged_data$rec_rate==0.05),]
sub_i0.1_r0.1<-merged_data[which(merged_data$inf_prob==0.1 & merged_data$rec_rate==0.1),]

sub_i0.25_r0.01<-merged_data[which(merged_data$inf_prob==0.25 & merged_data$rec_rate==0.01),]
sub_i0.25_r0.05<-merged_data[which(merged_data$inf_prob==0.25 & merged_data$rec_rate==0.05),]
sub_i0.25_r0.1<-merged_data[which(merged_data$inf_prob==0.25 & merged_data$rec_rate==0.1),]

sub_i0.5_r0.01<-merged_data[which(merged_data$inf_prob==0.5 & merged_data$rec_rate==0.01),]
sub_i0.5_r0.05<-merged_data[which(merged_data$inf_prob==0.5 & merged_data$rec_rate==0.05),]
sub_i0.5_r0.1<-merged_data[which(merged_data$inf_prob==0.5 & merged_data$rec_rate==0.1),]

sub_n125_i0.1_r0.01<-merged_data[which(merged_data$inf_prob==0.1 & merged_data$rec_rate==0.01 & merged_data$n.initial==125),]
sub_n125_i0.1_r0.05<-merged_data[which(merged_data$inf_prob==0.1 & merged_data$rec_rate==0.05 & merged_data$n.initial==125),]
sub_n125_i0.1_r0.1<-merged_data[which(merged_data$inf_prob==0.1 & merged_data$rec_rate==0.01 & merged_data$n.initial==125),]

sub_n125_i0.25_r0.05<-merged_data[which(merged_data$inf_prob==0.25 & merged_data$rec_rate==0.05 & merged_data$n.initial==125),]
sub_n125_i0.25_r0.1<-merged_data[which(merged_data$inf_prob==0.25 & merged_data$rec_rate==0.1 & merged_data$n.initial==125),]
sub_n125_i0.25_r0.01<-merged_data[which(merged_data$inf_prob==0.25 & merged_data$rec_rate==0.01 & merged_data$n.initial==125),]


sub_n125_i0.5_r0.05<-merged_data[which(merged_data$inf_prob==0.5 & merged_data$rec_rate==0.05 & merged_data$n.initial==125),]
sub_n125_i0.5_r0.1<-merged_data[which(merged_data$inf_prob==0.5 & merged_data$rec_rate==0.1 & merged_data$n.initial==125),]
sub_n125_i0.5_r0.01<-merged_data[which(merged_data$inf_prob==0.5 & merged_data$rec_rate==0.01 & merged_data$n.initial==125),]

sub_n250_i0.5_r0.05<-merged_data[which(merged_data$inf_prob==0.5 & merged_data$rec_rate==0.05 & merged_data$n.initial==250),]
sub_n250_i0.5_r0.1<-merged_data[which(merged_data$inf_prob==0.5 & merged_data$rec_rate==0.1 & merged_data$n.initial==250),]
sub_n250_i0.5_r0.01<-merged_data[which(merged_data$inf_prob==0.5 & merged_data$rec_rate==0.01 & merged_data$n.initial==250),]


# HEATMAPS ----------------------------------------------------------------

##Heatmap of max_prevelance for n.initial=125 (Decay rates)
sdf <- summaryBy(max_I~inf_prob + rec_rate+dur_scent +initial_load +scent_decay + inf_decay, data=size125, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(inf_prob), y=as.factor(rec_rate), fill=max_I.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Maximum number of \n infected individuals")
gg <- gg + coord_equal()
gg <- gg + facet_grid(scent_decay~ inf_decay, labeller = labeller(
  scent_decay = c("0.1" = "Slow scent decay", "0.5"= "Fast scent decay"),
  inf_decay = c("0.1" = "Slow pathogen decay", "0.5" = "Fast pathogen decay")
))
gg<- gg+ labs(x ="Infection probability", y="Recovery rate")
gg<- gg+theme(legend.position="bottom")
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg


##Heatmap of max_prevelance for n.initial=250 (Decay rates)
sdf <- summaryBy(max_I~inf_prob + rec_rate+dur_scent +initial_load +scent_decay + inf_decay, data=size250, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(inf_prob), y=as.factor(rec_rate), fill=max_I.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Maximum number of \n infected individuals")
gg <- gg + coord_equal()
gg <- gg + facet_grid(scent_decay~ inf_decay, labeller = labeller(
  scent_decay = c("0.1" = "Slow scent decay", "0.5"= "Fast scent decay"),
  inf_decay = c("0.1" = "Slow pathogen decay", "0.5" = "Fast pathogen decay")
))
gg<- gg+ labs(x ="Infection probability", y="Recovery rate")
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg

##Heatmap duration for n.initial=125 (Decay rates)
sdf <- summaryBy(duration~inf_prob + rec_rate+dur_scent +initial_load +scent_decay + inf_decay, data=size125, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(inf_prob), y=as.factor(rec_rate), fill=duration.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Duration")
gg <- gg + coord_equal()
gg <- gg + facet_grid(scent_decay~ inf_decay, labeller = labeller(
  scent_decay = c("0.1" = "Slow scent decay", "0.5"= "Fast scent decay"),
  inf_decay = c("0.1" = "Slow pathogen decay", "0.5" = "Fast pathogen decay")
))
gg<- gg+ labs(x ="Infection probability", y="Recovery rate")
gg<- gg+ labs(x ="Infection probability", y="Recovery rate")
gg<- gg+ theme(axis.title = element_text(face="bold", size=20))
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=14)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg


##Heatmap duration for n.initial=250 (Decay rates)
sdf <- summaryBy(duration~inf_prob + rec_rate+dur_scent +initial_load +scent_decay + inf_decay, data=size250, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(inf_prob), y=as.factor(rec_rate), fill=duration.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Duration")
gg <- gg + coord_equal()
gg <- gg + facet_grid( initial_load~ dur_scent)
gg<- gg+ labs(x ="Infection probability", y="Recovery rate")
gg<- gg+ theme(axis.title = element_text(face="bold", size=20))
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=14)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg

##Heatmap of max_prevelance for n.initial=125 (Initial deposits)
sdf <- summaryBy(max_I~inf_prob + rec_rate+dur_scent +initial_load +scent_decay + inf_decay, data=size125, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(inf_prob), y=as.factor(rec_rate), fill=max_I.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Maximum number of \n infected individuals")
gg <- gg + coord_equal()
gg <- gg + facet_grid(initial_load~dur_scent, labeller = labeller(
  initial_load = c("1" = "Low pathogen load", "10"= "High pathogen load"),
  dur_scent = c("1" = "Low scent load", "10" = "High scent load")
))
gg<- gg+ labs(x ="Infection probability", y="Recovery rate")
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg


##Heatmap of max_prevelance for n.initial=250 (Initial deposits)
sdf <- summaryBy(max_I~inf_prob + rec_rate+dur_scent +initial_load +scent_decay + inf_decay, data=size250, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(inf_prob), y=as.factor(rec_rate), fill=max_I.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Maximum number of \n infected individuals")
gg <- gg + coord_equal()
gg <- gg + facet_grid(initial_load~dur_scent, labeller = labeller(
  initial_load = c("1" = "Low pathogen load", "10"= "High pathogen load"),
  dur_scent = c("1" = "Low scent load", "10" = "High scent load")
))
gg<- gg+ labs(x ="Infection probability", y="Recovery rate")
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg

##Heatmap duration for n.initial=125 (Initial deposits)
sdf <- summaryBy(duration~inf_prob + rec_rate+dur_scent +initial_load +scent_decay + inf_decay, data=size125, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(inf_prob), y=as.factor(rec_rate), fill=duration.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Duration")
gg <- gg + coord_equal()
gg <- gg + facet_grid(scent_decay~ inf_decay, labeller = labeller(
  scent_decay = c("0.1" = "Slow scent decay", "0.5"= "Fast scent decay"),
  inf_decay = c("0.1" = "Slow pathogen decay", "0.5" = "Fast pathogen decay")
))
gg<- gg+ labs(x ="Infection probability", y="Recovery rate")
gg<- gg+ labs(x ="Infection probability", y="Recovery rate")
gg<- gg+ theme(axis.title = element_text(face="bold", size=20))
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=14)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg


##Heatmap duration for n.initial=250 (Initial deposits)
sdf <- summaryBy(duration~inf_prob + rec_rate+dur_scent +initial_load +scent_decay + inf_decay, data=size250, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(inf_prob), y=as.factor(rec_rate), fill=duration.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Duration")
gg <- gg + coord_equal()
gg <- gg + facet_grid(initial_load~dur_scent, labeller = labeller(
  initial_load = c("1" = "Low pathogen load", "10"= "High pathogen load"),
  dur_scent = c("1" = "Low scent load", "10" = "High scent load")
))
gg<- gg+ labs(x ="Infection probability", y="Recovery rate")
gg<- gg+ theme(axis.title = element_text(face="bold", size=20))
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=14)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg

##Heatmap duration for n.initial=125, inf_prob= 0.1 and recovery rate= 0.05
sdf <- summaryBy(duration~n.initial +inf_prob + rec_rate+dur_scent +initial_load +scent_decay + inf_decay, data=sub_n250_i0.5_r0.01, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(initial_load), y=as.factor(dur_scent), fill=duration.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Duration")
gg <- gg + coord_equal()
gg <- gg + facet_grid(scent_decay~ inf_decay, labeller = labeller(
  scent_decay = c("0.1" = "slow scent decay", "0.5"= "fast scent decay"),
  inf_decay = c("0.1" = "slow pathogen decay", "0.5" = "fast pathogen decay")
))
gg<- gg+ labs(x ="Initial Pathogen Load", y="Initial Scent Load")
gg<- gg+ theme(axis.title = element_text(face="bold", size=20))
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=14)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg

sdf <- summaryBy(max_prevalence~n.initial +inf_prob + rec_rate+dur_scent +initial_load +scent_decay + inf_decay, data=sub_n250_i0.5_r0.05, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(initial_load), y=as.factor(dur_scent), fill=max_prevalence.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Max\n Prevalence")
gg <- gg + coord_equal()
gg <- gg + facet_grid(scent_decay~ inf_decay, labeller = labeller(
  scent_decay = c("0.1" = "slow scent decay", "0.5"= "fast scent decay"),
  inf_decay = c("0.1" = "slow pathogen decay", "0.5" = "fast pathogen decay")
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
gg<- ggplot(data=size250, aes(y=max_prevalence, x=as.factor(rec_rate)))
#gg<- gg+ geom_bar(stat="identity", position=position_dodge(), colour="black")+ facet_grid(Rec.prob. ~ percep, scales= "free_x", labeller=label_both)
gg<-gg+ geom_boxplot()+ facet_grid( ~ inf_prob, scales= "free_x", labeller=label_value)
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

ii<- ggplot(data=sub_n250_i0.5_r0.05, aes(x=SL, y=max_prevalence, fill=as.factor(PL)))
ii<- ii+ geom_boxplot() #+ geom_jitter(alpha=0.5)
ii <- ii + facet_grid(scent_decay~ inf_decay, labeller = labeller(
  scent_decay = c("0.1" = "slow scent decay", "0.5"= "fast scent decay"),
  inf_decay = c("0.1" = "slow pathogen decay", "0.5" = "fast pathogen decay")
))
ii<- ii+ labs(x = "Pathogen Load", y = "Maximum prevalence", fill= "Scent Load")
ii<- ii+ theme(plot.title = element_text(size=12))
ii<- ii+ theme(axis.title = element_text(size=11))
ii<- ii + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=0, size=10))
ii<- ii + theme(axis.text.y = element_text(vjust=0.5, size=10)) 
ii<- ii + theme(strip.text = element_text(vjust=0.5, size=11)) 
ii<- ii + theme(legend.text = element_text(vjust=0.5, size=11), legend.title=element_text(vjust=0.5, size=11))
ii<- ii + geom_vline(xintercept=c(6.5,12.5), linetype="dashed")
ii<- ii + theme(plot.margin=unit(c(0.1,0.1,0,0.1), "cm"))
ii


sub_n250_i0.5_r0.01lPLlSL<-merged_data[which(merged_data$inf_prob==0.5 & merged_data$rec_rate==0.01 & merged_data$n.initial==250 & merged_data$initial_load==10 & merged_data$dur_scent==1),]
ii<- ggplot(data=sub_n250_i0.5_r0.1, aes(x=decay, y=max_I))
ii<- ii+ geom_boxplot() #+ geom_jitter(alpha=0.5)
ii<- ii+ stat_summary(fun.y=mean, geom="point", col="blue")
ii <- ii + facet_grid(dur_scent~ initial_load, labeller = labeller(
  dur_scent = c("1" = "Low Scent Load", "10"= "High Scent Load"),
  initial_load = c("1" = "Low Pathogen Load", "10" = "High Pathogen Load")))
ii<- ii+ labs(x = "Decay rate/time (Scent & Pathogen)", y = "Maximum prevalence", fill= "Scent Load")
ii<- ii+ theme(plot.title = element_text(size=12))
ii<- ii+ theme(axis.title = element_text(size=11))
ii<- ii + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=0, size=10))
ii<- ii + theme(axis.text.y = element_text(vjust=0.5, size=10)) 
ii<- ii + theme(strip.text = element_text(vjust=0.5, size=11)) 
ii<- ii + theme(legend.text = element_text(vjust=0.5, size=11), legend.title=element_text(vjust=0.5, size=11))
ii<- ii + geom_vline(xintercept=c(6.5,12.5), linetype="dashed")
ii<- ii + theme(plot.margin=unit(c(0.1,0.1,0,0.1), "cm"))
ii

