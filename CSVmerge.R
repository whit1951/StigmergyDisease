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
filenames <- list.files(path = "~/StigmergyDisease")
summaries<-filenames[grep("summary", filenames)] #summary data
infecteds<-filenames[grep("infected", filenames)] #I data


##Make a merged dataset from all available .csv files
merged_data<- do.call("rbind", lapply(summaries, read.csv, header = TRUE))
merged_data$X<-NULL
merged_data$inf_prob[which(merged_data$inf_prob==0.1)]<-0.10 #Same number of precision points
#merged_data$betas <- paste(merged_data$beta1,merged_data$beta2)
#merged_data$betas <- ordered(merged_data$betas, levels = c("0 -6", "0 -3", "0 0", "3 -6", "3 -3", "3 0", "6 -6", "6 -3", "6 0"))
#merged_data$pbyH <- paste(merged_data$p,merged_data$H) #, sep=',')
merged_data[which(is.na(merged_data$duration)),]
merged_data$duration[which(is.na(merged_data$duration))]<-1000


##Make some sub data sets by density and recovery rate
size25<-merged_data[which(merged_data$n.initial==25),]
size50<-merged_data[which(merged_data$n.initial==50),]

sub_n25_i0.1_r0.01<-merged_data[which(merged_data$inf_prob==0.1 & merged_data$rec_rate==0.01 & merged_data$n.initial==25),]
sub_i0.1_r0.05<-merged_data[which(merged_data$inf_prob==0.1 & merged_data$rec_rate==0.05),]
sub_i0.1_r0.1<-merged_data[which(merged_data$inf_prob==0.1 & merged_data$rec_rate==0.1),]

sub_i0.25_r0.01<-merged_data[which(merged_data$inf_prob==0.25 & merged_data$rec_rate==0.01),]
sub_i0.25_r0.05<-merged_data[which(merged_data$inf_prob==0.25 & merged_data$rec_rate==0.05),]
sub_i0.25_r0.1<-merged_data[which(merged_data$inf_prob==0.25 & merged_data$rec_rate==0.1),]

sub_i0.5_r0.01<-merged_data[which(merged_data$inf_prob==0.5 & merged_data$rec_rate==0.01),]
sub_i0.5_r0.05<-merged_data[which(merged_data$inf_prob==0.5 & merged_data$rec_rate==0.05),]
sub_i0.5_r0.1<-merged_data[which(merged_data$inf_prob==0.5 & merged_data$rec_rate==0.1),]


# HEATMAPS ----------------------------------------------------------------

##Heatmap of max_prevelance for n.initial=25
sdf <- summaryBy(max_I~inf_prob + rec_rate+dur_scent +initial_load +scent_decay + inf_decay, data=size25, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(inf_prob), y=as.factor(rec_rate), fill=max_I.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Maximum number of \n infected individuals")
gg <- gg + coord_equal()
gg <- gg + facet_grid( scent_decay~ inf_decay)
#gg <- gg + geom_text(aes(label=round(max_prevalence.mean,2)))
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg


##Heatmap of max_prevelance for n.initial=50
sdf <- summaryBy(max_I~inf_prob + rec_rate+dur_scent +initial_load +scent_decay + inf_decay, data=size50, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(inf_prob), y=as.factor(rec_rate), fill=max_I.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Maximum number of \n infected individuals")
gg <- gg + coord_equal()
gg <- gg + facet_grid( initial_load~ dur_scent)
#gg <- gg + geom_text(aes(label=round(max_prevalence.mean,2)))
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg

##Heatmap duration for n.initial=25
sdf <- summaryBy(duration~inf_prob + rec_rate+dur_scent +initial_load +scent_decay + inf_decay, data=size25, FUN=mean)
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


##Heatmap duration for n.initial=50
sdf <- summaryBy(duration~inf_prob + rec_rate+dur_scent +initial_load +scent_decay + inf_decay, data=size50, FUN=mean)
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

##Heatmap duration for n.initial=25, inf_prob= 0.1 and recovery rate= 0.01
sdf <- summaryBy(duration~n.initial +inf_prob + rec_rate+dur_scent +initial_load +scent_decay + inf_decay, data=sub_n25_i0.1_r0.01, FUN=mean)
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

sdf <- summaryBy(max_prevalence~n.initial +inf_prob + rec_rate+dur_scent +initial_load +scent_decay + inf_decay, data=sub_i0.1_r0.01, FUN=mean)
gg <- ggplot(sdf, aes(x=as.factor(initial_load), y=as.factor(dur_scent), fill=max_prevalence.mean))
gg <- gg + geom_tile(color="white", size=0.1)
gg <- gg + scale_fill_viridis(name="Max\n Prevalence")
gg <- gg + coord_equal()
gg <- gg + facet_grid( scent_decay~ inf_decay)
gg<- gg+ labs(x ="Initial Pathogen Load", y="Initial Scent Load")
gg<- gg+ theme(axis.title = element_text(face="bold", size=20))
gg<- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16))
gg<- gg + theme(axis.text.y = element_text(vjust=0.5, size=16)) 
gg<- gg + theme(strip.text = element_text(vjust=0.5, size=14)) 
gg<- gg + theme(legend.text = element_text(vjust=0.5, size=16), legend.title=element_text(vjust=0.5, size=20))
gg
