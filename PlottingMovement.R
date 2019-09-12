#' Plotting movement data
#' @author Lauren White
#' @date May 7, 2019


library(magrittr)
library(ggplot2)
library(tidyr)
library(gganimate)
# devtools::install_github('thomasp85/transformr')
library(transformr)
library(magick)

source('~/StigmergyDisease/ToyModelStigmergy.R')

movedat<-t(movedat)
rownames(movedat)<-1:nrow(movedat)
movedat<-as.data.frame(movedat)
movedat$time<-1:nrow(movedat)

movedf<-gather(movedat, key="AnimalID", value= "Vec",colnames(movedat)[1:n.initial])
test2<-Rmatrix(movedf$Vec-1, lsize)
movedf$xloc<-test2$col
movedf$yloc<-test2$row
movedf$dist<-rep(0, times=nrow(movedf))

linear_dist<-function(x1,x2, y1, y2)
{
  sqrt((x1-x2)^2+(y1-y2)^2)
}

for(i in 1:(nrow(movedf)-1)){
  movedf$dist[i+1]<-linear_dist(movedf$xloc[i], movedf$xloc[i+1], movedf$yloc[i], movedf$yloc[i+1])
}
movedf$yloc[which(movedf$dist>1.5)]<-NA
colnames(movedf)[colnames(movedf)=="time"] <- "time_step"

## Fig. 2 static movement trajectories
 # m <- ggplot(movedf, aes(xloc, yloc, color=AnimalID, alpha=time_step*0.01))+
    m <- ggplot(movedf, aes(xloc, yloc, color=AnimalID))+
  geom_path(show.legend=FALSE)+xlab(NULL)+ylab(NULL)+
  theme(legend.position = "none") 
m


#animate with gganimate
m <- ggplot(movedf, aes(xloc, yloc, color=AnimalID))+ geom_point(show.legend=FALSE)+xlab(NULL)+ylab(NULL)
m + transition_time(time_step) +
  shadow_trail(distance=0.01, size=0.75, alpha=0.75, max_frames=5) +
  labs(title = "Time step: {frame_time}")
anim_save("stigmergy.gif")

V2<-movedf[movedf$AnimalID=="V2",]
V2$yloc[which(V2$dist>1.5)]<-NA
m <- ggplot(V2, aes(xloc, yloc, alpha=time))
m + geom_path(na.rm=TRUE)


########################
#EXAMPLE
library(gapminder)
head(gapminder)
p <- ggplot(
  gapminder, 
  aes(x = gdpPercap, y=lifeExp, size = pop, colour = country)
) +
  geom_point(show.legend = FALSE, alpha = 0.7) +
  scale_color_viridis_d() +
  scale_size(range = c(2, 12)) +
  scale_x_log10() +
  labs(x = "GDP per capita", y = "Life expectancy")
p