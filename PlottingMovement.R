#' Plotting movement data
#' @author Lauren White
#' @date May 7, 2019

#install.packages('moveVis')
library(moveVis)
library(move)
library(magrittr)
library(ggplot2)
library(tidyr)

#colnames(movedat)<-1:100
#movedat$AnimalID<-1:20
#test<-matrix(movedat, dimnames=list(t(outer(colnames(movedat), rownames(movedat), FUN=paste)), NULL))
#movedf<-data.frame(Time= rep(1:100, each=20), AnimalID= rep(1:20, each=100), vec=NA, xloc=NA, yloc=NA)
movedat<-t(movedat)
rownames(movedat)<-1:nrow(movedat)
movedat<-as.data.frame(movedat)
movedat$time<-1:nrow(movedat)


movedf<-gather(movedat, key="AnimalID", value= "Vec",colnames(movedat)[1:20])
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

m <- ggplot(movedf, aes(xloc, yloc, col=AnimalID))
m + geom_path()

V2<-movedf[movedf$AnimalID=="V2",]
V2$yloc[which(V2$dist>1.5)]<-NA
m <- ggplot(V2, aes(xloc, yloc))
m + geom_path(na.rm=TRUE)



###################################################
#NEED To convert time to POSIXct in order to animate
#######################################################
movedf$time<-as.POSIXct(movedf$time, format="%Y", origin="0000")

test<-df2move(movedf, proj=NULL, x="xloc", y="yloc", time="time", track_id="AnimalID")
test<-move(movedf$xloc, movedf$yloc, time=movedf$time, proj= "+init=epsg:4326", animal=movedf$AnimalID)

m <- align_move(test, res = 5, digit = 0, unit = "secs")

# create spatial frames with a OpenStreetMap watercolour map
frames <- frames_spatial(m, path_colours = 1:n.initial,
                         map_service = "osm", map_type = "watercolor", alpha = 0.5) %>% 
  add_labels(x = "Longitude", y = "Latitude") %>% # add some customizations, such as axis labels
  add_northarrow() %>% 
  add_scalebar() %>% 
  add_timestamps(m, type = "label") %>% 
  add_progress()

frames[[1]] # preview one of the frames, e.g. the 100th frame
animate_frames(frames, out_file = "/research-home/lwhite/Puma-movement/test.gif")

