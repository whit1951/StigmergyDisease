#' Plot RF results for stigmergy simulations

## Define multiplot function

#' Multiplot function
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols), byrow=TRUE)
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



## Figure 1
#Random forest variable importance results

library(randomForest)
#setwd()
#Load importance values (not scaled by SD)
imp_logit<-read.csv("imp_logit.csv")
imp_logitprev<-read.csv("imp_prev.csv")
imp_logitdur<-read.csv("imp_dur.csv")

#Order according to increasing node purity
imp_logit<-imp_logit[order(imp_logit$IncNodePurity),]
imp_logitdur<-imp_logitdur[order(imp_logitdur$IncNodePurity),] 
imp_logitprev<-imp_logitprev[order(imp_logitprev$IncNodePurity),] 

#In ggplot
imp_logit<-imp_logit[order(imp_logit$X.IncMSE),] #order according to %increase in MSE
#imp_logit$SD<-rep(0.005, times=nrow(imp_logit))
A<-ggplot(imp_logit, aes(x=reorder(X, X.IncMSE), y=X.IncMSE)) + 
  geom_bar(position=position_dodge(), stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  # geom_errorbar(aes(ymin=X.IncMSE-SD, ymax=X.IncMSE+SD),
  #               size=.3,    # Thinner lines
  #               width=.2,
  #               position=position_dodge(.9)) +
  scale_x_discrete(element_blank(), labels = c("scent_decay" = expression(delta), "rec_rate"=expression(gamma),
                                               "dur_scent" =  expression(eta[0]),"n.initial" = "N", "inf_prob" = expression(beta),
                                               "initial_load" = expression(kappa[0]), "inf_decay"=expression(alpha))) +
  #ylab("Mean Decrease in Accuracy (%IncMSE)") +
  ylab("")+
  theme_bw()+
  ggtitle("(A)") + #" Variable Importance for Epidemic \nSuccess") +
  theme(axis.text=element_text(size=14), plot.title = element_text(size = 16)) +
  theme(plot.margin=unit(c(0.1,0.1,0,0), "cm"))+
  coord_flip()
# theme_bw(axis.text=element_text(size=14),    axis.title=element_text(size=14,face="bold"), plot.title = element_text(size = 40))
A
#imp_logitprev$SD<-rep(0.0000001, times=nrow(imp_logitprev))
B<-ggplot(imp_logitprev, aes(x=reorder(X, X.IncMSE), y=X.IncMSE)) + 
  geom_bar(position=position_dodge(), stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  # geom_errorbar(aes(ymin=X.IncMSE-SD, ymax=X.IncMSE+SD),
  #               size=.3,    # Thinner lines
  #               width=.2,
  #               position=position_dodge(.9)) +
  scale_x_discrete(element_blank(), labels = c("scent_decay" = expression(delta), "rec_rate"=expression(gamma),
                                               "dur_scent" =  expression(eta[0]),"n.initial" = "N", "inf_prob" = expression(beta),
                                               "initial_load" = expression(kappa[0]), "inf_decay"=expression(alpha))) +
  #ylab("Mean Decrease in Accuracy (%IncMSE)") +
  ylab("")+
  theme_bw()+
  ggtitle("(B)") + #" Variable Importance for Maximum \nPrevalence|Success") +
  theme(axis.text=element_text(size=14), plot.title = element_text(size = 16)) +
  theme(plot.margin=unit(c(-0.25,0.1,0,0), "cm"))+
  coord_flip()
B
#imp_logitdur$SD<-rep(100, times=nrow(imp_logitdur))
C<-ggplot(imp_logitdur, aes(x=reorder(X, X.IncMSE), y=X.IncMSE)) + 
  geom_bar(position=position_dodge(), stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  # geom_errorbar(aes(ymin=X.IncMSE-SD, ymax=X.IncMSE+SD),
  #               size=.3,    # Thinner lines
  #               width=.2,
  #               position=position_dodge(.9)) +
  scale_x_discrete(element_blank(), labels = c("scent_decay" = expression(delta), "rec_rate"=expression(gamma),
                                               "dur_scent" =  expression(eta[0]),"n.initial" = "N", "inf_prob" = expression(beta),
                                               "initial_load" = expression(kappa[0]), "inf_decay"=expression(alpha))) +
  ylab("Mean Decrease in Accuracy (%IncMSE)") +
  theme_bw()+
  ggtitle("(C)") + #" Variable Importance for Epidemic \nDuration|Success") +
  
  theme(axis.text=element_text(size=14), plot.title = element_text(size = 16)) +
  theme(plot.margin=unit(c(-0.25,0.1,0,0), "cm"))+
  coord_flip()
C
#setwd()
# tiff("Figure1_RF.tiff", height =10 , width =8.7, units = "cm", compression = "lzw", res = 1200)
multiplot(A, B, C, cols=1)
# dev.off()