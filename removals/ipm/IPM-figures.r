
# call from removal_analysis_wrapper.r

# box plots
if(!exists("baseline")) baseline <- read.csv("ipm/baselineCover.csv")
if(!exists("baseline.noARTR")) baseline.noARTR <- read.csv("ipm/baselineCover-noARTR.csv")
if(!exists("removal.noARTR")) removal.noARTR <- read.csv("ipm/removalCover-noARTR.csv")
if(!exists("removal.noARTR.maxCI")) removal.noARTR.maxCI <- read.csv("ipm/removalCover-noARTR-maxCI.csv")

myCols<-c("darkgrey","blue3","red3")
sppNames <-c("A. tripartita","H. comata","P. secunda","P. spicata")

png("ipm/boxplots.png",height=3.5, width=8, units="in",res=400)
  par(tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,1,1),cex.lab=1.2)
  plot(c(1:17),rep(NA,17),ylim=c(0,10),ylab="Cover (%)",xlab="",xaxt="n")
  for(doSpp in 1:4){
    tmp <- cbind(baseline[,doSpp],baseline.noARTR[,doSpp],removal.noARTR[,doSpp])
    boxplot(100*tmp,col=myCols,add=T, at=((2:4)+(doSpp-1)*4),names=rep("",3),xaxt="n")
  }
  axis(side=1,at=c(3,7,11,15),labels=sppNames,font=4)
  legend("topright",c("Baseline model","Baseline model, no ARTR","Removal model, no ARTR"),
         fill=myCols,bty="n",cex=0.9)
dev.off()

# same figure as previous, but add max treatment effect simulation
png("ipm/boxplots-maxCI.png",height=3.5, width=8, units="in",res=400)
  par(tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,1,1),cex.lab=1.2)
  plot(c(1:17),rep(NA,17),ylim=c(0,10),ylab="Cover (%)",xlab="",xaxt="n")
  for(doSpp in 1:4){
    tmp <- cbind(baseline[,doSpp],baseline.noARTR[,doSpp],removal.noARTR[,doSpp],removal.noARTR.maxCI[,doSpp])
    boxplot(100*tmp,col=c(myCols,"orange"),add=T, at=(seq(1:4)+(doSpp-1)*4),names=rep("",4),xaxt="n")
  }
  axis(side=1,at=c(3,7,11,15),labels=sppNames,font=4)
  abline(v=4.5);abline(v=8.5);abline(v=12.5)
  legend("topright",c("Baseline model","Baseline model, no ARTR","Removal model, no ARTR","Removal model, no ARTR, max CI"),
         fill=c(myCols,"orange"),bg="white",cex=0.9)
dev.off()
