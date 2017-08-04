
# call from validate wrapper
setwd("ibm/")

sppList=c("ARTR")
Nspp=length(sppList)
sppNames=c("A. tripartita")

# get observed cover
obsCov<-read.csv("../QuadYearCover.csv")
obsCov<-subset(obsCov,species=="Artemisia tripartita" & year>=2001)

# function to plot all simulated trajectories in one quad
plotSims <- function(mydata,maxCov,qName,trtName){
  mydata<-mydata/100 # convert to % cover
  myLims<-c(0,maxCov/100)
  meanCov<-colMeans(mydata)
  cutoffs<-quantile(meanCov,c(0.1,0.9))
  low<-which(meanCov<cutoffs[1])
  high<-which(meanCov>cutoffs[2])
  med<-which(meanCov>=cutoffs[1] & meanCov<=cutoffs[2])
  matplot(mydata,ylim=myLims,type="n",xlab="Year",ylab="Cover",main=paste(qName,trtName))
  par(new=T); matplot(mydata[,med],ylim=myLims,ylab="",xlab="",axes=F,type="l",lty=1,col="grey")
  par(new=T); matplot(mydata[,low],ylim=myLims,ylab="",xlab="",axes=F,type="l",lty=1,col="red")
  par(new=T); matplot(mydata[,high],ylim=myLims,ylab="",xlab="",axes=F,type="l",lty=1,col="blue")
  lines(1:NROW(mydata),rowMeans(mydata),col="black",lty="dashed")
  # add observed cover
  tmp<-subset(obsCov,quad==qName)
  lines(1:6,tmp$cover,col="black",lwd=2,lty=1)
}

# grass removal plots
qList <- c("Q48","Q49","Q51","Q55","Q57","Q58","Q60","Q62")  # no grass
pdf("ARTRtrajectories.pdf",height=10,width=7)
par(mfrow=c(4,2),tcl=-0.1,mgp=c(2,0.5,0),mar=c(3,3,2,1))
for(i in 1:length(qList)){
   infile=paste("simulations/",qList[i],"_ARTRcov_removals_noTrt.csv",sep="")
   tmpD1=read.csv(infile)
   infile=paste("simulations/",qList[i],"_ARTRcov_removals_Trt.csv",sep="")
   tmpD2=read.csv(infile)
   maxCov= max(c(max(tmpD1),max(tmpD2)))
   plotSims(tmpD1,maxCov,qList[i],"Baseline")
   plotSims(tmpD2,maxCov,qList[i],"Removal")
}
dev.off()

setwd("..")

