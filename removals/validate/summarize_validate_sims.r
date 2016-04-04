
# call from validate wrapper
setwd("validate/")

sppList=c("ARTR","HECO","POSE","PSSP")
Nspp=length(sppList)
sppNames=c("A. tripartita","H. comata","Poa secunda","P. spicata")
myCol=c("black","forestgreen","blue","red")

# control plots
qList <- c("Q1","Q2","Q3","Q4","Q5","Q6" ) 
covD=NULL
for(i in 1:length(qList)){
   infile=paste(qList[i],"_validation_cov_removals.csv",sep="")
   tmpD=read.csv(infile)
   covD=rbind(covD,tmpD)
}
control.mean=aggregate(covD[,2:NCOL(covD)],by=list(year=covD$year),FUN=mean,na.rm=T)
control.sd=aggregate(covD[,2:NCOL(covD)],by=list(year=covD$year),FUN=sd,na.rm=T)

# no shrub plots
qList <- c("Q47","Q50","Q52","Q53","Q54","Q56","Q59","Q61")
covD=NULL
for(i in 1:length(qList)){
   infile=paste(qList[i],"_validation_cov_removals.csv",sep="")
   tmpD=read.csv(infile)
   covD=rbind(covD,tmpD)
}
covD$ARTR=NA ; covD$ARTRpred = NA # get rid of ARTR so grass trends are easier to see
noshrub.mean=aggregate(covD[,2:NCOL(covD)],by=list(year=covD$year),FUN=mean,na.rm=T)
noshrub.sd=aggregate(covD[,2:NCOL(covD)],by=list(year=covD$year),FUN=sd,na.rm=T)

# no grass plots
qList <- c("Q48","Q49","Q51","Q55","Q57","Q58","Q60","Q62")  # no grass
covD=NULL
for(i in 1:length(qList)){
   infile=paste(qList[i],"_validation_cov_removals.csv",sep="")
   tmpD=read.csv(infile)
   covD=rbind(covD,tmpD)
}
nograss.mean=aggregate(covD[,2:NCOL(covD)],by=list(year=covD$year),FUN=mean,na.rm=T)
nograss.sd=aggregate(covD[,2:NCOL(covD)],by=list(year=covD$year),FUN=sd,na.rm=T)


#set up plotting function
plotObsPred<-function(mydata,mytitle){
  matplot(mydata[,1],mydata[,2:NCOL(mydata)]/100,type="n",
    xaxt="n",yaxt="n",xlab="",ylab="")
  par(new=T)
  matplot(mydata[,1],mydata[,2:NCOL(mydata)]/100,type="o",
    col=myCol,lty=c(rep("solid",Nspp),c(rep("dashed",Nspp))),
    pch=c(rep(16,Nspp),rep(1,Nspp)),xlab="",ylab="",yaxt="n")
  axis(2,las=1)
  title(main=mytitle,adj=0,font.main=1)  
}


png("obsVSpred_v2.png",units="in",height=4.5,width=9.5,res=600)
  
  par(mfrow=c(1,3),tcl=-0.2,mgp=c(2,0.5,0),mar=c(2,2,2,1),oma=c(2,2,0,0))
  
  plotObsPred(control.mean,"Control")
  par(font=3)
  legend(2013,15,sppNames,col=myCol,lty="solid",pch=16,bty="n")
  par(font=1)
  
  plotObsPred(noshrub.mean,"No ARTR")
  
  plotObsPred(nograss.mean,"No grass")
  
  mtext(side=1,"Year",line=0.5, outer=T)
  mtext(side=2,"Mean cover (%)",line=0.5, outer=T)

dev.off()

 