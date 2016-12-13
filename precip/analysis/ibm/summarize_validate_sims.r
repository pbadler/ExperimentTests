
# call from validate wrapper
setwd("ibm/")

sppList=c("ARTR","HECO","POSE","PSSP")
Nspp=length(sppList)
sppNames=c("A. tripartita","H. comata","Poa secunda","P. spicata")
myCol=c("black","forestgreen","blue","red")

# control plots
qList <- paste0("Q",c(1:6,19:26))  #c("Q1","Q2","Q3","Q4","Q5","Q6" ) 
covD=NULL
for(i in 1:length(qList)){
   infile=paste("simulations1step/",qList[i],"_validation_cov_removals_noTrt.csv",sep="")
   tmpD=read.csv(infile)
   covD=rbind(covD,tmpD)
}
control.mean=aggregate(covD[,2:NCOL(covD)],by=list(year=covD$year),FUN=mean,na.rm=T)
control.sd=aggregate(covD[,2:NCOL(covD)],by=list(year=covD$year),FUN=sd,na.rm=T)

# no shrub plots, no treatment effects
qList <- c("Q47","Q50","Q52","Q53","Q54","Q56","Q59","Q61")
covD=NULL
for(i in 1:length(qList)){
   infile=paste("simulations1step/",qList[i],"_validation_cov_removals_noTrt.csv",sep="")
   tmpD=read.csv(infile)
   covD=rbind(covD,tmpD)
}
covD$ARTR=NA ; covD$ARTRpred = NA # get rid of ARTR so grass trends are easier to see
noshrub.mean=aggregate(covD[,2:NCOL(covD)],by=list(year=covD$year),FUN=mean,na.rm=T)
noshrub.sd=aggregate(covD[,2:NCOL(covD)],by=list(year=covD$year),FUN=sd,na.rm=T)

# no shrub plots, WITH treatment effects
qList <- c("Q47","Q50","Q52","Q53","Q54","Q56","Q59","Q61")
covD=NULL
for(i in 1:length(qList)){
   infile=paste("simulations1step/",qList[i],"_validation_cov_removals_Trt.csv",sep="")
   tmpD=read.csv(infile)
   covD=rbind(covD,tmpD)
}
covD$ARTR=NA ; covD$ARTRpred = NA # get rid of ARTR so grass trends are easier to see
noshrubTRT.mean=aggregate(covD[,2:NCOL(covD)],by=list(year=covD$year),FUN=mean,na.rm=T)
noshrubTRT.sd=aggregate(covD[,2:NCOL(covD)],by=list(year=covD$year),FUN=sd,na.rm=T)

# no grass plots, no treatment effects
qList <- c("Q48","Q49","Q51","Q55","Q57","Q58","Q60","Q62")  # no grass
covD=NULL
for(i in 1:length(qList)){
   infile=paste("simulations1step/",qList[i],"_validation_cov_removals_noTrt.csv",sep="")
   tmpD=read.csv(infile)
   covD=rbind(covD,tmpD)
}
nograss.mean=aggregate(covD[,2:NCOL(covD)],by=list(year=covD$year),FUN=mean,na.rm=T)
nograss.sd=aggregate(covD[,2:NCOL(covD)],by=list(year=covD$year),FUN=sd,na.rm=T)

# no grass plots WITH treatment effects
qList <- c("Q48","Q49","Q51","Q55","Q57","Q58","Q60","Q62")  # no grass
covD=NULL
for(i in 1:length(qList)){
   infile=paste("simulations1step/",qList[i],"_validation_cov_removals_Trt.csv",sep="")
   tmpD=read.csv(infile)
   covD=rbind(covD,tmpD)
}
nograssTRT.mean=aggregate(covD[,2:NCOL(covD)],by=list(year=covD$year),FUN=mean,na.rm=T)
nograssTRT.sd=aggregate(covD[,2:NCOL(covD)],by=list(year=covD$year),FUN=sd,na.rm=T)


#set up plotting function
plotObsPred<-function(doSpp,mydata1,mydata2,mydata3,mytitle){
  # format data
  newD=data.frame(mydata1$year,mydata1[,1+doSpp],mydata1[,5+doSpp],  # control obs and pred
                  mydata2[,1+doSpp],mydata2[,5+doSpp], # removal obs and pred (no TRT effect)
                  mydata3[,5+doSpp])                               # removal pred (with TRT effect)
  names(newD)=c("year","control.obs","control.pred","remove.obs","remove.pred","remove.predTRT")
  matplot(newD$year,newD[,2:6]/100,type="o",xlab="",ylab="",
    col=c(rep("black",2),rep("blue",3)),
    pch=c(16,1,16,1,2),
    lty=c("solid","dashed","solid","dashed","dashed"))
  title(main=mytitle,adj=0,font.main=1)  
}

png("obsVSpred_project1step.png",units="in",height=3.5,width=8.5,res=600)
  
  par(mfrow=c(1,4),tcl=-0.2,mgp=c(2,0.5,0),mar=c(2,2,2,1),oma=c(2,2,0,0))
  
  plotObsPred(1,control.mean,nograss.mean,nograssTRT.mean,"ARTR")
  plotObsPred(2,control.mean,noshrub.mean,noshrubTRT.mean,"HECO")
  plotObsPred(3,control.mean,noshrub.mean,noshrubTRT.mean,"POSE")
  plotObsPred(4,control.mean,noshrub.mean,noshrubTRT.mean,"PSSP")
  
  mtext(side=1,"Year",line=0.5, outer=T)
  mtext(side=2,"Mean cover (%)",line=0.5, outer=T)

dev.off()

setwd("..")
 