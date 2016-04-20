
# call from validate wrapper
setwd("validate/")

sppList=c("ARTR","HECO","POSE","PSSP")
Nspp=length(sppList)
sppNames=c("A. tripartita","H. comata","Poa secunda","P. spicata")
myCol=c("black","forestgreen","blue","red")

# function to load observed and predicted cover
getSims<-function(qList,trtEffects){
   out=NULL
   for(i in qList){
     if(trtEffects==F){
       infile=paste("simulations1step/",i,"_validation_cov_removals_noTrt.csv",sep="")
     }else{
       infile=paste("simulations1step/",i,"_validation_cov_removals_Trt.csv",sep="")
     }
     tmpD=read.csv(infile)
     out <- rbind(out,tmpD)
   }
   return(out)
}

# function to load data and calculate log(cover_t+1 / cover_t)
getPopGrowth<-function(qList,trtEffects){
   out=NULL
   for(i in qList){
     if(trtEffects==F){
       infile=paste("simulations1step/",i,"_validation_cov_removals_noTrt.csv",sep="")
     }else{
       infile=paste("simulations1step/",i,"_validation_cov_removals_Trt.csv",sep="")
     }
     tmpD=read.csv(infile)
     # make sure that the denominator (t0) is always the observed cover, not predicted cover
     change=data.frame(year=tmpD$year[1:(NROW(tmpD)-1)],
            log(tmpD[2:NROW(tmpD),2:NCOL(tmpD)]/cbind(tmpD[1:(NROW(tmpD)-1),2:5],tmpD[1:(NROW(tmpD)-1),2:5])))
     change[change==Inf]<-NA; change[change==-Inf]<-NA
     out=rbind(out,change)
   }
   return(out)
}

# control plots
qList <- paste0("Q",c(1:6,19:26))  # all contemporary control plots#
#qList <- c("Q1","Q2","Q3","Q4","Q5","Q6" ) # just Group 1 
quadD <- getSims(qList,trtEffects=F)
control.cov <- aggregate(quadD[,2:NCOL(quadD)],by=list(year=quadD$year),FUN=mean,na.rm=T)
quadD <- getPopGrowth(qList,trtEffects=F)
control.pgr <- aggregate(quadD[,2:NCOL(quadD)],by=list(year=quadD$year),FUN=mean,na.rm=T)
#control.sd=aggregate(quadD[,2:NCOL(quadD)],by=list(year=quadD$year),FUN=sd,na.rm=T)

# no shrub plots, no treatment effects
qList <- c("Q47","Q50","Q52","Q53","Q54","Q56","Q59","Q61")
quadD <- getSims(qList,trtEffects=F)
noshrub.cov <- aggregate(quadD[,2:NCOL(quadD)],by=list(year=quadD$year),FUN=mean,na.rm=T)
quadD <- getPopGrowth(qList,trtEffects=F)
noshrub.pgr <- aggregate(quadD[,2:NCOL(quadD)],by=list(year=quadD$year),FUN=mean,na.rm=T)

# no shrub plots, WITH treatment effects
qList <- c("Q47","Q50","Q52","Q53","Q54","Q56","Q59","Q61")
quadD <- getSims(qList,trtEffects=T)
noshrubTRT.cov <- aggregate(quadD[,2:NCOL(quadD)],by=list(year=quadD$year),FUN=mean,na.rm=T)
quadD <- getPopGrowth(qList,trtEffects=T)
noshrubTRT.pgr <- aggregate(quadD[,2:NCOL(quadD)],by=list(year=quadD$year),FUN=mean,na.rm=T)

# no grass plots, no treatment effects
qList <- c("Q48","Q49","Q51","Q55","Q57","Q58","Q60","Q62")  # no grass
quadD <- getSims(qList,trtEffects=F)
nograss.cov <- aggregate(quadD[,2:NCOL(quadD)],by=list(year=quadD$year),FUN=mean,na.rm=T)
quadD <- getPopGrowth(qList,trtEffects=F)
nograss.pgr <- aggregate(quadD[,2:NCOL(quadD)],by=list(year=quadD$year),FUN=mean,na.rm=T)

# no grass plots, no treatment effects
qList <- c("Q48","Q49","Q51","Q55","Q57","Q58","Q60","Q62")  # no grass
quadD <- getSims(qList,trtEffects=T)
nograssTRT.cov <- aggregate(quadD[,2:NCOL(quadD)],by=list(year=quadD$year),FUN=mean,na.rm=T)
quadD <- getPopGrowth(qList,trtEffects=T)
nograssTRT.pgr <- aggregate(quadD[,2:NCOL(quadD)],by=list(year=quadD$year),FUN=mean,na.rm=T)

###
### plot observed and predicted cover chronologically
###
#set up plotting function

plotTimeSeries<-function(doSpp,mydata1,mydata2,mydata3,mytitle){
  # choose colors
  col1="dodgerblue3"
  col2="red3"
  # format data
  newD=data.frame(mydata1$year,mydata1[,1+doSpp],mydata1[,5+doSpp],  # control obs and pred
                  mydata2[,1+doSpp],mydata2[,5+doSpp], # removal obs and pred (no TRT effect)
                  mydata3[,5+doSpp])                               # removal pred (with TRT effect)
  names(newD)=c("year","control.obs","control.pred","remove.obs","remove.pred","remove.predTRT")
  matplot(newD$year,newD[,2:6]/100,type="o",xlab="",ylab="",
    col=c(rep(col1,2),rep(col2,3)),
    pch=c(16,1,16,1,2),
    lty=c("solid","blank","solid","blank","blank"))
  #line segments for controls
  segments(x0=newD$year[1:4],y0=newD$control.obs[1:4]/100,
           x1=newD$year[2:5],y1=newD$control.pred[2:5]/100,col=col1,lty="dashed")
  #line segments for removals, no TRT
  segments(x0=newD$year[1:4],y0=newD$remove.obs[1:4]/100,
           x1=newD$year[2:5],y1=newD$remove.pred[2:5]/100,col=col2,lty="dashed")
  #line segments for removals, no TRT
  segments(x0=newD$year[1:4],y0=newD$remove.obs[1:4]/100,
           x1=newD$year[2:5],y1=newD$remove.predTRT[2:5]/100,col=col2,lty="dashed")
  title(main=mytitle,adj=0,font.main=1)  
}

png("obsVSpred_project1step.png",units="in",height=3.5,width=8.5,res=600)
  
  par(mfrow=c(1,4),tcl=-0.2,mgp=c(2,0.5,0),mar=c(2,2,2,1),oma=c(2,2,0,0))
  
  plotTimeSeries(1,control.cov,nograss.cov,nograssTRT.cov,"ARTR")
  plotTimeSeries(2,control.cov,noshrub.cov,noshrubTRT.cov,"HECO")
  plotTimeSeries(3,control.cov,noshrub.cov,noshrubTRT.cov,"POSE")
  plotTimeSeries(4,control.cov,noshrub.cov,noshrubTRT.cov,"PSSP")
  
  mtext(side=1,"Year",line=0.5, outer=T)
  mtext(side=2,"Mean cover (%)",line=0.5, outer=T)

dev.off()

###
### plot observations vs predictions 1:1
###

plot1to1 <- function(doSpp,controlD,removalD,removalTRTD,mytitle,doLegend=F){
  #tmp <- rbind(controlD[,2:NCOL(controlD)],removalD[,2:NCOL(removalD)],removalTRTD[,2:NCOL(removalTRTD)])
  #mylims <- c(0.95*min(tmp,na.rm=T),1.05*max(tmp,na.rm=T))
  color1=rgb(0,100,255,alpha=175,maxColorValue = 255)
  color2=rgb(153,0,0,alpha=175,maxColorValue = 255)
  mylims <- c(-1.2,0.9)
  plot(controlD[,1+doSpp],controlD[,5+doSpp],xlim=mylims,ylim=mylims,
       xlab="",ylab="",type="n")
  abline(0,1)
  abline(h=0,lty="dotted",col="grey")
  abline(v=0,lty="dotted",col="grey")
  points(removalD[,1+doSpp],removalD[,5+doSpp],pch=22,bg=color2,cex=1.5)
  points(removalTRTD[,1+doSpp],removalTRTD[,5+doSpp],pch=24,bg=color2,cex=1.5)
  points(controlD[,1+doSpp],controlD[,5+doSpp],pch=21,bg=color1,cex=1.5)
  if(doLegend==T){
    legend("topleft",c("Control plots","Removal plots","Removal plots + effects"),pch=c(21,22,24),pt.cex=1.5,
        pt.bg=c(color1,color2,color2),bty="n")
  }
  title(main=mytitle,adj=0,font.main=1) 
}

png("obsVSpred1to1.png",units="in",height=2.75,width=8.5,res=600)
  par(mfrow=c(1,4),tcl=-0.2,mgp=c(2,0.5,0),mar=c(2,2,2,1),oma=c(2,2,0,0))
  
  plot1to1(1,control.pgr,nograss.pgr,nograssTRT.pgr,"ARTR",doLegend=T)
  plot1to1(2,control.pgr,noshrub.pgr,noshrubTRT.pgr,"HECO")
  plot1to1(3,control.pgr,noshrub.pgr,noshrubTRT.pgr,"POSE")
  plot1to1(4,control.pgr,noshrub.pgr,noshrubTRT.pgr,"PSSP")
  
  mtext(side=1,"Observed",line=0.5, outer=T)
  mtext(side=2,"Predicted",line=0.5, outer=T)

dev.off()


###
### plot growth rates chronologically 
###

plotObsPred<-function(doSpp,mydata1,mydata2,mydata3,mytitle){
  # format data
  newD=data.frame(mydata1$year,mydata1[,1+doSpp],mydata1[,5+doSpp],  # control obs and pred
                  mydata2[,1+doSpp],mydata2[,5+doSpp], # removal obs and pred (no TRT effect)
                  mydata3[,5+doSpp])                               # removal pred (with TRT effect)
  names(newD)=c("year","control.obs","control.pred","remove.obs","remove.pred","remove.predTRT")
  matplot(newD$year,newD[,2:6],type="o",xlab="",ylab="",
    col=c(rep("black",2),rep("blue",3)),
    pch=c(16,1,16,1,2),
    lty=c("solid","blank","solid","blank","blank"))
  title(main=mytitle,adj=0,font.main=1)  
}

png("cover_change_chrono.png",units="in",height=3.5,width=8.5,res=600)
  
  par(mfrow=c(1,4),tcl=-0.2,mgp=c(2,0.5,0),mar=c(2,2,2,1),oma=c(2,2,0,0))
  
  plotObsPred(1,control.pgr,nograss.pgr,nograssTRT.pgr,"ARTR")
  plotObsPred(2,control.pgr,noshrub.pgr,noshrubTRT.pgr,"HECO")
  plotObsPred(3,control.pgr,noshrub.pgr,noshrubTRT.pgr,"POSE")
  plotObsPred(4,control.pgr,noshrub.pgr,noshrubTRT.pgr,"PSSP")
  
  mtext(side=1,"Year",line=0.5, outer=T)
  mtext(side=2,"log Cover change",line=0.5, outer=T)

dev.off()

setwd("..")
 