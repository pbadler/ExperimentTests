
# call from validate wrapper
setwd("validate/")

sppList=c("ARTR","HECO","POSE","PSSP")
Nspp=length(sppList)
sppNames=c("A. tripartita","H. comata","Poa secunda","P. spicata")
myCol=c("black","forestgreen","blue","red")

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
     change=data.frame(year=tmpD$year[1:(NROW(tmpD)-1)],log(tmpD[2:NROW(tmpD),2:NCOL(tmpD)]/tmpD[1:(NROW(tmpD)-1),2:NCOL(tmpD)]))
     change[change==Inf]<-NA; change[change==-Inf]<-NA
     out=rbind(out,change)
   }
   return(out)
}

# control plots
qList <- c("Q1","Q2","Q3","Q4","Q5","Q6" ) 
quadD <- getPopGrowth(qList,trtEffects==F)
control.mean=aggregate(quadD[,2:NCOL(quadD)],by=list(year=quadD$year),FUN=mean,na.rm=T)
control.sd=aggregate(quadD[,2:NCOL(quadD)],by=list(year=quadD$year),FUN=sd,na.rm=T)

# no shrub plots, no treatment effects
qList <- c("Q47","Q50","Q52","Q53","Q54","Q56","Q59","Q61")
quadD <- getPopGrowth(qList,trtEffects==F)
noshrub.mean=aggregate(quadD[,2:NCOL(quadD)],by=list(year=quadD$year),FUN=mean,na.rm=T)
noshrub.sd=aggregate(quadD[,2:NCOL(quadD)],by=list(year=quadD$year),FUN=sd,na.rm=T)

# no shrub plots, WITH treatment effects
qList <- c("Q47","Q50","Q52","Q53","Q54","Q56","Q59","Q61")
quadD <- getPopGrowth(qList,trtEffects==T)
noshrubTRT.mean=aggregate(quadD[,2:NCOL(quadD)],by=list(year=quadD$year),FUN=mean,na.rm=T)
noshrubTRT.sd=aggregate(quadD[,2:NCOL(quadD)],by=list(year=quadD$year),FUN=sd,na.rm=T)

# no grass plots, no treatment effects
qList <- c("Q48","Q49","Q51","Q55","Q57","Q58","Q60","Q62")  # no grass
quadD <- getPopGrowth(qList,trtEffects==F)
nograss.mean=aggregate(quadD[,2:NCOL(quadD)],by=list(year=quadD$year),FUN=mean,na.rm=T)
nograss.sd=aggregate(quadD[,2:NCOL(quadD)],by=list(year=quadD$year),FUN=sd,na.rm=T)

# no grass plots, no treatment effects
qList <- c("Q48","Q49","Q51","Q55","Q57","Q58","Q60","Q62")  # no grass
quadD <- getPopGrowth(qList,trtEffects==T)
nograssTRT.mean=aggregate(quadD[,2:NCOL(quadD)],by=list(year=quadD$year),FUN=mean,na.rm=T)
nograssTRT.sd=aggregate(quadD[,2:NCOL(quadD)],by=list(year=quadD$year),FUN=sd,na.rm=T)


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
    lty=c("solid","blank","solid","blank","blank"))
  title(main=mytitle,adj=0,font.main=1)  
}


png("obsVSpred1step.png",units="in",height=3.5,width=8.5,res=600)
  
  par(mfrow=c(1,4),tcl=-0.2,mgp=c(2,0.5,0),mar=c(2,2,2,1),oma=c(2,2,0,0))
  
  plotObsPred(1,control.mean,nograss.mean,nograssTRT.mean,"ARTR")
  plotObsPred(2,control.mean,noshrub.mean,noshrubTRT.mean,"HECO")
  plotObsPred(3,control.mean,noshrub.mean,noshrubTRT.mean,"POSE")
  plotObsPred(4,control.mean,noshrub.mean,noshrubTRT.mean,"PSSP")
  
  mtext(side=1,"Year",line=0.5, outer=T)
  mtext(side=2,"Mean cover (%)",line=0.5, outer=T)

dev.off()

setwd("..")
 