
# call from validate wrapper
setwd("ibm/")

sppList=c("ARTR","HECO","POSE","PSSP")
Nspp=length(sppList)
sppNames=c("A. tripartita","H. comata","Poa secunda","P. spicata")
myCol=c("black","forestgreen","blue","red")

simD <- read.csv("simulations1step/ObsPred_1step.csv")

quad.info <- unique(simD[,c("quad","Treatment","Group")],MARGIN=2)

### calculate mean cover by treatment

covMeans <- aggregate(simD[,3:14],by=list(Treatment=simD$Treatment,year=simD$year),
    na.rm=T,FUN=mean)
covMeans[1:3,7:14] <- covMeans[1:3,3:6] # copy initial values across
covMeans <- covMeans[order(covMeans$Treatment),]

### calculate per capita growth rates

get.pgr <- function(myname){
  tmp<- simD[,c("quad","year",paste0(myname,sppList))]
  tmp$year <- tmp$year-1
  names(tmp)[3:6] <-paste0("next.",sppList)
  tmp2 <- merge(simD[,c(1:6)],tmp)
  out <- cbind(tmp2[,c("quad","year")],log(tmp2[,7:10]/tmp2[,3:6]))
  out[out==Inf | out==-Inf] <- NA
  out
}

obs.pgr <- get.pgr("obs.")
pred.pgr <- get.pgr("pred.")
pred.trt.pgr <- get.pgr("pred.trt.")

# aggregate by treatment

get.trt.means<-function(mydat){
  mydat<-merge(mydat,quad.info)
  out <- aggregate(mydat[,3:6],by=list(Treatment=mydat$Treatment,year=mydat$year),
                        FUN=mean,na.rm=T)
  out <- out[order(out$Treatment,out$year),]
  #consolidate removal treatments
  out[6:10,4:6] <- out[11:15,4:6]
  out <- out[-c(11:15),]
  out$Treatment <- as.character(out$Treatment)
  out$Treatment[6:10] <- "Removal"
  out
}

obs.pgr.mean <- get.trt.means(obs.pgr)
pred.pgr.mean <- get.trt.means(pred.pgr)
pred.trt.pgr.mean <- get.trt.means(pred.trt.pgr)

###
### plot observed and predicted cover chronologically
###

png("cover_projections_1step.png",res=400,width=8.5,height=3,units="in")

par(mfrow=c(1,4),tcl=0.2,mgp=c(2,0.5,0),mar=c(2,2,2,1),oma=c(2,2,0,0))

for(i in 1:4){
  if(i==1){
    doRows <- which(covMeans$Treatment=="No_grass")
  }else{
    doRows <- which(covMeans$Treatment=="No_shrub")
  }  
  matplot(covMeans$year[1:6],cbind(covMeans[1:6,2+i],covMeans[1:6,6+i], # control plots
          covMeans[doRows,2+i],covMeans[doRows,6+i],covMeans[doRows,10+i]),
          xlab="",ylab="",type="l",
          lty=c("solid","dashed","solid","dashed","dotted"),
          col=c(rep("blue3",2),rep("red3",3)))   # removal plots
  title(main=sppNames[i],adj=0,font.main=4) 
  if(i==1){
    legend("topleft",c("Control (observed)","Removal (observed)","Baseline (predicted)","Removal (predicted)"),
    col=c("blue3","red3","darkgray","darkgray"), 
    lty=c("solid","solid","dashed","dotted"),bty="n")
  }
}

mtext("Year",side=1,line=0.5,outer=T,cex=1)
mtext("Cover (%)",side=2,line=0.5,outer=T,cex=1)

dev.off()

###
### plot growth rates chronologically 
###

plotObsPred<-function(doSpp,mytitle,doLegend=F){
  
  # format data
  newD=data.frame(year=2011:2015,obs.pgr.mean[obs.pgr.mean$Treatment=="Control",2 + doSpp],
                 pred.pgr.mean[pred.pgr.mean$Treatment=="Control",2+doSpp], 
                 obs.pgr.mean[obs.pgr.mean$Treatment=="Removal",2+doSpp],
                 pred.pgr.mean[pred.pgr.mean$Treatment=="Removal",2+doSpp],
                 pred.trt.pgr.mean[pred.trt.pgr.mean$Treatment=="Removal",2+doSpp])                               # removal pred (with TRT effect)
  names(newD)=c("year","control.obs","control.pred","remove.obs","remove.pred","remove.predTRT")
  
  color1=rgb(0,100,255,alpha=175,maxColorValue = 255)
  color2=rgb(153,0,0,alpha=175,maxColorValue = 255)
  my.y <- c(-1.2,1.1) # hard wire ylims
  matplot(newD$year,newD[,2:6],type="o",xlab="",ylab="",ylim=my.y,
    col=c(rep(color1,2),rep(color2,3)),xaxt="n",
    pch=c(16,21,16,21,24),bg="white",
    lty=c("solid","dotted","solid","dotted","dotted"))
  axis(1,at=c(2011:2015))
  abline(h=0,lty="solid",col="darkgray")
  # add standard error bars to observed means
#   arrows(x0=mysd1$year,y0=c(mydata1[,1+doSpp]-mysd1[,1+doSpp]/sqrt(14)), # 14 = number of control plots
#          x1=mysd1$year,y1=c(mydata1[,1+doSpp]+mysd1[,1+doSpp]/sqrt(14)),length=0.05,angle=90,code=3,col=color1)  
#   arrows(x0=mysd2$year,y0=c(mydata2[,1+doSpp]-mysd2[,1+doSpp]/sqrt(8)), # 8 = number of control plots
#          x1=mysd2$year,y1=c(mydata2[,1+doSpp]+mysd2[,1+doSpp]/sqrt(8)),length=0.05,angle=90,code=3,col=color2)  
  title(main=mytitle,adj=0,font.main=4)  
  if(doLegend==T){
    legend("bottomleft",c("Control","Baseline model","Removal","Baseline model","Removal model"),
    col=c(rep(color1,2),rep(color2,3)), pch=c(16,21,16,21,24),pt.bg = "white",
    lty=c("solid","dotted","solid","dotted","dotted"),bty="n")
  }
}

png("cover_change_chrono.png",units="in",height=3,width=8.5,res=600)
  
  par(mfrow=c(1,4),tcl=-0.2,mgp=c(2,0.5,0),mar=c(2,2,2,1),oma=c(2,2,0,0))
  
  plotObsPred(1,sppNames[1],doLegend=T)
  plotObsPred(2,sppNames[2])
  plotObsPred(3,sppNames[3])
  plotObsPred(4,sppNames[4])
  
  mtext(side=1,"Year",line=0.5, outer=T)
  mtext(side=2,expression(paste("Mean ",log(Cover[t+1]/Cover[t]))),line=0, outer=T)

dev.off()

setwd("..")
 