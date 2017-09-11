
# call from validate wrapper
setwd("ibm/")

sppList=c("ARTR","HECO","POSE","PSSP")
Nspp=length(sppList)
sppNames=c("A. tripartita","H. comata","P. secunda","P. spicata")
myCol=c("black","forestgreen","blue","red")

infile <- ifelse(max.CI==F,"simulations1step/ObsPred_1step.csv","simulations1step/ObsPred_1step_maxCI.csv")
simD <- read.csv(infile)

quad.info <- unique(simD[,c("quad","Treatment","Group")],MARGIN=2)

### calculate mean cover by treatment

covMeans <- aggregate(simD[,3:14],by=list(Treatment=simD$Treatment,year=simD$year),
    na.rm=T,FUN=mean)
covMeans[1:3,7:14] <- covMeans[1:3,3:6] # copy initial values across
covMeans <- covMeans[order(covMeans$Treatment),]

### calculate standard deviations by treatment

covSDs <- aggregate(simD[,3:14],by=list(Treatment=simD$Treatment,year=simD$year),
    na.rm=T,FUN=sd)
covSDs[1:3,7:14] <- covSDs[1:3,3:6] # copy initial values across
covSDs <- covSDs[order(covSDs$Treatment),]

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

# aggregate growth rate (means) by treatment and year

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
color1="black"
color2="blue2" #rgb(0,100,255,alpha=175,maxColorValue = 255)
color3="red" #rgb(153,0,0,alpha=175,maxColorValue = 255)

figName <- ifelse(max.CI==F,"cover_projections_1step.png","cover_projections_1step_maxCI.png" )
png(figName,res=400,width=8.5,height=6,units="in")

par(mfrow=c(2,2),tcl=0.2,mgp=c(2,0.5,0),mar=c(2,2,2,1),oma=c(2,2,0,0))

for(i in 1:4){
  if(i==1){
    doRows <- which(covMeans$Treatment=="No_grass")
  }else{
    doRows <- which(covMeans$Treatment=="No_shrub")
  }  
  matplot(covMeans$year[1:6],cbind(covMeans[1:6,2+i],covMeans[1:6,6+i], # control plots
          covMeans[doRows,2+i],covMeans[doRows,6+i],covMeans[doRows,10+i]),
          xlab="",ylab="",type="o",
          col=c(color1,color2,color1,color2,color3),
          pch=c(15,15,21,21,21),bg="white", cex=1.1, 
          lty=c("solid","dotted","solid","dotted","dotted"))   # removal plots
  title(main=sppNames[i],adj=0,font.main=4) 
  if(i==1){
    legend("top",c("Control (obs.)","Control (baseline pred.)","Removal (obs.)","Removal (baseline pred.)","Removal (removal pred.)"),
    col=c(color1,color2,color1,color2,color3),
          pch=c(15,15,21,21,21),pt.bg="white",cex=0.9,
          lty=c("solid","dotted","solid","dotted","dotted"),bty="n")
  }
}

mtext("Year",side=1,line=0.5,outer=T,cex=1.1)
mtext("Cover (%)",side=2,line=0.5,outer=T,cex=1.1)

dev.off()

###
### plot observed and predicted cover chronologically WITH ERROR BARS
###
color1="black"
color2="blue2" #rgb(0,100,255,alpha=175,maxColorValue = 255)
color3="red" #rgb(153,0,0,alpha=175,maxColorValue = 255)

figName <- ifelse(max.CI==F,"cover_projections_1stepBARS.png","cover_projections_1stepBARS_maxCI.png" )
png(figName,res=400,width=8.5,height=6,units="in")

par(mfrow=c(2,2),tcl=0.2,mgp=c(2,0.5,0),mar=c(2,2,2,1),oma=c(2,2,0,0))

for(i in 1:4){
  if(i==1){
    doRows <- which(covMeans$Treatment=="No_grass")
  }else{
    doRows <- which(covMeans$Treatment=="No_shrub")
  }  
  mycol=c(color1,color2,color1,color2,color3)
  mypch=c(15,15,21,21,21)
  mylty=c("solid","dotted","solid","dotted","dotted")
  colindex=c(2+i,6+i,2+i,6+i,10+i)
  my.ylims=c(min(covMeans[,colindex]-covSDs[,colindex]),max(covMeans[,colindex]+covSDs[,colindex]))
  matplot(covMeans$year[1:6],cbind(covMeans[1:6,2+i],covMeans[1:6,6+i], # control plots
          covMeans[doRows,2+i],covMeans[doRows,6+i],covMeans[doRows,10+i]),
          ylim=my.ylims,xlab="",ylab="",type="o",
          col=mycol,
          pch=mypch,bg="white", cex=1.1, 
          lty=mylty)   # removal plots
  # add error bars

  for(j in 1:5){
    if(j<3) {rowindex = 1:6} else {rowindex=doRows}
    arrows(x0=covMeans$year[1:6],y0=covMeans[rowindex,colindex[j]]-covSDs[rowindex,colindex[j]],
           x1=covMeans$year[1:6],y1=covMeans[rowindex,colindex[j]]+covSDs[rowindex,colindex[j]],
           code=3,length=0.03,angle=90,col=mycol[j],lty=mylty[j])
  }
  title(main=sppNames[i],adj=0,font.main=4) 
  if(i==1){
    legend("top",c("Control (obs.)","Control (baseline pred.)","Removal (obs.)","Removal (baseline pred.)","Removal (removal pred.)"),
    col=c(color1,color2,color1,color2,color3),
          pch=c(15,15,21,21,21),pt.bg="white",cex=0.9,
          lty=c("solid","dotted","solid","dotted","dotted"),bty="n")
  }
}

mtext("Year",side=1,line=0.5,outer=T,cex=1.1)
mtext("Cover (%)",side=2,line=0.5,outer=T,cex=1.1)

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
  
  my.y <- c(-1.2,1.1) # hard wire ylims
  matplot(newD$year,newD[,2:6],type="o",xlab="",ylab="",ylim=my.y,
    col=c(color1,color2,color1,color2,color3),xaxt="n",
    pch=c(15,15,21,21,21), cex=1.1, bg="white",
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
    legend("topright",c("Control (obs.)","Control (baseline pred.)","Removal (obs.)","Removal (baseline pred.)","Removal (removal pred.)"),
    col=c(color1,color2,color1,color2,color3),
    pch=c(15,15,21,21,21), pt.bg="white",
    lty=c("solid","dotted","solid","dotted","dotted"),bty="n")
  }
}

figName <- ifelse(max.CI==F,"cover_change_chrono.png","cover_change_chrono_maxCI.png" )
png(figName,units="in",height=6,width=8.5,res=600)
  
  par(mfrow=c(2,2),tcl=-0.2,mgp=c(2,0.5,0),mar=c(2,2,2,1),oma=c(2,2,0,0))
  
  plotObsPred(1,sppNames[1])
  plotObsPred(2,sppNames[2])
  plotObsPred(3,sppNames[3],doLegend=T)
  plotObsPred(4,sppNames[4])
  
  mtext(side=1,"Year",line=0.5, outer=T,cex=1.1)
  mtext(side=2,expression(paste("Mean ",log(Cover[t+1]/Cover[t]))),line=0, outer=T,cex=1.1)

dev.off()


###
### plot quadrat level observed vs predicted cover
###

figName <- ifelse(max.CI==F,"obsVpred_quad_year.png","obsVpred_quad_year_maxCI.png" )
png(figName,height=7,width=7,units="in",res=450)

par(mfrow=c(2,2),tcl=-0.2,mgp=c(2,0.5,0),mar=c(2,2,2,1),oma=c(2,2,0,0))

for(i in 1:4){
  if(i==1) {myTrt="No_grass"}else{myTrt="No_shrub"}
  control=cbind(simD[simD$Treatment=="Control",(6+i)],simD[simD$Treatment=="Control",(2+i)])
  removal.base=cbind(simD[simD$Treatment==myTrt,(6+i)],simD[simD$Treatment==myTrt,(2+i)])
  removal.trt=cbind(simD[simD$Treatment==myTrt,(10+i)],simD[simD$Treatment==myTrt,(2+i)])
  maxCov=1.05*max(c(control,removal.base,removal.trt),na.rm=T)
  plot(control,ylim=c(0,maxCov),xlim=c(0,maxCov),xlab="",ylab="")
  title(sppNames[i],font.main=4)
  abline(0,1,lty="dashed")
  points(removal.base,pch=1,col="blue2")
  points(removal.trt,pch=1,col="red")
  abline(lm(control[,2]~0+control[,1]),col="black")
  abline(lm(removal.base[,2]~0+removal.base[,1]),col="blue2")
  abline(lm(removal.trt[,2]~0+removal.trt[,1]),col="red")
  
  if(i==1){
    legend("topleft",c("Control","Removal (baseline)","Removal (treatment)"),pch=1,
      col=c("black","blue2","red"),lty="solid",bty="n")
  }

}

mtext("Observed cover (%)",2,outer=T,line=0.5,cex=1.2)
mtext("Predicted cover (%)",1,outer=T,line=0.5,cex=1.2)

dev.off()


setwd("..")
 