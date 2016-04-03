# Individually-based model for a multiple species,
# with density-dependence and explicit space (random
# spatial pattern)

# Simulate one quadrat through time

# PBA  8-28-09

# qName = "Q1"
totSims=50
Nyrs=22
Ngroups=6
startYr=2011
# doGroup=1  # E1 exclosure
restartYrs=c(2,15,17,20)
L=100 # dimension of square quadrat (cm)
expand=1  # 1 = 1x1 m^2, 2 = 2x2m^2, etc
#sppList=c("ARTR","HECO","POSE","PSSP")
myCol=c("black","gold1","blue","red")
minSize=0.25
maxSize=c(8000,500,500,500)

outfile1=paste(qName,"_validation_cov_removals.csv",sep="")
outfile2=paste(qName,"_validation_den_removals.csv",sep="")


#GET OBSERVED DATA AND INITIAL CONDITIONS -------------------------------------------
Nspp=length(sppList)
obsA=data.frame(year=startYr:2015)
obsN=data.frame(year=startYr:2015)
init.plants=data.frame(NULL)
lastID=0
for(i in 1:length(sppList)){
  infile=paste("c:\\repos\\driversdata\\data\\idaho_modern\\speciesData\\",sppList[i],"\\",sppList[i],"_genet_xy.csv",sep="")
  tmpD=read.csv(infile)
  tmpD=subset(tmpD,quad==qName & tmpD$year>=startYr)
    if(dim(tmpD)[1]>0){
      obsCov=aggregate(tmpD$area,by=list(tmpD$year),FUN=sum)
      names(obsCov)=c("year",sppList[i])
      obsA=merge(obsA,obsCov,all.x=T)
      obsDen=aggregate(rep(1,dim(tmpD)[1]),by=list(tmpD$year),FUN=sum)
      names(obsDen)=c("year",sppList[i])
      obsN=merge(obsN,obsDen,all.x=T)
      tmpD2=subset(tmpD,year==startYr)
      if(dim(tmpD2)[1]>0){
          spp=rep(i,dim(tmpD2)[1])
          id=(lastID+1:length(spp))
          tmpD2=data.frame(cbind(spp,tmpD2[,c("area","x","y")],id))
          names(tmpD2)=c("spp","size","x","y","id")
          init.plants=rbind(init.plants,tmpD2)
          lastID=max(init.plants[,5])
      }    
    }
}
rm(tmpD)

# if one species missing, add in zeros
tmp=which(!is.element(sppList,names(obsA)))
if(length(tmp)==1){
  missing=rep(0,dim(obsA)[1])
  missing[rowSums(!is.na(obsA[,2:NCOL(obsA)]))==0]=NA
  obsA$missing=missing
  obsN$missing=missing
  names(obsA)[NCOL(obsA)]=sppList[tmp]
  names(obsN)[NCOL(obsN)]=sppList[tmp]
  obsA=obsA[,c("year",sppList)]
  obsN=obsN[,c("year",sppList)]
}

# find true zeros
tmp1=rowSums(is.na(obsA))<Nspp
tmp1=matrix(tmp1,length(tmp1),Nspp)
tmp2=is.na(obsA[,2:NCOL(obsA)])
tmp3=which(tmp1==T & tmp2==T,arr.ind=T)
tmp4=matrix(FALSE,dim(tmp2)[1],dim(tmp2)[2])
tmp4[tmp3]=T
tmp5=obsA[,2:NCOL(obsA)]
tmp5[tmp4]=0
obsA[,2:NCOL(obsA)]=tmp5
tmp5=obsN[,2:NCOL(obsA)]
tmp5[tmp4]=0
obsN[,2:NCOL(obsN)]=tmp5
rm(tmp1,tmp2,tmp3,tmp4,tmp5)


# FORMAT PARAMETERS ------------------------------------------------
Nspp=length(sppList)

# use model from Adler et al. 2010

# set up survival parameters and function
source("H:/idahochart/ipm/multispp_glm_v3/survival/import2ibm.r")
# set up growth parameters and function
source("H:/idahochart/ipm/multispp_glm_v3/growth/import2ibm.r")
# set up recruitment parameters and function
source("H:/idahochart/ipm/multispp_glm_v3/recruitment/import2ibm.r")

# model spatial group variation (or not)
if(!is.na(doGroup)){
  Spars$intcpt=Spars$intcpt+Spars$intcpt.gr[doGroup,]
  Gpars$intcpt=Gpars$intcpt+Gpars$intcpt.gr[doGroup,]
  Rpars$intcpt.yr=Rpars$intcpt.yr+matrix(Rpars$intcpt.gr[doGroup,],Nyrs,Nspp,byrow=T)
}

# FUNCTIONS---------------------------------------------------------
library(boot)
library(mvtnorm)
library(msm)

# crowding function, assumes toroidal landscape
getCrowding=function(plants,alpha,L,expand){
 # plants is a matrix: sizes in column 1; x,y coords in columns 2 and 3
 # d is the distance weighting parameter
 # functions returns a vector of length = rows in plants
 if(dim(plants)[1]>1){
   xdiff=abs(outer(plants[,3],plants[,3],FUN="-"))
   #tmp=which(xdiff>((L*expand)/2))
   #xdiff[tmp]=(L*expand)-xdiff[tmp]
   ydiff=abs(outer(plants[,4],plants[,4],FUN="-"))
   #tmp=which(ydiff>((L*expand)/2))
   #ydiff[tmp]=(L*expand)-ydiff[tmp]
   distMat=sqrt(xdiff^2+ydiff^2) 
   distMat[distMat==0]=NA
   #distMat[distMat<1]=1  # CLUDGE TO PREVENT GROWTH EXPLOSIONS
   distMat=exp(-1*alpha[plants[,1]]*distMat^2)
   sizeMat=matrix(plants[,2],dim(plants)[1],dim(plants)[1])
   out=aggregate(distMat*sizeMat,by=list("spp"=plants[,1]),FUN=sum,na.rm=T)
   # put in missing zeros
   tmp=data.frame("spp"=c(1:length(sppList)))
   out=merge(out,tmp,all.y=T)
   out[is.na(out)]=0
   out=out[order(out$spp),]
   out=as.matrix(out[,c(2:NCOL(out))])  # drop spp column
 }else{
   out=rep(0,Nspp)
 }
 out
}

# MAIN LOOP -------------------------------------------------------
if(qName=="Q11" |  qName=="Q13" | qName=="Q14")  yrList=yrList[yrList<50 | yrList>55]
doYrList=yrList[yrList>=startYr]
Nyrs=length(doYrList)
bigA=array(NA,dim=c((Nyrs+1),Nspp,totSims))
bigN=array(NA,dim=c((Nyrs+1),Nspp,totSims))
for(iSim in 1:totSims){
  # initialize with N.init plants of size.init for each species
  plants=init.plants[[1]]
  lastID=max(plants[,5])
  
  # arrays to store results
  N=matrix(0,(Nyrs+1),Nspp)
  N[1,1:4]=as.numeric(obsN[1,2:5])
  A=matrix(0,(Nyrs+1),Nspp)
  A[1,1:4]=as.numeric(obsA[1,2:5])
  
  # plot initial conditions
  par(mgp=c(2,0.5,0),tcl=-0.2)
  symbols(x = plants[,3], y = plants[,4], circles = sqrt(plants[,2]/pi),fg=myCol[plants[,1]],
    xlim=c(0,L*expand),ylim=c(0,L*expand),main ="Time=1",xlab="x",ylab="y",inches=F,lwd=2)
  
  for(tt in 1:Nyrs){
     # draw year effects
     doYr=which(yrList==doYrList[tt])
     
     if(doYr==restartYrs[1]) plants=init.plants[[2]]
     if(doYr==restartYrs[2]) plants=init.plants[[3]]
     if(doYr==restartYrs[3]) plants=init.plants[[4]]
     if(doYr==restartYrs[4]) plants=init.plants[[5]]
     
     nextplants=plants
     
     # recruitment
     newplants=recruit(Rpars,sizes=plants[,2],spp=plants[,1],doYear=doYr,lastID=lastID,L,expand)
     
     for(ss in 1:Nspp){
      if(N[tt,ss]>0){ # make sure spp ss is not extinct
  
        # growth
        W=getCrowding(plants,Gpars$alpha[ss,],L,expand)
        newsizes=grow(Gpars,doSpp=ss,doYear=doYr,sizes=plants[,2],crowding=W)
  
        if(sum(newsizes==Inf)>0) browser()
        if(is.na(sum(newsizes))) browser()
        
        # survival
        # uses same W as growth            
        live=survive(Spars,doSpp=ss,doYear=doYr,sizes=plants[,2],crowding=W)
  
        # put it all together
        tmp=which(plants[,1]==ss)  # only alter plants of focal spp        
        nextplants[tmp,2]=newsizes[tmp]*live[tmp]   #update with G and S
   

        
       } # end if no plants
     } # next ss  
     
     nextplants=nextplants[nextplants[,2]>0,]    # remove dead plants 
     nextplants=rbind(nextplants,newplants)     # add recruits

     # output cover and density
     tmp=aggregate(nextplants[,2],by=list(nextplants[,1]),FUN=sum)
     A[tt+1,tmp[,1]]=tmp[,2]
     tmp=aggregate(rep(1,dim(nextplants)[1]),by=list(nextplants[,1]),FUN=sum)
     N[tt+1,tmp[,1]]=tmp[,2]
     
     plants=nextplants
     lastID=max(plants[,5])  

     # plot
     if(sum(N[tt,])>0){
      symbols(x = plants[,3], y = plants[,4], circles = sqrt(plants[,2]/pi),fg=myCol[plants[,1]],
        xlim=c(0,L*expand),ylim=c(0,L*expand),main =paste("Time=",tt,sep=""),
        xlab="x",ylab="y",inches=F,lwd=2)
     }else{
       break
     }
     print(tt);flush.console() 
  } # next tt 
  bigN[,,iSim]=N ; bigA[,,iSim]=A
  print(paste("Sim ",iSim," complete",sep=""))
  flush.console()
} # next iSim

# average across simulations
predA=apply(bigA,MARGIN=c(1,2),FUN=mean,na.rm=T)
predN=apply(bigN,MARGIN=c(1,2),FUN=mean)

year=c(doYrList[1],1+doYrList)

predA=data.frame(cbind(year,predA))
names(predA)[2:dim(predA)[2]]=paste(sppList,"pred",sep="")
output1=merge(obsA,predA,all.x=T)
tmp=which(is.element(output1$year,init.yrs[2:5]))
output1[tmp,(2+Nspp):NCOL(output1)]=output1[tmp,2:(Nspp+1)]

predN=data.frame(cbind(year,predN))
names(predN)[2:dim(predN)[2]]=paste(sppList,"pred",sep="")
output2=merge(obsN,predN,all.x=T)
tmp=which(is.element(output2$year,init.yrs[2:5]))
output2[tmp,(2+Nspp):NCOL(output2)]=output2[tmp,2:(Nspp+1)]

par(mfrow=c(1,2),tcl=-0.2,mgp=c(2,0.5,0))
matplot(output1[,1],output1[,2:NCOL(output1)],type="o",
  col=myCol,lty=c(rep("solid",Nspp),c(rep("dashed",Nspp))),
  pch=c(rep(16,Nspp),rep(1,Nspp)),xlab="Year",ylab="Cover")
matplot(output2[,1],output2[,2:NCOL(output2)],type="o",
  col=myCol,lty=c(rep("solid",Nspp),c(rep("dashed",Nspp))),
  pch=c(rep(16,Nspp),rep(1,Nspp)),xlab="Year",ylab="Density")

setwd("H:/idahochart/ipm/multispp_glm_v3/validate")
write.table(output1,outfile1,row.names=F,sep=",")
write.table(output2,outfile2,row.names=F,sep=",")

setwd("H:/idahochart/ipm/multispp_glm_v3")
