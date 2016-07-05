
# import parameters
# recruitment parameters
Rpars=list(intcpt.mu=rep(0,Nspp),intcpt.yr=matrix(0,Nyrs,Nspp),intcpt.trt=matrix(0,3,Nspp),intcpt.tau=rep(100,Nspp),
  intcpt.gr=matrix(NA,6,Nspp),g.tau=rep(NA,Nspp),
  dd=matrix(NA,Nspp,Nspp),theta=rep(NA,Nspp),sizeMean=rep(NA,Nspp),sizeVar=rep(NA,Nspp),
  recSizes=list(1))

#   if(trtEffects==T){
   infile=paste0("recruitment/recruit_params_m1.csv")
#   }else{
#     infile=paste0("recruitment/recruit_params_m0.csv")
#   } 

 Rdata=read.csv(infile)
 
 # subset out non-essential parameters
 tmp=c(grep("lambda",row.names(Rdata)),grep("deviance",row.names(Rdata)),
  grep("DIC",row.names(Rdata)))   
 Rdata=Rdata[-tmp,]
 
 # format parameters
 tmp=paste("Rpars$",row.names(Rdata),"<-",Rdata[,1],sep="")
 eval(parse(n=dim(Rdata)[1],text=tmp))
 for(i in 1:Nspp){
   infile=paste(root,"/driversdata/data/idaho/speciesData/",sppList[i],"/recSize.csv",sep="")
   recSize=read.csv(infile)
   Rpars$sizeMean[i]=mean(log(recSize$area))
   Rpars$sizeVar[i]=var(log(recSize$area))
   #Rpars$recSizes[[i]]=recSize$area
 }
Rpars$dd=t(Rpars$dd) # c[i,j] = effect of j on i

# get edge of confidence interval for treatment effects
if(max.CI==T){
  keep <- grep("intcpt.trt",row.names(Rdata))
  tmp <- Rdata[keep,c("X2.5.","X97.5.")]
  tmp$max <- ifelse(abs(tmp$X2.5.)>abs(tmp$X97.5.),tmp$X2.5.,tmp$X97.5.)
  Rpars$intcpt.trt[2,2] <- tmp$max[which(row.names(tmp)=="intcpt.trt[2,2]")]
  Rpars$intcpt.trt[2,3] <- tmp$max[which(row.names(tmp)=="intcpt.trt[2,3]")]
  Rpars$intcpt.trt[2,4] <- tmp$max[which(row.names(tmp)=="intcpt.trt[2,4]")]
  Rpars$intcpt.trt[3,1] <- tmp$max[which(row.names(tmp)=="intcpt.trt[3,1]")]
}

# reformat treatment effects
Rpars$intcpt.trt <- c(Rpars$intcpt.trt[3,1], # ARTR in no_grass
                      Rpars$intcpt.trt[2,2:4])   # grasses in no_shrub

rm(Rdata)

# recruitment
recruit=function(Rpars,totArea,doGroup,doYear,Treatment){
    # dd is a matrix of dim Nspp X Nspp
    # sizes and spp are vectors of the same length (=N plants)
    trt.Effect = Rpars$intcpt.trt
    if(Treatment=="Control") trt.Effect[] = 0
    # calculate recruits
    out=lambda=rep(NA,Nspp) # seed production
    for(i in 1:Nspp){
       lambda[i]=totArea[i]*exp(Rpars$intcpt.yr[doYear,i]+Rpars$intcpt.gr[doGroup,i]+
                   trt.Effect[i]+sum(sqrt(totArea)*Rpars$dd[i,]))
       out[i]=lambda[i]*exp(Rpars$sizeMean[i])
    }
   out=out/100 # convert to % cover
   return(out)
}
