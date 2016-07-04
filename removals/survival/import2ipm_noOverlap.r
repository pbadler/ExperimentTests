# import & format survival parameters
# then define survival function

# survival parameters
Spars=list(intcpt=rep(NA,Nspp),intcpt.yr=matrix(0,Nyrs,Nspp),intcpt.gr=matrix(0,6,Nspp),
  intcpt.trt=rep(0,Nspp),slope=rep(NA,Nspp),slope.yr=matrix(0,Nyrs,Nspp),
  nb=matrix(0,Nspp,Nspp))

for(i in 1:Nspp){
  
  infile=paste0("survival/",sppList[i],"_surv.csv")
  
  Sdata=read.csv(infile)
  # main intercept
  Spars$intcpt[i]=Sdata$X.Intercept.[1]
  # Treatment effect
  if(trtEffects==T){
    tmp=grep("Treatment",names(Sdata))
    if(max.CI==T) {
      Spars$intcpt[i]=Spars$intcpt[i] + Sdata$TreatMaxCI[1] # add max treatment intercept
    }else{
      Spars$intcpt[i]=Spars$intcpt[i] + Sdata[1,tmp[1]] # add treatment intercept
    }
  }
  # Group effects
  tmp=which(names(Sdata)=="Group")
  if(length(tmp)>0) Spars$intcpt.gr[,i]=Sdata$Group[!is.na(Sdata$Group)] # get spatial average
  # random year effects
  Spars$intcpt.yr[,i]=Sdata$Intercept.yr
  Spars$slope[i]=Sdata$logarea[1]
  # random effects on slope
  tmp=which(names(Sdata)=="logarea.yr")
  if(length(tmp)>0) Spars$slope.yr[,i]=Sdata[,tmp]
  # get competition coefficients
  tmp=paste("W",sppList,sep=".")
  tmp=which(is.element(names(Sdata),tmp))
  if(length(tmp)>0) Spars$nb[i,]=as.numeric(Sdata[1,tmp])
} # next i
Spars$yrList=Sdata$year
rm(Sdata)


## survival function: probability an individual of size u survives  (u is on log scale)
S=function(u,W,Spars,doYear,doSpp){
  mu=Spars$intcpt[doSpp]+Spars$intcpt.yr[doYear,doSpp]+Spars$intcpt.trt[doSpp]+
    (Spars$slope[doSpp]+Spars$slope.yr[doYear,doSpp])*u+
    W%*%Spars$nb[doSpp,]
  return(inv.logit(mu))
}



