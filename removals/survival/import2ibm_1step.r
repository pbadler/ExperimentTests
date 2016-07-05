# import & format survival parameters
# then define survival function

# survival parameters
Spars=list(intcpt=rep(NA,Nspp),intcpt.yr=matrix(0,Nyrs,Nspp),intcpt.gr=matrix(0,6,Nspp),
  intcpt.trt=rep(0,Nspp),slope=rep(NA,Nspp),slope.yr=matrix(0,Nyrs,Nspp),
  nb=matrix(0,Nspp,Nspp+2))

for(i in 1:Nspp){

  infile=paste0("survival/",sppList[i],"_surv.csv")

  Sdata=read.csv(infile)
  # main intercept
  Spars$intcpt[i]=Sdata$X.Intercept.[1]
  # Treatment effect
  tmp=grep("Treatment",names(Sdata))
  if(max.CI!=T){
    Spars$intcpt.trt[i] <- Sdata[1,tmp[1]]
  }else{
    Spars$intcpt.trt[i] <- Sdata$TreatMaxCI[1] # add max treatment intercept
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
  tmp=grep("W.",names(Sdata))  # don't filter out W.allcov and W.allpts
  Spars$nb[i,]=as.numeric(Sdata[1,tmp])
} # next i
Spars$yrList=Sdata$year
rm(Sdata)

# survival function
survive=function(Spars,doSpp,doGroup,doYear,sizes,crowding,Treatment){
  trt.Effect=Spars$intcpt.trt[doSpp]
  if(Treatment=="Control") trt.Effect=0
   mu=Spars$intcpt[doSpp]+Spars$intcpt.gr[doGroup,doSpp]+Spars$intcpt.yr[doYear,doSpp]+trt.Effect+
    Spars$slope[doSpp]*sizes+sum(Spars$nb[doSpp,]*crowding)
   out=inv.logit(mu)
   #out=rbinom(length(sizes),1,out)
   out
}
