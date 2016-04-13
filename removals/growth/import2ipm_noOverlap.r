# Import and format growth parameters
# then define growth function

# growth parameters
Gpars=list(intcpt=rep(NA,Nspp),intcpt.yr=matrix(0,Nyrs,Nspp),intcpt.gr=matrix(0,6,Nspp),intcpt.trt=rep(0,Nspp),
  slope=rep(NA,Nspp),slope.yr=matrix(0,Nyrs,Nspp),
  nb=matrix(0,Nspp,Nspp),sigma2.a=rep(NA,Nspp),sigma2.b=rep(NA,Nspp))
for(i in 1:Nspp){
  if(trtEffects==T){
    infile=paste0("growth/",sppList[i],"_growth_Trt.csv")
  }else{
    infile=paste0("growth/",sppList[i],"_growth_noTrt.csv")
  }
  Gdata=read.csv(infile)
  # main intercept
  Gpars$intcpt[i]=Gdata$X.Intercept.[1]
  # random year effects on intercept
  Gpars$intcpt.yr[,i]=Gdata$Intercept.yr
  # Treatment effect
  tmp=grep("Treatment",names(Gdata))
  if(length(tmp)>0) Gpars$intcpt.trt[i]=Gdata[1,tmp]
  # group effects
  tmp=which(names(Gdata)=="Group")
  if(length(tmp)>0) Gpars$intcpt.gr[,i]=Gdata$Group[!is.na(Gdata$Group)] 
  # size slope
  Gpars$slope[i]=Gdata$logarea.t0[1]
  # random effects on slope
  tmp=which(names(Gdata)=="logarea.yr")
  if(length(tmp)>0) Gpars$slope.yr[,i]=Gdata[,tmp]
  # get competition coefficients
  tmp=paste("W",sppList,sep=".")
  tmp=which(is.element(names(Gdata),tmp))
  if(length(tmp)>0) Gpars$nb[i,]=as.numeric(Gdata[1,tmp])
  # variance parameters
  Gpars$sigma2.a[i]=Gdata$sigma.a[1]
  Gpars$sigma2.b[i]=Gdata$sigma.b[1]
} # next i
rm(Gdata)

# growth function
G=function(v,u,W,Gpars,doYear,doSpp){
  mu=Gpars$intcpt[doSpp]+Gpars$intcpt.yr[doYear,doSpp]+Gpars$intcpt.trt[doSpp]+
    (Gpars$slope[doSpp]+Gpars$slope.yr[doYear,doSpp])*u+
    W%*%Gpars$nb[doSpp,]
  sigma2=Gpars$sigma2.a[doSpp]*exp(Gpars$sigma2.b[doSpp]*mu)
  out=dnorm(v,mu,sqrt(sigma2))
  out
}

