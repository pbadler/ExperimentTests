# Import and format growth parameters
# then define growth function

# growth parameters
Gpars=list(intcpt=rep(NA,Nspp),intcpt.yr=matrix(0,Nyrs,Nspp),intcpt.gr=matrix(0,6,Nspp),
  slope=rep(NA,Nspp),slope.yr=matrix(0,Nyrs,Nspp),
  nb=matrix(0,Nspp,Nspp),sigma2.a=rep(NA,Nspp),sigma2.b=rep(NA,Nspp))
for(i in 1:Nspp){

  infile=paste0("growth/",sppList[i],"_growth.csv")

  Gdata=read.csv(infile)
  # main intercept
  Gpars$intcpt[i]=Gdata$X.Intercept.[1]
  # random year effects on intercept
  Gpars$intcpt.yr[,i]=Gdata$Intercept.yr
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
  # Treatment effects
  if(trtEffects==T){
    tmp=grep("Treatment",names(Gdata))
    if(length(tmp)>0) Gpars$intcpt[i]=Gpars$intcpt[i] + Gdata[1,tmp[1]] # add treatment intercept
    if(length(tmp)>1) Gpars$nb[i,i]=Gpars$nb[i,i] + Gdata[1,tmp[2]] # add treatment density dependence
    if(length(tmp)>2) stop("too many treatment effects")
  }
  # variance parameters
  Gpars$sigma2.a[i]=Gdata$sigma.a[1]
  Gpars$sigma2.b[i]=Gdata$sigma.b[1]
} # next i
rm(Gdata)

# growth
grow=function(Gpars,doSpp,doGroup,doYear,sizes,crowding){
  # crowding and nb are vectors of dim Nspp
  logsizes=log(sizes)
  mu=Gpars$intcpt[doSpp]+Gpars$intcpt.yr[doYear,doSpp]+Gpars$intcpt.gr[doGroup,doSpp]+
    (Gpars$slope[doSpp]+Gpars$slope.yr[doYear,doSpp])*logsizes+ Gpars$nb[doSpp,]%*%crowding
  tmp=which(mu<log(minSize)*1.5)  # we will kill vanishingly small plants to avoid numerical overflow...below
  mu[tmp]=log(minSize) # truncate tiny sizes (they cause problems in sigma2)
  sigma2=Gpars$sigma2.a[doSpp]*exp(Gpars$sigma2.b[doSpp]*mu)
  #out=exp(rnorm(length(sizes),mu,sqrt(sigma2)))
  out=exp(mu)
  if(sum(is.na(out))>0) browser()
  out[tmp]=0   # here's the killing
  #out[out>maxSize[doSpp]]=maxSize[doSpp] #truncate big plants
  out
}
