
# function to format and output parameters of growth models

formatGrowthPars <-function(model,outfile){
  
  # fit variance
  x=model$summary.fitted.values$mean
  y=(x-allD$logarea.t1)^2  # squared residuals
  #plot(x,y)
  outVar=nls(y~a*exp(b*x),start=list(a=1,b=0))
  
  # calculate scaled residuals and write to file
  V.pred=predict(outVar)
  scaledResiduals=(allD$logarea.t1-x)/sqrt(V.pred) # use observed - predicted
  scaledResiduals = (scaledResiduals - mean(scaledResiduals))/sd(scaledResiduals) #make sure mean = 0 and variance =1
  tmp=paste0(substring(outfile,1,4),"_scaled_residuals.csv")
  write.csv(scaledResiduals,tmp,row.names=F) # write to file
  
  #year effects
  Intercept.yr=model$summary.random$yearID$mean
  year=model$summary.random$year$ID
  logarea.yr=model$summary.random$year$mean
  params=data.frame(year,Intercept.yr,logarea.yr)
  #group effects
  tmp=as.data.frame(model$summary.random$GroupID[,1:2])
  if(dim(tmp)[1]>0){
    names(tmp)=c("GroupID","Group")
    tmp[(nrow(tmp)+1):dim(params)[1],]=NA
    params=cbind(params,tmp)
  }
  #fixed effects
  tmp=matrix(NA,dim(params)[1],length(model$names.fixed))
  colnames(tmp)=model$names.fixed
  tmp[1,]=model$summary.fixed$mean
  params=cbind(params,tmp)
  
  # record the edge of the 95% CI furthest from zero
  j <- grep("Treatment",row.names(model$summary.fixed))
  tmp <- model$summary.fixed[j,c("0.025quant","0.975quant")]
  k <- which(abs(tmp)==max(abs(tmp)))
  params$TreatMaxCI <- NA
  params$TreatMaxCI[1] <- as.numeric(tmp[k])
  
  #variance 
  params$sigma.a=NA; params$sigma.a[1]=coef(outVar)[1] 
  params$sigma.b=NA; params$sigma.b[1]=coef(outVar)[2]
  
  write.table(params,outfile,row.names=F,sep=",")
    
}