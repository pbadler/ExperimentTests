
# function to format and output parameters of survival models

formatSurvPars <-function(model,outfile){

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

  write.table(params,outfile,row.names=F,sep=",")
    
}