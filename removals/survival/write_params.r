
# function to format and output parameters of growth models

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

  write.table(params,outfile,row.names=F,sep=",")
    
}