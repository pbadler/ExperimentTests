
# call all scripts for removal experiment analysis

rm(list=ls(all=TRUE))
graphics.off();

root=ifelse(.Platform$OS.type=="windows","c:/Repos","~/repos"); # modify as needed
setwd(paste(root,"/ExperimentTests/removals/",sep="")); # modify as needed 

# required packages
library(texreg) # to save output
library(lme4)
library(xtable)
library(INLA)
library(boot)
library(R2WinBUGS)
library("TeachingDemos") # for inset plots
library("quantreg")

###
### 1. get treatment trends #################################
###

statsOutput <- paste0(getwd(),"/stats_tables.tex")
source("treatment_trends_removals.r")

# make climate figure
source("climate_fig.r")

# clean up
tmp=ls() ; tmp=tmp[tmp!="root" & tmp!="statsOutput"]
rm(list=tmp)

###
### 2. fit vital rate regressions ###########################
###

# table to store Treatment effects
trtTests <- data.frame("species"="c","stage"="c","effect"=1,"CI.02.5"=1,"CI.97.5"=1,stringsAsFactors = F)

# read in distance weights
#dists <- read.csv(paste0(root,"/ExperimentTests/data/idaho_modern/speciesdata/IdahoModDistanceWeights_noExptl.csv"))
dists <- read.csv(paste0(root,"/ExperimentTests/data/idaho/speciesdata/IdahoDistanceWeights.csv"))


###
### fit survival models (takes ~ 10 minutes)
###
setwd("survival")
source("write_params.r") # get function to format and output parameters
for(iSpp in c("ARTR","HECO","POSE","PSSP")){
  
  source(paste0(iSpp,"survival.r"))
  
  # save fixed effects summary to file
  cat("",file=statsOutput,sep="\n",append=T)
  cat(capture.output(print(xtable(m1$summary.fixed,digits=4,caption=paste("Summary of fixed effects for the",iSpp,"survival model"),
        label=paste0(iSpp,"survival")),caption.placement="top")),file=statsOutput,sep="\n",append=T)
  
  # save treatment test
  irow <- dim(trtTests)[1]
  trtTests[irow+1,] <- NA
  trtTests[irow+1,1:2] <- c(iSpp,"survival")
  tmp <- grep("Treatment",row.names(m1$summary.fixed))
  trtTests[irow+1,3:5] <- m1$summary.fixed[tmp,c("mean","0.025quant","0.975quant")]
  
  # write parameters
  formatSurvPars(m1,paste0(iSpp,"_surv.csv")) 
  
}
setwd("..")

###
### fit growth models (takes < 5 mins)
###

setwd("growth")
source("write_params.r") # get function to format and output parameters
growth_residuals <- list() # place to store growth residuals needed for plotting
sppList <- c("ARTR","HECO","POSE","PSSP")
for(i in 1:length(sppList)){
  
  iSpp <- sppList[i]
  source(paste0(iSpp,"growth.r"))
  
  # save fixed effects summary to file
  cat("",file=statsOutput,sep="\n",append=T)
  cat(capture.output(print(xtable(m1$summary.fixed,digits=4,caption=paste("Summary of fixed effects for the",iSpp,"growth model"),
        label=paste0(iSpp,"growth")),caption.placement="top")),file=statsOutput,sep="\n",append=T)

  # save treatment test
  irow <- dim(trtTests)[1]
  trtTests[irow+1,] <- NA
  trtTests[irow+1,1:2] <- c(iSpp,"growth")
  tmp <- grep("Treatment",row.names(m1$summary.fixed))
  trtTests[irow+1,3:5] <- m1$summary.fixed[tmp,c("mean","0.025quant","0.975quant")]
  
  # write parameters for best model
  formatGrowthPars(m1,paste0(iSpp,"_growth.csv")) 
  
  # save growth residuals
  allD$resids <- allD$logarea.t1-m1$summary.fitted.values$mean
  growth_residuals[[i]] <- allD[,c("quad","year","Treatment","logarea.t1","resids")]
  
}

# make growth residuals figure
source("growth_residuals_fig.R")

setwd("..")

###
### fit recruitment models (this can take a few hours)
###

setwd("recruitment")

source("call_recruit_m1.r")

# add treatment test data for ARTR 
pars.summary <- read.csv("recruit_params_m1.csv")
irow <- dim(trtTests)[1]
trtTests[irow+1,] <- NA
trtTests[irow+1,1:2] <- c("ARTR","recruitment")
tmp <- grep("intcpt.trt",row.names(pars.summary))
trtTests[irow+1,3:5] <- pars.summary[tmp[5],c("mean","X2.5.","X97.5.")]

# add treatment test data for the three grasses
irow <- dim(trtTests)[1]
trtTests[(irow+1):(irow+3),] <- NA
trtTests[(irow+1):(irow+3),1:2] <- cbind(c("HECO","POSE","PSSP"),rep("recruitment",3))
trtTests[(irow+1):(irow+3),3:5] <- pars.summary[tmp[2:4],c("mean","X2.5.","X97.5.")]

# create summary stats table
fixeffs <- c("intcpt.mu","intcpt.trt","dd","theta","u")
keep <- NULL
for(i in fixeffs){
  keep <- c(keep,grep(i,row.names(pars.summary),fixed=TRUE))
}
keep <- keep[-c(5,10:12,33:44)]  # some editing
output <- pars.summary[keep,c(1:3,7:9)]
cat("",file=statsOutput,sep="\n",append=T)
cat(capture.output(print(xtable(output,digits=4,caption=paste("Summary of fixed effects for the recruitment model"),
        label="table:recruitment"),caption.placement="top")),file=statsOutput,sep="\n",append=T)
  
setwd("..")


###
### write results to file
###

trtTests <- trtTests[-1,] # throw away first line junk
write.csv(trtTests,"treatment_test_results.csv",row.names=F)

# make treatment effect figure
source("treatment_test_figure.r")

# clean up
tmp=ls() ; tmp=tmp[tmp!="root" & tmp!="statsOutput"]
rm(list=tmp)

###
### 3. explore neighborhood composition ###################################
###

source("Wdistrib/exploreSurvivalWs.R")

# clean up
tmp=ls() ; tmp=tmp[tmp!="root" & tmp!="statsOutput"]
rm(list=tmp)

###
### 4. get IBM predictions for quadrat cover ###############################
###


sppList <-  c("ARTR","HECO","POSE","PSSP")

# read in distance weights
#dists <- read.csv(paste0(root,"/ExperimentTests/data/idaho_modern/speciesdata/IdahoModDistanceWeights_noExptl.csv"))
dists <- read.csv(paste0(root,"/ExperimentTests/data/idaho/speciesdata/IdahoDistanceWeights.csv"))

max.CI <- F  # TRUE means use maximum removal effect
source("ibm/ibm_removal_1step.r")
source("ibm/summarize_sims1step.r")

max.CI <- T  # TRUE means use maximum removal effect
source("ibm/ibm_removal_1step.r")
source("ibm/summarize_sims1step.r")


# clean up
tmp=ls() ; tmp=tmp[tmp!="root" & tmp!="statsOutput"]
rm(list=tmp)


###
### 5. get equilibrium cover from an IPM ###############################
###

max.CI <- F  # TRUE means use maximum removal effect

sppList <-  c("ARTR","HECO","POSE","PSSP")

source("ipm/get_W_functions.r")  # get neighbor distance decay functions

#no treatment effects, all species
init.species <- c(1:4)
tlimit <- 2500
burn.in <- 500
trtEffects=F
source("ipm/IPM-setup.r")
source("ipm/IPM-getEquilibrium.r")
write.csv(covSave[(burn.in+1):tlimit,],"ipm/baselineCover.csv",row.names=F)
meanCover1 <- meanCover

#no treatment effects, ARTR removal
init.species <- c(2:4)
source("ipm/IPM-getEquilibrium.r")
# SPE 
par(mfrow=c(2,1)); 
HECO.base = sizeSave[[2]]; 
plot(v[[2]],apply(HECO.base,1,mean));
#
write.csv(covSave[(burn.in+1):tlimit,],"ipm/baselineCover-noARTR.csv",row.names=F)
meanCover2 <- meanCover

# removal treatment effects, ARTR removal
init.species <- c(2:4)
trtEffects=T
max.CI=F
source("ipm/IPM-setup.r")
source("ipm/IPM-getEquilibrium.r")
# SPE 
HECO.trt = sizeSave[[2]]; 
plot(v[[2]],apply(HECO.trt,1,mean));
#
write.csv(covSave[(burn.in+1):tlimit,],"ipm/removalCover-noARTR.csv",row.names=F)
meanCover3 <- meanCover

# removal treatment effects, ARTR removal, max treat effect
init.species <- c(2:4)
trtEffects=T
max.CI=T
source("ipm/IPM-setup.r")
source("ipm/IPM-getEquilibrium.r")
write.csv(covSave[(burn.in+1):tlimit,],"ipm/removalCover-noARTR-maxCI.csv",row.names=F)
meanCover4 <- meanCover

# # removal treatment effects, ARTR removal, no PSSP, what happens to HECO and POSE?
# init.species <- c(2:3)
# trtEffects=T
# max.CI=F
# source("ipm/IPM-setup.r")
# source("ipm/IPM-getEquilibrium.r")
# write.csv(covSave[(burn.in+1):tlimit,],"ipm/removalCover-noARTRnoPSSP.csv",row.names=F)

# baseline model (no treatment effects), remove all grasses
init.species <- c(1)
trtEffects=F
max.CI=F
source("ipm/IPM-setup.r")
source("ipm/IPM-getEquilibrium.r")
write.csv(covSave[(burn.in+1):tlimit,],"ipm/baselineCover-noGrass.csv",row.names=F)

# baseline model (with treatment effects), remove all grasses
init.species <- c(1)
trtEffects=T
max.CI=F
source("ipm/IPM-setup.r")
source("ipm/IPM-getEquilibrium.r")
write.csv(covSave[(burn.in+1):tlimit,],"ipm/removalCover-noGrass.csv",row.names=F)

simResults <- rbind(meanCover1,meanCover2,meanCover3) 
colnames(simResults) <- sppList
write.csv(simResults,"ipm/simResults-meanCover.csv",row.names=F)

simResults <- rbind(meanCover1,meanCover2,meanCover4) 
colnames(simResults) <- sppList
write.csv(simResults,"ipm/simResults-meanCover-maxCI.csv",row.names=F)

#no treatment effects, all species, boost ARTR cover
init.species <- c(1:4)
tlimit <- 2500
burn.in <- 500
trtEffects=F
max.CI=F
source("ipm/IPM-setup.r")
Rpars$intcpt.yr[,1] <- Rpars$intcpt.yr[,1]+1
source("ipm/IPM-getEquilibrium.r")
print(rbind(meanCover1,meanCover)) # compare baseline run with this one

source("ipm/IPM-figures.r") # plot results of ARTR removals

# eyeball results for ARTR cover following grass removals
baseline<-mean(read.csv("ipm/baselineCover.csv")[,1])
nograss1<-mean(read.csv("ipm/baselineCover-noGrass.csv")[,1])
nograss2<-mean(read.csv("ipm/removalCover-noGrass.csv")[,1])
ARTRcover<-c(baseline,nograss1,nograss2)
names(ARTRcover)<-c("Baseline","Baseline, grass removal","Treatment effects, grass removal")
print(ARTRcover)

