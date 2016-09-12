
# call from PSSP or POSE growth models

# ARTR cover
# historical ARTR quad cover
tmp1 <- read.csv(paste0(dataDir1, "/speciesData/ARTR/quadratCover.csv"))
tmp1$year<-tmp1$year + 1900
# modern ARTR quad cover
tmp2 <- read.csv(paste0(dataDir2, "/speciesData/ARTR/quadratCover.csv"))
covD <- rbind(tmp1,tmp2)

# get list of all quad years
#  load quad inventories
QInv1 <- read.csv(paste0(dataDir1, "../../../driversdata/data/idaho/quad_inventory.csv")) # not included in drivers data 
QInv1$year <- QInv1$year + 1900  
QInv1 <- reshape(QInv1, idvar = "year", times=names(QInv1)[2:NCOL(QInv1)],
                 varying = list(2:NCOL(QInv1)),timevar = "quad", v.names="present",
                 direction = "long")
QInv2 <- read.csv(paste0("../../../driversdata/data/idaho_modern/quad_inventory.csv"))
QInv2$year <- c(2007:2015)
QInv2 <- QInv2[,c("year",names(QInv2)[1:(NCOL(QInv2)-1)])]
QInv2 <- reshape(QInv2, idvar = "year", times=names(QInv2)[2:NCOL(QInv2)],
                 varying = list(2:NCOL(QInv2)),timevar = "quad", v.names="present",
                 direction = "long")
QInv <-rbind(QInv1,QInv2)
QInv <- subset(QInv,!is.na(present))
QInv <- QInv[,c("quad","year")]

# merge all quad years with ARTR data
covD<-merge(covD,QInv,all.y=T)

# add in true zeros
covD[is.na(covD)] <- 0

covD <- reshape(covD,direction="wide",idvar="year",timevar="quad")
covD <- covD[order(covD$year),]

# get rid of new plots
keep <- which(colSums(!is.na(covD))>5)
covD <- covD[,keep]

exclude.quads <- names(sort(colMeans(covD[,2:NCOL(covD)],na.rm=T))[1:8])
exclude.quads <- unlist(strsplit(exclude.quads,"\\."))
exclude.quads <- exclude.quads[(grep("Q",exclude.quads))]
