
rm(list=ls(all=TRUE))

#set file dir in control code *or*
D<-read.csv("../../../driversdata/data/idaho/climate/Dubois_daily_alldata.txt")

data<-data.frame(D$year, D$month, D$day, D$Tmax, D$Tmin, D$Prcp)
names(data)<-c("year","month","day","Tmax", "Tmin", "Ppt")

data$Tmax<-ifelse(data$Tmax>9999, NA, data$Tmax)
data$Tmin<-ifelse(data$Tmin>9999, NA, data$Tmin)
data$Ppt<-as.numeric(levels(data$Ppt))[data$Ppt] 
data$Ppt<-ifelse(data$Ppt>99, NA, data$Ppt)
data$seqDay<-1:nrow(data)

#replace NA's with means of row before and after
library(zoo)
datan<-data
datan$Tmax<-na.approx(data$Tmax)
datan$Tmin<-na.approx(data$Tmin)
datan$Ppt[is.na(datan$Ppt)]<-0


#See what the function is doing:
par(mfrow=c(2,1))
#plot(Tmax~seqDay, data=data[1:1000,], type="l")
#plot(Tmax~seqDay, data=data[200:400,], type="l")
plot(Tmax~seqDay, data=data[290:350,], type="l")
points(Tmax~seqDay, data=datan[290:350,], type="l",lty=2, col="blue")

#plot(Ppt~seqDay, data=data[2000:3000,], type="l")
plot(Ppt~seqDay, data=data[2650:2780,], type="l")
points(Ppt~seqDay, data=datan[2650:2780,], type="l",lty=2, col="blue")


#flag data that is interpolated
naData<-is.na(data)*1
interpYN<-ifelse(rowSums(naData)<1, 0, 1)
datan<-cbind(datan,interpYN)
head(datan)
datan<-datan[,-7]


#Attempt to reduce variable number
#mean of Tmax and Tmin
datan$Tmid<-rowMeans(datan[,4:5], dims = 1)


#export interpolated file
write.csv(datan, "interpClim.csv")
