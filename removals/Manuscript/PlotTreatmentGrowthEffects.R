par(mfrow=c(2,2)); 
# ARTR
z=seq(0,8,length=100); 
z1 = 0.717910833 + 0.868422492*z; 
z1t = z1+0.135455557; 
sigma2 = 1.41632672*exp(-0.187017078*z1); 
matplot(z,cbind(z1,z1t,z1+1.6*sqrt(sigma2),z1-1.6*sqrt(sigma2)),xlab="Logarea t", ylab="Logarea t+1",type="l",
col=c("blue","red","blue","blue"),lty=c(1,1,2,2), lwd=c(2,2,1,1)); 
title(main="ARTR"); 

# HECO
z=seq(-1,5,length=100); 
z1 = 0.401102802 +0.820463327*z; 
z1t = z1+0.064933303
sigma2 = 0.697290865*exp(-0.062086902*z1); 
matplot(z,cbind(z1,z1t,z1+1.6*sqrt(sigma2),z1-1.6*sqrt(sigma2)),xlab="Logarea t", ylab="Logarea t+1",type="l",
col=c("blue","red","blue","blue"),lty=c(1,1,2,2), lwd=c(2,2,1,1)); 
title(main="HECO");

# POSE
z=seq(-1,4,length=100); 
z1 = 0.508775295	+ 0.668209953*z; 
z1t = z1 + 0.200615309
sigma2=1.065517749*exp(-0.05906554*z1)
matplot(z,cbind(z1,z1t,z1+1.6*sqrt(sigma2),z1-1.6*sqrt(sigma2)),xlab="Logarea t", ylab="Logarea t+1",type="l",
col=c("blue","red","blue","blue"),lty=c(1,1,2,2), lwd=c(2,2,1,1)); 
title(main="POSE");

#PSSP 
z=seq(-1,5,length=100); 
z1 = 0.416608292	+ 0.827328598*z; 
z1t = z1 + 0.233883335
sigma2 = 0.895560448*exp(-0.145014501*z1) 
matplot(z,cbind(z1,z1t,z1+1.6*sqrt(sigma2),z1-1.6*sqrt(sigma2)),xlab="Logarea t", ylab="Logarea t+1",type="l",
col=c("blue","red","blue","blue"),lty=c(1,1,2,2), lwd=c(2,2,1,1)); 
title(main="PSSP");


