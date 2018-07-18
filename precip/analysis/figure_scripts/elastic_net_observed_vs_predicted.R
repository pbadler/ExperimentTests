
# custom plotting function to show output of elastic net fits and predictions

obs_pred_fig <- function(x, trts, vital_rate, effect = "Intercept", legend_location){
  # x is a list containing observations and predictions named y, y_hat, y_new, and y_hat_new
  # trts is a vector of treatments for the out of sample data points
  # vital_rate is a string: either "growth" or "survival"
  # effect is a string with default "Intercept" and alternative value "slope"
  # legend_location is a string, such as "topleft"

  plot(c(x$y,x$y_new),c(x$y_hat,x$y_hat_new),type="n",xlab="Observed",ylab="Predicted",
       ylim=c(min(c(x$y,x$y_new,x$y_hat,x$y_hat_new)),max(c(x$y,x$y_new,x$y_hat,x$y_hat_new))),
       xlim=c(min(c(x$y,x$y_new,x$y_hat,x$y_hat_new)),max(c(x$y,x$y_new,x$y_hat,x$y_hat_new))),
       main=paste(sppList[i], vital_rate, effect))
  abline(0,1)
  points(x$y,x$y_hat)
  points(x$y_new[which(trts =="Control")],x$y_hat_new[which(trts =="Control")],pch=16)
  points(x$y_new[which(trts =="Drought")],x$y_hat_new[which(trts =="Drought")],pch=16,col="red")
  points(x$y_new[which(trts =="Irrigation")],x$y_hat_new[which(trts =="Irrigation")],pch=16,col="blue")
  legend(legend_location,c("Control (training)","Control (out-of-sample)","Drought (out-of-sample)",
                     "Irrigation (out-of-sample)"),pch=c(1,16,16,16),
                      col=c("black","black","red","blue"),bty="n",cex=0.8)

}