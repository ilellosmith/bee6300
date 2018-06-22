#doc: this script explores collinearity in flow data from two gages in upstate NY
#data is standardized, projected onto 1:1 line and Euclidean norms are calculated
setwd("~/r/bee6300/hw1")
#read in files
flow.data1 <- read.csv('gage1.csv',header=T)
flow.data2 <- read.csv('gage2.csv',header=T)
#separate flow data by gage
x1 <- flow.data1$flow
x2 <- flow.data2$flow
#Standardize time series - subtract mean and divide by sd
x1_scale <- (x1 - mean(x1))/sd(x1)
x2_scale <- (x2 - mean(x2))/sd(x2)
#Combine the data into a single n×2 matrix, X
X <- cbind(x1_scale, x2_scale)
# plot gage data
plot(X, pch = 20, yaxs ="i", xaxs = "i", xlim = c(-5,5), ylim = c(-5, 5), 
     xlab = "gage1 standardized", ylab = "gage2 standardized")
legend("bottomright", legend=c("gage data standardized", "projected data"),
       col=c("black", "red"), pch=20:20, cex=0.8, bty = "n")
abline(0,1, h = 0, v = 0, untf = FALSE, col = "black")
#create projection matrix
u <- c(2,2)
norm_u <- norm(u,type = "2")
P <- u%*%t(u)/norm_u^2
#project data onto 1:1 line
Xp <- X%*%P 
"
P <- X%*%(solve(t(X)%*%X))%*%t(X) 
#got computational singular (not invertible)
#P <- X%*%t(X)
#you don't make a projection matrix out of the data set 
#the projection matrix is formed by the variable you are trying to project onto 
#create projection matrix from any vector along the 1:1 line 
#define u as c(1,1) or c(5,5) P= u%*%t(u(/norm_u^2))
#Px 2X2 by 2X1 is stream flow for one year (scalar)
#Xp is how to multiply 82X2 by 2X2 for all recorded years 
"
dim(X)
dim(Xp)
dim(P)
points(Xp,type = "p", pch = 20, col = "red")
#calculate Euclidean norm of each new data point 
mynorm <- numeric(length = 80)
for (i in 1:(length(Xp)/2)) {
  mynorm[i] <- sqrt(2)*Xp[i,1]
}
mynorm 
#indicates which of the projected data are in the 3rd quadrant
#neg <- which(Xp[,1]<0 & Xp[,2]<0)
#multiplies the norm for those projected data points by -1
#mynorm[neg] <- mynorm[neg]*-1
#points(mynorm,type = "p", pch = 20, col = "green")
#mynorm over gage1 standardized
plot(x1_scale, pch = 20, yaxs ="i", xaxs = "i", xlim = c(0,90), ylim = c(-4,4),
     ylab = "annual flow", xlab = "observation") 
abline(h = 0, v = 0, untf = FALSE, col = "black")
#xlim = c(0,90), ylim = c(0, 400)
legend("bottomright", legend=c("gage1 data standardized", "mynorm"),
       col=c("black", "orange"), pch=20:20, cex=0.8, bty = "n")
points(mynorm,type = "p", pch = 20, col = "orange")
#mynorm over gage2 standardized
plot(x2_scale, pch = 20, yaxs ="i", xaxs = "i", xlim = c(0,90),ylim = c(-4,4),
     ylab = "annual flow", xlab = "observation")
abline(h = 0, v = 0, untf = FALSE, col = "black")
#xlim = c(0,90), ylim = c(0, 2000)
legend("bottomright", legend=c("gage2 data standardized", "mynorm"),
       col=c("black", "orange"), pch=20:20, cex=0.8,  bty = "n")
points(mynorm,type = "p", pch = 20, col = "orange")
#mynorm as a function of gage1
plot(x1, mynorm, pch = 20, col = "purple", yaxs ="i", xaxs = "i", xlim = c(0,350), ylim = c(-4, 4)
   , ylab = "mynorm" , xlab = "gage1")
abline(h = 0, v = 0, untf = FALSE, col = "black")
#mynorm as a function of gage1 standardized
plot(x1_scale, mynorm, pch = 20, col = "purple", yaxs ="i", xaxs = "i", xlim = c(-3,3), ylim = c(-4,4),
     ylab = "mynorm" , xlab = "gage1 standardized")
abline(h = 0, v = 0, untf = FALSE, col = "black")
#mynorm as a function of gage2
plot(x2, mynorm, pch = 20, col = "purple", yaxs ="i", xaxs = "i", xlim = c(0,1750), ylim = c(-4, 4)
  , ylab = "mynorm", xlab = "gage2")
abline(h = 0, v = 0, untf = FALSE, col = "black")
#mynorm as a function of gage2 standardized
plot(x2_scale, mynorm, pch = 20, col = "purple", yaxs ="i", xaxs = "i", xlim = c(-3,3), ylim = c(-4,4),
     ylab = "mynorm" , xlab = "gage2 standardized")
abline(h = 0, v = 0, untf = FALSE, col = "black")
