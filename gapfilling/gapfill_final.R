#doc: this R script plots probability plots and quantile plots of stream gage data, as well as Manahalobis distance
#from the mean vector to assess the normality of the data. The script then imputes missing data using properties of the 
#multivariate normal distribution. Finally, it calculates the Hoteling-T statistic to assess whether the man precipitation
#for the gages has changed between the first 15 years of data, and the subsequent 16 years of data. 
#NB: much of this code could be made more efficient with loops
setwd("~/r/bee6300/hw3")
#read in files
prcp <- read.csv('Ontario.prcp.csv',header=F)
#standardize and center data without disrupting dates
prcp.std <- apply(prcp[,c('V2', 'V3', 'V4', 'V5')], 2, function(x){scale(x)[,1]})
dates <-prcp[,1]
#prcp.std <- cbind(dates,prcp.std)
#plot probability plot of non-exceedence probability for each gage against
#normal distribution cdf
par(mfrow = c(2,2))
for (i in 1:4) {
  n <- length(prcp.std[,i][!is.na(prcp.std[,i])])
  probs <- sort(pnorm(prcp.std[,i],mean=0, sd=1))
  non_exceed <- c(1:n)/(n+1)
  plot(probs, non_exceed, xlab = paste("pnorm gage",as.character(i)), ylab = "non exceedance probability")
  abline(lm(non_exceed~probs))
}
#plot quantile plot of each gage against normal distribution inverse cdf
for (j in 1:4) {
  n <- length(prcp.std[,j][!is.na(prcp.std[,j])])
  non_exceed <- c(1:n)/(n+1)
  gage <- sort(prcp.std[,j][!is.na(prcp.std[,j])])
  quantile <- sort(qnorm(non_exceed,mean=0, sd=1))
  plot(quantile, gage, ylab = paste("gage",as.character(j)), xlab = "inverse cdf")
  abline(lm(gage~quantile))
} 
#create correlation matrix for gage data
S <- cor(prcp.std, use = "pairwise.complete")
mean_vect <- c(mean(prcp.std[,1], na.rm=TRUE),mean(prcp.std[,2], na.rm=TRUE)
               ,mean(prcp.std[,3], na.rm=TRUE) ,mean(prcp.std[,4], na.rm=TRUE))
#calculate Manahalobis distances between each observation and mean vector
D.sq <- vector (length = 31)
for (i in 1:31) {
  D.sq[i] <- t(prcp.std[i,]-mean_vect)%*%(solve(S))%*%(prcp.std[i,]-mean_vect)
}
#plot probability plot of Manahalobis against Chi-Squared distribution (df=4)
par(mfrow=c(1,1))
n <-26
non_exceed <- c(1:n)/(n+1)
chi.prob<- sort(pchisq(D.sq,4))
plot(chi.prob, non_exceed, ylab = "non exceedance probability", xlab = "Chi-Squared Distribution CDF, Df = 4")
abline(lm(non_exceed~chi.prob))
#plot quantile plot of Manahalobis against Chi-Squared distribution (df=4)
chi.quant<- sort(qchisq(non_exceed,4))
D.sq.nona <- sort(D.sq[!is.na(D.sq)])
plot(chi.quant, D.sq.nona, ylab = "Ordered Manahalobis Distances", xlab = "Chi-Squared Distribution CDF^-1, Df = 4")
abline(lm(D.sq.nona~chi.quant))
#imputate data
#find NA indices 
na.index <- which(is.na(prcp.std))
na.ind <- (which(is.na(prcp.std), arr.ind = T))
na.ind <- data.frame(na.ind)
unknowns <- na.ind[order(na.ind$row),]
#row 1, two unknowns at columns 3,4
S11 <- S[c(1,2),c(1,2)]
S11
S21 <- S[c(3,4),c(1,2)]
S21
est <- S21%*%solve(S11)%*%(c(prcp.std[1,1],prcp.std[1,2]))
#descale 
sd3 <- sd(prcp$V4, na.rm = T)
mean3 <-mean(prcp$V4, na.rm = T )
sd4 <- sd(prcp$V5, na.rm = T)
mean4 <-mean(prcp$V5, na.rm = T )
descaled <- c(est[1]*sd3+mean3 , est[2]*sd4+mean4)
descaled
#row 4, two unkowns at columns 2,3
S11 <- S[c(1,4),c(1,4)]
S11
S21 <- S[c(2,3),c(1,4)]
S21
est <- S21%*%solve(S11)%*%(c(prcp.std[4,1],prcp.std[4,4]))
#descale 
sd2 <- sd(prcp$V3, na.rm = T)
mean2 <-mean(prcp$V3, na.rm = T )
sd3 <- sd(prcp$V4, na.rm = T)
mean3 <-mean(prcp$V4, na.rm = T )
descaled <- c(est[1]*sd2+mean2 , est[2]*sd3+mean3)
descaled
#row 9, one unknown at column 4
S11 <- S[c(1,2,3),c(1,2,3)]
S11
S21 <- S[c(4),c(1,2,3)]
S21
est <- S21%*%solve(S11)%*%(c(prcp.std[9,1:3]))
#descale 
sd4 <- sd(prcp$V5, na.rm = T)
mean4 <-mean(prcp$V5, na.rm = T )
descaled <- c(est[1]*sd4+mean4)
descaled
#row 15, one unknown at column 3
S11 <- S[c(1,2,4),c(1,2,4)]
S11
S21 <- S[c(3),c(1,2,4)]
S21
est <- S21%*%solve(S11)%*%(c(prcp.std[15,1:2],prcp.std[15,4]))
#descale 
sd3 <- sd(prcp$V4, na.rm = T)
mean3 <-mean(prcp$V4, na.rm = T )
descaled <- c(est[1]*sd3+mean3)
descaled
#row 24, two unknowns at columns 1,4
S11 <- S[c(2,3),c(2,3)]
S11
S21 <- S[c(1,4),c(2,3)]
S21
est <- S21%*%solve(S11)%*%(c(prcp.std[24,2],prcp.std[24,3]))
#descale 
sd1 <- sd(prcp$V2, na.rm = T)
mean1 <-mean(prcp$V2, na.rm = T )
sd4 <- sd(prcp$V5, na.rm = T)
mean4 <-mean(prcp$V5, na.rm = T )
descaled <- c(est[1]*sd1+mean1 , est[2]*sd4+mean4)
descaled
#plot raw data for each gage over time
par(mfrow = c(2,2))
for (i in 2:5) {
  gage <- prcp[,i]
  plot(prcp[,1],gage, ylab = paste("flow at gage",as.character(i-1)), xlab = "year")
  abline(lm(prcp[,i] ~ prcp[,1]))
}
#run linear regressions for each gage against time
lin_reg1 <- lm(prcp[,2] ~ prcp[,1])
lin_reg1
summary(lin_reg1)
lin_reg2 <- lm(prcp[,3] ~ prcp[,1])
lin_reg2
summary(lin_reg2)
lin_reg3 <- lm(prcp[,4] ~ prcp[,1])
lin_reg3
summary(lin_reg3)
lin_reg4 <- lm(prcp[,5] ~ prcp[,1])
lin_reg4
summary(lin_reg4)
#partition data into 1980-1994 and 1995-2010 segments
pre94 <- prcp[1:15,2:5]
post95 <-prcp[16:31,2:5]
#calculate multivariate means
pre94_means <- colMeans(pre94[,], na.rm = T, dims = 1)
post95_means <- colMeans(post95[,], na.rm = T, dims = 1)
#calculate covariance matrix
pre94_S <- cov(pre94, use = "pairwise.complete")
post95_S <- cov(post95, use = "pairwise.complete")
#estimate pooled covariance matrix 
n1 <- nrow(pre94)
n2 <- nrow(post95)
n <- n1+n2
S_pool <- (((n1-1)/(n-2))*pre94_S+((n2-1)/(n-2))*post95_S)/(n)
#Calculate Hoteling T-statistic 
ht_stat <- t(pre94_means-post95_means)%*%solve(S_pool)%*%(pre94_means-post95_means)
ht_stat
#Calculate F-stat
k <- 4
df2 <- n-k-1
alpha <- 0.05 
f_stat <- df2*ht_stat/((n-2)*k)
f_stat
#Calculate critical level in F distribution
crit <- qf(1-alpha,k,df2)
crit
