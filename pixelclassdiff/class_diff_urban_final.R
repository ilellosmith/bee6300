#This script takes pixel and superobject data from http://archive.ics.uci.edu/ml/datasets/Urban+Land+Cover# 
#and applies Fischer's linear discrimination to classify testing data into each of 7 land type categories: 
#building, grass, pavement (includes asphalt, concrete, car), pool, shadow, soil, and tree. 
#The script applies Fischer's linear discrimination twice to compare model accuracy using just spectral information 
#to model accuracy using spectral information and super object information 
#NB: This code need stylistic changes. It would benefit greatly from function definition with some of the repetitive
#calls to for loops. 
setwd("~/r/bee6300/f_proj")
# load dependencies, installing if necessary
REQUIRED_PACKAGES <- c("GGally", "ggplot2", "plyr")
package.check <- lapply(REQUIRED_PACKAGES, FUN = function(x) {
  if (! require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
covariates <- read.csv("training.csv")
select_covariates <- data.frame(cbind(covariates$Area, covariates$Round, covariates$Bright, covariates$Compact, covariates$LW, covariates$Assym
                                      , covariates$Dens, covariates$Mean_G, covariates$Mean_R, covariates$Mean_NIR, covariates$SD_G, covariates$SD_R
                                      , covariates$SD_NIR, covariates$GLCM1, covariates$NDVI, covariates$class))
colnames(select_covariates) <- c('area', 'round', 'bright', 'compact', 'lw', 'assym', 'dens', 'mG', 'mR', 'mNIR', 'sdG', 'sdR',
                                 'sdNIR', 'GLCM', 'NDVI', 'class')
#scale covariates
#shape covariates
scale.area <- (covariates$Area-mean(covariates$Area))/sd(covariates$Area)
scale.round <- (covariates$Round-mean(covariates$Round))/sd(covariates$Round)
scale.bright <- (covariates$Bright-mean(covariates$Bright))/sd(covariates$Bright)
scale.compact <- (covariates$Compact-mean(covariates$Compact))/sd(covariates$Compact)
scale.lw <- (covariates$LW-mean(covariates$LW))/sd(covariates$LW)
scale.assym <- (covariates$Assym-mean(covariates$Assym))/sd(covariates$Assym)
#fill covariates
scale.dens <- (covariates$Dens-mean(covariates$Dens))/sd(covariates$Dens)
#spectral covariates
scale.mG <- (covariates$Mean_G-mean(covariates$Mean_G))/sd(covariates$Mean_G)
scale.mR <- (covariates$Mean_R-mean(covariates$Mean_R))/sd(covariates$Mean_R)
scale.mNIR <- (covariates$Mean_NIR-mean(covariates$Mean_NIR))/sd(covariates$Mean_NIR)
scale.sdG <- (covariates$SD_G-mean(covariates$SD_G))/sd(covariates$SD_G)
scale.sdR <- (covariates$SD_R-mean(covariates$SD_R))/sd(covariates$SD_R)
scale.sdNIR <- (covariates$SD_NIR-mean(covariates$SD_NIR))/sd(covariates$SD_NIR)
scale.GLCM <- (covariates$GLCM1-mean(covariates$GLCM1))/sd(covariates$GLCM1)
scale.NDVI <- (covariates$NDVI-mean(covariates$NDVI))/sd(covariates$NDVI)
#PART 1: run analysis with just spectral covariates
#combine spectral covariates
k <- 8 
scale.spect <- data.frame(matrix(c(select_covariates$class, scale.mG, scale.mR, scale.mNIR
                                   , scale.sdG, scale.sdR, scale.sdNIR, scale.GLCM, scale.NDVI), ncol = k+1, nrow = 168))
colnames(scale.spect) <- c('class','mG', 'mR', 'mNIR', 'sdG', 'sdR',
                           'sdNIR', 'GLCM', 'NDVI')
class <- (split(scale.spect[,2:(k+1)], scale.spect[,1]))
class.even <- do.call(rbind.fill, class)
#plot covariance matrixes for each variable and each group
cols <- c('aquamarine', 'firebrick4', 'darkorange3', 'goldenrod3', 'forestgreen', 'blue','blueviolet')
cols.all <- class
cols.all$`1`[!is.na(cols.all$`1`)] <- cols[1]
final.cols <- c(cols.all$`1`[,1])
cols.all$`2`[!is.na(cols.all$`2`)] <- cols[2]
final.cols <- c(final.cols, cols.all$`2`[,1])
cols.all$`3`[!is.na(cols.all$`3`)] <- cols[3]
final.cols <- c(final.cols, cols.all$`3`[,1])
cols.all$`4`[!is.na(cols.all$`4`)] <- cols[4]
final.cols <- c(final.cols, cols.all$`4`[,1])
cols.all$`5`[!is.na(cols.all$`5`)] <- cols[5]
final.cols <- c(final.cols, cols.all$`5`[,1])
cols.all$`6`[!is.na(cols.all$`6`)] <- cols[6]
final.cols <- c(final.cols, cols.all$`6`[,1])
cols.all$`7`[!is.na(cols.all$`7`)] <- cols[7]
final.cols <- c(final.cols, cols.all$`7`[,1])
p <- ggpairs(class.even, aes(color = final.cols))+ theme_bw()
# Change color manually.
# Loop through each plot changing relevant scales
for(i in 1:p$nrow) {
  for(j in 1:p$ncol){
    p[i,j] <- p[i,j] + 
      scale_fill_manual(values= cols) +
      scale_color_manual(values= cols)  
  }
}
dev.new()
p
#multiple discriminant analysis
grp1 <- class$`1`
n1 <- dim(grp1)[1]
grp2 <- class$`2`
n2 <- dim(grp2)[1]
grp3 <- class$`3`
n3 <- dim(grp3)[1]
grp4 <- class$`4`
n4 <- dim(grp4)[1]
grp5 <- class$`5`
n5 <- dim(grp5)[1]
grp6 <- class$`6`
n6 <- dim(grp6)[1]
grp7 <- class$`7`
n7 <- dim(grp7)[1]
N <- sum(n1,n2,n3,n4,n5,n6,n7)
G <- 7
#calculate pooled covariance matrix
sum_cov <- (n1-1)*cov(grp1)+(n2-1)*cov(grp2)+(n3-1)*cov(grp3)+(n4-1)*cov(grp4)+(n5-1)*cov(grp5)+
  (n6-1)*cov(grp6)+(n7-1)*cov(grp7)
spool <- 1/(N-G)*sum_cov
#calculate scaled grand means
grandmeans <- vector()
for (i in 1:(dim(class.even)[2])) {
  grandmeans <- c(grandmeans, mean(class.even[,i]))
}
#calculate scaled group means
grp1_m <- vector()
for (i in 1: length(grp1)) {
  grp1_m <- c(grp1_m,mean(grp1[,i]))
}
grp2_m <- vector()
for (i in 1: length(grp2)) {
  grp2_m <- c(grp2_m,mean(grp2[,i]))
}
grp3_m <- vector()
for (i in 1: length(grp3)) {
  grp3_m <- c(grp3_m,mean(grp3[,i]))
}
grp4_m <- vector()
for (i in 1: length(grp4)) {
  grp4_m <- c(grp4_m,mean(grp4[,i]))
}
grp5_m <- vector()
for (i in 1: length(grp5)) {
  grp5_m <- c(grp5_m,mean(grp5[,i]))
}
grp6_m <- vector()
for (i in 1: length(grp6)) {
  grp6_m <- c(grp6_m,mean(grp6[,i]))
}
grp7_m <- vector()
for (i in 1: length(grp7)) {
  grp7_m <- c(grp7_m,mean(grp7[,i]))
}
groupmeans <- rbind(grp1_m,grp2_m,grp3_m,grp4_m,grp5_m,grp6_m,grp7_m)
groupmeans <- t(groupmeans)
#calculate between groups covariance matrix
Sb <- 1/(G-1)*(groupmeans-grandmeans)%*%t((groupmeans-grandmeans))
sig_noise <- solve(spool)%*%Sb
J <- G-1
#calculate alphas 
alpha <- matrix(nrow=k,ncol=6)
for (i in 1:J){
  eigenvect <- as.vector(sig_noise[,i])
  a <- eigenvect/c((t(eigenvect)%*%spool%*%eigenvect)^(1/2))
  alpha[,i] <- t(a)
}
write.csv(alpha, 'alpha1.csv')
#classify test data 
test <- read.csv("testing.csv")
select_test <- data.frame(cbind(test$Area, test$Round, test$Bright, test$Compact, test$LW, test$Assym
                                , test$Dens, test$Mean_G, test$Mean_R, test$Mean_NIR, test$SD_G, test$SD_R
                                , test$SD_NIR, test$GLCM1, test$NDVI, test$class))
colnames(select_test) <- c('area', 'round', 'bright', 'compact', 'lw', 'assym', 'dens', 'mG', 'mR', 'mNIR', 'sdG', 'sdR',
                           'sdNIR', 'GLCM', 'NDVI', 'class')
#scale test data
scale.mG <- (test$Mean_G-mean(test$Mean_G))/sd(test$Mean_G)
scale.mR <- (test$Mean_R-mean(test$Mean_R))/sd(test$Mean_R)
scale.mNIR <- (test$Mean_NIR-mean(test$Mean_NIR))/sd(test$Mean_NIR)
scale.sdG <- (test$SD_G-mean(test$SD_G))/sd(test$SD_G)
scale.sdR <- (test$SD_R-mean(test$SD_R))/sd(test$SD_R)
scale.sdNIR <- (test$SD_NIR-mean(test$SD_NIR))/sd(test$SD_NIR)
scale.GLCM <- (test$GLCM1-mean(test$GLCM1))/sd(test$GLCM1)
scale.NDVI <- (test$NDVI-mean(test$NDVI))/sd(test$NDVI)
test.spect <- (matrix(c(test$class, scale.mG, scale.mR, scale.mNIR
                        , scale.sdG, scale.sdR, scale.sdNIR, scale.GLCM, scale.NDVI), ncol = k+1, nrow = 507))
colnames(test.spect) <- c('class','mG', 'mR', 'mNIR', 'sdG', 'sdR',
                          'sdNIR', 'GLCM', 'NDVI')
test.spect.dat <- data.frame(test.spect)
class.test <- (split(test.spect.dat[,2:(k+1)], test.spect.dat[,1]))
grp1 <- class.test$`1`
n1 <- dim(grp1)[1]
grp2 <- class.test$`2`
n2 <- dim(grp2)[1]
grp3 <- class.test$`3`
n3 <- dim(grp3)[1]
grp4 <- class.test$`4`
n4 <- dim(grp4)[1]
grp5 <- class.test$`5`
n5 <- dim(grp5)[1]
grp6 <- class.test$`6`
n6 <- dim(grp6)[1]
grp7 <- class.test$`7`
n7 <- dim(grp7)[1]
N <- sum(n1,n2,n3,n4,n5,n6,n7)
#calculate distance of each entry from groups 
classified <- c(1:dim(test.spect)[1])
for (i in 1:dim(test.spect)[1]) {
  dist.g <- matrix(nrow = G, ncol = J)
  for(g in 1:G) {
    for (j in 1:J) {
      dist.g[g,j] <- sum(c((t(alpha[,j])*(((test.spect[i,2:(k+1)]))-(groupmeans[,g])))^2))
    }
  }
  final.dist <- rowSums(dist.g)
  classified[i] <-  which.min(final.dist)
}
confusion.test <- table(test.spect[,1], classified, dnn = c('Actual Group', 'Predicted Group'))
rownames(confusion.test) <- c('building','grass', 'pavement','pool','shadow','soil','tree')
colnames(confusion.test) <- c('building','grass', 'pavement','pool','shadow','soil','tree')
test.m <- as.matrix(confusion.test)
correct <- sum(diag(test.m))
msclass <- sum(test.m)-correct
err_rate <- msclass/(correct+msclass)
msclass
correct+msclass
err_rate
T_count <- c(n1,n2,n3,n4,n5,n6,n7)
correct_g <- vector()
for (g in 1:G) {
  correct_g[g] <- test.m[g,g]/T_count[g]
}
test.m <- rbind(confusion.test, correct_g, T_count)
write.csv(test.m, "confusion_test.csv")
#classify training data to compare prediction accuracy 
classified <- c(1:dim(scale.spect)[1])
for (i in 1:dim(scale.spect)[1]) {
  dist.g <- matrix(nrow = G, ncol = J)
  for(g in 1:G) {
    for (j in 1:J) {
      dist.g[g,j] <- sum(c((t(alpha[,j])*(((scale.spect[i,2:(k+1)]))-(groupmeans[,g])))^2))
    }
  }
  final.dist <- rowSums(dist.g)
  classified[i] <-  which.min(final.dist)
}
confusion.train <- table(scale.spect[,1], classified, dnn = c('Actual Group', 'Predicted Group'))
rownames(confusion.train) <- c('building','grass', 'pavement','pool','shadow','soil','tree')
colnames(confusion.train) <- c('building','grass', 'pavement','pool','shadow','soil','tree')
grp1 <- class$`1`
n1 <- dim(grp1)[1]
grp2 <- class$`2`
n2 <- dim(grp2)[1]
grp3 <- class$`3`
n3 <- dim(grp3)[1]
grp4 <- class$`4`
n4 <- dim(grp4)[1]
grp5 <- class$`5`
n5 <- dim(grp5)[1]
grp6 <- class$`6`
n6 <- dim(grp6)[1]
grp7 <- class$`7`
n7 <- dim(grp7)[1]
N <- sum(n1,n2,n3,n4,n5,n6,n7)
test.m <- as.matrix(confusion.train)
correct <- sum(diag(test.m))
msclass <- sum(test.m)-correct
err_rate <- msclass/(correct+msclass)
msclass
correct+msclass
err_rate
T_count <- c(n1,n2,n3,n4,n5,n6,n7)
correct_g <- vector()
for (g in 1:G) {
  correct_g[g] <- test.m[g,g]/T_count[g]
}
test.m <- rbind(confusion.train, correct_g, T_count)
write.csv(test.m, "confusion_train.csv")
#PART 2 run analysis with spectral and super-object covariates
covariates <- read.csv("training.csv")
select_covariates <- data.frame(cbind(covariates$Area, covariates$Round, covariates$Bright, covariates$Compact, covariates$LW, covariates$Assym
                                      , covariates$Dens, covariates$Mean_G, covariates$Mean_R, covariates$Mean_NIR, covariates$SD_G, covariates$SD_R
                                      , covariates$SD_NIR, covariates$GLCM1, covariates$NDVI, covariates$class))
colnames(select_covariates) <- c('area', 'round', 'bright', 'compact', 'lw', 'assym', 'dens', 'mG', 'mR', 'mNIR', 'sdG', 'sdR',
                                 'sdNIR', 'GLCM', 'NDVI', 'class')
#scale all covariates
scale.area <- (covariates$Area-mean(covariates$Area))/sd(covariates$Area)
scale.round <- (covariates$Round-mean(covariates$Round))/sd(covariates$Round)
scale.bright <- (covariates$Bright-mean(covariates$Bright))/sd(covariates$Bright)
scale.compact <- (covariates$Compact-mean(covariates$Compact))/sd(covariates$Compact)
scale.lw <- (covariates$LW-mean(covariates$LW))/sd(covariates$LW)
scale.assym <- (covariates$Assym-mean(covariates$Assym))/sd(covariates$Assym)
scale.dens <- (covariates$Dens-mean(covariates$Dens))/sd(covariates$Dens)
scale.mG <- (covariates$Mean_G-mean(covariates$Mean_G))/sd(covariates$Mean_G)
scale.mR <- (covariates$Mean_R-mean(covariates$Mean_R))/sd(covariates$Mean_R)
scale.mNIR <- (covariates$Mean_NIR-mean(covariates$Mean_NIR))/sd(covariates$Mean_NIR)
scale.sdG <- (covariates$SD_G-mean(covariates$SD_G))/sd(covariates$SD_G)
scale.sdR <- (covariates$SD_R-mean(covariates$SD_R))/sd(covariates$SD_R)
scale.sdNIR <- (covariates$SD_NIR-mean(covariates$SD_NIR))/sd(covariates$SD_NIR)
scale.GLCM <- (covariates$GLCM1-mean(covariates$GLCM1))/sd(covariates$GLCM1)
scale.NDVI <- (covariates$NDVI-mean(covariates$NDVI))/sd(covariates$NDVI)
#combine covariates
k <- 15
scale.spect.super <- data.frame(matrix(c(select_covariates$class, scale.mG, scale.mR, scale.mNIR
                                         , scale.sdG, scale.sdR, scale.sdNIR, scale.GLCM, scale.NDVI, scale.area,scale.round,scale.bright,scale.compact,scale.lw
                                         ,scale.assym, scale.dens), ncol = k+1, nrow = 168))
colnames(scale.spect.super) <- c('class','mG', 'mR', 'mNIR', 'sdG', 'sdR',
                                 'sdNIR', 'GLCM', 'NDVI', 'area', 'round', 'bright', 'compact', 'lw', 'assym', 'dens')
class2 <- (split(scale.spect.super[,2:(k+1)], scale.spect.super[,1]))
class.even2 <- do.call(rbind.fill, class2)
#plot covariance matrixes for each variable and each group
cols.all <- class2
cols.all$`1`[!is.na(cols.all$`1`)] <- cols[1]
final.cols <- c(cols.all$`1`[,1])
cols.all$`2`[!is.na(cols.all$`2`)] <- cols[2]
final.cols <- c(final.cols, cols.all$`2`[,1])
cols.all$`3`[!is.na(cols.all$`3`)] <- cols[3]
final.cols <- c(final.cols, cols.all$`3`[,1])
cols.all$`4`[!is.na(cols.all$`4`)] <- cols[4]
final.cols <- c(final.cols, cols.all$`4`[,1])
cols.all$`5`[!is.na(cols.all$`5`)] <- cols[5]
final.cols <- c(final.cols, cols.all$`5`[,1])
cols.all$`6`[!is.na(cols.all$`6`)] <- cols[6]
final.cols <- c(final.cols, cols.all$`6`[,1])
cols.all$`7`[!is.na(cols.all$`7`)] <- cols[7]
final.cols <- c(final.cols, cols.all$`7`[,1])
p <- ggpairs(class.even2, aes(color = final.cols))+ theme_bw()
# Change color manually.
# Loop through each plot changing relevant scales
for(i in 1:p$nrow) {
  for(j in 1:p$ncol){
    p[i,j] <- p[i,j] + 
      scale_fill_manual(values= cols) +
      scale_color_manual(values= cols)  
  }
}
dev.new()
p
#multiple discriminant analysis
grp1 <- class2$`1`
n1 <- dim(grp1)[1]
grp2 <- class2$`2`
n2 <- dim(grp2)[1]
grp3 <- class2$`3`
n3 <- dim(grp3)[1]
grp4 <- class2$`4`
n4 <- dim(grp4)[1]
grp5 <- class2$`5`
n5 <- dim(grp5)[1]
grp6 <- class2$`6`
n6 <- dim(grp6)[1]
grp7 <- class2$`7`
n7 <- dim(grp7)[1]
N <- sum(n1,n2,n3,n4,n5,n6,n7)
G <- 7
#calculate pooled covariance matrix
sum_cov <- (n1-1)*cov(grp1)+(n2-1)*cov(grp2)+(n3-1)*cov(grp3)+(n4-1)*cov(grp4)+(n5-1)*cov(grp5)+
  (n6-1)*cov(grp6)+(n7-1)*cov(grp7)
spool <- 1/(N-G)*sum_cov
#calculate scaled grand means
grandmeans <- vector()
for (i in 1:(dim(class.even2)[2])) {
  grandmeans <- c(grandmeans, mean(class.even2[,i]))
}
#calculate scaled group means
grp1_m <- vector()
for (i in 1: length(grp1)) {
  grp1_m <- c(grp1_m,mean(grp1[,i]))
}
grp2_m <- vector()
for (i in 1: length(grp2)) {
  grp2_m <- c(grp2_m,mean(grp2[,i]))
}
grp3_m <- vector()
for (i in 1: length(grp3)) {
  grp3_m <- c(grp3_m,mean(grp3[,i]))
}
grp4_m <- vector()
for (i in 1: length(grp4)) {
  grp4_m <- c(grp4_m,mean(grp4[,i]))
}
grp5_m <- vector()
for (i in 1: length(grp5)) {
  grp5_m <- c(grp5_m,mean(grp5[,i]))
}
grp6_m <- vector()
for (i in 1: length(grp6)) {
  grp6_m <- c(grp6_m,mean(grp6[,i]))
}
grp7_m <- vector()
for (i in 1: length(grp7)) {
  grp7_m <- c(grp7_m,mean(grp7[,i]))
}
groupmeans <- rbind(grp1_m,grp2_m,grp3_m,grp4_m,grp5_m,grp6_m,grp7_m)
groupmeans <- t(groupmeans)
#calculate between groups covariance matrix
Sb <- 1/(G-1)*(groupmeans-grandmeans)%*%t((groupmeans-grandmeans))
sig_noise <- solve(spool)%*%Sb
J <- G-1
#calculate alphas 
alpha <- matrix(nrow=k,ncol=J)
for (i in 1:J){
  eigenvect <- as.vector(sig_noise[,i])
  a <- eigenvect/c((t(eigenvect)%*%spool%*%eigenvect)^(1/2))
  alpha[,i] <- t(a)
}
write.csv(alpha, 'alpha2.csv')
#classify test data 
select_test2 <- data.frame(cbind(test$Area, test$Round, test$Bright, test$Compact, test$LW, test$Assym
                                 , test$Dens, test$Mean_G, test$Mean_R, test$Mean_NIR, test$SD_G, test$SD_R
                                 , test$SD_NIR, test$GLCM1, test$NDVI, test$class))
colnames(select_test2) <- c('area', 'round', 'bright', 'compact', 'lw', 'assym', 'dens', 'mG', 'mR', 'mNIR', 'sdG', 'sdR',
                            'sdNIR', 'GLCM', 'NDVI', 'class')
#scale test data
scale.mG <- (test$Mean_G-mean(test$Mean_G))/sd(test$Mean_G)
scale.mR <- (test$Mean_R-mean(test$Mean_R))/sd(test$Mean_R)
scale.mNIR <- (test$Mean_NIR-mean(test$Mean_NIR))/sd(test$Mean_NIR)
scale.sdG <- (test$SD_G-mean(test$SD_G))/sd(test$SD_G)
scale.sdR <- (test$SD_R-mean(test$SD_R))/sd(test$SD_R)
scale.sdNIR <- (test$SD_NIR-mean(test$SD_NIR))/sd(test$SD_NIR)
scale.GLCM <- (test$GLCM1-mean(test$GLCM1))/sd(test$GLCM1)
scale.NDVI <- (test$NDVI-mean(test$NDVI))/sd(test$NDVI)
scale.area <- (test$Area-mean(test$Area))/sd(test$Area)
scale.round <- (test$Round-mean(test$Round))/sd(test$Round)
scale.bright <- (test$Bright-mean(test$Bright))/sd(test$Bright)
scale.compact <- (test$Compact-mean(test$Compact))/sd(test$Compact)
scale.lw <- (test$LW-mean(test$LW))/sd(test$LW)
scale.assym <- (test$Assym-mean(test$Assym))/sd(test$Assym)
scale.dens <- (test$Dens-mean(test$Dens))/sd(test$Dens)
test.spect2 <- (matrix(c(test$class, scale.mG, scale.mR, scale.mNIR
                         , scale.sdG, scale.sdR, scale.sdNIR, scale.GLCM, scale.NDVI, scale.area,scale.round,scale.bright,scale.compact,scale.lw
                         ,scale.assym, scale.dens), ncol = k+1, nrow = 507))
colnames(test.spect2) <- c('class','mG', 'mR', 'mNIR', 'sdG', 'sdR',
                           'sdNIR', 'GLCM', 'NDVI', 'area', 'round', 'bright', 'compact', 'lw', 'assym', 'dens')
test.spect2.dat <- data.frame(test.spect)
class.test <- (split(test.spect2.dat[,2:(k+1)], test.spect2.dat[,1]))
grp1 <- class.test$`1`
n1 <- dim(grp1)[1]
grp2 <- class.test$`2`
n2 <- dim(grp2)[1]
grp3 <- class.test$`3`
n3 <- dim(grp3)[1]
grp4 <- class.test$`4`
n4 <- dim(grp4)[1]
grp5 <- class.test$`5`
n5 <- dim(grp5)[1]
grp6 <- class.test$`6`
n6 <- dim(grp6)[1]
grp7 <- class.test$`7`
n7 <- dim(grp7)[1]
N <- sum(n1,n2,n3,n4,n5,n6,n7)
#calculate distance of each entry from groups 
classified <- c(1:dim(test.spect2)[1])
for (i in 1:dim(test.spect2)[1]) {
  dist.g <- matrix(nrow = G, ncol = J)
  for(g in 1:G) {
    for (j in 1:J) {
      dist.g[g,j] <- sum(c((t(alpha[,j])*(((test.spect2[i,2:(k+1)]))-(groupmeans[,g])))^2))
    }
  }
  final.dist <- rowSums(dist.g)
  classified[i] <-  which.min(final.dist)
}
confusion.test <- table(test.spect2[,1], classified, dnn = c('Actual Group', 'Predicted Group'))
rownames(confusion.test) <- c('building','grass', 'pavement','pool','shadow','soil','tree')
colnames(confusion.test) <- c('building','grass', 'pavement','pool','shadow','soil','tree')
test.m <- as.matrix(confusion.test)
correct <- sum(diag(test.m))
msclass <- sum(test.m)-correct
err_rate <- msclass/(correct+msclass)
msclass
correct+msclass
err_rate
T_count <- c(n1,n2,n3,n4,n5,n6,n7)
correct_g <- vector()
for (g in 1:G) {
  correct_g[g] <- test.m[g,g]/T_count[g]
}
test.m <- rbind(confusion.test, correct_g, T_count)
write.csv(test.m, "confusion_test2.csv")
#classify training data to compare prediction accuracy 
classified <- c(1:dim(scale.spect.super)[1])
for (i in 1:dim(scale.spect.super)[1]) {
  dist.g <- matrix(nrow = G, ncol = J)
  for(g in 1:G) {
    for (j in 1:J) {
      dist.g[g,j] <- sum(c((t(alpha[,j])*(((scale.spect.super[i,2:(k+1)]))-(groupmeans[,g])))^2))
    }
  }
  final.dist <- rowSums(dist.g)
  classified[i] <-  which.min(final.dist)
}
confusion.train <- table(scale.spect.super[,1], classified, dnn = c('Actual Group', 'Predicted Group'))
rownames(confusion.train) <- c('building','grass', 'pavement','pool','shadow','soil','tree')
colnames(confusion.train) <- c('building','grass', 'pavement','pool','shadow','soil','tree')
grp1 <- class$`1`
n1 <- dim(grp1)[1]
grp2 <- class$`2`
n2 <- dim(grp2)[1]
grp3 <- class$`3`
n3 <- dim(grp3)[1]
grp4 <- class$`4`
n4 <- dim(grp4)[1]
grp5 <- class$`5`
n5 <- dim(grp5)[1]
grp6 <- class$`6`
n6 <- dim(grp6)[1]
grp7 <- class$`7`
n7 <- dim(grp7)[1]
N <- sum(n1,n2,n3,n4,n5,n6,n7)
test.m <- as.matrix(confusion.train)
correct <- sum(diag(test.m))
msclass <- sum(test.m)-correct
err_rate <- msclass/(correct+msclass)
msclass
correct+msclass
err_rate
T_count <- c(n1,n2,n3,n4,n5,n6,n7)
correct_g <- vector()
for (g in 1:G) {
  correct_g[g] <- test.m[g,g]/T_count[g]
}
test.m <- rbind(confusion.train, correct_g, T_count)
write.csv(test.m, "confusion_train2.csv")
