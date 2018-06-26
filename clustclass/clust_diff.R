#doc: this R script clusters stream catchments using single linkage, complete linkage and centroid methods. 
#It compares the mean properties of classified groups and whether or not they contain brook trout. 
#Using K-means clustering, the script assesses how many groups are appropriate to parsimoniously explain as much variance as 
#possible. It then uses Fischer's linear discrimination to classify new stream catchments into two groups - likely to contain
#brook trout or unlikely to contain brook trout. 
setwd("~/r/bee6300/hw5")
covariates <- read.csv("Maryland.Trout.csv")
#scale covariates
scale.bkt <- (covariates$bkt-mean(covariates$bkt))/sd(covariates$bkt)
scale.X.Ag <- (covariates$X.Ag-mean(covariates$X.Ag))/sd(covariates$X.Ag)
scale.RIFFQUAL <- (covariates$RIFFQUAL-mean(covariates$RIFFQUAL))/sd(covariates$RIFFQUAL)
scale.LOG_RD <- (covariates$LOG_RD-mean(covariates$LOG_RD))/sd(covariates$LOG_RD)
scale.TEMP_FLD <- (covariates$TEMP_FLD-mean(covariates$TEMP_FLD))/sd(covariates$TEMP_FLD)
scale.DO_FLD <- (covariates$DO_FLD-mean(covariates$DO_FLD))/sd(covariates$DO_FLD)
scale.covariates <- matrix(c(scale.X.Ag, scale.RIFFQUAL, scale.LOG_RD, scale.TEMP_FLD,
                             scale.DO_FLD), ncol = 4, nrow = 84)
#calculate distance matrix
d <- dist(scale.covariates)
#plot dendrograms 
myclust <- hclust(d,method= "single")
plot(myclust, main = "Cluster Dendrogram Single Linkage")
myclust <- hclust(d,method="complete")
plot(myclust,  main = "Cluster Dendrogram Complete Linkage")
myclust <- hclust(d,method="centroid")
plot(myclust,  main = "Cluster Dendrogram Centroid")
#cut complete linkage tree
myclust <- hclust(d,method="complete")
complete_twotree <- cutree(myclust, k = 2, h = NULL)
#compare two resulting groups' covariates
grp1 <- which(complete_twotree == 1)
grp2 <- which(complete_twotree == 2)
mean(covariates$X.Ag[grp1])
mean(covariates$X.Ag[grp2])
mean(covariates$RIFFQUAL[grp1])
mean(covariates$RIFFQUAL[grp2])
mean(covariates$LOG_RD[grp1])
mean(covariates$LOG_RD[grp2])
mean(covariates$TEMP_FLD[grp1])
mean(covariates$TEMP_FLD[grp2])
mean(covariates$DO_FLD[grp1])
mean(covariates$DO_FLD[grp2])
#compare brook trout
bkt_grp1_sum <- sum(covariates$bkt[grp1])
bkt_grp1_sum
bkt_grp2_sum <- sum(covariates$bkt[grp2])
bkt_grp2_sum
bkt_grp1_rate <- sum(covariates$bkt[grp1])/length(grp1)
bkt_grp1_rate 
bkt_grp2_rate <- sum(covariates$bkt[grp2])/length(grp2)
bkt_grp2_rate
covariates$bkt[grp1]
covariates$bkt[grp2]
#cluster with K-means clustering
k <- c(2:10)
var_explained <- vector()
var_explained[1] = 0
for (i in k) {
  clust.k <- kmeans(scale.covariates, centers = i, nstart=10)
  var_explained [i] <- clust.k$betweenss/clust.k$totss
}
clusters <- c(1,k)
plot(clusters, var_explained, xlab ="number of clusters", ylab = "variance explained")
#cluster with k=2
clust.k <- kmeans(scale.covariates, centers = 2, nstart=10)
#compare two resulting groups' covariates
grp1 <- which(clust.k$cluster == 1)
grp2 <- which(clust.k$cluster == 2)
mean(covariates$X.Ag[grp1])
mean(covariates$X.Ag[grp2])
mean(covariates$RIFFQUAL[grp1])
mean(covariates$RIFFQUAL[grp2])
mean(covariates$LOG_RD[grp1])
mean(covariates$LOG_RD[grp2])
mean(covariates$TEMP_FLD[grp1])
mean(covariates$TEMP_FLD[grp2])
mean(covariates$DO_FLD[grp1])
mean(covariates$DO_FLD[grp2])
#compare two resulting groups' brook trout
bkt_grp1_sum <- sum(covariates$bkt[grp1])
bkt_grp1_sum
bkt_grp2_sum <- sum(covariates$bkt[grp2])
bkt_grp2_sum
bkt_grp1_rate <- sum(covariates$bkt[grp1])/length(grp1)
bkt_grp1_rate 
bkt_grp2_rate <- sum(covariates$bkt[grp2])/length(grp2)
bkt_grp2_rate
covariates$bkt[grp1]
covariates$bkt[grp2]
#build a predictive trout model using classification/discrimination
#split data set into training and test data
train <- covariates[1:60,]
test <- covariates[61:84,]
#split training data into bkt = T and bkt = F
bkt1 <- train[which(train$bkt == 1),]
bkt1 <- bkt1[,3:7]
bkt0 <- train[which(train$bkt == 0),]
bkt0 <- bkt0[,3:7]
bkt1_obs <- length(bkt1$X.Ag)
bkt0_obs <- length(bkt0$X.Ag)
#calculate pooled covariance matrix
spool <- ((bkt1_obs-1)/(bkt1_obs+bkt0_obs-2))*cov(bkt1) + ((bkt0_obs-1)/(bkt0_obs+bkt1_obs-2))*cov(bkt0)
#estimate alpha 
grp_mean_bkt1 <- c(mean(bkt1$X.Ag),mean(bkt1$RIFFQUAL),mean(bkt1$LOG_RD),mean(bkt1$TEMP_FLD),mean(bkt1$DO_FLD))
grp_mean_bkt0 <- c(mean(bkt0$X.Ag),mean(bkt0$RIFFQUAL),mean(bkt0$LOG_RD),mean(bkt0$TEMP_FLD),mean(bkt0$DO_FLD))
alpha <- solve(spool)%*%(grp_mean_bkt1-grp_mean_bkt0)
#calculate m discrimination value
m <- 0.5*(t(grp_mean_bkt1-grp_mean_bkt0)%*%solve(spool)%*%((grp_mean_bkt1+grp_mean_bkt0)))
#classify testing data 
test <- test[,3:7]
classified <- vector()
for (i in 1:length(test[,1])){
  delta <- t(alpha)%*%t(test[i,])
  if (delta > m) {
    classified[i] <- 1
  }
  else {classified[i] <- 0}
}
#compare estimates with observations
obs <- covariates[61:84,]
obs$bkt
classified == obs$bkt
#trout present, estimate no trout present 
correct_notrout = 0
incorrect_trout = 0
incorrect_notrout = 0 
correct_trout = 0 
for (i in 1:length(obs$bkt)){
  if (obs$bkt[i] == 0 && classified[i] == 0){
    correct_notrout = correct_notrout + 1}
  else if (obs$bkt[i] == 0 && classified[i] == 1){
    incorrect_trout = incorrect_trout + 1}
  else if (obs$bkt[i] == 1 && classified[i] == 0){
    incorrect_notrout = incorrect_notrout + 1}
  else {
    correct_trout = correct_trout + 1}
}
correct_notrout/length(obs$bkt)*100
incorrect_trout/length(obs$bkt)*100
incorrect_notrout/length(obs$bkt)*100
correct_trout/length(obs$bkt)*100
