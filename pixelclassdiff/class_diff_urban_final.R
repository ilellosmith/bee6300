"This script takes pixel and superobject data from an aerial image of Deerfield Beach, FL
source(http://archive.ics.uci.edu/ml/datasets/Urban+Land+Cover#)
and applies Fischer's linear discrimination to classify testing data into each of 7 land type categories: 
building, grass, pavement (pavement includes asphalt, concrete, and car), pool, shadow, soil, and tree. 
The script applies Fischer's linear discrimination twice to compare model accuracy using just spectral information 
to model accuracy using spectral information and super object information"

##########################################################################
# globals used to establish workspace, data sets, covariates
#
WORKING_DIRECTORY <- "~/r/bee6300/f_proj"
#
TRAINING_DATA <- "training.csv"
#
TEST_DATA <- "testing.csv"
#
SPECTRAL_COVARIATES <- c('Mean_G', 'Mean_R', 'Mean_NIR', 'SD_G', 'SD_R','SD_NIR', 'GLCM1', 'NDVI')
#
FILL_COVARIATES <- 'Dens'
#
SHAPE_COVARIATES <- c('Area', 'Round', 'Bright', 'Compact', 'LW', 'Assym')
#
CLASSIFICATION <- 'class'
##########################################################################

# load dependencies, installing if necessary
REQUIRED_PACKAGES <- c("GGally", "ggplot2", "dplyr", "tidyr", 'PerformanceAnalytics', 'corrplot')
package.check <- lapply(REQUIRED_PACKAGES, FUN = function(x) {
  if (! require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

# read in and scale covariates
setwd(WORKING_DIRECTORY)
covariates <- as_tibble(read.csv(TRAINING_DATA)) %>%
  select(one_of(c(SHAPE_COVARIATES, FILL_COVARIATES, SPECTRAL_COVARIATES, CLASSIFICATION)))
scale_covariates <- as_tibble(cbind(scale(covariates %>% 
                                          select(one_of(SHAPE_COVARIATES, FILL_COVARIATES, SPECTRAL_COVARIATES)),
                                          center = T, scale =T),
                                          covariates %>% select(CLASSIFICATION)))

####### PART 1: run classification with just spectral covariates
spect <-  scale_covariates %>% select(one_of(c(SPECTRAL_COVARIATES, CLASSIFICATION))) %>% group_by(class)

# check for similar covariance structures 
grp_count <- spect %>% summarize(n = n())
N <- sum(grp_count$n)
G <- nrow(grp_count)
split_tibble <- function(tibble, col = 'col') tibble %>% split(.[,col]) # function to split tibble into list of classes
cov_tbl <- spect %>% do(data.frame((cov(.[,1:(ncol(spect)-1)]))))  # calculate covariance matrix for each group
cov_list <- cov_tbl %>% split_tibble(col = 'class')  # create list of covariance matrices by group
cov_list_no_class <- lapply(cov_list, function(x) x[,-1]) # drop the first column (class label) from each group

# plot covariance matrices for each group
dev.new()
par(mfrow = c(2,4))
mapply(x = cov_list_no_class, n = names(cov_list_no_class), FUN = function(x, n, i) 
{a <- n[i]
m <- cov2cor(as.matrix(as.data.frame(x), row.names = SPECTRAL_COVARIATES))
corrplot(m, type = "upper", method = "color", title = a, mar=c(0,0,2,0))
}) 

# calculate pooled covariance, group and grand means
cov_list_n_no_class <- lapply(1:nrow(grp_count), function(i, d, g) d[i]*(g[i,2]-1), d = cov_list_no_class, g = grp_count)
sum_cov <- Reduce('+', cov_list_n_no_class) # sum (n-1)*covariance matrix for all groups
spool <- 1/(N-G)*sum_cov # calculate pooled covariance matrix
group_means <- t(data.matrix((spect %>% summarize_all(.funs = mean)))[,-1]) # scaled group means
grand_means <- t((as.matrix(spect %>% ungroup() %>% summarize_at(SPECTRAL_COVARIATES, mean)))) # scaled grand means, should be ~ 0

# calculate between-groups covariance matrix Sb
Sb <- 1/(G-1)*((group_means)%*%t(group_means)) 
sig_noise <- solve(spool)%*%Sb # signal to noise matrix
J <- G-1 # number of principal directions (eigenvectors) in sig_noise

# calculate alphas 
alpha <- matrix(nrow=nrow(sig_noise),ncol=J)
for (i in 1:J){
  eigenvect <- as.vector(sig_noise[,i])
  a <- eigenvect/c((t(eigenvect)%*%data.matrix(spool)%*%eigenvect)^(1/2)) 
  alpha[,i] <- t(a)
}
write.csv(alpha, 'alpha1.csv')

# prepare test data 
covariates_test <- as_tibble(read.csv(TEST_DATA)) %>%
  select(one_of(c(SHAPE_COVARIATES, FILL_COVARIATES, SPECTRAL_COVARIATES, CLASSIFICATION)))
scale_covariates_test <- as_tibble(cbind(scale(covariates_test %>% 
                                            select(one_of(SHAPE_COVARIATES, FILL_COVARIATES, SPECTRAL_COVARIATES)),
                                          center = T, scale =T),
                                    covariates_test %>% select(CLASSIFICATION)))
spect_test <- scale_covariates_test %>% select(one_of(c(SPECTRAL_COVARIATES, CLASSIFICATION))) %>% group_by(class)
grp_count_test <- spect_test %>% summarize(n = n())
N_test <- sum(grp_count_test$n)
spect_test_matrix <- data.matrix(spect_test[,-9])

# calculate distance of each spectral measurement from groups, classify measurement into closest group  
classified <- c(1:nrow(spect_test))
for (i in 1:nrow(spect_test)) {
  dist_g <- matrix(nrow = G, ncol = J)
  for(g in 1:G) {
    for (j in 1:J) {
      dist_g[g,j] <- sum(c((t(alpha[,j])*((spect_test_matrix[i,])-(group_means[,g]))))^2)
    }
  }
  final_dist <- rowSums(dist_g)
  classified[i] <-  which.min(final_dist)
}

# assess error rate in classification of test data
confusion_test <- table(spect_test$class, classified, dnn = c('Actual Group', 'Predicted Group'))
colnames(confusion_test) <- c('building','grass', 'pavement','pool','shadow','soil','tree')
confusion_test <- as.matrix(confusion_test) # matrix of predicted group vs. actual
correct <- sum(diag(confusion_test)) # total correct classifications
msclass <- sum(confusion_test)-correct # total misclassifications
err_rate <- msclass/(correct+msclass) 
print(paste('missclassified count:', msclass, ', total:', correct+msclass, ', error rate:', err_rate)) 

# calculate correct prediction rate by group
correct_g <- unlist(sapply(1:G, function(g, diag, count) diag[g,g]/count[g,2], diag = confusion_test, count = grp_count_test))
rbind(confusion_test, correct_g, t(grp_count_test[2])) %>% write.csv("confusion_test.csv") #output results in csv

# classify training data to compare prediction accuracy 
no_group_spect <- data.matrix(ungroup(spect))[,-9]
classified <- c(1:nrow(spect))
for (i in 1:nrow(spect)) {
  dist_g <- matrix(nrow = G, ncol = J)
  for(g in 1:G) {
    for (j in 1:J) {
      dist_g[g,j] <- sum(c((t(alpha[,j])*((no_group_spect[i,])-(group_means[,g]))))^2)
    }
  }
  final_dist <- rowSums(dist_g)
  classified[i] <-  which.min(final_dist)
}

# assess error rate in classification of training data
confusion_train <- table(spect$class, classified, dnn = c('Actual Group', 'Predicted Group'))
colnames(confusion_train) <- c('building','grass', 'pavement','pool','shadow','soil','tree')
confusion_train <- as.matrix(confusion_train) # matrix of predicted group vs. actual
correct <- sum(diag(confusion_train)) # total correct classifications
msclass <- sum(confusion_train)-correct # total misclassifications
err_rate <- msclass/(correct+msclass) 
print(paste('missclassified count:', msclass, ', total:', correct+msclass, ', error rate:', err_rate)) 

# calculate correct prediction rate by group
correct_g <- unlist(sapply(1:G, function(g, diag, count) diag[g,g]/count[g,2], diag = confusion_train, count = grp_count))
rbind(confusion_train, correct_g, t(grp_count[2])) %>% write.csv("confusion_train.csv") #output results in csv

####### PART 2: run analysis with spectral and super-object covariates
spect_so <- scale_covariates %>% 
  select(one_of(c(SPECTRAL_COVARIATES, FILL_COVARIATES, SHAPE_COVARIATES, CLASSIFICATION))) %>% 
  group_by(class)

# check for similar covariance structures 
grp_count <- spect_so %>% summarize(n = n())
N <- sum(grp_count$n)
G <- nrow(grp_count)
split_tibble <- function(tibble, col = 'col') tibble %>% split(.[,col]) # function to split tibble into list of classes
cov_tbl <- spect_so %>% do(data.frame((cov(.[,1:(ncol(spect_so)-1)]))))  # calculate covariance matrix for each group
cov_list <- cov_tbl %>% split_tibble(col = 'class')  # create list of covariance matrices by group
cov_list_no_class <- lapply(cov_list, function(x) x[,-1]) # drop the first column (class label) from each group

# plot covariance matrices for each group
dev.new()
par(mfrow = c(2,4))
mapply(x = cov_list_no_class, n = names(cov_list_no_class), FUN = function(x, n, i) 
{a <- n[i]
m <- cov2cor(as.matrix(as.data.frame(x), row.names = c(SPECTRAL_COVARIATES,FILL_COVARIATES, SHAPE_COVARIATES)))
corrplot(m, type = "upper", method = "color", title = a, mar=c(0,0,2,0))
}) 

# calculate pooled covariance, group and grand means
cov_list_n_no_class <- lapply(1:nrow(grp_count), function(i, d, g) d[i]*(g[i,2]-1), d = cov_list_no_class, g = grp_count)
sum_cov <- Reduce('+', cov_list_n_no_class) # sum (n-1)*covariance matrix for all groups
spool <- 1/(N-G)*sum_cov # calculate pooled covariance matrix
group_means <- t(data.matrix((spect_so %>% summarize_all(.funs = mean)))[,-1]) # scaled group means
grand_means <- t((as.matrix(spect_so %>% ungroup() %>% summarize_at(c(SPECTRAL_COVARIATES, FILL_COVARIATES, SHAPE_COVARIATES), mean)))) # scaled grand means, should be ~ 0

# calculate between-groups covariance matrix Sb
Sb <- 1/(G-1)*((group_means)%*%t(group_means)) 
sig_noise <- solve(spool)%*%Sb # signal to noise matrix
J <- G-1 # number of principal directions (eigenvectors) in sig_noise

# calculate alphas 
alpha <- matrix(nrow=nrow(sig_noise),ncol=J)
for (i in 1:J){
  eigenvect <- as.vector(sig_noise[,i])
  a <- eigenvect/c((t(eigenvect)%*%data.matrix(spool)%*%eigenvect)^(1/2)) 
  alpha[,i] <- t(a)
}
write.csv(alpha, 'alpha2.csv')

# prepare test data 
covariates_test <- as_tibble(read.csv(TEST_DATA)) %>%
  select(one_of(c(SHAPE_COVARIATES, FILL_COVARIATES, SPECTRAL_COVARIATES, CLASSIFICATION)))
scale_covariates_test <- as_tibble(cbind(scale(covariates_test %>% 
                                                 select(one_of(SHAPE_COVARIATES, FILL_COVARIATES, SPECTRAL_COVARIATES)),
                                               center = T, scale =T),
                                         covariates_test %>% select(CLASSIFICATION)))
spect_so_test <- scale_covariates_test %>% select(one_of(c(SPECTRAL_COVARIATES, FILL_COVARIATES, SHAPE_COVARIATES, CLASSIFICATION))) %>% group_by(class)
grp_count_test <- spect_so_test %>% summarize(n = n())
N_test <- sum(grp_count_test$n)
spect_so_test_matrix <- data.matrix(spect_so_test[,-ncol(spect_so_test)])

# calculate distance of each spectral & super object measurement from groups, classify measurement into closest group
classified <- c(1:nrow(spect_so_test))
for (i in 1:nrow(spect_so_test)) {
  dist_g <- matrix(nrow = G, ncol = J)
  for(g in 1:G) {
    for (j in 1:J) {
      dist_g[g,j] <- sum(c((t(alpha[,j])*((spect_so_test_matrix[i,])-(group_means[,g]))))^2)
    }
  }
  final_dist <- rowSums(dist_g)
  classified[i] <-  which.min(final_dist)
}

# assess error rate in classification of test data
confusion_test <- table(spect_so_test$class, classified, dnn = c('Actual Group', 'Predicted Group'))
colnames(confusion_test) <- c('building','grass', 'pavement','pool','shadow','soil','tree')
confusion_test <- as.matrix(confusion_test) # matrix of predicted group vs. actual
correct <- sum(diag(confusion_test)) # total correct classifications
msclass <- sum(confusion_test)-correct # total misclassifications
err_rate <- msclass/(correct+msclass) 
print(paste('missclassified count:', msclass, ', total:', correct+msclass, ', error rate:', err_rate)) 

# calculate correct prediction rate by group
correct_g <- unlist(sapply(1:G, function(g, diag, count) diag[g,g]/count[g,2], diag = confusion_test, count = grp_count_test))
rbind(confusion_test, correct_g, t(grp_count_test[2])) %>% write.csv("confusion_test2.csv") #output results in csv

# classify training data to compare prediction accuracy 
no_group_spect_so <- data.matrix(ungroup(spect_so))[,-ncol(spect_so_test)]
classified <- c(1:nrow(spect_so))
for (i in 1:nrow(spect_so)) {
  dist_g <- matrix(nrow = G, ncol = J)
  for(g in 1:G) {
    for (j in 1:J) {
      dist_g[g,j] <- sum(c((t(alpha[,j])*((no_group_spect_so[i,])-(group_means[,g]))))^2)
    }
  }
  final_dist <- rowSums(dist_g)
  classified[i] <-  which.min(final_dist)
}

# assess error rate in classification of training data
confusion_train <- table(spect_so$class, classified, dnn = c('Actual Group', 'Predicted Group'))
colnames(confusion_train) <- c('building','grass', 'pavement','pool','shadow','soil','tree')
confusion_train <- as.matrix(confusion_train) # matrix of predicted group vs. actual
correct <- sum(diag(confusion_train)) # total correct classifications
msclass <- sum(confusion_train)-correct # total misclassifications
err_rate <- msclass/(correct+msclass) 
print(paste('missclassified count:', msclass, ', total:', correct+msclass, ', error rate:', err_rate)) 

# calculate correct prediction rate by group
correct_g <- unlist(sapply(1:G, function(g, diag, count) diag[g,g]/count[g,2], diag = confusion_train, count = grp_count))
rbind(confusion_train, correct_g, t(grp_count[2])) %>% write.csv("confusion_train2.csv") #output results in csv
