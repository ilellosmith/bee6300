setwd("~/r/bee6300/hw2")
#read in files
gages_data <- read.csv('gages2.data.csv',header=T)
#center and scale each column
X.final <- apply(gages_data,2,function(x){scale(x)[,1]})
X.final <- data.frame(X.final)
#fit standard linear regression for average annual runoff vs. other covariates 
regular.lm <- lm(runoff~.,data=X.final)
table <- summary(regular.lm)
#install tools to export tables 
install.packages("devtools")
library("devtools", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
install_github("dgrtwo/broom")
library("broom", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
regtable <- tidy(regular.lm)
#regtable
regtable2 <-tidy(table)
write.csv(regtable, "tidy_lmfit.csv")
write.csv(regtable2, "regtable2.csv")
#calculate pearson correlation and VIF
pcor <- cor(X.final)
pcortable <- tidy(pcor)
#pcortable
write.csv(pcortable, "tidy_pcortable.csv")
library("car", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
VIF <- vif(regular.lm)
viftable <- tidy(VIF)
#viftable
write.csv(viftable, "tidy_viftable.csv")
library("glmnet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
#run ridge and lasso regressions
x <- data.matrix(X.final[,2:ncol(X.final)]) #all of the covariates
y <- as.vector(X.final[,1]) #annual flow
"
for (i in 1:20) {
#cv.lasso <- cv.glmnet(x,y, nfolds = 10, alpha = 1) 
#cv.lasso
#print(min(cv.lasso$lambda))
#plot(cv.lasso)
cv.lasso <- cv.glmnet(x,y, alpha = 1) # nfolds=10,
print(min(cv.lasso$lambda))
#plot(cv.lasso)
}
"
#cross validation lasso
cv.lasso <- cv.glmnet(x,y, alpha = 1) # nfolds=10,
lambda.lasso <-cv.lasso$lambda.min
print(lambda.lasso)
plot(cv.lasso, main = "lasso")
#cross validation ridge
cv.ridge <- cv.glmnet(x,y, alpha = 0)
#cv.ridge
lambda.ridge <- cv.ridge$lambda.min
print(lambda.ridge)
plot(cv.ridge, main = "ridge")
#lasso regression
lasso_reg <- glmnet(x,y,lambda = lambda.lasso, alpha = 1)
lasso_table <- tidy(lasso_reg)
coef_lasso <- tidy(coef(lasso_reg))
write.csv(lasso_table, "lasso_table.csv")
write.csv(coef_lasso, "coef_lasso.csv")
#ridge regression
ridge_reg <- glmnet(x,y,lambda = lambda.ridge, alpha = 0)
ridge_table <-tidy(ridge_reg)
coef_ridge <- tidy(coef(ridge_reg))
write.csv(ridge_table, "ridge_table.csv")
write.csv(coef_ridge, "coef_ridge.csv")
#run ridge, lasso and OLS and calc rmse n times
n <- 100
ols_rmse <- vector (length = n)
lasso_rmse <- vector (length = n)
ridge_rmse <- vector (length = n)
for (i in 1:n){
#split the data
all <- 1:dim(X.final)[1]
subset1 <- sample(all,size=round(length(all)/2),replace=F)
subset2 <- all[-subset1]
#OLS with subset1, calculate root mean squared error
regular.lm <- lm(runoff~.,data=X.final[subset1,])
lm.pred <- predict(regular.lm,newdata=X.final[subset2,])
rmse.lm <- sqrt(mean((lm.pred-X.final$runoff[subset2])^2))
ols_rmse[i] <-rmse.lm
#Select lambdas
x <- data.matrix(X.final[subset1,2:ncol(X.final)])
y <- as.vector(X.final[subset1,1])
#cross validation lasso
cv.lasso <- cv.glmnet(x,y, alpha = 1) 
lambda.lasso <-cv.lasso$lambda.min
#print(lambda.lasso)
#cross validation ridge
cv.ridge <- cv.glmnet(x,y, alpha = 0)
lambda.ridge <- cv.ridge$lambda.min
#print(lambda.ridge)
#run lasso and ridge regressions, calc RMSE
lasso.lm <- glmnet(x,y,lambda = lambda.lasso, alpha = 1)
ridge.lm <- glmnet(x,y,lambda = lambda.ridge, alpha = 0)
newx <- data.matrix(X.final[subset2,2:ncol(X.final)])
lasso.pred <- predict(lasso.lm,newx=newx)
rmse.lasso <- sqrt(mean((lasso.pred-X.final$runoff[subset2])^2))
lasso_rmse[i] <-rmse.lasso
ridge.pred <- predict(ridge.lm,newx=newx)
rmse.ridge <- sqrt(mean((ridge.pred-X.final$runoff[subset2])^2))
ridge_rmse[i] <-rmse.ridge
}
#library("ggplot2", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
par(mfrow = c(1,3))
std_y <- c(0.5,1.2)
boxplot(ols_rmse, ylim = std_y, main = "OLS RMSE")
boxplot(lasso_rmse, ylim = std_y, main = "Lasso RMSE")
boxplot(ridge_rmse, ylim = std_y, main = "Ridge RMSE")




 
