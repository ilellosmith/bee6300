setwd("~/r/bee6300/hw4")
load("PDSI.Final.RData")
load("mylon.RData")
load("mylat.RData")
#remove columns with no data 
dim(PDSI.final)
#track indices of columns with complete data
not.missing <- !is.na(colSums(PDSI.final))
PDSI.final.nonan <- PDSI.final[,not.missing == TRUE]
#get lat/long values for columns with complete data
lon.lat <- expand.grid('lon'=mylon,'lat'=mylat) 
lon.lat.nonan <- lon.lat[not.missing,]
#run PCA
pca.PDSI <- prcomp(PDSI.final.nonan, center = TRUE, scale. = TRUE, rank. = 10, retx = T)
rotate <- data.frame(pca.PDSI$rotation)
sdev <- data.frame(pca.PDSI$sdev)
variance <- data.frame((sdev[1:10,])^2)
write.csv(variance, file = "variance.csv")
n <- length(t(sdev))
s.error <- variance*sqrt(2/n)
#create scree plot with error bars
require(ggplot2)
ggplot(variance, aes(c(1:10),variance)) + geom_point() + geom_line() +
  geom_errorbar(aes(ymin = (variance-s.error), ymax = (variance +s.error)), width=0.3,position=position_dodge(.9)) + ggtitle("scree plot") 
require(maps)
mybreaks <- c(seq(-0.04,-0.001,by=0.004),seq(0.001,0.04,by=0.004))
#assign each element of a particular EOF to one of the categories based on the intervals above
eof.retain <- rotate[,1:4]
mycat1 <- as.numeric(cut(eof.retain$PC1,breaks=mybreaks))
mycat2 <- as.numeric(cut(eof.retain$PC2,breaks=mybreaks))
mycat3 <- as.numeric(cut(eof.retain$PC3,breaks=mybreaks))
mycat4 <- as.numeric(cut(eof.retain$PC4,breaks=mybreaks))
#create a color palette that spans from red to blue
mycolpal <- colorRampPalette(c("red","white","blue"))
#use the color palette and the assigned categories to assign a color to each loading value in an EOF
mycol1 <- mycolpal(length(mybreaks)-1)[mycat1]
mycol2 <- mycolpal(length(mybreaks)-1)[mycat2]
mycol3 <- mycolpal(length(mybreaks)-1)[mycat3]
mycol4 <- mycolpal(length(mybreaks)-1)[mycat4]
#plot categories at their lat/long with associated color  
plot(lon.lat.nonan$lon,lon.lat.nonan$lat,col = mycol1, main = "EOF1", xlab = "longitude", ylab = "latitude") + map("world",add=T)
plot(lon.lat.nonan$lon,lon.lat.nonan$lat,col = mycol2, main = "EOF2", xlab = "longitude", ylab = "latitude") + map("world",add=T)
plot(lon.lat.nonan$lon,lon.lat.nonan$lat,col = mycol3, main = "EOF3", xlab = "longitude", ylab = "latitude") + map("world",add=T)
plot(lon.lat.nonan$lon,lon.lat.nonan$lat,col = mycol4, main = "EOF4", xlab = "longitude", ylab = "latitude") + map("world",add=T)
#forecast next year's PDSI
#create simple autoregressive lag-1 time series for each PC
rot.data <- data.frame(pca.PDSI$x[,1:4])
n <- length(rot.data$PC1)
PC1 <-rot.data$PC1
PC2 <-rot.data$PC2
PC3 <-rot.data$PC2
PC4 <-rot.data$PC2
lmPC1 <- lm(PC1[2:n] ~ PC1[1:(n-1)])
slope1 <- lmPC1$coefficients[2]
intercept1 <- lmPC1$coefficients[1]
lmPC2 <- lm(PC2[2:n] ~ PC2[1:(n-1)])
slope2 <- lmPC2$coefficients[2]
intercept2 <- lmPC2$coefficients[1]
lmPC3 <- lm(PC3[2:n] ~ PC3[1:(n-1)])
slope3 <- lmPC3$coefficients[2]
intercept3 <- lmPC3$coefficients[1]
lmPC4 <- lm(PC4[2:n] ~ PC4[1:(n-1)])
slope4 <- lmPC4$coefficients[2]
intercept4 <- lmPC4$coefficients[1]
#project U for year 2006
est1 <- slope1*PC1[n] + intercept1
est2 <- slope2*PC2[n] + intercept2
est3 <- slope3*PC3[n] + intercept3
est4 <- slope4*PC4[n] + intercept4
PDSI.predict <- c(est1,est2,est3,est4)
U <- PDSI.predict
W <- rotate[,1:4]
#Synthesis formula 
X <- U%*%t(W)
#Plot PDSI forecast 
#reset mybreaks
mybreaks <- c(seq(-1,-0.1,by=0.1),seq(0.1,1,by=0.1))
#use the color palette and the assigned categories to assign a color to each loading value in an EOF
mycat <- as.numeric(cut(X,breaks=mybreaks))
mycol <- mycolpal(length(mybreaks)-1)[mycat]
#use the color palette and the assigned categories to assign a color to each loading value in an EOF
plot(lon.lat.nonan$lon,lon.lat.nonan$lat, col = mycol, main = "PDSI Prediction for 2006", xlab = "longitude", ylab = "latitude") + map("world",add=T)
#Compare EOFs and Rotated EOFs
#generate rotation matrix
my.varimax <-varimax(pca.PDSI$rotation[,1:4])
#multiply data by reofs
dim(PDSI.final.nonan)
dim(my.varimax$rotmat)
reofs <- (pca.PDSI$rotation[,1:4]%*%my.varimax$rotmat)
reofs.data <- PDSI.final.nonan%*%reofs
#correlation matrix of original PCs
pc.cor.orig <- cor(pca.PDSI$x[,1:4])
pc.cor.orig
write.csv(pc.cor.orig, file = "correlation_original_PCs.csv")
#correlation matrix of rotated PCs
pc.cor.rot <- cor(reofs.data)
pc.cor.rot
write.csv(pc.cor.rot, file = "correlation_rotated_PCs.csv")
#plot rotated EOFs
reofs.std <- (reofs - mean(reofs))/sd(reofs)
low <- round(range(reofs.std)[1],1)
high <- round(range(reofs.std)[2],1)
intvl <- abs((high - low)/20)
mybreaks <- c(seq(low,(0-intvl),by=intvl),seq(0,high,by=intvl))
#assign each element of a particular EOF to one of the categories based on the intervals above
#keep the rotated EOFS on the same scale as the originals!
#^tried this, but found faded versions of patterns with modified scale. Code is commented below
"mybreaks <- c(seq(-0.04,-0.001,by=0.004),seq(0.001,0.04,by=0.004))
mycat1 <- as.numeric(cut(reofs[,1],breaks=mybreaks))
mycat2 <- as.numeric(cut(reofs[,2],breaks=mybreaks))
mycat3 <- as.numeric(cut(reofs[,3],breaks=mybreaks))
mycat4 <- as.numeric(cut(reofs[,4],breaks=mybreaks))
"
mycat1 <- as.numeric(cut(reofs.std[,1],breaks=mybreaks))
mycat2 <- as.numeric(cut(reofs.std[,2],breaks=mybreaks))
mycat3 <- as.numeric(cut(reofs.std[,3],breaks=mybreaks))
mycat4 <- as.numeric(cut(reofs.std[,4],breaks=mybreaks))
#use the color palette and the assigned categories to assign a color to each loading value in an EOF
mycol1 <- mycolpal(length(mybreaks)-1)[mycat1]
mycol2 <- mycolpal(length(mybreaks)-1)[mycat2]
mycol3 <- mycolpal(length(mybreaks)-1)[mycat3]
mycol4 <- mycolpal(length(mybreaks)-1)[mycat4]
#plot categories at their lat/long with associated color  
plot(lon.lat.nonan$lon,lon.lat.nonan$lat,col = mycol1, main = "rotated EOF1", xlab = "longitude", ylab = "latitude") + map("world",add=T)
plot(lon.lat.nonan$lon,lon.lat.nonan$lat,col = mycol2, main = "rotated EOF2", xlab = "longitude", ylab = "latitude") + map("world",add=T)
plot(lon.lat.nonan$lon,lon.lat.nonan$lat,col = mycol3, main = "rotated EOF3", xlab = "longitude", ylab = "latitude") + map("world",add=T)
plot(lon.lat.nonan$lon,lon.lat.nonan$lat,col = mycol4, main = "rotated EOF4", xlab = "longitude", ylab = "latitude") + map("world",add=T)