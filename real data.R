if (!requireNamespace("BiocManager", quietly=TRUE))
       install.packages("BiocManager")
BiocManager::install("transcriptR")

BiocManager::install("impute")

library(impute)
library(rgl)
library(misc3d)
library(ggpattern)
require(gridExtra)


### import source code from other files
source("ACS_final.r")
source("ACS_prepare.r")
source("kernel_matrix_prepare.r")

### import data
trainx<-t(data.matrix(GSE24417_Training_DataMatrix[,-1]))
trainx<-impute.knn(trainx ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)$data
trainy<-trainy$x
trainx<-scale(trainx)
trainx<-atan(trainx)

### ACS estimate
kernel<-dr_kernel(trainx,trainy,2)
kernel<-cbind(kernel,sir_kernel(trainx,trainy,2)[,1])
result<-ACS_SDR(X,kernel,5,15) ### Better to manually choose moderate lambda between 0.06-0.12

###Take the index set of the selected columns
index<-result$index
#### code for drawing figure for ACS-method
outcome<-svd(result)$u[,1:2]
data1=data.frame(Class=factor(trainy,levels=c(0,1),labels=c("Negative Samples","Positive Samples")),x1=trainx%*%outcome[,1],x2=trainx%*%outcome[,2])
par(mar=c(4,4.5,1,1)+0.1)
gridRange <- apply(data1[c('x1', 'x2')], 2, range)
x1 <- seq(from = -2.07, to = 2, length = 150)
x2 <- seq(from = -2.27, to = 2, length = 150)
grid <- expand.grid(x1 = x1, x2 = x2)
fit<-svm(Class~x1+x2,data1,type="nu-classification")
grid$class <- predict(fit,grid)
decisionValues <- predict(fit, grid, decision.values = TRUE)
grid$z <- as.vector(attributes(decisionValues)$decision)
plot(data1[c('x1', 'x2')], pch = ifelse(data1$Class == "Positive Samples", 1, 19),xlab='Direction 1',ylab='Direction2',xlim=c(-1.9,2.3),ylim=c(-2.1,2.1),cex.axis=1.45,cex.lab=1.45)
contour(x1, x2, matrix(grid$z, length(x1), length(x2)), level=0, lwd = 1.5, drawlabels = FALSE, add=TRUE)
legend("topright",legend=c("Positive Samples", "Negative Samples"),pch=c(1,19), cex=1.2)
plot.new()


### use the selected columns to do the leave-one-out cross validation to reduce the computation time
index<-c()
countDR=0
countDRQDA=0
for(i in 1:216)
{
   kernel<-cbind(dr_kernel(trainx[-i,],trainy[-i],2),sir_kernel(trainx[-i,],trainy[-i],2))[,index]
   solution_equal<-coef(glmnet(x=trainx[-i,],y=kernel[,index],family="mgaussian",parallel=TRUE,lambda=0.08))
   temp<-do.call(cbind,(lapply(solution_equal,deal_list)))
   weight<-apply(temp,MARGIN=1,choose_weight)
   solution_equal<-Matrix_estimation(data=t(t(trainx[-i,])/weight),ytilde=kernel)
   result1<-do.call(cbind,(lapply(solution_equal,deal_list)))/weight
   outcome1<-svd(result1)$u[,1:2]
   train_data<-data.frame(Class=factor(trainy[-i],levels=c(0,1),labels=c("0","1")),x1=(trainx[-i,]%*%outcome)[,1],x2=(trainx[-i,]%*%outcome)[,2])
   valid_data<-data.frame(x1=(trainx[i,]%*%outcome)[,1],x2=(trainx[i,]%*%outcome)[,2])
   if(predict(qda(Class~.,train_data),valid_data)$class==trainy[i])
   {
      print("TRUE")
      countDRQDA=countDRQDA+1
   }else
   {
      print("FALSE")
   }
}

countSIRQDA=0
### since SC-SIR uses matlab, resultSCSIR is the file recording the leave-one-out estimate of central subspace from SC-SIR
for(i in 1:216)
{
   index_train<-c(1:216)[-i]
   index_test<-c(1:216)[i]
   outcome<-resultSCSIR[i,]
   train_data<-data.frame(Class=factor(trainy[index_train],levels=c(0,1),labels=c("0","1")),x1=(trainx[index_train,]%*%t(outcome)))
   valid_data<-data.frame(x1=(trainx[i,]%*%t(outcome)))
   if(predict(lda(Class~.,train_data),valid_data)$class==trainy[i])
   {
      print("TRUE")
      countSIRQDA=countSIRQDA+1
   }else
   {
      print("FALSE")
   }
}

resultSCSIR<-as.matrix(resultSCSIR)

#### code for drawing figure for SC-SIR
outcome1<-as.matrix(SC-SIRestimate)
data2<-data.frame(Class=factor(trainy,levels=c(0,1),labels=c("Negative Samples","Positive Samples")),x1=trainx%*%outcome1)

X1<-as.numeric(data2$x1[which(data2$Class=='negative')])
X2<-as.numeric(data2$x1[which(data2$Class=='positive')])
p1<-hist(X1,plot=FALSE)
p2<-hist(X2,plot=FALSE)
par(mar = c(4,4,1,7),xpd=TRUE)
plot(0,0,type="n",xlim=c(-0.7,0.8),ylim=c(0,2.5),ylab="Density",xlab="Central direction by SC-SIR",cex.lab=1.2)
plot(p1,freq=FALSE,density=15,angle=135,add=TRUE)
plot(p2,freq=FALSE,density=15,angle=45,add=TRUE)
s1<-seq(min(X1),max(X1),length=50)
fun1<-dnorm(s1,mean=mean(X1),sd=sd(X1))
lines(s1,fun1,col='black',lwd=2,lty=1)
s2<-seq(min(X2),max(X2),length=50)
fun2<-dnorm(s2,mean=mean(X2),sd=sd(X2))
lines(s2,fun2,col='black',lwd=2,lty=2)

#legend(x="topright",legend = c("positive","negative"),
#       density = c(15,15),
#       angle = c(45,135),inset = c(-0.25, 0.05),xpd = TRUE,
#       horiz = FALSE,cex=1)
#legend(x="bottomright",legend = c("positive","negative"),
#       lwd=2,lty=c(2,1),inset = c(-0.25, 0.65),xpd = TRUE,horiz = FALSE,cex=0.8438)


