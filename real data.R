if (!requireNamespace("BiocManager", quietly=TRUE))
       install.packages("BiocManager")
BiocManager::install("transcriptR")

BiocManager::install("impute")

library(impute)
library(rgl)
library(misc3d)
library(ggpattern)
require(gridExtra)

current_file_path<-rstudioapi::getSourceEditorContext()$path
current_dir<-dirname(current_file_path)
setwd(current_dir)
### import source code from other files
source("ACS_final.r")
source("ACS_prepare.r")
source("kernel_matrix_prepare.r")

### import data
trainx<-as.matrix(trainx)
trainy<-trainy$x


###SEAS estimate
result_SEAS=cv.seas(x=trainx, y=trainy, yclass=trainy, category=TRUE)
outcome1<-result_SEAS$beta

outcome1<-SEAS.SIR_estimate[[1]]

data2<-data.frame(Class=factor(trainy,levels=c(0,1),labels=c("Negative Samples","Positive Samples")),x1=trainx%*%outcome1)

X1 <- as.numeric(data2$x1[which(data2$Class == 'Negative Samples')])
X2 <- as.numeric(data2$x1[which(data2$Class == 'Positive Samples')])

# 直方图
p1 <- hist(X1, plot = FALSE, breaks = 15)
p2 <- hist(X2, plot = FALSE, breaks = 15)

# 设置图形参数
par(mar = c(5, 4.5, 4.5,1), xpd = TRUE)

# 绘制空白图
plot(0, 0, type = "n", xlim = c(-8, 8), ylim = c(0, 0.23), ylab = "", xlab = "", cex.axis = 1.5, cex.lab = 1.25)

mtext("Density", side = 2, line = 2.75, cex = 1.5)
mtext("The reduced predictor from sparse SIR by [46]", side = 1, line = 2.75, cex = 1.5)

# 绘制直方图
plot(p1, freq = FALSE, density = 15, angle = 135, col = "black", add = TRUE)
plot(p2, freq = FALSE, density = 15, angle = 45, col = "black", add = TRUE)

# 绘制密度曲线
s1 <- seq(min(X1), max(X1), length = 50)
fun1 <- dnorm(s1, mean = mean(X1), sd = sd(X1))
lines(s1, fun1, col = 'black', lwd = 2, lty = 1)

s2 <- seq(min(X2), max(X2), length = 50)
fun2 <- dnorm(s2, mean = mean(X2), sd = sd(X2))
lines(s2, fun2, col = 'black', lwd = 2, lty = 1)

par(mar = c(5, 4.5, 4.5,1), xpd = TRUE)
mtext(" ", side = 3, line = 3)

countSIRQDA=0
for(i in 1:216)
{
  index_train<-c(1:216)[-i]
  index_test<-c(1:216)[i]
  result_SEAS=cv.seas(x=trainx[index_train,], y=trainy[index_train], yclass=trainy[index_train], category=TRUE)
  outcome1<-result_SEAS$beta
  train_data<-data.frame(Class=factor(trainy[index_train],levels=c(0,1),labels=c("0","1")),x1=(trainx[index_train,]%*%outcome1))
  valid_data<-data.frame(x1=(trainx[i,]%*%outcome1))
  if(predict(qda(Class~.,train_data),valid_data)$class==trainy[i])
  {
    print("TRUE")
    countSIRQDA=countSIRQDA+1
  }else
  {
    print("FALSE")
  }
}


dr_kernel<-function(x,y,H)
{
  n <- dim(x)[1]
  p <-dim(x)[2]
  dy <- slicing(y,H)
  H <- dy$H
  ytilde <- dy$ytilde
  prop <- dy$prop
  xxbar<-numeric(0)
  freq<-rep(0,H)
  if(H>=3)
  {
    for (j in 2:(H-1))
    {
      for(i in 1:(j-1))
      {
        indj<-((y >= ytilde[j])&(y < ytilde[j+1]))
        indi<-((y >= ytilde[i])&(y < ytilde[i+1]))
        freq[j]=length(which(indj==TRUE))
        freq[i]=length(which(indi==TRUE))
        Sj<-matrix(0,ncol=n,nrow=freq[j])
        Sj[,indj]<-diag(freq[j])
        Si<-matrix(0,ncol=n,nrow=freq[i])
        Si[,indi]<-diag(freq[i])
        xxbar<-cbind(xxbar,sqrt(prop[j]*prop[i])*(2*x-1/prop[i]*t(Si)%*%Si%*%x-1/prop[j]*t(Sj)%*%Sj%*%x
                                                  +n*t(Si)%*%rep(1/freq[i],freq[i])%*%t(rep(1/freq[j],freq[j]))%*%Sj%*%x
                                                  +n*t(Sj)%*%rep(1/freq[j],freq[j])%*%t(rep(1/freq[i],freq[i]))%*%Si%*%x))
      }
    }
  }
  for(i in 1:(H-1))
  {
    indj<-(y >= ytilde[H])
    indi<-((y >= ytilde[i])&(y < ytilde[i+1]))
    freq[H]=length(which(indj==TRUE))
    freq[i]=length(which(indi==TRUE))
    Sj<-matrix(0,ncol=n,nrow=freq[H])
    Sj[,indj]<-diag(freq[H])
    Si<-matrix(0,ncol=n,nrow=freq[i])
    Si[,indi]<-diag(freq[i])
    xxbar<-cbind(xxbar,sqrt(prop[H]*prop[i])*(2*x-1/prop[i]*t(Si)%*%Si%*%x-1/prop[H]*t(Sj)%*%Sj%*%x
                                              +t(Si)%*%rep(1/freq[i],freq[i])%*%t(rep(1/freq[H],freq[H]))%*%Sj%*%x
                                              +t(Sj)%*%rep(1/freq[H],freq[H])%*%t(rep(1/freq[i],freq[i]))%*%Si%*%x))
  }
  return(xxbar)
}

Sparse_vector<-function(Y,design_matrix,nfolds=5,lambda=NULL)
{
  if(is.null(lambda))
  {
    fit=cv.glmnet(x=design_matrix,y=Y,nfolds=nfolds,standardize=FALSE, type.measure="mae",nlambda=100)
  }else
  {
    fit=cv.glmnet(x=design_matrix,y=Y,nfolds=nfolds,lambda=lambda,standardize=FALSE,type.measure="mae",nlambda=100)
  }
  return(as.numeric(coef(fit))[-1])
}


### ACS estimate
kernel<-dr_kernel(trainx,trainy,2)
kernel<-cbind(kernel,trainy)
sfInit(parallel = TRUE,cpus=15)
sfLibrary(glmnet)
sparse_estimate<-sfApply(kernel,margin=2,fun=Sparse_vector,design_matrix=trainx,nfolds=10,lambda=NULL)
sfStop()
index<-c()

for(i in 1:dim(sparse_estimate)[2])
{
  if(all(sparse_estimate[,i]==0))
  {
    index=append(index,i)
  }
}
index<-c(1:ncol(kernel))[-index]

sfInit(parallel = TRUE,cpus=15)
sfLibrary(glmnet)
sparse_estimate[,index]<-sfApply(kernel[,index],margin=2,fun=Sparse_vector,design_matrix=trainx,nfolds=5,lambda=NULL)
sfStop()
index<-c()

for(i in 1:dim(sparse_estimate)[2])
{
  if(all(sparse_estimate[,i]==0))
  {
    index=append(index,i)
  }
}
index<-c(1:ncol(kernel))[-index]
print(index)

selected_column_index<-forward_column_selection(sparse_estimate[,index],threshold=1e-12,D=15)
selected_index<-index[selected_column_index]

selected_index<-c(639,7417,8807,6273,2035,9971,9624,8900,10753,2194,6884,6734,7911,779,4679)

Matrix_estimation<-function(data,ytilde,nfolds=5)
{
  fit<-cv.glmnet(x=data,y=ytilde,family="mgaussian",nfolds=nfolds,lambda.min.ratio=1e-4,nlambda=1000,type.measure = "mae")
  print(fit$lambda.min)
  return(coef(fit,s=fit$lambda.min))
}  

solution_equal<-Matrix_estimation(data=trainx,ytilde=kernel[,selected_index])
temp=do.call(cbind,(lapply(solution_equal,deal_list)))
#solution_equal<-coef(glmnet(x=trainx,y=kernel[,selected_index],lambda=0.06,family="mgaussian"))
#temp=do.call(cbind,(lapply(solution_equal,deal_list)))
weight<-apply(temp,MARGIN=1,choose_weight)
solution_equal<-Matrix_estimation(data=t(t(trainx)/weight),ytilde=kernel[,selected_index])
result<-do.call(cbind,(lapply(solution_equal,deal_list)))/weight

svd(result)$d
outcome<-svd(result)$u[,1:5]

data1=data.frame(Class=factor(trainy,levels=c(0,1),labels=c("Negative Samples","Positive Samples")),x1=trainx%*%outcome[,1],x2=trainx%*%outcome[,2])
par(mar=c(5,4.5,4.5,1)+0.1)
gridRange <- apply(data1[c('x1', 'x2')], 2, range)
x1 <- seq(from = range(data1$x1)[1]-0.24, to = range(data1$x1)[2]+0.24, length = 150)
x2 <- seq(from = range(data1$x2)[1]-0.24, to = range(data1$x2)[2]+0.24, length = 150)
grid <- expand.grid(x1 = x1, x2 = x2)
fit<-svm(Class~x1+x2,data1,type="nu-classification")
grid$class <- predict(fit,grid)
decisionValues <- predict(fit, grid, decision.values = TRUE)
grid$z <- as.vector(attributes(decisionValues)$decision)
plot(data1[c('x1', 'x2')], pch = ifelse(data1$Class == "Positive Samples", 1, 19),xlab='The first proposed reduced predictor',ylab='The second proposed reduced predictor',xlim=c(range(data1$x1)[1]-0.1,range(data1$x1)[2]+0.1),ylim=c(range(data1$x2)[1]-0.1,range(data1$x2)[2]+0.1),cex.axis=1.45,cex.lab=1.45)
contour(x1, x2, matrix(grid$z, length(x1), length(x2)), level=0, lwd = 1.5, drawlabels = FALSE, add=TRUE)
which(abs(outcome[,1])>1e-7)


p=ncol(trainx)
n=nrow(trainx)


M=cov(trainx)
U=t(trainx)%*%kernel[,selected_index]/215

foldid <- c()
for(cnt in c(143,73)){
  foldid <- c(foldid, sample(rep(seq(5), length = cnt)))}
U_train=vector("list", length = 5)
for(k in 1:5)
{
  x_val <- trainx[foldid==k,,drop=FALSE]
  y_val <- trainy[foldid==k]
  x_train <- trainx[foldid!=k,,drop=FALSE]
  y_train <- trainy[foldid!=k]
  U_train[[k]]=t(x_train)%*%cbind(dr_kernel(x_train,y_train,2),y_train)[,selected_index]/215
}

lam2=0
lam1=seq(0.02,0.04,0.00025)
gam=c(0.01)

eval_all<-rep(0,50)
for(i in 1:50)
{
lam1_train=lam1[i]
eval_fold=0
for(k in 1:5)
{ # Cross-validation
  fit_fold <- admm(M, U_train[[k]], p,p, lam1=lam1_train, lam2=lam2, gam=gam, eps=1e-4, maxit=1e+5)
  B_l<-fit_fold$B
  rank_l <- fit_fold$rank
  step_l <- fit_fold$step
  time_l <- fit_fold$time
  eval_fold <- eval_dc(B_l, x_val, y_val)+eval_fold  # The evaluation: distance correlation.
  ind <- which(sapply(beta_l, is.null))
  rank_l[ind] <- -1
  eval_fold[ind] <- min(eval_fold, na.rm = TRUE)
}
eval_all[i]=eval_fold/5
print(eval_all[i])
}

cvm <- eval_all
gam_max=gam
cvm_max=which.max(cvm)
lam1_max=lam1[cvm_max]
fit<-admm(M, U, p,p, lam1_max, lam2_max, gam_max, eps=1e-10, maxit=1e+8)$B[[1]]
which(abs(fit[,1])>1e-7)
weight<-1/apply(as.matrix(fit),MARGIN=1, choose_weight)
which(weight!=0)

lam2=0
lam1=seq(0.02,0.04,0.00025)*max(weight)
gam=c(0.01)
eval_all<-rep(0,50)
M1=M*outer(weight,weight)
for(i in 1:50)
{
  lam1_train=lam1[i]
  eval_fold=0
  for(k in 1:5)
  { # Cross-validation
    fit_fold <- admm(M1, weight*U_train[[k]], p,p, lam1=lam1_train, lam2=lam2, gam=gam, eps=1e-3, maxit=1e+3)
    beta_l <- fit_fold$beta
    B_l<-fit_fold$B
    rank_l <- fit_fold$rank
    step_l <- fit_fold$step
    time_l <- fit_fold$time
    eval_fold <- eval_dc(beta_l, x_val, y_val)+eval_fold  # The evaluation: distance correlation.
    ind <- which(sapply(beta_l, is.null))
    rank_l[ind] <- -1
    eval_fold[ind] <- min(eval_fold, na.rm = TRUE)
  }
  eval_all[i]=eval_fold/5
  print(eval_fold/5)
}

cvm <- eval_all
gam_max=gam
cvm_max=which.max(cvm)
lam1_max=lam1[cvm_max]
lam2_max <- 0
M1=M*outer(weight,weight)
U1=weight*U
weight[weight==0]=Inf
fit1<-admm(M1, U1, p,p, lam1_max, lam2_max, gam_max, eps=1e-10, maxit=1e+8)$B[[1]]/weight
which(fit1[,1]!=0)

svd(fit1)$d
outcome<-svd(fit1)$u[,1:2]

data1=data.frame(Class=factor(trainy,levels=c(0,1),labels=c("Negative Samples","Positive Samples")),x1=trainx%*%outcome[,1],x2=trainx%*%outcome[,2])
par(mar=c(5,4.5,4.5,1)+0.1)
gridRange <- apply(data1[c('x1', 'x2')], 2, range)
x1 <- seq(from = range(data1$x1)[1]-0.5, to = range(data1$x1)[2]+0.5, length = 150)
x2 <- seq(from = range(data1$x2)[1]-0.5, to = range(data1$x2)[2]+0.5, length = 150)
grid <- expand.grid(x1 = x1, x2 = x2)
fit<-svm(Class~x1+x2,data1,type="nu-classification")
grid$class <- predict(fit,grid)
decisionValues <- predict(fit, grid, decision.values = TRUE)
grid$z <- as.vector(attributes(decisionValues)$decision)
plot(data1[c('x1', 'x2')], pch = ifelse(data1$Class == "Positive Samples", 1, 19),xlab='The first proposed reduced predictor',ylab='The second proposed reduced predictor',xlim=c(range(data1$x1)[1]-0.5,range(data1$x1)[2]+0.5),ylim=c(range(data1$x2)[1]-0.5,range(data1$x2)[2]+0.5),cex.axis=1.45,cex.lab=1.45)
contour(x1, x2, matrix(grid$z, length(x1), length(x2)), level=0, lwd = 1.5, drawlabels = FALSE, add=TRUE)
which(abs(outcome[,1])>1e-7)
### use the selected columns and glmnet to do the leave-one-out cross validation to reduce the computation time
index<-selected_index
countDR=0
countDRQDA=0

for(i in 1:216)
{
   kernel<-cbind(dr_kernel(trainx[-i,],trainy[-i],2),sir_kernel(trainx[-i,],trainy[-i],2))[,index]
   solution_equal<-Matrix_estimation(data=trainx[-i,],ytilde=kernel)
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

#### code for drawing figure for SC-SIR
resultSCSIR<-as.matrix(resultSCSIR)
outcome1<-as.matrix(SC.SIRestimate)
data2<-data.frame(Class=factor(trainy,levels=c(0,1),labels=c("Negative Samples","Positive Samples")),x=trainx%*%outcome1)

X1<-as.numeric(data2$x[which(data2$Class=='Negative Samples')])
X2<-as.numeric(data2$x[which(data2$Class=='Positive Samples')])
p1<-hist(X1,plot=FALSE)
p2<-hist(X2,plot=FALSE)
par(mar = c(5, 4.5, 4.5,1), xpd = TRUE)
plot(0,0,type="n",xlim=c(-0.7,0.8),ylim=c(0,2.5),ylab=" ",xlab=" ",cex.axis=1.5, cex.lab=1.25)
plot(p1,freq=FALSE,density=15,angle=135,add=TRUE,col="black")
plot(p2,freq=FALSE,density=15,angle=45,add=TRUE,col="black")
s1<-seq(min(X1),max(X1),length=50)
fun1<-dnorm(s1,mean=mean(X1),sd=sd(X1))
lines(s1,fun1,col='black',lwd=2,lty=1)
s2<-seq(min(X2),max(X2),length=50)
fun2<-dnorm(s2,mean=mean(X2),sd=sd(X2))
lines(s2,fun2,col='black',lwd=2,lty=1)

mtext("The reduced predictor from sparse SIR by [31]", side = 1, line = 2.75, cex = 1.5)
mtext("Density", side = 2, line = 2.75, cex = 1.5)
#legend(x="topright",legend = c("positive","negative"),
#       density = c(15,15),
#       angle = c(45,135),inset = c(-0.25, 0.05),xpd = TRUE,
#       horiz = FALSE,cex=1)
#legend(x="bottomright",legend = c("positive","negative"),
#       lwd=2,lty=c(2,1),inset = c(-0.25, 0.65),xpd = TRUE,horiz = FALSE,cex=0.8438)


### TCSAVE estimate
gamma<-as.matrix(tcsavegamma)
gamma=gamma/1000
data1=data.frame(Class=factor(trainy,levels=c(0,1),labels=c("Negative Samples","Positive Samples")),x1=trainx%*%gamma[,1],x2=trainx%*%gamma[,2])
x1<-seq(from = min(data1$x1)-0.38, to = max(data1$x1)+0.38, length = 200)
x2<-seq(from= min(data1$x2)-0.38, to = max(data1$x2)+0.38, length = 200)
par(mar = c(5, 4.5, 4.5,1), xpd = TRUE)
grid <- expand.grid(x1 = x1, x2 = x2)
fit<-svm(Class~x1+x2,data1,type="nu-classification")
grid$class <- predict(fit,grid)
decisionValues <- predict(fit, grid, decision.values = TRUE)
grid$z <- as.vector(attributes(decisionValues)$decision)
plot(data1[c('x1', 'x2')], pch = ifelse(data1$Class == "Positive Samples", 1, 19),xlab='',ylab='',xlim=c(min(data1$x1)-0.1,max(data1$x1)+0.1),ylim=c(min(data1$x2)-0.1,max(data1$x2)+0.1),cex.axis=1.5,cex.lab=1.25)
contour(x1, x2, matrix(grid$z, length(x1), length(x2)), level=0, lwd = 1.5, drawlabels = FALSE, add=TRUE)
mtext("The first reduced predictor from TC-SAVE", side = 1, line = 2.75, cex = 1.5)
mtext("The second reduced predictor from TC-SAVE", side = 2, line = 2.75, cex = 1.5)
plot.new()


