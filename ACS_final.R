library(Matrix)
library(MASS)
require(snowfall)
sfLibrary(glmnet)
library(e1071)
library(lattice)
library(parallel)
library(doParallel)
library(LassoSIR)
library(snowfall)
library(stringr)
library(lattice)
library(ggplot2)
library(glmnet)
library(doParallel)


current_file_path<-rstudioapi::getSourceEditorContext()$path
current_dir<-dirname(current_file_path)
setwd(current_dir)

source("ACS_prepare.r")
source("kernel_matrix_prepare.r")


orth<-function(data,D,threshold=1e-8)
{
  data<-as.matrix(data)
  colnum<-ncol(data)
  count=1
  if(colnum==1) # deal with only one column sceanrio
  {
    return(1)
  }else
  {
  column_norm<-apply(data,2,norm,type="2")
  index_chosen<-order(column_norm,decreasing = TRUE)[1]
  working_index<-index_chosen
  }
  count=count+1
  while(count<=colnum)
  {
      for(i in c(1:colnum)[-index_chosen])
      {
        data[,i]=data[,i]-1/(t(data[,working_index])%*%data[,working_index])*(t(data[,working_index]%*%data[,i]))*data[,working_index]
      }
      column_norm[-index_chosen]<-apply(as.matrix(data[,-index_chosen]),2,norm,type="2")
      working_index<-order(column_norm,decreasing = TRUE)[count]
      if(column_norm[working_index]<threshold)
      {
        break
      }
      index_chosen<-c(index_chosen,working_index)
    count=count+1
    if((count>D)||(length(index_chosen)==colnum))
    {
      break
    }
  }
  return(index_chosen)
}

Sparse_vector<-function(Y,design_matrix,nfolds=5)
{
  fit=cv.glmnet(x=design_matrix,y=Y,nfolds=nfolds,type.measure = "mse",lambda.min.ratio=1e-6)
  return(as.numeric(coef(fit))[-1])
}

Matrix_estimation<-function(data,ytilde,nfolds=5)
{
  fit<-cv.glmnet(x=data,y=ytilde,family="mgaussian",nfolds=nfolds,lambda.min.ratio=1e-5,nlambda=10000,type.measure = "mse")
  return(coef(fit))
}

choose_weight<-function(x,rho=0.25)                
{
  x<-as.numeric(x)
  return(sum(x^2)^(-rho))
}

ACS_SDR<-function(X, kernel_matrix, nfolds=5, D=15)
{
  p<-ncol(as.matrix(X))
  n<-nrow(as.matrix(Y))
  ##column estimation #######
  sfInit(parallel = TRUE,cpus=15)
  sfLibrary(glmnet)
  sparse_estimate<-sfApply(kernel_matrix,margin=2,fun=Sparse_vector,design_matrix=X,nfolds=nfolds)
  sfStop()
  #sparse_estimate<-apply(kernel_matrix,MARGIN=2,FUN=Sparse_vector,design_matrix=X,nfolds=nfolds)
  index=vector()
  ## equivalent to take repeated column estimations
  num=ncol(sparse_estimate)
  for(i in 1:dim(sparse_estimate)[2])
  {
    if(all(sparse_estimate[,i]==0))
    {
    index=append(index,i)
    }
  }
  sparse_estimate<-sparse_estimate[,-index]
  kernel_matrix<-kernel_matrix[,-index]
  print((c(1:num)[-index]))
  sfInit(parallel = TRUE,cpus=15)
  sfLibrary(glmnet)
  sparse_estimate<-apply(kernel_matrix,MARGIN=2,FUN=Sparse_vector,design_matrix=X,nfolds=nfolds)
  sfStop()
  index1=vector()
  for(i in 1:dim(sparse_estimate)[2])
  {
    if(all(sparse_estimate[,i]==0))
    {
      index1=append(index1,i)
    }
  }
  if((length(index1)!=0)&&(is.null(ncol(kernel_matrix))==FALSE))
  {
  kernel_matrix<-kernel_matrix[,-index1]
  sparse_estimate<-sparse_estimate[,-index1]
  }
  ###column_selection
  if(is.null(ncol(kernel_matrix))==FALSE)
  {
  index_chosen<-orth(sparse_estimate,D=D)
  kernel_matrix<-kernel_matrix[,index_chosen]
  print(((c(1:num)[-index])[-index1])[index_chosen])
  }
  ###final estimate
  if(ncol(as.matrix(kernel_matrix))==1)
  {
  result<-as.numeric(coef(cv.glmnet(x=X,y=kernel_matrix,nfolds=5)))[-1]
  weight<-apply(as.matrix(result),MARGIN=1,choose_weight)
  if(all(1/weight==0))
  {
    result<-rep(0,p)
  }else
  {
  result<-as.numeric(coef(cv.glmnet(x=t(t(X)/weight),y=kernel_matrix,nfolds=5)))[-1]/weight
  }
  }else
  {
  solution_equal<-Matrix_estimation(data=X,ytilde=kernel_matrix)
  temp=do.call(cbind,(lapply(solution_equal,deal_list)))
  if(allzero(temp))
  {
    result<-temp
  }else
  {
    weight<-apply(temp,MARGIN=1,choose_weight)
    solution_equal<-Matrix_estimation(data=t(t(X)/weight),ytilde=kernel_matrix)
    result<-do.call(cbind,(lapply(solution_equal,deal_list)))/weight
  }
  }
  return(list(estimate=sparse_estimate,kernel=kernel_matrix,result=result,X=X))
}

order_determination<-function(X,Y,H=2,r=20)
{
  p<-ncol(as.matrix(X))
  n<-nrow(as.matrix(X))
  s<-5
  if(p==1)
  {
    return(1)
  }
  hnzero<-rep(0,p)
  lambda<-eigen(mhat_dr(x=X,y=Y,H=H))$values
  lambda[p+1]=0
  for(i in 1:r)
  {
    augX<-cbind(X,matrix(rnorm(s*n,sd=0.5),n,s))
    V=eigen(mhat_dr(x=augX,y=Y,H=H))
    for(j in 2:p)
    {
      hnzero[j]=hnzero[j]+sum(V$vectors[(p+1):(p+s),j-1]^2)/r
    }
  }
  phi=rep(0,p)
  phi[1]=lambda[1]/(1+lambda[1])
  phi[2:p]=lambda[2:p]/sum(lambda)+hnzero[2:p] 
  order_stat=which.min(phi)-1
  return(order_stat)
}





