rm(list=ls())
current_file_path<-rstudioapi::getSourceEditorContext()$path
current_dir<-dirname(current_file_path)
setwd(current_dir)
sfLibrary(glmnet)

source("ACS_final.r")
source("ACS_prepare.r")
source("kernel_matrix_prepare.r")
source("SEAS.r")

p=1000
n=200
H=2
nfolds=5
D=15

rho=0.5

cov=diag(p)
cov_decay<-function(p,rho)
{
  cov1<-diag(p)
  for(i in 1:p)
  {
    for(j in 1:p)
    {
      cov1[i,j]=rho^(abs(i-j))
    }
  }
  return(cov1)
}

cov_dense_biway<-function(p,rho,a)
{
  cov<-diag(p)
  for(i in 1:p)
  {
    for(j in 1:p)
    {
      if((i<=a)||(j<=a)||(i>=(p-a+1))||(j>=(p-a+1)))
      {
        cov[i,j]=rho^abs(i-j)
      }else
        if(i==j)
        {
          cov[i,j]=1
        }
      else
      {
        cov[i,j]=rho
      }
    }
  }
  return(cov)
}

cov_dense_oneway<-function(p,rho,a)
{
  cov<-diag(p)
  for(i in 1:p)
    {
      for(j in 1:p)
        {
          if((i<=a)||(j<=a))
            {
              cov[i,j]=rho^abs(i-j)
            }else
          if(i==j)
           {
             cov[i,j]=1
           }
          else
          {
            cov[i,j]=rho
          }
        }
     }
  return(cov)
}

Model5<-function(X)
{
  p=ncol(X)
  n=nrow(X)
  beta<-rep(0,p)
  beta[1:5]<-rep(1,5)
  prob<-1/(1+exp((X%*%beta)^2-rep(t(beta)%*%cov(X)%*%beta*pchisq(0.5,1),n)))
  Y=rbinom(n,1,prob)
  return(Y)
}

Model6<-function(X)
{
  p=ncol(X)
  beta1<-rep(0,p)
  beta2<-rep(0,p)
  beta1[1:3]<-1
  beta2[(p-2):p]<-1
  Y=0.4*(X%*%beta1)^2+3*abs((X%*%beta2))^0.5+0.2*rnorm(n)
  return(Y)
}

Model7<-function(X)
{
  p=ncol(X)
  beta=rep(0,p)
  beta[1:2]=1
  Y=exp(X%*%beta)+0.2*rnorm(n)
  return(Y)
}
  
Model8<-function(X)
{
  p=ncol(X)
  beta1<-rep(0,p)
  beta2<-rep(0,p)
  beta1[1:2]<-1
  beta2[(p-1):p]<-1
  Y=(2*sign(X%*%beta2)-1)*(sign(X%*%beta1)*abs(X%*%beta1)^(1/3)+0.5)+0.2*rnorm(n)
  return(Y)
}
  
Exp_model<-function(X,idx)
{
  if(idx==5)
  {
    return(Model5(X))
  }
  if(idx==6)
  {
    return(Model6(X))
  }
  if(idx==7)
  {
    return(Model7(X))
  }
  if(idx==8)
  {
    return(Model8(X))
  }
}

CvDR=0
ICvDR=0
CvSAVE=0
ICvSAVE=0
loss2DR<-rep(0,200)
loss2SAVE<-rep(0,200)
colDR=0
colSAVE=0
for(k in 1:200)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=cov)
  Y=Exp_model(X,7)
  X=scale(X,scale=FALSE)
  kernel<-cbind(save_kernel(X,Y,2),sir_kernel(X,Y,5))
  result<-ACS_SDR(X,kernel,nfolds,D)
  if(is.null(dim(result_DR$result)))
  {
    suppDR=which(result_DR$result!=0)
  } else
  {
    suppDR=which(result_DR$result[,1]!=0)
  }
  if(is.null(dim(result_SAVE$result)))
  {
    suppSAVE=which(result_SAVE$result!=0)
  } else
  {
    suppSAVE=which(result_SAVE$result[,1]!=0)
  }
  d1=order_determination(X=X[,suppDR],Y=Y,H=5)
  d2=order_determination(X=X[,suppSAVE],Y=Y,H=5)
  if(is.null(dim(result_DR$result)))
  {
    outcomeDR<-result_DR$result
  }else
  {
    outcomeDR<-svd(result_DR$result)$u[,1:d]
  }
  if(is.null(dim(result_SAVE$result)))
  {
    outcomeSAVE<-result_SAVE$result
  }else
  {
    outcomeSAVE<-svd(result_SAVE$result)$u[,1:d]
  }
  #outcome<-svd(result$result)$u[,1]
  CvDR=CvDR+length(which(suppDR<4))+length(which(suppDR>197))
  ICvDR=ICvDR+(length(which(suppDR>3))+length(which(suppDR<198))-(length(which(suppDR<4))-length(which(suppDR>197))))/2
  CvSAVE=CvSAVE+length(which(suppSAVE<4))+length(which(suppSAVE>197))
  ICvSAVE=ICvSAVE+(length(which(suppSAVE>3))+length(which(suppSAVE<198))-(length(which(suppSAVE<4))-length(which(suppSAVE>197))))/2
  loss2DR[k]<-norm(projection(cbind(beta1,beta2))-projection(outcomeDR),type="2")
  loss2SAVE[k]<-norm(projection(cbind(beta1,beta2))-projection(outcomeSAVE),type="2")
  colDR=colDR+dim(result_DR$result)[2]
  colSAVE=colSAVE+dim(result_SAVE$result)[2]
}

p=1000
CvSEAS=0
ICvSEAS=0
col=0
loss2=rep(0,200)
cov=cov_decay(p,rho=0.5)
#cov=cov_dense_biway(p,rho=0.5,a=2)
for(k in 1:200)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=cov)
  Y=Exp_model(X,8)
  X=scale(X,scale=FALSE)
  kernel<-cbind(save_kernel(X,Y,2),sir_kernel(X,Y,5))
  kernel<-save_kernel(X,Y,2)
  result<-ACS_SDR(X,kernel,nfolds,D)
  if(is.null(dim(result$result)))
  {
    supp=which(result$result!=0)
  } else
  {
    supp=which(result$result[,1]!=0)
  }
  d=order_determination(X=X[,supp],Y=Y,H=5)
  print(d)
  #d=2
  if(is.null(dim(result$result)))
  {
    outcome<-result$result
    col=col+1
  }else
  {
    outcome<-svd(result$result)$u[,1:d]
    col=col+ncol(result$result)
  }
  #outcome<-svd(result$result)$u[,1]
  #CvSEAS=CvSEAS+2-length(which(supp<=2))
  #ICvSEAS=ICvSEAS+length(which(supp>2))
  CvSEAS=CvSEAS+4-length(which(supp<3))-length(which(supp>998))
  ICvSEAS=ICvSEAS+(length(which(supp>2))+length(which(supp<999))-(length(which(supp<3))-length(which(supp>998))))/2 
  beta1<-rep(0,p)
  beta1[1:2]<-1
  beta2<-rep(0,p)
  beta2[999:1000]<-1
  beta<-rep(0,p)
  beta[1:2]<-1
  #loss2[k]<-norm(projection(beta)-projection(outcome),type="2")
  loss2[k]<-norm(projection(cbind(beta1,beta2))-projection(outcome),type="2")
  print(loss2[k])
}

print(mean(loss2))
print(sd(loss2)/sqrt(200))
print(CvSEAS/200)
print(ICvSEAS/200)
col/200

supp=c(1,2,999,1000,498,499,500,501,502)
for(k in 1:200)
{
X<-mvrnorm(n=n,rep(0,p),Sigma=cov)
Y=Exp_model(X,8)
d=order_determination(X=X[,supp],Y=Y,H=2)
print(d)
}

p=200
n=200
CvSEAS=0
ICvSEAS=0
#cov=cov_decay(p,rho=0.5)
cov=cov_dense_oneway(p,rho=0.5,a=2)
loss=rep(0,200)

for(k in 1:200)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=cov)
  Y=Exp_model(X,7)
  X=scale(X,scale=FALSE)
  result=cv.seas(x=X,y=Y,d=NULL)
  supp=which(result$beta[,1]!=0)
  CvSEAS=CvSEAS+2-length(which(supp<=2))
  ICvSEAS=ICvSEAS+length(which(supp>2))
  #CvSEAS=CvSEAS+4-length(which(supp<3))-length(which(supp>198))
  #ICvSEAS=ICvSEAS+(length(which(supp>2))+length(which(supp<199))-(length(which(supp<3))-length(which(supp>198))))/2                                     
  beta1<-rep(0,p)
  beta2<-rep(0,p)
  beta1[1:2]<-1
  beta2[(p-1):p]<-1
  #loss[k]=norm(projection(cbind(beta1, beta2))-projection(result$beta),type="2")
  beta<-rep(0,p)
  beta[1:2]<-1
  loss[k]=norm(projection(beta)-projection(result$beta),type="2")
}

print(mean(loss))
print(sd(loss)/sqrt(200))
print(CvSEAS/200)
print(ICvSEAS/200)


print(mean(loss2DR))
print(mean(loss2SAVE))
print(colDR/200)
print(colSAVE/200)
print(CvDR/200)
print(CvSAVE/200)
print(ICvDR/200)
print(ICvSAVE/200)

write.csv(loss2SAVE,"D:/Simulation_result/ACS-SAVE/High-dimensional/Model2Sigma2p1000n200save.csv")
write.csv(loss2DR,"D:/Simulation_result/ACS-SAVE/High-dimensional/Model2Sigma2p1000n200Hybrid.csv")


