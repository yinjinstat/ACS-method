library(MASS)
result_table.8<-matrix(0,16,4)
sdcol<-rep(0,16)
sdTSR<-rep(0,16)
col=0


save<-function(x,y,H)
{
  n<-nrow(x)
  p<-ncol(x)
  sig<-cov(x)
  xc<- t(t(x)-apply(x,2,mean))
  isg<-matpower(sig,-1)
  dy <- slicing(y,H)
  H <- dy$H
  ytilde <- dy$ytilde
  prop <- dy$prop
  xxbar<-numeric(0)
  for (j in 1:(H-1))
  {
    ind<-((y >= ytilde[j])&(y < ytilde[j+1]))
    xxbar<-cbind(xxbar,isg%*%(sig-cov(xc[ind,])))
  } 
  ind<-(y >= ytilde[H])
  xxbar<-cbind(xxbar,isg%*%(sig-cov(xc[ind,])))
  return(eigen((xxbar)%*%t(xxbar)))
}

save_partial<-function(x,y,H)
{
  n<-nrow(x)
  p<-ncol(x)
  sig<-cov(x)
  xc<- t(t(x)-apply(x,2,mean))
  isg<-matpower(sig,-1)
  dy <- slicing(y,H)
  H <- dy$H
  ytilde <- dy$ytilde
  prop <- dy$prop
  xxbar<-numeric(0)
  for (j in 1:(H-1))
  {
    ind<-((y >= ytilde[j])&(y < ytilde[j+1]))
    xxbar<-cbind(xxbar,(isg%*%(sig-cov(xc[ind,])))[,1])
  } 
  ind<-(y >= ytilde[H])
  xxbar<-cbind(xxbar,(isg%*%(sig-cov(xc[ind,])))[,1])
  return(eigen((xxbar)%*%t(xxbar)))
}


ctsave_oracle<-function(x,y,H,d,b0)
{
  n<-nrow(x)
  p<-ncol(x)
  sig<-cov(x)
  xc<- t(t(x)-apply(x,2,mean))
  isg<-matpower(sig,-1)
  dy <- slicing(y,H)
  H <- dy$H
  ytilde <- dy$ytilde
  prop <- dy$prop
  xxbar<-numeric(0)
  for (j in 1:(H-1))
  {
    ind<-((y >= ytilde[j])&(y < ytilde[j+1]))
    xxbar<-cbind(xxbar,isg%*%(sig-cov(xc[ind,]))/sum(ind))
  } 
  ind<-(y >= ytilde[H])
  xxbar<-cbind(xxbar,isg%*%(sig-cov(xc[ind,]))/sum(ind))
  ssall<-numeric(0)
  for (k in 1:p)
  {ssall<-cbind(ssall,rep(c(rep(FALSE,each=2^(k-1)),rep(TRUE,each=2^(k-1))),
                          2^(p-k)))}
  ssall<-ssall[-1,]
  ssalls<-numeric(0)
  for (h in 1:H) {ssalls<-cbind(ssalls,ssall)}
  
  K<-nrow(ssalls) 
  s<-rep(0,K)
  for (i in 1:K) 
  {
    u<-svd(xxbar[,ssalls[i,]==1])$u[,1:d]
    s[i]<-dist(u,b0)
  }
  res<-list()
  res$s<-s
  res$subset<-ssalls[which.min(s),]
  res$b1<-svd(xxbar[,ssalls[1,]==1])$u[,1:d]
  res$bf<-svd(xxbar[,ssalls[2^p-1,]==1])$u[,1:d]
  return(res)
}

############################
### ctsave for a single sample, useful in the official function below
############################
ctsave_single<-function(x,y,H,d)
{
  n<-nrow(x)
  p<-ncol(x)
  sig<-cov(x)
  xc<- t(t(x)-apply(x,2,mean))
  isg<-matpower(sig,-1)
  dy <- slicing(y,H)
  H <- dy$H
  ytilde <- dy$ytilde
  prop <- dy$prop
  xxbar<-numeric(0)
  for (j in 1:(H-1))
  {
    ind<-((y >= ytilde[j])&(y < ytilde[j+1]))
    xxbar<-cbind(xxbar,isg%*%(sig-cov(xc[ind,])))
  } 
  ind<-(y >= ytilde[H])
    xxbar<-cbind(xxbar,isg%*%(sig-cov(xc[ind,])))
  
  ssall<-numeric(0)  ###family of all the index sets of columns
  for (k in 1:p)
  {ssall<-cbind(ssall,rep(c(rep(FALSE,each=2^(k-1)),rep(TRUE,each=2^(k-1))),
                          2^(p-k)))}
  ssall<-ssall[-1,]
  ssalls<-numeric(0)
  for (h in 1:H) {ssalls<-cbind(ssalls,ssall)}
  
  K<-nrow(ssalls) 
  s<-array(0,dim=c(p,d,K))
  r=3
  order_stat=rep(0,K)
  for (i in 1:K) 
  {
    u<-svd(xxbar[,ssalls[i,]==1])$u[,1:d]
    s[,,i]<-u
  }
  for(i in 1:K)
  {
    phi_final=rep(0,length(which(ssalls[i,]==1)))
    for(q in 1:10)
    {
    augX<-cbind(x,matrix(rnorm(r*n,sd=0.5),n,r))
    xxbar_aug<-numeric(0)
    sig_aug<-cov(augX)
    isg_aug<-matpower(sig_aug,-1)
    xc_aug<- t(t(augX)-apply(augX,2,mean))
    for (j in 1:(H-1))
    {
        ind<-((y >= ytilde[j])&(y < ytilde[j+1]))
        xxbar_aug<-cbind(xxbar_aug,sqrt(prop[j])*isg_aug%*%(sig_aug-cov(xc_aug[ind,])))
    } 
    ind<-(y >= ytilde[H])
    xxbar_aug<-cbind(xxbar_aug,sqrt(prop[H])*isg_aug%*%(sig_aug-cov(xc_aug[ind,])))
    eigens=svd(xxbar[,ssalls[i,]==1])$d
    V=svd(xxbar_aug)
    t=length(eigens)
    hnzero<-rep(0,t)
    phi=rep(0,t)
    phi[1]=eigens[1]/(1+eigens[1])
    hnzero[1]=0
    for(j in 2:t)
    {
        hnzero[j]=hnzero[j]+sum(V$u[(t+1):(t+r),j-1]^2)
    }
    phi[2:t]=eigens[2:t]/sum(eigens)+hnzero[1:(t-1)]   
    phi_final=phi_final+phi
    }
    order_stat[i]=which.min(phi_final)-1
  }
  order_max=order_stat[length(order_stat)]
  return(list(s=s,keep_index=which(order_stat==order_max)))
}

ctsave_single_change<-function(x,y,H,d,keep_index)
{
  n<-nrow(x)
  p<-ncol(x)
  sig<-cov(x)
  xc<- t(t(x)-apply(x,2,mean))
  isg<-matpower(sig,-1)
  dy <- slicing(y,H)
  H <- dy$H
  ytilde <- dy$ytilde
  prop <- dy$prop
  xxbar<-numeric(0)
  for (j in 1:(H-1))
  {
    ind<-((y >= ytilde[j])&(y < ytilde[j+1]))
    xxbar<-cbind(xxbar,isg%*%(sig-cov(xc[ind,])))
  } 
  ind<-(y >= ytilde[H])
  xxbar<-cbind(xxbar,isg%*%(sig-cov(xc[ind,])))
  
  ssall<-numeric(0)  ###family of all the index sets of columns
  for (k in 1:p)
  {ssall<-cbind(ssall,rep(c(rep(FALSE,each=2^(k-1)),rep(TRUE,each=2^(k-1))),
                          2^(p-k)))}
  ssall<-ssall[-1,]
  ssalls<-numeric(0)
  for (h in 1:H) {ssalls<-cbind(ssalls,ssall)}
  
  K<-nrow(ssalls) 
  s<-array(0,dim=c(p,d,K))

  for (i in keep_index) 
  {
    u<-svd(xxbar[,ssalls[i,]==1])$u[,1:d]
    s[,,i]<-u
  }
  return(s)
}

#########################################
###
### the working ctsave: run all the subsets exhaustively and pick the optimal
### 
### we measure the goodness of fit by 
###   (bootstrap variance of $S(\hat \beta)$) 
### + n^(-1/4)*distance between $S(\hat \beta)$ & $\hat {full column space}$
###  i.e. estimation variance + n^(-1/4)*estimation bias
###
### note: the first term is for asymptotic relative efficiency,
###     the second term is to make sure $S(\hat \beta)$ spans the full space
#########################################
### x:n*p, y:n-dim vector, H=#of slices, d=dim of central subspace, 
### nboot = #of bootstap samples
#########################################

ctsave_optimal<-function(x,y,H,d,beta0)
{
  n<-nrow(x)
  p<-ncol(x)
  ssall<-numeric(0)
  for (k in 1:p)
  {ssall<-cbind(ssall,rep(c(rep(FALSE,each=2^(k-1)),rep(TRUE,each=2^(k-1))),
                          2^(p-k)))}
  ssall<-ssall[-1,]
  K<-nrow(ssall)
  ssalls<-numeric(0)
  for (h in 1:H) {ssalls<-cbind(ssalls,ssall)}
  
  sf<-ctsave_single_change(x,y,H,d,1:nrow(ssall)) ###set of $\hat \beta$ for the full sample
  res<-rep(0,K)
  for(i in 1:K)
  {
    res[i]<-sqrt(sum((proj(sf[,,i])-proj(beta0))^2))
  }
  return(res)
}

ctsave<-function(x,y,H,d,nboot,Fbest)
{
  n<-nrow(x)
  p<-ncol(x)
  if (missing(nboot)) {nboot=n}
  ssall<-numeric(0)
  for (k in 1:p)
  {ssall<-cbind(ssall,rep(c(rep(FALSE,each=2^(k-1)),rep(TRUE,each=2^(k-1))),
                          2^(p-k)))}
  ssall<-ssall[-1,]
  ssalls<-numeric(0)
  for (h in 1:H) {ssalls<-cbind(ssalls,ssall)}
  
  result<-ctsave_single(x,y,H,d)
  sf<-result$s ###set of $\hat \beta$ for the full sample
  keep_index<-result$keep_index
  print(result$keep_index)
  
  obs<-cbind(x,y)
  K<-nrow(ssall) 
  s<-array(0,dim=c(p,d,K,nboot)) ###set of $\hat \beta$ for the bs samples
  for (i in 1:nboot)
  {
    sami<-sample(1:n,size=n,replace=TRUE)
    xi<-x[sami,]
    yi<-y[sami]
    s[,,,i]<-ctsave_single_change(xi,yi,H,d,keep_index)}
  
  rsb<-rep(0,K)  
  rsv<-rep(0,K)
  rsv[-keep_index]<-100
  len<-rep(0,K)
  for(i in 1:K)
  {
    len[i]=length((1:p)[ssall[i,]==1])
  }
  for (i in keep_index) 
  {
    ui<-sf[,,i]
    rsb[i]<-norm(proj(ui)-proj(sf[,,K]),type="F")*n^(-0.3)
    for (j in 1:nboot)
    {
      uij<-s[,,i,j]
      #rsv[i]<-rsv[i]+norm(proj(uij)-proj(ui),type="F")/n
      rsv[i]<-rsv[i]+norm(proj(uij)-proj(ui),type="F")^2/sqrt(n)
    }
  }
  rs<-rsv
  rs=rs+n^(-0.7)*len
  ropt<-which.min(rs) ###optimal subset
  I<-(1:p)[ssall[ropt,]==1]
  res<-list()
  res$vb<-rsv[ropt] ###variance of optimal $\hat \beta$
  res$beta<-sf[,,ropt] ###optimal $\hat \beta$
  res$I<-I  ###corresponding optimal subset
  res$bfull<-sf[,,K]  ###$\hat \beta$ for original SAVE
  res$var<-rsv ###variances of all the choices
  res$bias<-rsb ###bias of all the choices
  res$optimal<-sf[,,Fbest]
  return(res)
}

p=6
n=200
lossboot<-rep(0,200)
loss<-rep(0,200)
lossbest<-rep(0,200)
rho<-0.2
covariance<-diag(p)

for(i in 1:p)
{
  for(j in 1:p)
  {
    if(i!=j)
    {
      covariance[i,j]=rho
    }
  }
}

H=5
beta1<-rep(0,p)
beta2<-rep(0,p)
col_num=0
ICv=0
Cv=0
Cvper=rep(0,500)
col=rep(0,500)
for(k in 1:500)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=covariance)
  #beta<-rep(0,p)
  #beta[1:3]<-rep(1,3)
  #beta1<-c(1,0,0,0,0,0)
  beta1<-c(1,0,0,0,0,0)
  beta2<-c(0,1,0,0,0,0)
  #beta2<-c(0,0,0,0,0,1,1,1,1,1)
  #beta2[2]<-1
  error<-rnorm(n)
  Y=exp((X%*%beta1)^2)+0.2*error
  #Y=log((X%*%beta1)^2+5)+0.2*error
  #Y=0.4*(X%*%beta1)^2+3*abs((X%*%beta2))^0.5+0.2*error
  #Y=3*sin(X%*%beta1/4)+0.4*(X%*%beta2)^2+0.2*error
  #Y=exp((X%*%beta1)^2)+2*exp((X%*%beta2)^2)+0.2*error
  #supp=which(result$result!=0)
  resultboot<-ctsave(X,Y,H,1,200,1)
  col[k]=length(resultboot$I)
  col_num=col_num+length(resultboot$I)
  result<-save(X,Y,H)
  lossboot[k]<-norm(proj(resultboot$beta)-proj(beta1),"2")
  loss[k]<-norm(proj(result$vectors[,1])-proj(beta1),"2")
  lossbest[k]<-norm(proj(resultboot$optimal)-proj(beta1),"2")
  print(resultboot$I)
  print(lossboot[k]/loss[k])
  #lossboot[k]<-norm(projection(resultboot$beta)-projection(cbind(beta1,beta2)),"2")
  #loss[k]<-norm(projection(result$vectors[,1:2])-projection(cbind(beta1,beta2)),"2")
  Cv=Cv+length(which(resultboot$I<=1))
  Cvper[k]=length(which(resultboot$I<=1))
  ICv=ICv+length(which(resultboot$I>1))
}

result_table.8[1,1]<-mean(lossboot)/mean(loss)
result_table.8[1,2]<-col_num/100
result_table.8[1,3]<-Cv/500

bootvar<-rep(0,500)
colboot<-rep(0,500)
for(i in 1:500)
{
  index=sample(1:500,500,replace=TRUE)
  bootboot=lossboot[index]
  bootloss=loss[index]
  bootvar[i]=mean(bootboot)/mean(bootloss)
}
result_table.8[1,4]<-sd(bootvar)
write.csv(lossboot,"C:/Simulation_result/ACS-SAVE/Low-dimensional/Model1rho0.2p6ctsave.csv")
write.csv(loss,"C:/Simulation_result/ACS-SAVE/Low-dimensional/Model1rho0.2p6save.csv")
write.csv(lossbest,"C:/Simulation_result/ACS-SAVE/Low-dimensional/Model1rho0.2p6oracle.csv")


p=6
n=200
lossboot<-rep(0,200)
loss<-rep(0,200)
lossbest<-rep(0,200)
rho<-0.8
covariance<-diag(p)
for(i in 1:p)
{
  for(j in 1:p)
  {
    if(i!=j)
    {
      covariance[i,j]=rho
    }
  }
}

H=5
beta1<-rep(0,p)
beta2<-rep(0,p)
col_num=0
ICv=0
Cv=0
Cvper=rep(0,500)
col=rep(0,500)
for(k in 1:500)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=covariance)
  #beta<-rep(0,p)
  #beta[1:3]<-rep(1,3)
  #beta1<-c(1,0,0,0,0,0)
  beta1<-c(1,0,0,0,0,0)
  beta2<-c(0,1,0,0,0,0)
  #beta2<-c(0,0,0,0,0,1,1,1,1,1)
  #beta2[2]<-1
  error<-rnorm(n)
  Y=exp((X%*%beta1)^2)+0.2*error
  #Y=log((X%*%beta1)^2+5)+0.2*error
  #Y=0.4*(X%*%beta1)^2+3*abs((X%*%beta2))^0.5+0.2*error
  #Y=3*sin(X%*%beta1/4)+0.4*(X%*%beta2)^2+0.2*error
  #Y=exp((X%*%beta1)^2)+2*exp((X%*%beta2)^2)+0.2*error
  #supp=which(result$result!=0)
  resultboot<-ctsave(X,Y,H,1,200,1)
  col_num=col_num+length(resultboot$I)
  col[k]=length(resultboot$I)
  result<-save(X,Y,H)
  lossboot[k]<-norm(proj(resultboot$beta)-proj(beta1),"2")
  loss[k]<-norm(proj(result$vectors[,1])-proj(beta1),"2")
  lossbest[k]<-norm(proj(resultboot$optimal)-proj(beta1),"2")
  #lossboot[k]<-norm(projection(resultboot$beta)-projection(cbind(beta1,beta2)),"2")
  #loss[k]<-norm(projection(result$vectors[,1:2])-projection(cbind(beta1,beta2)),"2")
  Cv=Cv+length(which(resultboot$I<=1))
  Cvper=length(which(resultboot$I<=1))
  print(resultboot$I)
  ICv=ICv+length(which(resultboot$I>1))
}

result_table.8[2,1]<-mean(lossboot)/mean(loss)
result_table.8[2,2]<-col_num/100
result_table.8[2,3]<-Cv/500
bootvar<-rep(0,500)
for(i in 1:500)
{
  index=sample(1:500,500,replace=TRUE)
  bootboot=lossboot[index]
  bootloss=loss[index]
  bootvar[i]=mean(bootboot)/mean(bootloss)
}
result_table.8[2,4]<-sd(bootvar)

res=rep(0,63)
for(k in 1:500)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=covariance)
  beta1<-c(1,0,0,0,0,0)
  error<-rnorm(n)
  Y=exp((X%*%beta1)^2)+0.2*error
  res=res+ctsave_optimal(X,Y,H,1,beta1)
  print(k)
}


write.csv(bootboot,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model1rho0.8p6ctsave.csv")
write.csv(bootloss,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model1rho0.8p6save.csv")
sdcol[2]=sd(col)/sqrt(500)
sdTSR[2]=sd(Cvper)/sqrt(500)

p=10
n=200
lossboot<-rep(0,200)
loss<-rep(0,200)
lossbest<-rep(0,200)
rho<-0.2
covariance<-diag(p)
for(i in 1:p)
{
  for(j in 1:p)
  {
    if(i!=j)
    {
      covariance[i,j]=rho
    }
  }
}

H=5
beta1<-rep(0,p)
beta2<-rep(0,p)
col_num=0
ICv=0
Cv=0
Cvper=rep(0,500)
col=rep(0,500)
for(k in 1:500)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=covariance)
  #beta<-rep(0,p)
  #beta[1:3]<-rep(1,3)
  #beta1<-c(1,0,0,0,0,0)
  beta1<-c(1,0,0,0,0,0,0,0,0,0)
  beta2<-c(0,1,0,0,0,0)
  #beta2<-c(0,0,0,0,0,1,1,1,1,1)
  #beta2[2]<-1
  error<-rnorm(n)
  Y=exp((X%*%beta1)^2)+0.2*error
  #Y=log((X%*%beta1)^2+5)+0.2*error
  #Y=0.4*(X%*%beta1)^2+3*abs((X%*%beta2))^0.5+0.2*error
  #Y=3*sin(X%*%beta1/4)+0.4*(X%*%beta2)^2+0.2*error
  #Y=exp((X%*%beta1)^2)+2*exp((X%*%beta2)^2)+0.2*error
  #supp=which(result$result!=0)
  resultboot<-ctsave(X,Y,H,1,200,1)
  col_num=col_num+length(resultboot$I)
  result<-save(X,Y,H)
  lossboot[k]<-norm(proj(resultboot$beta)-proj(beta1),"2")
  loss[k]<-norm(proj(result$vectors[,1])-proj(beta1),"2")
  lossbest[k]<-norm(proj(resultboot$optimal)-proj(beta1),"2")
  #lossboot[k]<-norm(projection(resultboot$beta)-projection(cbind(beta1,beta2)),"2")
  #loss[k]<-norm(projection(result$vectors[,1:2])-projection(cbind(beta1,beta2)),"2")
  Cv=Cv+length(which(resultboot$I<=1))
  Cvper=length(which(resultboot$I<=1))
  col[k]=length(resultboot$I)
  print(resultboot$I)
  ICv=ICv+length(which(resultboot$I>1))
}

result_table.8[3,1]<-mean(lossboot)/mean(loss)
result_table.8[3,2]<-col_num/100
result_table.8[3,3]<-Cv/500
bootvar<-rep(0,500)
for(i in 1:500)
{
  index=sample(1:500,500,replace=TRUE)
  bootboot=lossboot[index]
  bootloss=loss[index]
  bootvar[i]=mean(bootboot)/mean(bootloss)
}
result_table.8[3,4]<-sd(bootvar)
write.csv(bootboot,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model1rho0.2p10ctsave.csv")
write.csv(bootloss,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model1rho0.2p10save.csv")
sdcol[3]=sd(col)/sqrt(500)
sdTSR[3]=sd(Cvper)/sqrt(500)

p=10
n=200
lossboot<-rep(0,200)
loss<-rep(0,200)
lossbest<-rep(0,200)
rho<-0.8
covariance<-diag(p)
for(i in 1:p)
{
  for(j in 1:p)
  {
    if(i!=j)
    {
      covariance[i,j]=rho
    }
  }
}

H=5
beta1<-rep(0,p)
beta2<-rep(0,p)
col_num=0
ICv=0
Cv=0
for(k in 1:500)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=covariance)
  #beta<-rep(0,p)
  #beta[1:3]<-rep(1,3)
  #beta1<-c(1,0,0,0,0,0)
  beta1<-c(1,0,0,0,0,0,0,0,0,0)
  beta2<-c(0,1,0,0,0,0)
  #beta2<-c(0,0,0,0,0,1,1,1,1,1)
  #beta2[2]<-1
  error<-rnorm(n)
  Y=exp((X%*%beta1)^2)+0.2*error
  #Y=log((X%*%beta1)^2+5)+0.2*error
  #Y=0.4*(X%*%beta1)^2+3*abs((X%*%beta2))^0.5+0.2*error
  #Y=3*sin(X%*%beta1/4)+0.4*(X%*%beta2)^2+0.2*error
  #Y=exp((X%*%beta1)^2)+2*exp((X%*%beta2)^2)+0.2*error
  #supp=which(result$result!=0)
  resultboot<-ctsave(X,Y,H,1,200,1)
  col_num=col_num+length(resultboot$I)
  result<-save(X,Y,H)
  lossboot[k]<-norm(proj(resultboot$beta)-proj(beta1),"2")
  loss[k]<-norm(proj(result$vectors[,1])-proj(beta1),"2")
  lossbest[k]<-norm(proj(resultboot$optimal)-proj(beta1),"2")
  #lossboot[k]<-norm(projection(resultboot$beta)-projection(cbind(beta1,beta2)),"2")
  #loss[k]<-norm(projection(result$vectors[,1:2])-projection(cbind(beta1,beta2)),"2")
  Cv=Cv+length(which(resultboot$I<=1))
  Cvper=length(which(resultboot$I<=1))
  col[k]=length(resultboot$I)
  print(resultboot$I)
  ICv=ICv+length(which(resultboot$I>1))
}

result_table.8[4,1]<-mean(lossboot)/mean(loss)
result_table.8[4,2]<-col_num/100
result_table.8[4,3]<-Cv/500
bootvar<-rep(0,500)
for(i in 1:500)
{
  index=sample(1:500,500,replace=TRUE)
  bootboot=lossboot[index]
  bootloss=loss[index]
  bootvar[i]=mean(bootboot)/mean(bootloss)
}
result_table.8[4,4]<-sd(bootvar)
write.csv(bootboot,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model1rho0.8p10ctsave.csv")
write.csv(bootloss,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model1rho0.8p10save.csv")
sdcol[4]=sd(col)/sqrt(500)
sdTSR[4]=sd(Cvper)/sqrt(500)


p=6
n=200
lossboot<-rep(0,200)
loss<-rep(0,200)
lossbest<-rep(0,200)
rho<-0.2
covariance<-diag(p)
for(i in 1:p)
{
  for(j in 1:p)
  {
    if(i!=j)
    {
      covariance[i,j]=rho
    }
  }
}

H=5
beta1<-rep(0,p)
beta2<-rep(0,p)
col_num=0
ICv=0
Cv=0
for(k in 1:500)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=covariance)
  beta1<-c(1,1,1,1,1,1)
  beta2<-c(0,1,0,0,0,0)
  error<-rnorm(n)
  Y=log((X%*%beta1)^2+5)+0.2*error
  resultboot<-ctsave(X,Y,H,1,200,63)
  #resultboot<-ctsave_change(X,Y,H,1,500,63)
  col_num=col_num+length(resultboot$I)
  result<-save(X,Y,H)
  lossboot[k]<-norm(proj(resultboot$beta)-proj(beta1),"2")
  loss[k]<-norm(proj(result$vectors[,1])-proj(beta1),"2")
  lossbest[k]<-norm(proj(resultboot$optimal)-proj(beta1),"2")
  Cv=Cv+length(which(resultboot$I<=6))
  Cvper=length(which(resultboot$I<=6)/6)
  print(lossboot[k]/loss[k])
  col[k]=length(resultboot$I)
  print(resultboot$I)
  ICv=ICv+length(which(resultboot$I>6))
}

result_table.8[5,1]<-mean(lossboot)/mean(loss)
result_table.8[5,2]<-col_num/100
result_table.8[5,3]<-Cv/500/6
bootvar<-rep(0,500)
for(i in 1:500)
{
  index=sample(1:500,500,replace=TRUE)
  bootboot=lossboot[index]
  bootloss=loss[index]
  bootvar[i]=mean(bootboot)/mean(bootloss)
}
result_table.8[5,4]<-sd(bootvar)
write.csv(bootboot,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model2rho0.2p6n500ctsave.csv")
write.csv(bootloss,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model2rho0.2p6n500save.csv")
sdcol[5]=sd(col)/sqrt(500)
sdTSR[5]=sd(Cvper)/sqrt(500)

p=6
n=200
lossboot<-rep(0,200)
loss<-rep(0,200)
lossbest<-rep(0,200)
rho<-0.8
covariance<-diag(p)
for(i in 1:p)
{
  for(j in 1:p)
  {
    if(i!=j)
    {
      covariance[i,j]=rho
    }
  }
}

col=0
H=5
beta1<-rep(0,p)
beta2<-rep(0,p)
col_num=0
ICv=0
Cv=0
for(k in 1:500)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=covariance)
  #beta<-rep(0,p)
  #beta[1:3]<-rep(1,3)
  #beta1<-c(1,0,0,0,0,0)
  beta1<-c(1,1,1,1,1,1)
  beta2<-c(0,1,0,0,0,0)
  #beta2<-c(0,0,0,0,0,1,1,1,1,1)
  #beta2[2]<-1
  error<-rnorm(n)
  Y=log((X%*%beta1)^2+5)+0.2*error
  #Y=exp((X%*%beta1)^2)+0.2*error
  #Y=sin((X%*%beta1)^2)+0.2*error
  #Y=0.4*(X%*%beta1)^2+3*abs((X%*%beta2))^0.5+0.2*error
  #Y=3*sin(X%*%beta1/4)+0.4*(X%*%beta2)^2+0.2*error
  #Y=exp((X%*%beta1)^2)+2*exp((X%*%beta2)^2)+0.2*error
  #supp=which(result$result!=0)
  resultboot<-ctsave(X,Y,H,1,200,63)
  col_num=col_num+length(resultboot$I)
  result<-save(X,Y,H)
  lossboot[k]<-norm(proj(resultboot$beta)-proj(beta1),"2")
  loss[k]<-norm(proj(result$vectors[,1])-proj(beta1),"2")
  lossbest[k]<-norm(proj(resultboot$optimal)-proj(beta1),"2")
  print(lossboot[k]/loss[k])
  #lossboot[k]<-norm(projection(resultboot$beta)-projection(cbind(beta1,beta2)),"2")
  #loss[k]<-norm(projection(result$vectors[,1:2])-projection(cbind(beta1,beta2)),"2")
  Cv=Cv+length(which(resultboot$I<=6))
  Cvper=length(which(resultboot$I<=6)/6)
  col[k]=length(resultboot$I)
  print(resultboot$I)
  ICv=ICv+length(which(resultboot$I>6))
}

result_table.8[6,1]<-mean(lossboot)/mean(loss)
result_table.8[6,2]<-col_num/100
result_table.8[6,3]<-Cv/500/6
bootvar<-rep(0,500)
for(i in 1:500)
{
  index=sample(1:500,500,replace=TRUE)
  bootboot=lossboot[index]
  bootloss=loss[index]
  bootvar[i]=mean(bootboot)/mean(bootloss)
}
result_table.8[6,4]<-sd(bootvar)


res=rep(0,63)
for(k in 1:200)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=covariance)
  beta1<-c(1,1,1,1,1,1)
  beta2<-c(0,1,0,0,0,0)
  error<-rnorm(n)
  Y=log((X%*%beta1)^2+5)+0.2*error
  res=res+ctsave_optimal(X,Y,H,1,beta1)
  print(k)
}

write.csv(bootboot,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model2rho0.8p6n500ctsave.csv")
write.csv(bootloss,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model2rho0.8p6n500save.csv")
sdcol[6]=sd(col)/sqrt(500)
sdTSR[6]=sd(Cvper)/sqrt(500)

p=10
n=200
lossboot<-rep(0,200)
loss<-rep(0,200)
lossbest<-rep(0,200)
rho<-0.2
covariance<-diag(p)
for(i in 1:p)
{
  for(j in 1:p)
  {
    if(i!=j)
    {
      covariance[i,j]=rho
    }
  }
}

H=5
beta1<-rep(0,p)
beta2<-rep(0,p)
col_num=0
ICv=0
Cv=0
for(k in 1:500)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=covariance)
  #beta<-rep(0,p)
  #beta[1:3]<-rep(1,3)
  #beta1<-c(1,0,0,0,0,0)
  beta1<-c(1,1,1,1,1,1,1,1,1,1)
  beta2<-c(0,1,0,0,0,0)
  #beta2<-c(0,0,0,0,0,1,1,1,1,1)
  #beta2[2]<-1
  error<-rnorm(n)
  #Y=exp((X%*%beta1)^2)+0.2*error
  Y=abs(X%*%beta1)
  #Y=log((X%*%beta1)^2+5)+0.2*error
  #Y=0.4*(X%*%beta1)^2+3*abs((X%*%beta2))^0.5+0.2*error
  #Y=3*sin(X%*%beta1/4)+0.4*(X%*%beta2)^2+0.2*error
  #Y=exp((X%*%beta1)^2)+2*exp((X%*%beta2)^2)+0.2*error
  #supp=which(result$result!=0)
  resultboot<-ctsave(X,Y,H,1,500,1023)
  col_num=col_num+length(resultboot$I)
  result<-save(X,Y,H)
  lossboot[k]<-norm(proj(resultboot$beta)-proj(beta1),"2")
  loss[k]<-norm(proj(result$vectors[,1])-proj(beta1),"2")
  lossbest[k]<-norm(proj(resultboot$optimal)-proj(beta1),"2")
  #lossboot[k]<-norm(projection(resultboot$beta)-projection(cbind(beta1,beta2)),"2")
  #loss[k]<-norm(projection(result$vectors[,1:2])-projection(cbind(beta1,beta2)),"2")
  Cv=Cv+length(which(resultboot$I<=10))
  Cvper=length(which(resultboot$I<=10)/10)
  col[k]=length(resultboot$I)
  print(resultboot$I)
  print(lossboot[k]/loss[k])
  ICv=ICv+length(which(resultboot$I>10))
}

result_table.8[7,1]<-mean(lossboot)/mean(loss)
result_table.8[7,2]<-col_num/100
result_table.8[7,3]<-Cv/500/10
bootvar<-rep(0,500)
for(i in 1:500)
{
  index=sample(1:500,500,replace=TRUE)
  bootboot=lossboot[index]
  bootloss=loss[index]
  bootvar[i]=mean(bootboot)/mean(bootloss)
}
result_table.8[7,4]<-sd(bootvar)
write.csv(bootboot,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model2rho0.2p10n500ctsave.csv")
write.csv(bootloss,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model2rho0.2p10n500save.csv")
sdcol[7]=sd(col)/sqrt(500)
sdTSR[7]=sd(Cvper)/sqrt(500)

p=10
n=200
lossboot<-rep(0,200)
loss<-rep(0,200)
lossbest<-rep(0,200)
rho<-0.8
covariance<-diag(p)
for(i in 1:p)
{
  for(j in 1:p)
  {
    if(i!=j)
    {
      covariance[i,j]=rho
    }
  }
}

H=5
beta1<-rep(0,p)
beta2<-rep(0,p)
col_num=0
ICv=0
Cv=0
for(k in 1:500)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=covariance)
  #beta<-rep(0,p)
  #beta[1:3]<-rep(1,3)
  #beta1<-c(1,0,0,0,0,0)
  beta1<-c(1,1,1,1,1,1,1,1,1,1)
  beta2<-c(0,1,0,0,0,0)
  #beta2<-c(0,0,0,0,0,1,1,1,1,1)
  #beta2[2]<-1
  error<-rnorm(n)
  #Y=exp((X%*%beta1)^2)+0.2*error
  Y=log((X%*%beta1)^2+5)+0.2*error
  #Y=0.4*(X%*%beta1)^2+3*abs((X%*%beta2))^0.5+0.2*error
  #Y=3*sin(X%*%beta1/4)+0.4*(X%*%beta2)^2+0.2*error
  #Y=exp((X%*%beta1)^2)+2*exp((X%*%beta2)^2)+0.2*error
  #supp=which(result$result!=0)
  resultboot<-ctsave(X,Y,H,1,200,1023)
  col_num=col_num+length(resultboot$I)
  result<-save(X,Y,H)
  lossboot[k]<-norm(proj(resultboot$beta)-proj(beta1),"2")
  loss[k]<-norm(proj(result$vectors[,1])-proj(beta1),"2")
  lossbest[k]<-norm(proj(resultboot$optimal)-proj(beta1),"2")
  #lossboot[k]<-norm(projection(resultboot$beta)-projection(cbind(beta1,beta2)),"2")
  #loss[k]<-norm(projection(result$vectors[,1:2])-projection(cbind(beta1,beta2)),"2")
  Cv=Cv+length(which(resultboot$I<=10))
  Cvper=length(which(resultboot$I<=10)/10)
  col[k]=length(resultboot$I)
  print(resultboot$I)
  ICv=ICv+length(which(resultboot$I>10))
  print(lossboot[k]/loss[k])
}

result_table.8[8,1]<-mean(lossboot)/mean(loss)
result_table.8[8,2]<-col_num/100
result_table.8[8,3]<-Cv/500/10
bootvar<-rep(0,500)
for(i in 1:500)
{
  index=sample(1:500,500,replace=TRUE)
  bootboot=lossboot[index]
  bootloss=loss[index]
  bootvar[i]=mean(bootboot)/mean(bootloss)
}
result_table.8[8,4]<-sd(bootvar)

write.csv(bootboot,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model2rho0.8p10n500ctsave.csv")
write.csv(bootloss,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model2rho0.8p10n500save.csv")
sdcol[8]=sd(col)/sqrt(500)
sdTSR[8]=sd(Cvper)/sqrt(500)

p=6
n=500
lossboot<-rep(0,200)
loss<-rep(0,200)
lossbest<-rep(0,200)
rho<-0.2
covariance<-diag(p)
for(i in 1:p)
{
  for(j in 1:p)
  {
    if(i!=j)
    {
      covariance[i,j]=rho
    }
  }
}

H=5
beta1<-rep(0,p)
beta2<-rep(0,p)
col_num=0
ICv=0
Cv=0
for(k in 1:500)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=covariance)
  #beta<-rep(0,p)
  #beta[1:3]<-rep(1,3)
  #beta1<-c(1,0,0,0,0,0)
  beta1<-c(1,0,0,0,0,0)
  beta2<-c(0,1,0,0,0,0)
  #beta2<-c(0,0,0,0,0,1,1,1,1,1)
  #beta2[2]<-1
  error<-rnorm(n)
  #Y=exp((X%*%beta1)^2)+0.2*error
  #Y=log((X%*%beta1)^2+5)+0.2*error
  Y=0.4*(X%*%beta1)^2+3*abs((X%*%beta2))^0.5+0.2*error
  #Y=3*sin(X%*%beta1/4)+0.4*(X%*%beta2)^2+0.2*error
  #Y=exp((X%*%beta1)^2)+2*exp((X%*%beta2)^2)+0.2*error
  #supp=which(result$result!=0)
  resultboot<-ctsave(X,Y,H,2,500,3)
  col_num=col_num+length(resultboot$I)
  result<-save(X,Y,H)
  #lossboot[k]<-norm(proj(resultboot$beta)-proj(beta1),"2")
  #loss[k]<-norm(proj(result$vectors[,1])-proj(beta1),"2")
  #lossbest[k]<-norm(proj(resultboot$optimal)-proj(beta1),"2")
  lossboot[k]<-norm(projection(resultboot$beta)-projection(cbind(beta1,beta2)),"2")
  loss[k]<-norm(projection(result$vectors[,1:2])-projection(cbind(beta1,beta2)),"2")
  lossbest[k]<-norm(proj(resultboot$optimal)-proj(cbind(beta1,beta2)),"2")
  Cv=Cv+length(which(resultboot$I<=2))
  Cvper=length(which(resultboot$I<=2)/2)
  col[k]=length(resultboot$I)
  print(resultboot$I)
  ICv=ICv+length(which(resultboot$I>2))
}

result_table.8[9,1]<-mean(lossboot)/mean(loss)
result_table.8[9,2]<-col_num/100
result_table.8[9,3]<-Cv/500/2
bootvar<-rep(0,500)
for(i in 1:500)
{
  index=sample(1:500,500,replace=TRUE)
  bootboot=lossboot[index]
  bootloss=loss[index]
  bootvar[i]=mean(bootboot)/mean(bootloss)
}
result_table.8[9,4]<-sd(bootvar)

write.csv(bootboot,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model3rho0.2p6ctsave.csv")
write.csv(bootloss,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model3rho0.2p6save.csv")
sdcol[9]=sd(col)/sqrt(500)
sdTSR[9]=sd(Cvper)/sqrt(500)


p=6
n=500
lossboot<-rep(0,200)
loss<-rep(0,200)
lossbest<-rep(0,200)
rho<-0.8
covariance<-diag(p)
for(i in 1:p)
{
  for(j in 1:p)
  {
    if(i!=j)
    {
      covariance[i,j]=rho
    }
  }
}

H=5
beta1<-rep(0,p)
beta2<-rep(0,p)
col_num=0
ICv=0
Cv=0
for(k in 1:500)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=covariance)
  #beta<-rep(0,p)
  #beta[1:3]<-rep(1,3)
  #beta1<-c(1,0,0,0,0,0)
  beta1<-c(1,0,0,0,0,0)
  beta2<-c(0,1,0,0,0,0)
  #beta2<-c(0,0,0,0,0,1,1,1,1,1)
  #beta2[2]<-1
  error<-rnorm(n)
  #Y=exp((X%*%beta1)^2)+0.2*error
  #Y=log((X%*%beta1)^2+5)+0.2*error
  Y=0.4*(X%*%beta1)^2+3*abs((X%*%beta2))^0.5+0.2*error
  #Y=3*sin(X%*%beta1/4)+0.4*(X%*%beta2)^2+0.2*error
  #Y=exp((X%*%beta1)^2)+2*exp((X%*%beta2)^2)+0.2*error
  #supp=which(result$result!=0)
  resultboot<-ctsave(X,Y,H,2,500,2)
  col_num=col_num+length(resultboot$I)
  result<-save(X,Y,H)
  #lossboot[k]<-norm(proj(resultboot$beta)-proj(beta1),"2")
  #loss[k]<-norm(proj(result$vectors[,1])-proj(beta1),"2")
  #lossbest[k]<-norm(proj(resultboot$optimal)-proj(beta1),"2")
  lossboot[k]<-norm(projection(resultboot$beta)-projection(cbind(beta1,beta2)),"2")
  loss[k]<-norm(projection(result$vectors[,1:2])-projection(cbind(beta1,beta2)),"2")
  lossbest[k]<-norm(proj(resultboot$optimal)-proj(cbind(beta1,beta2)),"2")
  Cv=Cv+length(which(resultboot$I==2))
  Cvper=length(which(resultboot$I==2)/2)
  col[k]=length(resultboot$I)
  print(resultboot$I)
  print(lossboot[k]/loss[k])
  ICv=ICv+length(which(resultboot$I>1))
}

res=rep(0,31)
for(k in 1:500)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=covariance)
  beta1<-c(1,0,0,0,0,0)
  beta2<-c(0,1,0,0,0,0)
  error<-rnorm(n)
  Y=0.4*(X%*%beta1)^2+3*abs((X%*%beta2))^0.5+0.2*error
  res=res+ctsave_optimal(X,Y,H,2,cbind(beta1,beta2))
  print(k)
}

result_table.8[10,1]<-mean(lossboot)/mean(loss)
result_table.8[10,2]<-col_num/100
result_table.8[10,3]<-Cv/500
#bootvar<-rep(0,500)
#for(i in 1:500)
##{
#  index=sample(1:500,500,replace=TRUE)
#  bootboot=lossboot[index]
#  bootloss=loss[index]
#  bootvar[i]=mean(bootboot)/mean(bootloss)
#}
result_table.8[10,4]<-sd(bootvar)
result_table.8[10,5]<-mean(lossbest)/mean(loss)

write.csv(bootboot,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model3rho0.8p6ctsave.csv")
write.csv(bootloss,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model3rho0.8p6save.csv")
sdcol[10]=sd(col)/sqrt(500)
sdTSR[10]=sd(Cvper)/sqrt(500)

p=10
n=500
lossboot<-rep(0,200)
loss<-rep(0,200)
lossbest<-rep(0,200)
rho<-0.2
covariance<-diag(p)
for(i in 1:p)
{
  for(j in 1:p)
  {
    if(i!=j)
    {
      covariance[i,j]=rho
    }
  }
}

H=5
beta1<-rep(0,p)
beta2<-rep(0,p)
col_num=0
ICv=0
Cv=0
for(k in 1:500)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=covariance)
  #beta<-rep(0,p)
  #beta[1:3]<-rep(1,3)
  #beta1<-c(1,0,0,0,0,0)
  beta1<-c(1,0,0,0,0,0,0,0,0,0)
  beta2<-c(0,1,0,0,0,0,0,0,0,0)
  #beta2<-c(0,0,0,0,0,1,1,1,1,1)
  #beta2[2]<-1
  error<-rnorm(n)
  #Y=exp((X%*%beta1)^2)+0.2*error
  #Y=log((X%*%beta1)^2+5)+0.2*error
  Y=0.4*(X%*%beta1)^2+3*abs((X%*%beta2))^0.5+0.2*error
  #Y=3*sin(X%*%beta1/4)+0.4*(X%*%beta2)^2+0.2*error
  #Y=exp((X%*%beta1)^2)+2*exp((X%*%beta2)^2)+0.2*error
  #supp=which(result$result!=0)
  resultboot<-ctsave(X,Y,H,2,50,3)
  col_num=col_num+length(resultboot$I)
  result<-save(X,Y,H)
  #lossboot[k]<-norm(proj(resultboot$beta)-proj(beta1),"2")
  #loss[k]<-norm(proj(result$vectors[,1])-proj(beta1),"2")
  #lossbest[k]<-norm(proj(resultboot$optimal)-proj(beta1),"2")
  lossboot[k]<-norm(projection(resultboot$beta)-projection(cbind(beta1,beta2)),"2")
  loss[k]<-norm(projection(result$vectors[,1:2])-projection(cbind(beta1,beta2)),"2")
  lossbest[k]<-norm(proj(resultboot$optimal)-proj(cbind(beta1,beta2)),"2")
  Cv=Cv+length(which(resultboot$I<=2))
  Cvper=length(which(resultboot$I<=2)/2)
  col[k]=length(resultboot$I)
  print(resultboot$I)
  ICv=ICv+length(which(resultboot$I>2))
}

result_table.8[11,1]<-mean(lossboot)/mean(loss)
result_table.8[11,2]<-col_num/100
result_table.8[11,3]<-Cv/500/2
bootvar<-rep(0,500)
for(i in 1:500)
{
  index=sample(1:500,500,replace=TRUE)
  bootboot=lossboot[index]
  bootloss=loss[index]
  bootvar[i]=mean(bootboot)/mean(bootloss)
}
result_table.8[11,4]<-sd(bootvar)

write.csv(bootboot,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model3rho0.2p10ctsave.csv")
write.csv(bootloss,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model3rho0.2p10save.csv")
sdcol[11]=sd(col)/sqrt(500)
sdTSR[11]=sd(Cvper)/sqrt(500)

p=10
n=500
lossboot<-rep(0,200)
loss<-rep(0,200)
lossbest<-rep(0,200)
rho<-0.8
covariance<-diag(p)
for(i in 1:p)
{
  for(j in 1:p)
  {
    if(i!=j)
    {
      covariance[i,j]=rho
    }
  }
}

H=5
beta1<-rep(0,p)
beta2<-rep(0,p)
col_num=0
ICv=0
Cv=0
for(k in 52:500)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=covariance)
  #beta<-rep(0,p)
  #beta[1:3]<-rep(1,3)
  #beta1<-c(1,0,0,0,0,0)
  beta1<-c(1,0,0,0,0,0,0,0,0,0)
  beta2<-c(0,1,0,0,0,0,0,0,0,0)
  #beta2<-c(0,0,0,0,0,1,1,1,1,1)
  #beta2[2]<-1
  error<-rnorm(n)
  #Y=exp((X%*%beta1)^2)+0.2*error
  #Y=log((X%*%beta1)^2+5)+0.2*error
  Y=0.4*(X%*%beta1)^2+3*abs((X%*%beta2))^0.5+0.2*error
  #Y=3*sin(X%*%beta1/4)+0.4*(X%*%beta2)^2+0.2*error
  #Y=exp((X%*%beta1)^2)+2*exp((X%*%beta2)^2)+0.2*error
  #supp=which(result$result!=0)
  resultboot<-ctsave(X,Y,H,2,500,2)
  col_num=col_num+length(resultboot$I)
  result<-save(X,Y,H)
  #lossboot[k]<-norm(proj(resultboot$beta)-proj(beta1),"2")
  #loss[k]<-norm(proj(result$vectors[,1])-proj(beta1),"2")
  #lossbest[k]<-norm(proj(resultboot$optimal)-proj(beta1),"2")
  lossboot[k]<-norm(projection(resultboot$beta)-projection(cbind(beta1,beta2)),"2")
  loss[k]<-norm(projection(result$vectors[,1:2])-projection(cbind(beta1,beta2)),"2")
  lossbest[k]<-norm(proj(resultboot$optimal)-proj(cbind(beta1,beta2)),"2")
  Cv=Cv+length(which(resultboot$I==2))
  Cvper=length(which(resultboot$I==2))
  col[k]=length(resultboot$I)
  print(resultboot$I)
  ICv=ICv+length(which(resultboot$I>1))
}

res=rep(0,1023)
for(k in 1:500)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=covariance)
  #beta<-rep(0,p)
  #beta[1:3]<-rep(1,3)
  #beta1<-c(1,0,0,0,0,0)
  beta1<-c(1,0,0,0,0,0,0,0,0,0)
  beta2<-c(0,1,0,0,0,0,0,0,0,0)
  #beta2<-c(0,0,0,0,0,1,1,1,1,1)
  #beta2[2]<-1
  error<-rnorm(n)
  #Y=exp((X%*%beta1)^2)+0.2*error
  #Y=log((X%*%beta1)^2+5)+0.2*error
  Y=0.4*(X%*%beta1)^2+3*abs((X%*%beta2))^0.5+0.2*error
  res=res+ctsave_optimal(X,Y,H,2,cbind(beta1,beta2))
  print(k)
}

result_table.8[12,1]<-mean(lossboot)/mean(loss)
result_table.8[12,2]<-col_num/100
result_table.8[12,3]<-Cv/500
bootvar<-rep(0,500)
for(i in 1:500)
{
  index=sample(1:500,500,replace=TRUE)
  bootboot=lossboot[index]
  bootloss=loss[index]
  bootvar[i]=mean(bootboot)/mean(bootloss)
}
result_table.8[12,4]<-sd(bootvar)

write.csv(bootboot,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model3rho0.8p10ctsave.csv")
write.csv(bootloss,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model3rho0.8p10save.csv")

sdcol[12]=sd(col)/sqrt(500)
sdTSR[12]=sd(Cvper)/sqrt(500)


p=6
n=500
lossboot<-rep(0,200)
loss<-rep(0,200)
lossbest<-rep(0,200)
rho<-0.2
covariance<-diag(p)
for(i in 1:p)
{
  for(j in 1:p)
  {
    if(i!=j)
    {
      covariance[i,j]=rho
    }
  }
}

H=5
beta1<-rep(0,p)
beta2<-rep(0,p)
col_num=0
ICv=0
Cv=0
for(k in 1:500)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=covariance)
  #beta<-rep(0,p)
  #beta[1:3]<-rep(1,3)
  #beta1<-c(1,0,0,0,0,0)
  beta1<-c(1,0,0,0,0,0)
  beta2<-c(0,1,0,0,0,0)
  #beta2<-c(0,0,0,0,0,1,1,1,1,1)
  #beta2[2]<-1
  error<-rnorm(n)
  #Y=exp((X%*%beta1)^2)+0.2*error
  #Y=log((X%*%beta1)^2+5)+0.2*error
  #Y=0.4*(X%*%beta1)^2+3*abs((X%*%beta2))^0.5+0.2*error
  Y=3*sin(X%*%beta1/4)+0.4*(X%*%beta2)^2+0.2*error
  #Y=exp((X%*%beta1)^2)+2*exp((X%*%beta2)^2)+0.2*error
  #supp=which(result$result!=0)
  resultboot<-ctsave(X,Y,H,2,200,3)
  col_num=col_num+length(resultboot$I)
  result<-save(X,Y,H)
  #lossboot[k]<-norm(proj(resultboot$beta)-proj(beta1),"2")
  #loss[k]<-norm(proj(result$vectors[,1])-proj(beta1),"2")
  #lossbest[k]<-norm(proj(resultboot$optimal)-proj(beta1),"2")
  lossboot[k]<-norm(projection(resultboot$beta)-projection(cbind(beta1,beta2)),"2")
  loss[k]<-norm(projection(result$vectors[,1:2])-projection(cbind(beta1,beta2)),"2")
  lossbest[k]<-norm(proj(resultboot$optimal)-proj(cbind(beta1,beta2)),"2")
  Cv=Cv+length(which(resultboot$I<=2))
  Cvper=length(which(resultboot$I<=2)/2)
  col[k]=length(resultboot$I)
  print(resultboot$I)
  print(lossboot[k]/loss[k])
  ICv=ICv+length(which(resultboot$I>2))
}

result_table.8[13,1]<-mean(lossboot)/mean(loss)
result_table.8[13,2]<-col_num/100
result_table.8[13,3]<-Cv/500/2
bootvar<-rep(0,500)
for(i in 1:500)
{
  index=sample(1:500,500,replace=TRUE)
  bootboot=lossboot[index]
  bootloss=loss[index]
  bootvar[i]=mean(bootboot)/mean(bootloss)
}
result_table.8[13,4]<-sd(bootvar)

write.csv(bootboot,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model4rho0.2p6ctsave.csv")
write.csv(bootloss,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model4rho0.2p6save.csv")

sdcol[13]=sd(col)/sqrt(500)
sdTSR[13]=sd(Cvper)/sqrt(500)

p=6
n=500
lossboot<-rep(0,200)
loss<-rep(0,200)
lossbest<-rep(0,200)
rho<-0.8
covariance<-diag(p)
for(i in 1:p)
{
  for(j in 1:p)
  {
    if(i!=j)
    {
      covariance[i,j]=rho
    }
  }
}

H=5
beta1<-rep(0,p)
beta2<-rep(0,p)
col_num=0
ICv=0
Cv=0
col=rep(0,500)
for(k in 1:500)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=covariance)
  #beta<-rep(0,p)
  #beta[1:3]<-rep(1,3)
  #beta1<-c(1,0,0,0,0,0)
  beta1<-c(1,0,0,0,0,0)
  beta2<-c(0,1,0,0,0,0)
  #beta2<-c(0,0,0,0,0,1,1,1,1,1)
  #beta2[2]<-1
  error<-rnorm(n)
  #Y=exp((X%*%beta1)^2)+0.2*error
  #Y=log((X%*%beta1)^2+5)+0.2*error
  #Y=0.4*(X%*%beta1)^2+3*abs((X%*%beta2))^0.5+0.2*error
  Y=3*sin(X%*%beta1/4)+0.4*(X%*%beta2)^2+0.2*error
  #Y=exp((X%*%beta1)^2)+2*exp((X%*%beta2)^2)+0.2*error
  #supp=which(result$result!=0)
  resultboot<-ctsave(X,Y,H,2,500,3)
  col_num=col_num+length(resultboot$I)
  col[i]=length(resultboot$I)
  result<-save(X,Y,H)
  #lossboot[k]<-norm(proj(resultboot$beta)-proj(beta1),"2")
  #loss[k]<-norm(proj(result$vectors[,1])-proj(beta1),"2")
  #lossbest[k]<-norm(proj(resultboot$optimal)-proj(beta1),"2")
  lossboot[k]<-norm(projection(resultboot$beta)-projection(cbind(beta1,beta2)),"2")
  loss[k]<-norm(projection(result$vectors[,1:2])-projection(cbind(beta1,beta2)),"2")
  lossbest[k]<-norm(proj(resultboot$optimal)-proj(cbind(beta1,beta2)),"2")
  Cv=Cv+length(which(resultboot$I<=2))
  Cvper=length(which(resultboot$I<=2)/2)
  col[k]=length(resultboot$I)
  print(resultboot$I)
  print(lossboot[k]/loss[k])
  ICv=ICv+length(which(resultboot$I>2))
}
mean(lossbest)/mean(loss)
result_table.8[14,1]<-mean(lossboot)/mean(loss)
result_table.8[14,2]<-col_num/100
result_table.8[14,3]<-Cv/500/2
bootvar<-rep(0,500)
for(i in 1:500)
{
  index=sample(1:500,500,replace=TRUE)
  bootboot=lossboot[index]
  bootloss=loss[index]
  bootvar[i]=mean(bootboot)/mean(bootloss)
}
result_table.8[14,4]<-sd(bootvar)
result_table.8[14,5]<-mean(lossbest)/mean(loss)

write.csv(bootboot,"C:/Simulation_result/ACS-SAVE/Low-dimensional/Model4rho0.8p6ctsave.csv")
write.csv(bootloss,"C:/Simulation_result/ACS-SAVE/Low-dimensional/Model4rho0.8p6save.csv")
sdcol[14]=sd(col)/sqrt(500)
sdTSR[14]=sd(Cvper)/sqrt(500)

p=10
n=500
lossboot<-rep(0,500)
loss<-rep(0,500)
lossbest<-rep(0,500)
rho<-0.2
covariance<-diag(p)
for(i in 1:p)
{
  for(j in 1:p)
  {
    if(i!=j)
    {
      covariance[i,j]=rho
    }
  }
}

H=5
beta1<-rep(0,p)
beta2<-rep(0,p)
col_num=0
ICv=0
Cv=0
for(k in 1:500)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=covariance)
  #beta<-rep(0,p)
  #beta[1:3]<-rep(1,3)
  #beta1<-c(1,0,0,0,0,0)
  beta1<-c(1,0,0,0,0,0,0,0,0,0)
  beta2<-c(0,1,0,0,0,0,0,0,0,0)
  #beta2<-c(0,0,0,0,0,1,1,1,1,1)
  #beta2[2]<-1
  error<-rnorm(n)
  #Y=exp((X%*%beta1)^2)+0.2*error
  #Y=log((X%*%beta1)^2+5)+0.2*error
  #Y=0.4*(X%*%beta1)^2+3*abs((X%*%beta2))^0.5+0.2*error
  Y=3*sin(X%*%beta1/4)+0.4*(X%*%beta2)^2+0.2*error
  #Y=exp((X%*%beta1)^2)+2*exp((X%*%beta2)^2)+0.2*error
  #supp=which(result$result!=0)
  resultboot<-ctsave(X,Y,H,2,500,3)
  col_num=col_num+length(resultboot$I)
  result<-save(X,Y,H)
  #lossboot[k]<-norm(proj(resultboot$beta)-proj(beta1),"2")
  #loss[k]<-norm(proj(result$vectors[,1])-proj(beta1),"2")
  #lossbest[k]<-norm(proj(resultboot$optimal)-proj(beta1),"2")
  lossboot[k]<-norm(projection(resultboot$beta)-projection(cbind(beta1,beta2)),"2")
  loss[k]<-norm(projection(result$vectors[,1:2])-projection(cbind(beta1,beta2)),"2")
  lossbest[k]<-norm(proj(resultboot$optimal)-proj(cbind(beta1,beta2)),"2")
  Cv=Cv+length(which(resultboot$I<=2))
  Cvper=length(which(resultboot$I<=2)/2)
  col[k]=length(resultboot$I)
  print(resultboot$I)
  print(lossboot[k]/loss[k])
  ICv=ICv+length(which(resultboot$I>2))
}

result_table.8[15,1]<-mean(lossboot)/mean(loss)
result_table.8[15,2]<-col_num/100
result_table.8[15,3]<-Cv/500/2
bootvar<-rep(0,500)
for(i in 1:500)
{
  index=sample(1:500,500,replace=TRUE)
  bootboot=lossboot[index]
  bootloss=loss[index]
  bootvar[i]=mean(bootboot)/mean(bootloss)
}
result_table.8[15,4]<-sd(bootvar)

for(k in 1:500)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=covariance)
  #beta<-rep(0,p)
  #beta[1:3]<-rep(1,3)
  #beta1<-c(1,0,0,0,0,0)
  beta1<-c(1,0,0,0,0,0,0,0,0,0)
  beta2<-c(0,1,0,0,0,0,0,0,0,0)
  error<-rnorm(n)
  Y=3*sin(X%*%beta1/4)+0.4*(X%*%beta2)^2+0.2*error
  res=res+ctsave_optimal(X,Y,H,2,cbind(beta1,beta2))
  print(k)
}


write.csv(bootboot,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model4rho0.2p10ctsave.csv")
write.csv(bootloss,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model4rho0.2p10save.csv")

sdcol[15]=sd(col)/sqrt(500)
sdTSR[15]=sd(Cvper)/sqrt(500)

p=10
n=500
lossboot<-rep(0,200)
loss<-rep(0,200)
lossbest<-rep(0,200)
rho<-0.8
covariance<-diag(p)
for(i in 1:p)
{
  for(j in 1:p)
  {
    if(i!=j)
    {
      covariance[i,j]=rho
    }
  }
}

H=5
beta1<-rep(0,p)
beta2<-rep(0,p)
col_num=0
ICv=0
Cv=0
for(k in 158:500)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=covariance)
  #beta<-rep(0,p)
  #beta[1:3]<-rep(1,3)
  #beta1<-c(1,0,0,0,0,0)
  beta1<-c(1,0,0,0,0,0,0,0,0,0)
  beta2<-c(0,1,0,0,0,0,0,0,0,0)
  #beta2<-c(0,0,0,0,0,1,1,1,1,1)
  #beta2[2]<-1
  error<-rnorm(n)
  #Y=exp((X%*%beta1)^2)+0.2*error
  #Y=log((X%*%beta1)^2+5)+0.2*error
  #Y=0.4*(X%*%beta1)^2+3*abs((X%*%beta2))^0.5+0.2*error
  Y=3*sin(X%*%beta1/4)+0.4*(X%*%beta2)^2+0.2*error
  #Y=exp((X%*%beta1)^2)+2*exp((X%*%beta2)^2)+0.2*error
  #supp=which(result$result!=0)
  resultboot<-ctsave(X,Y,H,2,500,3)
  col_num=col_num+length(resultboot$I)
  result<-save(X,Y,H)
  #lossboot[k]<-norm(proj(resultboot$beta)-proj(beta1),"2")
  #loss[k]<-norm(proj(result$vectors[,1])-proj(beta1),"2")
  #lossbest[k]<-norm(proj(resultboot$optimal)-proj(beta1),"2")
  lossboot[k]<-norm(projection(resultboot$beta)-projection(cbind(beta1,beta2)),"2")
  loss[k]<-norm(projection(result$vectors[,1:2])-projection(cbind(beta1,beta2)),"2")
  lossbest[k]<-norm(proj(resultboot$optimal)-proj(cbind(beta1,beta2)),"2")
  Cv=Cv+length(which(resultboot$I<=2))
  Cvper=length(which(resultboot$I<=2)/2)
  col[k]=length(resultboot$I)
  print(resultboot$I)
  print(lossboot[k]/loss[k])
  ICv=ICv+length(which(resultboot$I>2))
}

  result_table.8[16,1]<-mean(lossboot)/mean(loss)
result_table.8[16,2]<-col_num/100
result_table.8[16,3]<-Cv/500/2
bootvar<-rep(0,500)
for(i in 1:500)
{
  index=sample(1:500,500,replace=TRUE)
  bootboot=lossboot[index]
  bootloss=loss[index]
  bootvar[i]=mean(bootboot)/mean(bootloss)
}
result_table.8[16,4]<-sd(bootvar)
write.csv(bootboot,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model4rho0.8p10ctsave.csv")
write.csv(bootloss,"D:/Simulation_result/ACS-SAVE/Low-dimensional/Model4rho0.8p10save.csv")
sdcol[16]=sd(col)/sqrt(500)
sdTSR[16]=sd(Cvper)/sqrt(500)

