##take power of a positive definite matrix
matpower <- function(a,alpha){
  small <- .00000001
  if (length(c(a))==1) {return(a^(alpha))} else {
    p1<-nrow(a)
    eva<-eigen(a)$values
    eve<-eigen(a)$vectors
    eve<-eve/t(matrix((diag(t(eve)%*%eve)^0.5),p1,p1))
    index<-(1:p1)[abs(eva)>small]
    evai<-eva
    evai[index]<-(eva[index])^(alpha)
    ai<-eve%*%diag(evai,length(evai))%*%t(eve)
    return(ai)}
}

##take projection matrix of v
projection<-function(X)
{
  return(X%*%ginv(t(X)%*%X)%*%t(X))
}

proj<-function(v)
{return(v%*%matpower(t(v)%*%v,-1)%*%t(v))}

##take distance of two column spaces
dist<-function(u,v)
{return(sqrt(sum((proj(u)-proj(v))^2)))}

##standarize x
stand <- function(x){
  n<-nrow(x)
  p<-ncol(x)
  xb <- apply(x, 2, mean)
  xb <- t(matrix(xb, p, n))
  x1 <- x - xb
  sigma <- t(x1) %*% (x1)/(n-1)
  eva <- eigen(sigma)$values
  eve <- eigen(sigma)$vectors
  sigmamrt <- eve %*% diag(1/sqrt(eva)) %*% t(eve)
  z <- sigmamrt %*% t(x1)
  return(t(z))
}

##take rsquare
rsquare<-function(u,v)
{
  k=ncol(as.matrix(u))
  if(all(u==0)||all(v==0))
  {
    return(0)
  }
  result<-matpower(var(u),-0.5)%*%cov(u,v)%*%matpower(var(v),-1)%*%cov(u,v)%*%matpower(var(u),-0.5)
  return(sum(sqrt(svd(result)$d))/k)
}

##deal with list
deal_list<-function(x)
{
  return(as.vector(x)[-1])
}

allzero<-function(X)  
{
  return(all(X==0))
}

sumzero<-function(X)
{
  return(length(which(X!=0)))
}

thresholding<-function(X,threshold)
{
  return(all(sum(X^2)<=threshold))
}

##slicing Y
slicing<-function(y,H) 
{
  n<-length(y)
  if (length(levels(as.factor(y)))>H)
  {
    ytilde<-rep(0,H+1)
    ytilde[1]<-min(y)
    for (h in 1:(H-1))
    {
      ytilde[h+1]<-quantile(y,h/H)
      ytilde[h+1]<-max(y[y<=ytilde[h+1]])
      if (ytilde[h+1]==ytilde[h]) ytilde[h+1]=min(y[y>ytilde[h]])
    }  
  }
  if (length(levels(as.factor(y)))<=H)
  {
    H <- length(levels(as.factor(y)))
    ytilde<-rep(0,H+1)
    ytilde[1]=min(y)
    for (h in 1:(H-1))
    {
      ytilde[h+1]<-min(y[y>ytilde[h]])
    }
  } 
  ytilde[H+1]=max(y)+1
  prop<-rep(1,H)
  for (i in 1:H)
  {
    prop[i] = sum((y >= ytilde[i])&(y < ytilde[i+1]))/n
  }
  res<-list()
  res$H<-H
  res$ytilde<-ytilde
  res$prop<-prop
  return(res)
}

mhat_dr<-function(x,y,H){
  n<-nrow(x)
  p<-ncol(x)
  z<-stand(x)
  dy <- slicing(y,H)
  H <- dy$H
  ytilde <- dy$ytilde
  prop <- dy$prop
  ind<-matrix(0,n,H)
  zbar<-matrix(0,p,H)
  for (j in 1:(H-1))
  {
    ind[,j]<-((y >= ytilde[j])&(y < ytilde[j+1]))
    zbar[,j]<- (t(z)%*%(ind[,j]))/sum(ind[,j])
  }
  ind[,H]<-(y >= ytilde[H])
  zbar[,H]<- (t(z)%*%(ind[,H]))/sum(ind[,H])
  A<-matrix(0,p,p)
  B<-matrix(0,p,p)
  C<-0
  for (q in 1:H)
  {
    Z<-(t(z))[,ind[,q]==1]-zbar[,q]  
    A<-A + prop[q]*((Z%*%t(Z)/(sum(ind[,q])-1)+zbar[,q]%*%t(zbar[,q]))%*%  
                      (Z%*%t(Z)/(sum(ind[,q])-1)+zbar[,q]%*%t(zbar[,q])) - diag(1,p))
    B<-B + sqrt(prop[j])*(zbar[,q]%*%t(zbar[,q]))
    C<-C + sqrt(prop[j])*(t(zbar[,q])%*%zbar[,q])
  }
  C<-as.vector(C)
  M<-2*A + 2*(B%*%B) + 2*B*C
  return(M)
}
