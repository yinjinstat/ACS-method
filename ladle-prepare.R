##########################################################
###  The file contains functions needed in "ladle.R":  ### 
###   1. standardization of x (stand)                  ###
###   2. power of matrix (matpower)                    ###
###   3. (negative) power of matrix (mppower)          ###
###   4. center x to have zero mean (center)           ###
###   4. Mhat for CCA (mhat_cca)                       ###
###   5. Mhat for ICA (mhat_ica)                       ###
###   6. vector transform (diffin & diffbw)            ###
###       -- used in sequential test for CCA           ###
###   7. discretize Y into slices for SDR (slicing)    ###
###   8. Mhat for directional regression (mhat_dr)     ###
###   9. Mhat for (quadratic) kernel SIR (mhat_ksir)   ###
###   10. action from a sequence of p-values (pick)    ###  
###  x is n*p matrix                                   ###
###  y is n-dim vector (n*q matrix in CCA)             ###
##########################################################

####################################
###       standardize x          ### 
###  returns standardized x (z)  ###              
####################################
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

###############################
###    power of a matrix    ###
###  a is original matrix   ###
###  alpha is the power     ###
###  returns new matrix     ###
###############################
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

#######################################
###     Moore-Penrose type power    ###
###    ignoring criterion=ignore    ###                 
#######################################
mppower = function(matrix,power,ignore){
eig = eigen(matrix)
eval = eig$values
evec = eig$vectors
m = length(eval[abs(eval)>ignore])
tmp = evec[,1:m]%*%diag(eval[1:m]^power,m)%*%
    t(evec[,1:m])
return(tmp)
}

################################
###  center X (n*p matrix)  ####
################################
center = function(x){
return(t(t(x)-apply(x,2,mean)))}


######################
###  Mhat for CCA  ###
######################
mhat_cca<-function(x,y)
{
  n<-nrow(x)
  p<-ncol(x)
  q<-ncol(y)
  nsx<-matpower(cov(x),-0.5)
  nvy<-matpower(cov(y),-1)
  cxy<-cov(x,y)
  M<-nsx%*%cxy%*%nvy%*%t(cxy)%*%nsx  
return(M)
}

###############################
#####    Mhat for ICA     #####
###############################
mhat_ica<-function(x)
{
n<-nrow(x)
p<-ncol(x)
z<-stand(x)
w<-apply(z^2,1,sum)
m<-t(z)%*%diag(w)%*%z/n - (p+2)*diag(p)
#m<-m%*%m
return(m%*%m)
}

##################################################
### vector transformation needed in ST for CCA ###
##################################################
diffin<-function(v)
{
 s<-length(v)
 if (s==1) {return(1)} else {m<-(v%*%t(rep(1,s))-rep(1,s)%*%t(v))
 r<-numeric(0)
 for (i in 1:(s-1)) {r<-c(r,m[i,(i+1):s])}
return(r)}
}

diffbw<-function(u,v)
{
 s<-length(u)
 r<-numeric(0)
 for (i in 1:s)
 {r<-c(r,u[i]-v)}
return(r)
}


#####################################
### discretizing y into H slices  ###
#####################################
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

#########################################
###  mhat for directional regression  ###
#########################################
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


#####################################################
###              mhat for kernel SIR              ###      
###              (quadratic kernel)               ###
###   x is the original predictor                 ###
###   phi is the transformed predictor            ###
###     -- in the kernel space                    ### 
###   it returns:                                 ###
###    1. mhat - Mhat, but not used in "ladle"    ###
###    2. nx - phi, kernel-transformed x          ###
###    3. evectors & evalues                      ###
###     -- eigenvectors and eigenvalues from      ###
###        singular value decomposition (fast)    ###
#####################################################
mhat_ksir<-function(x=x,phi=phi,y=y,nslices=nslices)
{
n<-length(y)
H<-nslices
if ((missing(phi))) 
 {
  p<-ncol(x)
  phi=cbind(1,x)
  for(i in 1:p) 
   {phi = cbind(phi,x[,i]*x[,1:i])}
 }
oy = order(y)   # ordered y
ophi = phi[oy,] # ordered phi
bslice = H        # number of slices; b means beween
wslice = n/bslice  # number of observations in each slice: w means within
ephiy=numeric()    # E(Y|slice)
for(i in 1:bslice){
ephiy=rbind(ephiy,apply(ophi[((i-1)*wslice+1):(i*wslice),],2,mean))}
ephi = apply(phi,2,mean)
m<-var(t(t(ephiy)-ephi))
vd<-svd(t(t(ephiy)-ephi))
v<-vd$v[,order((vd$d)^2,decreasing=TRUE)]
vc<-svd(diag(1,ncol(phi))-v%*%t(v))$u
lam<-(vd$d)^2/(H-1)
res<-list()
res$mhat<-m
res$nx<-phi
res$evectors<-cbind(v,vc)
res$evalues<-lam
return(res)
}


##################################
### pick d from p-values of ST ###
##################################
pick<-function(v,alpha){
 p<-length(v)
 m=0
 while ((m<p)&(v[m+1]<alpha))
 {m=m+1}
 return(m)
}

#######################################
###  projection matrix of matrix v  ###
#######################################
proj<-function(v)
{return(v%*%matpower(t(v)%*%v,-1)%*%t(v))}

############################################
###  distance between two linear spaces  ###
############################################
dist<-function(u,v)
{return(sqrt(sum((proj(u)-proj(v))^2)))}

###################
###     end     ###
###################





 