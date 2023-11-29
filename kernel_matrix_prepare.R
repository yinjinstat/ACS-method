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

save_kernel<-function(x,y,H)
{
  n <- dim(x)[1]
  p <-dim(x)[2]
  dy <- slicing(y,H)
  H <- dy$H
  ytilde <- dy$ytilde
  prop <- dy$prop
  freq<-prop*n
  xxbar<-numeric(0)
  for (j in 1:(H-1))
  {
    ind<-((y >= ytilde[j])&(y < ytilde[j+1]))
    freq[j]=length(which(ind==TRUE))
    S<-matrix(0,ncol=n,nrow=freq[j])
    S[,ind]<-diag(freq[j])-matrix(1/freq[j],freq[j],freq[j])
    xxbar<-cbind(xxbar,sqrt(prop[j])*(x-1/prop[j]*t(S)%*%S%*%x))
  } 
  ind<-(y >= ytilde[H])
  freq[H]=length(which(ind==TRUE))
  S<-matrix(0,ncol=n,nrow=freq[H])
  S[,ind]<-diag(freq[H])-matrix(1/freq[H],freq[H],freq[H])
  xxbar<-cbind(xxbar,sqrt(prop[H])*(x-1/prop[H]*t(S)%*%S%*%x))
  return(xxbar)
}

sir_kernel<-function(x,y,H)
{
  n <- dim(x)[1]
  p <-dim(x)[2]
  dy <- slicing(y,H)
  H <- dy$H
  ytilde <- dy$ytilde
  prop <- dy$prop
  freq<-prop*n
  xxbar<-numeric(0)
  for (j in 1:(H-1))
  {
    ind<-((y >= ytilde[j])&(y < ytilde[j+1]))
    freq[j]<-length(which(ind==TRUE))
    S<-matrix(0,ncol=n,nrow=freq[j])
    S[,ind]<-diag(freq[j])
    xxbar<-cbind(xxbar,sqrt(prop[j])*n*(rep(1/n,n)-t(S)%*%rep(1/freq[j],freq[j])))
  } 
  ind<-(y >= ytilde[H])
  freq[H]<-length(which(ind==TRUE))
  S<-matrix(0,ncol=n,nrow=freq[H])
  S[,ind]<-diag(freq[H])
  xxbar<-cbind(xxbar,sqrt(prop[H])*(rep(1/n,n)-t(S)%*%rep(1/freq[H],freq[H])))
  return(xxbar)
}
