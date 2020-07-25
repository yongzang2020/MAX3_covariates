####################################################################################################################################################
## Rtest can be used to detect the genetic association between a diallelic marker and a binary phenotype outcome (e.g., case-control status). ######
## Rtest will run the MAX3, LRT and score tests optimal for the recessive, additive and dominant models and reports the corresponding p-values. ####
## Arguments:                                                                                                                                   ####
## data: This is the data contains the genotype, phenotype and covariates information. Each row of data represents a different individual. #########
## The first column of data records the phenotype information, with 0 representing a case status and 1 representing a control status. ##############
## The second column of data records the genotype information, with 0, 1 and 2 representing the number of risk allele. #############################
## The remaining columns of data records all the covariate information. ############################################################################
####################################################################################################################################################


Rtest=function(data){
  library(mvtnorm)
  
  library(logistf)
  
  est=function(D,X){
    return( logistf(D~X)$coefficients    )
  }
  
  
  estf=function(esta,X){
    n=dim(X)[1]
    X=cbind(  rep(1,n),X )
    re=NULL
    for(i in 1:n){
      re[i]=1/(1+exp( -1*sum(X[i,]*esta  )  ))
    }
    return(re)
  }
  
  score=function(esta,z,D,X,G){
    n=length(G)
    for(i in 1:n){
      if(G[i]==1){
        G[i]=z
      }else{G[i]=G[i]/2}
    }
    int=rep(1,n)
    X=cbind(int,X)
    re=0
    for(i in 1:n){
      estf=1/(1+exp( -1*sum(X[i,]*esta  )  ))
      temp=G[i]*(D[i]-estf)
      re=re+temp  
    }
    return(re)
    
    
  }
  
  inforb=function(estpen,z,D,X,G){
    n=dim(X)[1]
    X=cbind(  rep(1,n),X )
    for(i in 1:n){
      if(G[i]==1){
        G[i]=z
      }else{G[i]=G[i]/2}
    }
    Ib=sum(G^2*(1-estpen)*estpen)
    return(Ib)
  }
  
  
  inforba=function(estpen,z,D,X,G){
    n=dim(X)[1]
    X=cbind(  rep(1,n),X )
    l=dim(X)[2]
    for(i in 1:n){
      if(G[i]==1){
        G[i]=z
      }else{G[i]=G[i]/2}
    }
    Iba=NULL
    for(i in 1:l){
      Iba[i]=sum(X[,i]*G*(1-estpen)*estpen)
    }
    return(Iba)
  }
  
  infora=function(estpen,D,X,G){
    n=dim(X)[1]
    X=cbind(  rep(1,n),X )
    l=dim(X)[2]
    Ia=matrix(0, nrow=l,ncol=l)
    for(i in 1:l){
      for(j in 1:l){
        Ia[i,j]=sum(X[,i]*X[,j]*(1-estpen)*estpen)
      }
    }
    return(Ia)
  }
  
  inforb1b2=function(estpen,z1,z2,D,X,G){
    n=dim(X)[1]
    G1=NULL
    G2=NULL
    X=cbind(  rep(1,n),X )
    for(i in 1:n){
      if(G[i]==1){
        G1[i]=z1
        G2[i]=z2
      }else{
        G1[i]=G[i]/2
        G2[i]=G[i]/2
      }
    }
    I1=sum(G1^2*(1-estpen)*estpen)
    I12=sum(G1*G2*(1-estpen)*estpen)
    I2=sum(G2^2*(1-estpen)*estpen)
    re=matrix(  c(I1,I12,I12,I2),nrow=2 )
    return(re)
  }
  
  v=function(Ib,Iba,Ia){
    return(Ib-t(as.matrix(Iba))%*%solve(Ia)%*%as.matrix(Iba))
  }
  
  
  
  stest=function(z,D,X,G,esta,Ib,Iba,Ia){
    num=score(esta,z,D,X,G)
    den=sqrt(v(Ib,Iba,Ia))
    re=num/den ## test statistic
    return(re)
  }
  
  
  chi=function(s0,s1,r){
    return(   (s0^2+s1^2-2*r*s0*s1)/(1-r^2)  )
  }
    D=as.matrix(data[,1])
    G=as.matrix(data[,2])
    X=as.matrix(data[, -c(1,2)])
    esta=est(D,X)
    estpen=estf(esta,X)
    Ia=infora(estpen,D,X,G)
    Ib0=inforb(estpen,0,D,X,G)
    Ib05=inforb(estpen,0.5,D,X,G)
    Ib1=inforb(estpen,1,D,X,G)
    Iba0=inforba(estpen,0,D,X,G)
    Iba05=inforba(estpen,0.5,D,X,G)
    Iba1=inforba(estpen,1,D,X,G)
    s0=stest(0,D,X,G,esta,Ib0,Iba0,Ia)
    s05=stest(0.5,D,X,G,esta,Ib05,Iba05,Ia)
    s1=stest(1,D,X,G,esta,Ib1,Iba1,Ia)
    p0=  as.numeric(2*(1-pnorm( abs( s0  ))))
    p05= as.numeric(2*(1-pnorm( abs( s05 ))))
    p1=  as.numeric(2*(1-pnorm( abs( s1  ))))
    
    Iba005=rbind(Iba0,Iba05)
    Iba01=rbind(Iba0,Iba1)
    Iba051=rbind(Iba05,Iba1)
    Ib1b2005=inforb1b2(estpen,0,0.5,D,X,G)
    Ib1b201=inforb1b2(estpen,0,1,D,X,G)
    Ib1b2051=inforb1b2(estpen,0.5,1,D,X,G)
    C005=Ib1b2005-as.matrix(Iba005)%*%solve(Ia)%*%t(as.matrix(Iba005))
    C01=Ib1b201-as.matrix(Iba01)%*%solve(Ia)%*%t(as.matrix(Iba01))
    C051=Ib1b2051-as.matrix(Iba051)%*%solve(Ia)%*%t(as.matrix(Iba051))
    rho005=C005[1,2]/sqrt(C005[1,1]*C005[2,2])
    rho01=C01[1,2]/sqrt(C01[1,1]*C01[2,2])
    rho051=C051[1,2]/sqrt(C051[1,1]*C051[2,2])
    
    chisq=chi(s0,s1,rho01)
    pchi=as.numeric(1-pchisq(chisq,df=2))
    
    vacov=matrix(c(1,rho005,rho01,rho005,1,rho051,rho01,rho051,1),ncol=3)
    zmax=max(abs(s0),abs(s05),abs(s1))
    if(zmax>10000){zmax=10000}
    pmax=1-pmvnorm(lower=rep(-1*zmax,3),upper=rep(zmax,3),sigma=vacov)[1]
  return(list("p-value for MAX3"=pmax, "p-value for LRT"=pchi, "p-value for S(0) (optimal for recessive model)"=p0, "p-value for S(1/2) (optimal for additive model)"=p05, "p-value for S(1) (optimal for dominant model)"=p1   )) 

  
}








