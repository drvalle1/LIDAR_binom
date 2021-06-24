LIDAR_binom_theta_given_phi=function(y,n,ncomm,gamma,ngibbs,nburn,
                                     phi.post){
  #useful stuff
  nloc=nrow(y)
  nspp=ncol(y)
  hi=0.999999
  lo=0.000001
  NminusY=n-y
  
  #initial values
  theta=matrix(1/ncomm,nloc,ncomm)
  ooo=1
  phi=matrix(phi.post[ooo,],ncomm,nspp)
  npost=nrow(phi)
  # ordem.phi=matrix(1:(ncomm*nspp),ncomm,nspp)
  
  array.LSKP=array(NA,dim=c(nloc,nspp,ncomm,2))
  for (i in 1:nloc){
    for (j in 1:nspp){
      for (oo in 1:2){
        if (oo==1) {                #=0 (not observed)
          n1=n[i,j]-y[i,j] 
          tmp=theta[i,]*(1-phi[,j])
        }
        if (oo==2) {                #=1 (observed)
          n1=y[i,j]         
          tmp=theta[i,]*phi[,j]
        }
        prob1=tmp/sum(tmp)
        array.LSKP[i,j,,oo]=rmultinom(1,size=n1,prob=prob1)
      }
    }
  }
  # k=apply(array.LSKP,c(1,3),sum)
  # k1=k/apply(k,1,sum)
  # plot(theta.init,k1)
  # plot(theta,k1)
  
  #to store outcomes from gibbs sampler
  theta.out=matrix(NA,ngibbs,ncomm*nloc)
  llk=rep(NA,ngibbs)
  
  #run gibbs sampler
  options(warn=2)
  zeroes=array(0,dim=c(nloc,nspp,ncomm))
  for (i in 1:ngibbs){
    print(i)   

    #sample z when y=0
    tmp0=samplez0(theta=theta, OneMinusPhi=1-phi, 
                  NminusY=NminusY, ncommun=ncomm, nloc=nloc, nspp=nspp,
                  zeroes=zeroes)
    array.LSKP[,,,1]=tmp0$ArrayLSK
    nlk0=tmp0$nlk
    # k=apply(array.LSKP[,,,1],c(1,3),sum)
    # unique(nlk0-k)
    
    #sample z when y=1
    tmp1=samplez1(theta=theta, phi=phi, 
                  y=y, ncommun=ncomm, nloc=nloc, nspp=nspp,
                  zeroes=zeroes)
    array.LSKP[,,,2]=tmp1$ArrayLSK
    nlk1=tmp1$nlk
    # k=apply(array.LSKP[,,,2],c(1,3),sum)
    # unique(nlk1-k)
    
    #get parameters  
    theta=get.theta(nlk=nlk0+nlk1,gamma,ncomm,nloc) #theta.true#
    # theta[theta>hi]=hi; theta[theta<lo]=lo
    
    ooo=ooo+1;
    if (ooo>npost) ooo=1
    phi=matrix(phi.post[ooo,],ncomm,nspp)
    
    prob=theta%*%phi
    prob[prob>hi]=hi; prob[prob<lo]=lo

    #calculate logl and store results  
    llk[i]=sum(dbinom(y,size=n,prob=prob,log=T))
    theta.out[i,]=theta
    
    #re-order groups
    # if (i%%50==0 & i<nburn){
    #   med=apply(theta,2,mean)
    #   ordem=order(med,decreasing=T)
    #   theta=theta[,ordem]
    #   phi=phi[ordem,]
    #   array.LSKP[,,,1]=array.LSKP[,,ordem,1]
    #   array.LSKP[,,,2]=array.LSKP[,,ordem,2]
    # }
  }
  seq1=nburn:ngibbs
  list(llk=llk,theta=theta.out[seq1,])  
}