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
  
  array.LSKP=array(NA,dim=c(nloc,nspp,ncomm,2))
  prob1=rep(1/ncomm,ncomm)
  for (i in 1:nloc){
    for (j in 1:nspp){
      for (oo in 1:2){
        if (oo==1) n1=n[i,j]-y[i,j]
        if (oo==2) n1=y[i,j]
        array.LSKP[i,j,,oo]=rmultinom(1,size=n1,prob=prob1)
      }
    }
  }

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

    #sample z when y=1
    tmp1=samplez1(theta=theta, phi=phi, 
                  y=y, ncommun=ncomm, nloc=nloc, nspp=nspp,
                  zeroes=zeroes)
    array.LSKP[,,,2]=tmp1$ArrayLSK
    nlk1=tmp1$nlk

    #get parameters  
    theta=get.theta(nlk=nlk0+nlk1,gamma,ncomm,nloc) #theta.true#
    theta[theta>hi]=hi; theta[theta<lo]=lo
    
    ooo=ooo+1;
    if (ooo>npost) ooo=1
    phi=matrix(phi.post[ooo,],ncomm,nspp)
    
    prob=theta%*%phi
    prob[prob>hi]=hi; prob[prob<lo]=lo

    #calculate logl and store results  
    llk[i]=sum(dbinom(y,size=n,prob=prob,log=T))
    theta.out[i,]=theta
  }
  seq1=nburn:ngibbs
  list(llk=llk,theta=theta.out[seq1,])  
}