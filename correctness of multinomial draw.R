rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(59)

#get functions
setwd('U:\\GIT_models\\LIDAR_binom')
sourceCpp('teste.cpp')

#basic settings
nobs=100000
ngrp=10
niter=100

#test to see if it works
res=numeric()
for (i in 1:niter){
  tmp=runif(ngrp)
  prob=tmp/sum(tmp)
  k=rmultinom_1(probs=prob, size=nobs)
  tmp=cbind(prob,k/sum(k))
  res=rbind(res,tmp)
}

#compare
plot(res[,2],res[,1])
hist(res[,2]-res[,1])
