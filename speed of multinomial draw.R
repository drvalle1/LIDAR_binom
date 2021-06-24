rm(list=ls(all=TRUE))
library('microbenchmark')
library('Rcpp')
set.seed(59)

#get functions
setwd('U:\\GIT_models\\LIDAR_binom')
sourceCpp('aux1.cpp')
sourceCpp('teste.cpp')

#basic settings
nobs=100000
ngrp=10
tmp=runif(ngrp)
prob=tmp/sum(tmp)
niter=400

#my customized function
f=function(){
  for (i in 1:niter){
    k=rmultinom1(runif1=runif(nobs), prob=prob, ncommun=ngrp)
  }
}
f()
res1 <- microbenchmark(NULL, f(), times=10)
res1

#built in function
f=function(){
  for (i in 1:niter){
    k=rmultinom_1(probs=prob, size=nobs)
  }
}
f()
res2 <- microbenchmark(NULL, f(), times=10)
res2

