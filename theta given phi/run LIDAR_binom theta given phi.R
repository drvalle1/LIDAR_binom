rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(59)

#get functions
setwd('U:\\GIT_models\\LIDAR_binom')
source('LIDAR_binom aux functions.R')
sourceCpp('aux1.cpp')

setwd('U:\\GIT_models\\LIDAR_binom\\theta given phi')
source('LIDAR_binom theta given phi.R')

#get data
setwd('U:\\GIT_models\\LIDAR_binom\\fake data')
dat=read.csv('fake data y7.csv',as.is=T)
y=data.matrix(dat)

dat=read.csv('fake data n7.csv',as.is=T)
n=data.matrix(dat)

#get phi
setwd('U:\\GIT_models\\LIDAR_binom\\theta given phi')
phi.post=data.matrix(read.csv('phi post.csv'))

#useful stuff
ncomm=10
ngibbs=3000
nburn=ngibbs/2

#priors
a.phi=1
b.phi=1
gamma=0.1

#run gibbs
mod=LIDAR_binom_theta_given_phi(y=y,n=n,ncomm=ncomm,
                gamma=gamma,ngibbs=ngibbs,nburn=nburn,
                phi.post=phi.post)

plot(mod$llk,type='l')
seq1=500:ngibbs
plot(mod$llk[seq1],type='l')

tmp=apply(mod$theta,2,mean)
nloc=nrow(y)
theta.estim=matrix(tmp,nloc,ncomm)
boxplot(theta.estim)

plot(NA,NA,xlim=c(0,nloc),ylim=c(0,1))
for (i in 1:ncomm) lines(1:nloc,theta.estim[,i],col=i)