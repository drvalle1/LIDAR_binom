rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(166)

#get functions
setwd('U:\\GIT_models\\LIDAR_binom')
source('LIDAR_binom aux functions.R')
source('LIDAR_binom main function.R')
sourceCpp('aux1.cpp')

#get data
setwd('U:\\GIT_models\\LIDAR_binom\\fake data')
dat=read.csv('fake data y7.csv',as.is=T)
y=data.matrix(dat)

dat=read.csv('fake data n7.csv',as.is=T)
n=data.matrix(dat)

#useful stuff
ncomm=10
ngibbs=10000
nburn=ngibbs/2

#priors
a.phi=0.1
b.phi=0.1
gamma=0.1

#run gibbs
mod=LIDAR_binom(y=y,n=n,ncomm=ncomm,a.phi=a.phi,b.phi=b.phi,
                gamma=gamma,ngibbs=ngibbs,nburn=nburn)

plot(mod$llk,type='l')
seq1=nburn:ngibbs
plot(mod$llk[seq1],type='l')

# seq1=(2500-nburn):(ngibbs-nburn)
tmp=apply(mod$theta,2,mean)
nloc=nrow(y)
theta.estim=matrix(tmp,nloc,ncomm)
boxplot(theta.estim)

plot(NA,NA,xlim=c(0,nloc),ylim=c(0,1))
for (i in 1:ncomm) lines(1:nloc,theta.estim[,i],col=i)

#export phi to test fold in operation
setwd('U:\\GIT_models\\LIDAR_binom\\theta given phi')
write.csv(mod$phi,'phi post.csv',row.names=F)