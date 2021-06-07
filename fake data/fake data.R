rm(list=ls(all=TRUE))
set.seed(421)

nloc=2000
nspp=50
ncommun=7
base=floor(nloc/(ncommun-2))

#generate thetas
x=seq(from=-1,to=1,length.out=base)
y=sqrt(1-(x^2))*0.1
min1=0.0001
y[y<min1]=min1
# plot(x,y)
  
init=floor(nloc/ncommun)
seq1=c(seq(from=1,to=nloc,by=init),nloc)
  
theta=matrix(min1,nloc,ncommun)
for (i in 1:ncommun){
  seq2=seq1[i]:(seq1[i]+base-1)
  seq3=seq2[seq2<=nloc]
  theta[seq3,i]=y[1:length(seq3)]
}
theta=theta/matrix(apply(theta,1,sum),nloc,ncommun)
theta.true=theta
  
plot(NA,NA,xlim=c(0,nloc),ylim=c(0,1))
for (i in 1:ncommun) lines(1:nloc,theta[,i],col=i)
  
#generate phi  
tmp=matrix(rbeta(ncommun*nspp,0.5,0.5),ncommun,nspp)
tmp[,1:ncommun]=diag(0.8,ncommun)
phi=tmp
round(phi[,1:10],2)
hist(phi)
phi.true=phi

#get n
setwd('U:\\independent studies\\LIDAR\\experim\\TAN y_n\\2012')
n2=matrix(round(runif(nloc*nspp,min=10,max=200)),nloc,nspp)

#generate actual observations y
y=matrix(0,nloc,nspp)
array1=array(NA,dim=c(nloc,nspp,ncommun,2))
for (i in 1:nloc){
  for (j in 1:nspp){
    tmp=rmultinom(1,size=n2[i,j],prob=theta[i,])  
    for (k in 1:ncommun){
      tmp1=rbinom(1,size=tmp[k],prob=phi[k,j])
      array1[i,j,k,1]=tmp1
      array1[i,j,k,2]=tmp[k]-tmp1
      y[i,j]=y[i,j]+tmp1
    }
  }
}
mean(y<=n2)
image(data.matrix(y/n2))

#make sure things make sense
teste=apply(array1[,,,1],c(1,2),sum)
unique(y-teste)
teste=apply(array1,c(1,2),sum)
unique(n2-teste)

#make nice output table
colnames(y)=paste('spp',1:nspp,sep='')

#export results
setwd('U:\\GIT_models\\LIDAR_binom\\fake data')
nome=paste('fake data y',ncommun,'.csv',sep='')    
write.csv(y,nome,row.names=F)    
nome=paste('fake data n',ncommun,'.csv',sep='')    
write.csv(n2,nome,row.names=F)    