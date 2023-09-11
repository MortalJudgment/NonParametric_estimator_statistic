library(MASS)
data(galaxies)
help(galaxies)
hist(galaxies,fre=F,nclass = 100)
plot(density(galaxies,kemel = 'n',bw = 500))
###---------------- Mix of 2 density
### n: the size of the sample
### alpha: the poportion between two model
### l0, l1: which kind of distribution
### p0, p1: parameter
###---------------- f_X(x) = alpha*f_0(x,p0)

rmelange.m=function(n,alpha,l0,l1,p0,p1)
{
  z=rbinom(n,1,alpha)
  f1=eval(parse(text=paste('r',l1,'(',paste(c(n,p1),collapse=','),')',sep='')))
  f0=eval(parse(text=paste('r',l0,'(',paste(c(n,p0),collapse=','),')',sep='')))
  x=z*f1+(1-z)*f0
  return(x=x)
}
###-----------------
dmelange=function(t,alpha,l0,l1,p0,p1)
{
  res=alpha*eval(parse(text=paste('d',l1,'(t,',paste(p1,collapse=','),')',sep='')))+(1-alpha)*eval(parse(text=paste('d',l0,'(t,',paste(p0,collapse=','),')',sep='')))
}


n=300
alpha=0.3
l0='norm'
p0=c(8,1)
l1='norm'
p1=c(0,2)
type = 'r'
s=seq(-10,10,0.001)
x=rmelange.m(n,alpha,l0,l1,p0,p1)
par(mfrow=c(2,3))
plot(density(x,bw=0.001,kernel=type))
lines(s,dmelange(s,alpha,l0,l1,p0,p1),col='red')
plot(density(x,bw=0.01,kernel=type))
lines(s,dmelange(s,alpha,l0,l1,p0,p1),col='red')
plot(density(x,bw=0.1,kernel=type))
lines(s,dmelange(s,alpha,l0,l1,p0,p1),col='red')
plot(density(x,kernel=type))
lines(s,dmelange(s,alpha,l0,l1,p0,p1),col='red')
plot(density(x,bw=1,kernel=type))
lines(s,dmelange(s,alpha,l0,l1,p0,p1),col='red')
plot(density(x,kernel=type,bw=10))
lines(s,dmelange(s,alpha,l0,l1,p0,p1),col='red')
######################
par(mfrow=c(1,1))
res=density(x,kernel = 'r')
plot(res,main='Compare different type of kernel')
lines(density(x,kernel = 'g',bw = res$bw),col='blue')
lines(density(x,kernel = 'e',bw = res$bw),col='green')
lines(s,dmelange(s,alpha,l0,l1,p0,p1),col='red')
##########
x = galaxies
par(mfrow=c(1,1))
res=density(x,kernel = 'rectangular')
plot(res,main='Compare different type of kernel')
lines(density(x,kernel = 'gaussian',bw = res$bw),col='blue')
lines(density(x,kernel = 'epanechnikov',bw = res$bw),col='green')
legend("topright", title="Type of the Kernel", legend=c("Data","Gaussian", "Epanechnikov"),col=c("black","blue","green"), lty=1,cex=0.8)
################
library(MASS)
data(galaxies)
n=2000
alpha=0.3
l0='norm'
p0=c(8,1)
l1='norm'
p1=c(0,2)
s=seq(-10,10,0.001)

x=rmelange.m(n,alpha,l0,l1,p0,p1)
plot(density(x,bw = 'nrd0'),main = 'Different type of Bandwidth')
lines(density(x,bw='nrd'),col='blue')
lines(density(x,bw='SJ'),col='green')
lines(density(x,bw='ucv'),col='yellow')
lines(s,dmelange(s,alpha,l0,l1,p0,p1),col='red')
legend("topleft", title="Type of the Bandwidth", legend=c("Gaussian","common variation","Sheather&Jones","unbiased cross-validation"),col=c("black","blue","green","yellow"), lty=1,cex=0.8)

######## Quadratic loss
library('sm')
# number of simulations
J=30
hs=(1:20)/20
s=seq(-7,13,0.01) 

h0=5
s0=1700


QUAD_LOSS=function(s,hs,J,n,alpha,l0,l1,p0,p1,h0,s0)
{
  ls= length(s)
  lh=length(hs)  
  EST=array(NA,c(J,ls,lh))
  for (j in 1:J)
  {
    x=rmelange.m(n,alpha,l0,l1,p0,p1)
    for (h in 1:lh) 
    {
      EST[j,,h]=sm.density(x,h=hs[h],display='none',eval.points=s)$estimate
    }
  }
  BIAS=apply(EST,c(2,3),mean)-dmelange(s,alpha,l0,l1,p0,p1)
  VAR=apply(EST,c(2,3),var)
  EQ=BIAS^2+VAR
  
  nl=2
  #if (!is.null(h0)) nl=nl+1
  if (!is.null(s0)) nl=nl+1
  
  layout(matrix(c(1:3,rep(4,3),5:(3*nl+1)),byrow=TRUE, ncol=3))
  plot(hs,abs(apply(BIAS,2,mean)),type='l',ylab='|BIAS|')
  plot(hs,apply(VAR,2,mean),type='l',ylab='VAR')
  plot(hs,apply(EQ,2,mean),type='l',ylab='EQ')
  abline(v=density(x,bw='nrd')$bw,col='blue')
  abline(v=density(x,bw='sj')$bw,col='green')
  abline(v=density(x,bw='ucv')$bw,col='red')
  
  
  hopt=which(apply(EQ,2,mean)==min(apply(EQ,2,mean)))
  if (is.null(h0)) h0=hopt
  EST2=EST[,,h0]#apply(EST,2,as.vector)
  
  plot(s,EST2[1,],type='l',ylab='Estimates',main=paste('h=',hs[h0],sep=''))
  for (j in 1:J) lines(s,EST2[j,])
  lines(s,dmelange(s,alpha,l0,l1,p0,p1),col='red')
  
  
  #plot(s,apply(BIAS,1,mean),type='l',ylab='BIAS')
  #plot(s,apply(VAR,1,mean),type='l',ylab='VAR')
  #plot(s,apply(EQ,1,mean),type='l',ylab='EQ')
  
  
  if (!is.null(h0)) 
  {
    plot(s,abs(BIAS[,h0]),type='l',ylab='|BIAS|',main=paste('h=',hs[h0],sep=''))
    plot(s,VAR[,h0],type='l',ylab='VAR',main=paste('h=',hs[h0],sep=''))
    plot(s,EQ[,h0],type='l',ylab='EQ',main=paste('h=',hs[h0],sep=''))
  }
  if (!is.null(s0))
  {
    plot(hs,abs(BIAS[s0,]),type='l',ylab='|BIAS|',main=paste('s=',s[s0],sep=''))
    plot(hs,VAR[s0,],type='l',ylab='VAR',main=paste('s=',s[s0],sep=''))
    plot(hs,EQ[s0,],type='l',ylab='EQ',main=paste('s=',s[s0],sep=''))
  }
  return(list(BIAS=BIAS,VAR=VAR,EQ=EQ,hopt=hs[hopt]))
}

RES=QUAD_LOSS(s,hs,J,n,alpha,l0,l1,p0,p1,NULL,NULL)


#####
de.mode=function(x,a,b,M,bw='ucv',kernel='e')
{
  disc=seq(a,b,length.out=M)
  dens=density(x,from=a,to=b,n=M,bw=bw,kernel=kernel)$y
  mod=disc[(order(dens))[M]]
  max=max(dens)
  plot(disc,dens,type='l')
  return(list(mode=mod,max=max))
}

sm.mode=function(x,a,b,M,method='cv')
{
  disc=seq(a,b,length.out=M)
  dens=sm.density(x,eval.points=disc,method=method)$estimate
  mod=disc[(order(dens))[M]]  
  max=max(dens)
  #min=min(dens)
  return(list(mode=mod,max=max))
}

library(sm)
par(mfrow=c(1,1))
sm.density(x,method='normal',col='red')
sm.density(x,method='sj',add=T,col='blue')
sm.density(x,method='cv',add=T,col='green')
sm.mode(x,-10,15,1000,method='cv')
