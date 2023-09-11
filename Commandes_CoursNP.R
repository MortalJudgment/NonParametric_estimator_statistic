######## Melange de densites
rmelange.m=function(n,alpha,l0,l1,p0,p1)
{
  z=rbinom(n,1,alpha)
  f1=eval(parse(text=paste('r',l1,'(',paste(c(n,p1),collapse=','),')',sep='')))
  f0=eval(parse(text=paste('r',l0,'(',paste(c(n,p0),collapse=','),')',sep='')))
  x=z*f1+(1-z)*f0
  return(x=x)
}

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
s=seq(-10,10,0.001)

x=rmelange.m(n,alpha,l0,l1,p0,p1)
par(mfrow=c(2,3))
plot(density(x,bw=0.001,kernel='rectangular'))
lines(s,dmelange(s,alpha,l0,l1,p0,p1),col='red')
plot(density(x,bw=0.01,kernel='rectangular'))
lines(s,dmelange(s,alpha,l0,l1,p0,p1),col='red')
plot(density(x,bw=0.1,kernel='rectangular'))
lines(s,dmelange(s,alpha,l0,l1,p0,p1),col='red')
plot(density(x,kernel='rectangular'))
lines(s,dmelange(s,alpha,l0,l1,p0,p1),col='red')
plot(density(x,bw=1,kernel='rectangular'))
lines(s,dmelange(s,alpha,l0,l1,p0,p1),col='red')
plot(density(x,kernel='rectangular',bw=10))
lines(s,dmelange(s,alpha,l0,l1,p0,p1),col='red')




x=rmelange.m(n,alpha,l0,l1,p0,p1)
plot(density(x,bw='nrd',kernel='g'))
lines(s,dmelange(s,alpha,l0,l1,p0,p1),col='red')
lines(density(x,kernel='g'),col='blue')
lines(density(x,bw='SJ',kernel='g'),col='green')





h=4
x=rmelange.m(n,alpha,l0,l1,p0,p1)
par(mfrow=c(2,2))
res=density(x,bw=h,kernel='rectangular')
plot(res)
lines(s,dmelange(s,alpha,l0,l1,p0,p1),col='red')
plot(density(x,bw=h,kernel='e'))
lines(s,dmelange(s,alpha,l0,l1,p0,p1),col='red')
plot(density(x,bw=h,kernel='g'))
lines(s,dmelange(s,alpha,l0,l1,p0,p1),col='red')
plot(density(x,bw=h,kernel='b'))
lines(s,dmelange(s,alpha,l0,l1,p0,p1),col='red')


par(mfrow=c(1,1))
x=rmelange.m(n,alpha,l0,l1,p0,p1)
res=density(x,bw='ucv',kernel='e')
plot(res)
lines(s,dmelange(s,alpha,l0,l1,p0,p1),col='red')
N=100
for (i in 1:100)
{
  x=rmelange.m(n,alpha,l0,l1,p0,p1)
  lines(density(x,bw=res$bw,kernel='e'))
}


x=rmelange.m(n,alpha,l0,l1,p0,p1)
res=density(x,bw='ucv',kernel='e',from=-5,to=15,n=10000)
RES=res$y/100
for (i in 1:99)
{
  x=rmelange.m(n,alpha,l0,l1,p0,p1)
  res=density(x,bw='ucv',kernel='e',from=-5,to=15,n=10000)
  RES=RES+res$y/100
}
z=seq(-5,15,len=10000)
plot(z,RES,type='l')
lines(z,dmelange(z,alpha,l0,l1,p0,p1),col='red')




######## Quadratic loss
library('sm')
# number of simulations
J=30
hs=(1:20)/20
s=seq(-10,10,0.01)

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
  
######## Mode estimation

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
  return(list(mode=mod,max=max))
}
#######
Compute_CV=function(x,y,method,weights=NA)
{
  hopt=sm.regression(x,y,method=method,display='none')$h
  n=length(x)
  if (is.na(weights)) weights=rep(1,n)
  CV=0
  for (i in 1:n)
  {
    CV=CV+(y[i]-sm.regression(x[-i],y[-i],h=hopt,eval.points=x[i],display='none')$estimate)^2*weights[i] 
  }
  return(CV)
}
####### ROC 
ROC=function(z,p1)
{
  p=sort(unique(p1))
  lp=length(p)
  ROC=matrix(NA,lp,2)
  for (s in 1:lp)
  {
    ROC[s,1]=mean(p1[z==0]>p[s])
    ROC[s,2]=mean(p1[z==1]>p[s])  
  }
  AUC=sum(ROC[-lp,2]*(ROC[-lp,1]-ROC[-1,1]))
  
  colnames(ROC)=c("FPR","TPR")
  plot(ROC,type='l',main=paste('AUC =', AUC))
  return(list(ROC=ROC,AUC=AUC))
  }

#######
Classif_NP=function(X,Y,X0=NA)
{
  X=as.matrix(X)
  n=length(Y)
  V=sort(unique(Y))
  n_V=length(V)
  Prob=matrix(NA,n,n_V)
  colnames(Prob)=V
  Class=rep(NA,n)
  if (!is.na(max(X0)))
  {
    X0=as.matrix(X0)
    P0=matrix(NA,nrow(X0),n_V)
    Class0=rep(NA,n)
  }
  for (v in 1:n_V)
  {
    z=as.numeric(Y==V[v])
    for (i in 1:n )
    {
      Prob[i,v]=sm.regression(X[-i,],z[-i], eval.points=X,eval.grid=FALSE,display='none')$estimate[i]
    }
    if (!is.na(max(X0))) {P0[,v]=sm.regression(X,z,eval.points=X0,eval.grid=FALSE,display='none')$estimate}
  }
  if (n_V==2) {ROC(Y==V[2],Prob[,2])}
  Class=V[apply(Prob,1,which.max)]
  if (!is.na(max(X0))) {Class0=V[apply(P0,1,which.max)]}
  if (!is.na(max(X0))) {return(list(Class=Class, Prob=Prob, M_table=table(Y,Class), Class0=Class0,Prob0=P0))}
  else {return(list(Class=Class, Prob=Prob, M_table=table(Y,Class)))}
}