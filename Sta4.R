data(pressure)
p = pressure$pressure
t = pressure$temperature
plot(t,p)
library(sm)
sm.regression(t,p,col='red',method='cv')


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
  
  colnames(ROC)=c("FPR","TPR") #True Positive Rate, False Positive Rate
  plot(ROC,type='l',main=paste('AUC =', AUC)) #Area under the curve
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
### Similated data
n = 1000
x = rexp(n,2)
y = 10*exp(-3*x) + rnorm(n,0,exp(-x/2))
plot(function(x) 10*exp(-3*x),0,3,lwd = 2)
lines(x,y,col='red',type='p')
sm.regression(x,y,method = 'cv',add= T,col='blue',lwd = 2)
sm.regression(x,y,method = 'aicc',add= T,col='green',lwd = 2)
sm.regression(x,y,method = 'df',add= T,col='yellow',lwd = 2)



e = rnorm(n,0,1)
x = rnorm(n,2,0.5)
y = cos(2*pi*x) + 30*exp(-3*x) + 0.1*e
#plot(function(x) exp(-3*x),0,3,lwd = 2)
plot(x,y,col='red',type='p')
res = sm.regression(x,y,method = 'cv',add= T,col='blue',lwd = 2,eval.points=sort(x))
sm.regression(x,y,method = 'aicc',add= T,col='green',lwd = 2)
sm.regression(x,y,method = 'df',add= T,col='yellow',lwd = 2)
names(res)
## stardard error
cv=0 # MSE
for (i in 1:n)
  {
    res1 = sm.regression(x[-i],y[-i],method = 'cv',add= T,display='none',col='black',lwd = 2,eval.points=x[i])
    cv = cv + 1/n*(y[i] - res1$estimate)^2
    print(cv)
}

res = sm.regression(x,y,method = 'cv',col='blue',lwd = 2,eval.points=sort(x))
plot(res$estimate, y[order(x)],col='red')
abline(0,1)

load("Dopage.RData")
Classif_NP(hema,test)

Classif_NP(hema,test,c(51,27,72))

data("faithful")
x = attach(faithful)
sm.density(faithful)
e0 = seq(min(eruptions),max(eruptions),len=10)
w0 = seq(min(waiting),max(waiting),len=10)
a = sm.density(cbind(eruptions,waiting),eval.points = cbind(e0,w0),eval.grid = T,structure = 'common')

### Prite data
load("PRITE.RDATA")
sm.regression(cbind(x,y),eval.points = cbind(x,y),p,structure = 'common', method = 'cv',eval.grid = T)


### Tecator data
Tec = read.table("E:/KHTN/a_PUF 2022/Pro_Desol/Data_in_R/Tecator.dat")
X = Tec[,1:100]
Y = Tec[,101]
par(mfrow=c(2,2))
X = as.matrix(X)
R = funopare.kernel.cv(Y,X,X,semimetric='deriv',0,5,c(0,1))
names(R)
plot(R$Estimated.values,Y,col='red')
abline(0,1)
R = funopare.kernel.cv(Y,X,X,semimetric='deriv',1,5,c(0,1))
plot(R$Estimated.values,Y,col='red')
abline(0,1)
R = funopare.kernel.cv(Y,X,X,semimetric='deriv',2,5,c(0,1))
plot(R$Estimated.values,Y,col='red')
abline(0,1)
R = funopare.kernel.cv(Y,X,X,semimetric='deriv',3,5,c(0,1))
plot(R$Estimated.values,Y,col='red')
abline(0,1)

###
par(mfrow=c(1,1))
R = funopare.kernel.cv(Y,X,X,semimetric='pca',100)
plot(R$Estimated.values,Y,col='red')
abline(0,1)
