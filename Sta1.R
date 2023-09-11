# Addition
2+3

# Seqence
7:104
21:1
seq(0,1,0.1)
seq(0,1,length=50)
seq(0,1,length.out=50)

# Vector
x = c(-pi,021,sqrt(2),-10)
y = c('A','B','C')
z = c(T,T,F,T,F)
mode(z)
length(x)
c(x,y)

#repeat
rep(x,3)
rep(x,c(10,2,3,7))
#elements of vector
x[3]
x[-3]
z[c(2,4,5)]
# all the elements that positive
x[x>0]

Names=c('A','B','C')
Scores=c(20,7,15)
Names[Scores>12]

sort(x)
sort(x,dec=T)
order(x)

data(trees)
attach(trees)
plot(Girth,Volume)
# hist(Girth)
hist(Girth,fre=F)
x = Girth
moda=sort(unique(x))
moda
e = table(x)
hist(Girth,fre=F,breaks=c(8,9,10,11,12,13,14,15,16,17,18,19,20,21,22))
hist(Girth,fre=F,nclass=10)

###
#Boy/Girl
#     B   |   G
#   ____________
#     10  |   10 
mode=c('B','G')
e=c(10,10)
x = rep(mode,e)
f = e/sum(e)
pie(f)
barplot(f)

###
#Do you like statistic?
#     0  |   1   |   2
#     ___________________
#     1  |   3   |   16

e= c(1,3,16)
mode=c(0,1,2)
x=rep(mode,e)
boxplot(x)
f = e/sum(e)
F = cumsum(f)

### 
# Number of brothers/sisters
#     0  |   1   |   2    |   3   |   4
# ________________________________________
#     6  |   9   |   1    |   3   |   1
# numerical discrete data

moda = c(0,1,2,3,4)
e = c(6,9,1,3,1)
x=rep(moda,e)
mean(x)
median(x)
summary(x)
quantile(x,0.1)
f= e/sum(e)
plot(moda,f,type='h')
F = cumsum(f)
tau=1
plot(c(min(x)-tau,moda,max(x)+tau),c(0,F,1),type='s')

var(x)
sd(x)
boxplot(x)

####
x = rnorm(100)
moda = sort(unique(x))
e=table(x)
f = e/sum(e)
plot(moda,f,type='h')
summary(x)
F = FALSE
hist(x,fre=F,nclass=10)
#### X ~ N(mu,sigma) 
# mu-estimate = bar(X)
# sigma-estimate = sqrt(bar(X^2) - bar(X)^2)

y = seq(-10,10,0.01)
lines(y,dnorm(y,mean(x),sd(x)),col='blue')

### Galaxies data
library(MASS)
data(galaxies)
help(galaxies)
hist(galaxies,fre=F,nclass = 100)
y = seq(10000,35000,1)
lines(y,dnorm(y,mean(galaxies),sd(galaxies)),col='blue')

### Pressure data
data(pressure)
p = pressure$pressure
t = pressure$temperature
plot(t,p)
