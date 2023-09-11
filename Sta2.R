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

### Doping(dopage) data
load("Dopage.RData")

### Tecator
Tec = read.table("Tecator.dat")
X = Tec[,1:100]
Y = Tec[,101]
wl= seq(750,948,2)
plot(wl, X[1,], type='l', ylim=range(X), xlab = 'wavelength', ylab='obsorbtion')
for (i in 2:213)
{
  lines(wl,X[i,])
}

### Phonemes
load("PHONEMES.Rdata")
par(mfrow = c(2,3))
plot(CURVES[1,],type='l', ylim=range(CURVES) ,main = 'SH',ylab = 'logarirm')
for (i in 2:200)
{
  lines(CURVES[i,])
}
plot(CURVES[201,],type='l', ylim=range(CURVES) ,main = 'IY',ylab = 'logarirm')
for (i in 202:400)
{
  lines(CURVES[i,])
}
plot(CURVES[401,],type='l', ylim=range(CURVES) ,main = 'DCL',ylab = 'logarirm')
for (i in 402:600)
{
  lines(CURVES[i,])
}
plot(CURVES[601,],type='l', ylim=range(CURVES) ,main = 'AA',ylab = 'logarirm')
for (i in 602:800)
{
  lines(CURVES[i,])
}
plot(CURVES[801,],type='l', ylim=range(CURVES) ,main = 'AO',ylab = 'logarirm')
for (i in 802:1000)
{
  lines(CURVES[i,])
}
