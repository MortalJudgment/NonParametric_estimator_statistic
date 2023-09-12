#Question 1: Treasure Hunt
PT=load('PIRATES.RData')
#1
PT
# t, x, y are real values data while p is integer valued
#sample size : 1000
length(y)
#When do the readings stop
t[length(t)]
#Since final value of t : measurement instants is 12.562903917
#2

#3
m1=min(t)
M1=max(t)
m=min(p)
M=max(p)
par(mfrow=c(2,2))
plot(density(t,kernel="e",bw=0.005,from =m1, to = M1))
plot(density(t,kernel="e",bw=0.05,from =m1, to = M1))
plot(density(t,kernel="e",bw=5,from =m1, to = M1))
plot(density(t,kernel="e",bw="ucv",from =m1, to = M1))

par(mfrow=c(2,2))
plot(density(p,kernel="e",bw=0.07,from =m, to = M))
plot(density(p,kernel="e",bw=0.7,from =m, to = M))
plot(density(p,kernel="e",bw=7,from =m, to = M))
plot(density(p,kernel="e",bw="ucv",from =m, to = M))

par(mfrow=c(1,1))
#6
plot(t, p, type = "p", xlab = "Time", ylab = "Depth")
#7
library(sm)
sm.regression(t,p,method = 'cv',eval.points = t,add = T, col = 'blue',xlab = "Time", ylab = "Depth",lwd=3)
sm.regression(t,p,method = 'df',eval.points = t,add = T, col = 'red',xlab = "Time", ylab = "Depth",lwd=3)
sm.regression(t,p,method = 'aicc',eval.points = t,add = T, col = 'yellow',xlab = "Time", ylab = "Depth",lwd=3)
# method cross-validation is best in three of them.
## Reason: Method df is obviously fits bad. 
## Even though aicc fits well, some regions are underestimate. Only cv fits good
#8
res = sm.regression(t, p, method = "cv", eval.points = c(3.228297), display = "none")
res$estimate
#1.1 PART A
#1
plot(t,y,type="p",xlab = "Time", ylab = "Longitude")
#Comment:Visible pattern/relationship between the longitude and time,
## i.e, longitude can be expressed as a function of time.
#2
y_e=sm.regression(t,y,method='cv',eval.points=t,col='blue',display="none")
y_e
##Bonus
ln_res = lm(y ~ sin(t)) # y ~ a + b*sin(t) LINEAR REGRESSION METHOD
intercept = ln_res$coefficients[1]
coeff = ln_res$coefficients[2]
y_e_regress = intercept + coeff * sin(t)

tplot = seq(min(t), max(t), by = 0.01)
plot(t, y, type = "p", xlab = "time", ylab = "Longitude")
lines(tplot, intercept + coeff * sin(tplot), lwd = 3, col = "red")
#1.2 PART B
#1
plot(t,x,type="p",xlab = "Time", ylab = "Latitude")
#Comment:Visible pattern/relationship between the longitude and time,
## i.e, longitude can be expressed as a function of time.
#2
x_e=sm.regression(t,x,method='cv',eval.points=t,col='blue',display="none")
x_e
#Bonus
ln_res1=lm(x~cos(t)+cos(t/2))
ln_res1
intercept=ln_res1$coefficients[1]
intercept
coeff_1=ln_res1$coefficients[2]
coeff_1
coeff_2=ln_res1$coefficients[3]
x_e_regress=intercept+coeff_1*cos(t)+coeff_2*cos(t/2)

tplot_1 = seq(min(t), max(t), by = 0.01)
plot(t, x, type = "p", cex = 0.75, xlab = "time", ylab = "Latitude")
lines(tplot, intercept + coeff_1 * cos(tplot_1) + coeff_2 * cos(tplot_1/2), 
      lwd = 3, col = "red")
#Part C
#1
plot(x,y,type="p",xlab = "Latitude", ylab = "Longitude")
#2
plot(x, y, type = "p",  xlab = "Latitude", ylab = "Longitude")
lines(x_e$estimate, y_e$estimate, col = "yellow", lwd = 3)
#BONUS
plot(x, y, type = "p", cex = 0.5, xlab = "Latitude", ylab = "Longitude")
lines(x_e_regress, y_e_regress, col = "green", lwd = 3)
#BONUS
xs = seq(from = -2, to = 4, by = 0.03)
ys = seq(from = -1, to = 1, by = 0.01)
xy_eval = cbind(xs, ys)
p_e = sm.regression(cbind(x, y), p, eval.points = xy_eval, eval.grid = F, structure.2d = "common")

##2
Data=load('PHONEMES.RData')
Data
#1
typeof(PHONEME)
typeof(CURVES)
#2
#sample size n=2000
length(PHONEME)
table(PHONEME)#How many observations do you have for the different phonemes?
#3
train_sample=CURVES[1:1000,]
train_sample
test_sample=CURVES[1001:length(PHONEME),]
test_sample
#5
t=matrix(CURVES)
X=log(t)
X
Y_SH=ifelse(PHONEME[1:1000]=="SH",1,0)
Y_SH
Y_IY=ifelse(PHONEME[1:1000]=="IY",1,0)
Y_IY
Y_DCL=ifelse(PHONEME[1:1000]=="DCL",1,0)
Y_DCL
Y_AA=ifelse(PHONEME[1:1000]=="AA",1,0)
Y_AA
Y_AO=ifelse(PHONEME[1:1000]=="AO",1,0)
Y_AO
