load("E:/KHTN/a_PUF 2022/Pro_Desol/Data_in_R/PIRATES.RData")

## exercise 3
par(mfrow=c(2,3))
plot(density(t,kernel='e',bw = 0.005))
lines(density(t,kernel='e',bw = "ucv"),col = 'red')
plot(density(t,kernel='e',bw = 0.05))
lines(density(t,kernel='e',bw = "ucv"),col = 'red')
plot(density(t,kernel='e',bw = 5))
lines(density(t,kernel='e',bw = "ucv"),col = 'red')

plot(density(p,kernel='e',bw = 0.007))
lines(density(p,kernel='e',bw = "ucv"),col = 'red')
plot(density(p,kernel='e',bw = 0.07))
lines(density(p,kernel='e',bw = "ucv"),col = 'red')
plot(density(p,kernel='e',bw = 7))
lines(density(p,kernel='e',bw = "ucv"),col = 'red')

####
library(sm)
par(mfrow=c(1,1))
plot(t,p,type='p')
sm.regression(t,p,method = 'cv',add= T,col='blue',lwd = 2)
sm.regression(t,p,method = 'aicc',add= T,col='green',lwd = 2)
sm.regression(t,p,method = 'df',add= T,col='yellow',lwd = 2)
## Exercise 8
r=sm.regression(t,p,method = 'cv',add= T,col='blue',lwd = 2,eval.points=3.228297)
abline(v = 3.228297,col="red")
abline(h = r$estimate,col="red")
r$estimate

###
plot(t,y,type='p')
y_e = sm.regression(t,y,method = 'cv',add= T,col='blue',lwd = 2)

par(mfrow=c(1,1))
res = lm(y ~ sin(t))
ts = seq(min(t),max(t), by = 0.01)
plot(t,y,type='p')
w_e = lines(ts, res$coefficients[1] + res$coefficients[2]*sin(ts),col = 'red')

###
#par(mfrow=c(1,2))
plot(t,x,type='p')
x_e = sm.regression(t,x,method = 'cv',add= T,col='blue',lwd = 2)
###
res1 = lm(x ~ cos(t)+cos(t/2))
ts = seq(min(t),max(t), by = 0.01)
plot(t,x,type='p')
z_e = lines(ts, res1$coefficients[1] + res1$coefficients[2]*cos(ts) + res1$coefficients[3]*cos(ts/2),col = 'red')

###
x0=seq(-2,4,by = 0.03)
y0=seq(-1,1,by = 0.01)
sm.density(cbind(x,y), eval.points=cbind(x0,y0), eval.grid = T,structure = 'common')

par(mfrow=c(1,1))
plot(x,y,type = 'p')
lines(x_e$estimate,y_e$estimate,col='blue',lwd=2)
#lines(z_e$estimate,w_e$estimate,col='red',lwd=2)

rx = sm.regression(t,x,method = 'cv',add= T,col='blue',lwd = 2,eval.points = pi)
ry = sm.regression(t,y,method = 'cv',add= T,col='green',lwd = 2,eval.points = pi)

rrx = sm.regression(t,x,method = 'cv',add= T,col='blue',lwd = 2,eval.points = 3*pi)
rry = sm.regression(t,y,method = 'cv',add= T,col='green',lwd = 2,eval.points = 3*pi)

abline(v = rx$estimate,col = 'green')
abline(h = ry$estimate,col = 'green')
abline(v = rrx$estimate,col = 'purple')
abline(h = rry$estimate,col = 'purple')


rrz = sm.regression(x,y,method = 'cv',add= T,col='red',lwd = 2,eval.points = 3*pi)

abline(v = rrz$estimate,col = 'purple')






load("E:/KHTN/a_PUF 2022/Pro_Desol/Data_in_R/PHONEMES.RData")

table(PHONEME)

n=length(PHONEME)
tranning = CURVES[1:1000,]
test = CURVES[1001:n,]

y_training = PHONEME[1:1000]
y_test = PHONEME[1001:n]
y_sh = ifelse(y_training == "SH",1,0)
y_iy = ifelse(y_training == "IY",1,0)
y_dcl = ifelse(y_training == "DCL",1,0)
y_aa = ifelse(y_training == "AA",1,0)
y_ao = ifelse(y_training == "AO",1,0)


funopare.kernel.cv <- function(Response, CURVES, PRED, ..., kind.of.kernel = "quadratic", semimetric = "deriv", h.range = NULL)
{
  ################################################################
  # Performs functional prediction (regression) of a scalar response 
  # from a sample of curves via the functional kernel estimator. 
  # A global bandwidth is automatically selected with a 
  # cross-validation procedure. 
  #    "Response" vector containing the observations of the scalar 
  #               response
  #    "CURVES" matrix containing the curves dataset (row by row) 
  #             used for the estimating stage
  #    "PRED" matrix containing new curves stored row by row
  #           used for computing predictions
  #    "..." arguments needed for the call of the function computing 
  #          the semi-metric between curves
  #    "kind.of.kernel" the kernel function used for computing of 
  #                     the kernel estimator; you can choose 
  #                     "indicator", "triangle" or "quadratic (default)
  #    "semimetric" character string allowing to choose the function 
  #                 computing the semimetric;  you can select 
  #                 "deriv" (default), "fourier", "hshift", "mplsr", 
  #                 and "pca"
  #    "h.range" a vector of length 2 giving the range for the bandwidth. 
  #              By default, the procedure defines a sequence of candidates
  #              of bandwidths according to the values of the matrix CURVES 
  # Returns a list containing:
  #    "Estimated.values" vector containing estimated reponses for 
  #                        each curve of "CURVES"
  #    "Predicted.values" if PRED different from CURVES, this vector 
  #                       contains predicted responses for each 
  #                       curve of PRED
  #    "hopt" value of the optimal bandwidth
  #    "hseq" the used sequence of possible bandwidths
  #    "Mse" mean squared error between estimated values 
  #          and observed values 
  ################################################################
  Response <- as.vector(Response)
  if(is.vector(PRED)) PRED <- as.matrix(t(PRED))
  testfordim <- sum(dim(CURVES)==dim(PRED))==2
  twodatasets <- T
  if(testfordim) twodatasets <- sum(CURVES==PRED)!=prod(dim(CURVES))
  sm <- get(paste("semimetric.", semimetric, sep = ""))
  if(semimetric == "mplsr")
    SEMIMETRIC1 <- sm(Response, CURVES, CURVES, ...)
  else SEMIMETRIC1 <- sm(CURVES, CURVES, ...)
  Semimetric1 <- SEMIMETRIC1[row(SEMIMETRIC1) > col(SEMIMETRIC1)]
  if(is.null(h.range)) {
    h.seq <- quantile(Semimetric1, seq(0.05, 0.5, length = 20))
  }else {
    h.seq <- seq(h.range[1], h.range[2], length = 20)
  }
  kernel <- get(kind.of.kernel)
  count1 <- 0
  count2 <- 0
  Mse <- 0
  h.seq.corrected <- 0
  h.seq.length <- length(h.seq)
  h <- h.seq[1]
  while(count1 < h.seq.length) {
    count1 <- count1 + 1
    h <- 1.1 * h
    KERNEL1 <- kernel(SEMIMETRIC1/h)
    diag(KERNEL1) <- 0
    KERNEL1[KERNEL1 < 0] <- 0
    KERNEL1[KERNEL1 > 1] <- 0
    Denom1 <- apply(KERNEL1, 2, sum)
    Logic <- (Denom1 == 0)
    while(sum(Logic) >= 1) {
      h <- 1.1 * h
      KERNEL1 <- kernel(SEMIMETRIC1/h)
      diag(KERNEL1) <- 0
      KERNEL1[KERNEL1 < 0] <- 0
      KERNEL1[KERNEL1 > 1] <- 0
      Denom1 <- apply(KERNEL1, 2, sum)
      Logic <- (Denom1 == 0)
    }
    RESPKERNEL1 <- KERNEL1 * Response
    Response.estimated <- apply(RESPKERNEL1, 2, sum)/Denom1
    count2 <- count2 + 1
    h.seq.corrected[count2] <- h
    Mse[count2] <- sum((Response.estimated - Response)^2)
  }
  index.opt <- order(Mse)[1]
  h.opt <- h.seq.corrected[index.opt]
  KERNEL1 <- kernel(SEMIMETRIC1/h.opt)
  KERNEL1[KERNEL1 < 0] <- 0
  KERNEL1[KERNEL1 > 1] <- 0
  diag(KERNEL1) <- 0
  Denom1 <- apply(KERNEL1, 2, sum)
  RESPKERNEL1 <- KERNEL1 * Response
  Response.estimated <- apply(RESPKERNEL1, 2, sum)/Denom1
  Mse.estimated <- sum((Response.estimated - Response)^2)/length(Response)
  if(twodatasets) {
    if(semimetric == "mplsr")
      SEMIMETRIC2 <- sm(Response, CURVES, PRED, ...)
    else SEMIMETRIC2 <- sm(CURVES, PRED, ...)
    KERNEL2 <- kernel(SEMIMETRIC2/h.opt)
    KERNEL2[KERNEL2 < 0] <- 0
    KERNEL2[KERNEL2 > 1] <- 0
    Denom2 <- apply(KERNEL2, 2, sum)
    Logic <- (Denom2 == 0)
    while(sum(Logic) >= 1) {
      h.opt <- 1.1 * h.opt
      KERNEL2 <- kernel(SEMIMETRIC2/h.opt)
      diag(KERNEL2) <- 0
      KERNEL2[KERNEL2 < 0] <- 0
      KERNEL2[KERNEL2 > 1] <- 0
      Denom2 <- apply(KERNEL2, 2, sum)
      Logic <- (Denom2 == 0)
    }
    RESPKERNEL2 <- KERNEL2 * Response
    Response.predicted <- apply(RESPKERNEL2, 2, sum)/Denom2
    return(list(Estimated.values = Response.estimated, 
                Predicted.values = Response.predicted, hopt = h.opt, 
                hseq = h.seq.corrected, Mse = Mse.estimated))
  }else {
    return(list(Estimated.values = Response.estimated, hopt = h.opt,
                hseq = h.seq.corrected, Mse = Mse.estimated))
  }
}


#pSH = funopare.kernel.cv()




