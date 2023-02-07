m = 1e3
k = 1e2
n = 1e3 # Binomial number of trials
p = 0.01 # probability of success, Bernoulli

lambda = n*p # Poisson rate = 10, Poisson approximation to binomial
mu = n*p # Normal mean = 10, normal approximation to binomial
sigma=sqrt(n*p*(1-p)) # Normal standard deviation = 3.14

#--------------------------------------------------------------
#--------------------------------------------------------------
c = 0.1 # result 0-1 constant
#--------------------------------------------------------------
#--------------------------------------------------------------
# Poisson prior
#--------------------------------------------------------------
# alpha = 10 # parameter for conjugate prior gamma
alpha = 5 # parameter for conjugate prior gamma
# beta = 10 # parameter for conjugate prior gamma
beta = 5 # parameter for conjugate prior gamma

epsilon.poi = c*sqrt(lambda) # result acceptance
#--------------------------------------------------------------
# Normal prior
#--------------------------------------------------------------
tau0.sq = 1 # parameter for conjugate prior normal
# mu0 = 0 # parameter for conjugate prior normal
mu0 = 10 # parameter for conjugate prior normal
tau.sq = 1/(sigma^2)
epsilon.norm = c*(sigma) # result acceptance
#--------------------------------------------------------------

reproducibility.rate.poisson.smalln.mle = matrix(NA, k, m)
reproducibility.rate.poisson.smalln.bayes = matrix(NA, k, m)
reproducibility.rate.poisson.largen.mle = matrix(NA, k, m)
reproducibility.rate.poisson.largen.bayes = matrix(NA, k, m)
reproducibility.rate.normal.smalln.mle = matrix(NA, k, m)
reproducibility.rate.normal.smalln.bayes = matrix(NA, k, m)
reproducibility.rate.normal.largen.mle = matrix(NA, k, m)
reproducibility.rate.normal.largen.bayes = matrix(NA, k, m)


for( i in 1:k){
#--------------------------------------------------------------
# Poisson
#--------------------------------------------------------------
# n = Small, M_A= Poisson, Spost = MLE
n = 30
x = matrix(rpois(n*m, lambda), m, n) # data
mean.x = rowMeans(x)
result.poisson.smalln.mle  = (mean.x<lambda+epsilon.poi & mean.x>lambda-epsilon.poi)
reproducibility.rate.poisson.smalln.mle[i,] = cumsum(result.poisson.smalln.mle)/(1:m)

# n = Small, M_A= Poisson, Spost = Bayes
n = 30
x = matrix(rpois(n*m, lambda), m, n) # data
mode.x = (rowSums(x)+alpha-1)/(n+beta)
result.poisson.smalln.bayes = (mode.x<lambda+epsilon.poi & mode.x>lambda-epsilon.poi)
reproducibility.rate.poisson.smalln.bayes[i,] = cumsum(result.poisson.smalln.bayes)/(1:m)

# n = Large, M_A= Poisson, Spost = MLE
n = 200
x = matrix(rpois(n*m, lambda), m, n) # data
mean.x = rowMeans(x)
result.poisson.largen.mle = (mean.x<lambda+epsilon.poi & mean.x>lambda-epsilon.poi)
reproducibility.rate.poisson.largen.mle[i,] = cumsum(result.poisson.largen.mle)/(1:m)

# n = Large, M_A= Poisson, Spost = Bayes
n = 200
x = matrix(rpois(n*m, lambda), m, n) # data
mode.x = (rowSums(x)+alpha-1)/(n+beta)
result.poisson.largen.bayes = (mode.x<lambda+epsilon.poi & mode.x>lambda-epsilon.poi)
reproducibility.rate.poisson.largen.bayes[i,] = cumsum(result.poisson.largen.bayes)/(1:m)
#--------------------------------------------------------------
#--------------------------------------------------------------
# Normal
#--------------------------------------------------------------
# n = Small, M_A= Normal, Spost = MLE
n = 30
x =  matrix(rnorm(m*n, mu, sigma), m,n) # data
mean.x = rowMeans(x)
result.normal.smalln.mle = (mean.x<mu+epsilon.norm & mean.x>mu-epsilon.norm)
reproducibility.rate.normal.smalln.mle[i,] = cumsum(result.normal.smalln.mle)/(1:m)

# n = Small, M_A= Normal, Spost = bayes
n = 30  
x =  matrix(rnorm(m*n, mu, sigma), m,n) # data
mode.x = rowSums(x)*(tau.sq/(tau0.sq+n*tau.sq))+mu0*(tau0.sq/(tau0.sq+n*tau.sq)) # posterior mode
result.normal.smalln.bayes = (mode.x<mu+epsilon.norm & mode.x>mu-epsilon.norm)
reproducibility.rate.normal.smalln.bayes[i,] = cumsum(result.normal.smalln.bayes)/(1:m)

# n = Large, M_A= Normal, Spost = MLE
n = 200
x =  matrix(rnorm(m*n, mu, sigma), m,n) # data
mean.x = rowMeans(x) # posterior mode
result.normal.largen.mle = (mean.x<mu+epsilon.norm & mean.x>mu-epsilon.norm)
reproducibility.rate.normal.largen.mle[i,] = cumsum(result.normal.largen.mle)/(1:m)

# n = Large, M_A= Normal, Spost = bayes
n = 200
x =  matrix(rnorm(m*n, mu, sigma), m,n) # data
mode.x = rowSums(x)*(tau.sq/(tau0.sq+n*tau.sq))+mu0*(tau0.sq/(tau0.sq+n*tau.sq))
result.normal.largen.bayes = (mode.x<mu+epsilon.norm & mode.x>mu-epsilon.norm)
reproducibility.rate.normal.largen.bayes[i,] = cumsum(result.normal.largen.bayes)/(1:m)
}

mean.reproducibility.rate.poisson.smalln.mle = mean(reproducibility.rate.poisson.smalln.mle[,m])
mean.reproducibility.rate.poisson.smalln.bayes = mean(reproducibility.rate.poisson.smalln.bayes[,m])
mean.reproducibility.rate.poisson.largen.mle = mean(reproducibility.rate.poisson.largen.mle[,m])
mean.reproducibility.rate.poisson.largen.bayes = mean(reproducibility.rate.poisson.largen.bayes[,m])

mean.reproducibility.rate.normal.smalln.mle = mean(reproducibility.rate.normal.smalln.mle[,m])
mean.reproducibility.rate.normal.smalln.bayes = mean(reproducibility.rate.normal.smalln.bayes[,m])
mean.reproducibility.rate.normal.largen.mle = mean(reproducibility.rate.normal.largen.mle[,m])
mean.reproducibility.rate.normal.largen.bayes = mean(reproducibility.rate.normal.largen.bayes[,m])

# This is the order of models
# 1: poisson, small n, mle
# 2: poisson, small n, bayes
# 3: poisson, large n, mle
# 4: poisson, large n, bayes
# 5: normal, small n, mle
# 6: normal, small n, bayes
# 7: normal, large n, mle
# 8: normal, large n, bayes

vec.ind = c(mean.reproducibility.rate.poisson.smalln.mle,
                                         mean.reproducibility.rate.poisson.smalln.bayes,
                                         mean.reproducibility.rate.poisson.largen.mle,
                                         mean.reproducibility.rate.poisson.largen.bayes,
                                         mean.reproducibility.rate.normal.smalln.mle,
                                         mean.reproducibility.rate.normal.smalln.bayes,
                                         mean.reproducibility.rate.normal.largen.mle,
                                         mean.reproducibility.rate.normal.largen.bayes)
ind.sample.low = order(vec.ind)[1:4]
ind.sample.high = order(vec.ind)[5:8]

#--------------------------------------------------------------
# Randomly chosen scenario, for three cases:
# 1. Each experiment in the sequence is uniformly randomly chosen
# with equal probability 1/8 among all experiments
#  2. Each experiment in the sequence is uniformly randomly chosen
# with equal probability 1/4 among 4 lowest reproducibility rate experiments
#  3. Each experiment in the sequence is uniformly randomly chosen
# with equal probability 1/4 among 4 highest reproducibility rate experiments

#--------------------------------------------------------------
# Case 1:
# Mixture of 8 random variables with equal probability
# Reminder: This is the order of models
# 1: poisson, small n, mle
# 2: poisson, small n, bayes
# 3: poisson, large n, mle
# 4: poisson, large n, bayes
# 5: normal, small n, mle
# 6: normal, small n, bayes
# 7: normal, large n, mle
# 8: normal, large n, bayes
#--------------------------------------------------------------
result.matrix = rbind(result.poisson.smalln.mle, result.poisson.smalln.bayes, result.poisson.largen.mle,
                      result.poisson.largen.bayes,
                      result.normal.smalln.mle, result.normal.smalln.bayes, result.normal.largen.mle,
                      result.normal.largen.bayes)

reproducibility.rate.random = matrix(NA, k, m) 
result.random = rep(NA, m)

for( i in 1:k){
  # data
  ind = sample(8, m, replace=TRUE)
  result.random = rep(NA,m)
  for( j in 1:m){
    if(ind[j] == 1){result.random[j] =  result.matrix[1,j]}
    if(ind[j] == 2){result.random[j] =  result.matrix[2,j]}
    if(ind[j] == 3){result.random[j] =  result.matrix[3,j]}
    if(ind[j] == 4){result.random[j] =  result.matrix[4,j]}
    if(ind[j] == 5){result.random[j] =  result.matrix[5,j]}
    if(ind[j] == 6){result.random[j] =  result.matrix[6,j]}
    if(ind[j] == 7){result.random[j] =  result.matrix[7,j]}
    if(ind[j] == 8){result.random[j] =  result.matrix[8,j]}
  }
  reproducibility.rate.random[i,] = cumsum(result.random)/(1:m)
}
mean.reproducibility.rate.random = mean(reproducibility.rate.random[,m])
#--------------------------------------------------------------
# Case 2:
# Mixture of 4 random variables: Each experiment in the sequence is uniformly randomly chosen
# with equal probability 1/4 among 4 lowest reproducibility rate experiments
#--------------------------------------------------------------
reproducibility.rate.random.low = matrix(NA, k, m) 
result.random.low = rep(NA, m)


for( i in 1:k){
  # data
  ind = sample(ind.sample.low, m, replace=TRUE)
  result.random.low = rep(NA,m)
  for( j in 1:m){
    if(ind[j] == ind.sample.low[1]){result.random.low[j] =  result.matrix[ind.sample.low[1],j]}
    if(ind[j] == ind.sample.low[2]){result.random.low[j] =  result.matrix[ind.sample.low[2],j]}
    if(ind[j] == ind.sample.low[3]){result.random.low[j] =  result.matrix[ind.sample.low[3],j]}
    if(ind[j] == ind.sample.low[4]){result.random.low[j] =  result.matrix[ind.sample.low[4],j]}
  }
  reproducibility.rate.random.low[i,] = cumsum(result.random.low)/(1:m)
}
mean.reproducibility.rate.random.low = mean(reproducibility.rate.random.low[,m])
#--------------------------------------------------------------
# Case 3:
# Mixture of 4 random variables: Each experiment in the sequence is uniformly randomly chosen
# with equal probability 1/4 among 4 highest reproducibility rate experiments
#--------------------------------------------------------------
reproducibility.rate.random.high = matrix(NA, k, m) 
result.random.high = rep(NA, m)

for( i in 1:k){
  # data
  ind = sample(ind.sample.high, m, replace=TRUE)
  result.random.high = rep(NA,m)
  for( j in 1:m){
    if(ind[j] == ind.sample.high[1]){result.random.high[j] =  result.matrix[ind.sample.high[1],j]}
    if(ind[j] == ind.sample.high[2]){result.random.high[j] =  result.matrix[ind.sample.high[2],j]}
    if(ind[j] == ind.sample.high[3]){result.random.high[j] =  result.matrix[ind.sample.high[3],j]}
    if(ind[j] == ind.sample.high[4]){result.random.high[j] =  result.matrix[ind.sample.high[4],j]}
  }
  reproducibility.rate.random.high[i,] = cumsum(result.random.high)/(1:m)
}
mean.reproducibility.rate.random.high = mean(reproducibility.rate.random.high[,m])
#--------------------------------------------------------------

#--------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# PLOTS
#---------------------------------------------------------------------------
library(scales)

# Global plot parameters
offset.x = 50
alpha.best = 0.1
lwd.path = 0.1
#--------------------------------------------------------------
# colors:
# poisson
# rgb(0.7,1,0.4):green, 
# rgb(1,0.58,0):yellow, 
# rgb(0.21,0.57,0.90):blue, 
# rgb(0.68,0.5,1):violet, 

# normal
# rgb(0.80, 0.40, 0.46) : salmon
# rgb(0.53, 0.80, 0.93) : robin blue
# rgb(0.20, 0.13, 0.53) : purple
# rgb(0, 0.62, 0.45) : green

# mean
# rgb(1, 0, 0.8): pink
# mean low
# rgb(0, 0.78, 1): aqua
# mean high
# rgb(1,0.78,0):yellow

#--------------------------------------------------------------
pdf('/home/erkanb/Desktop/lab/Reproducibility/TheoryOfReproducibility/plots/Figure.ConvergenceOfReproducibilityRate.NonExact.pdf',
    width=13, height=8)
# Plot Layout
#--------------------------------------------
# my colorblind friendly colors
mygrass = "#B2FF66"
myblue = "#3691E6"
myorange = "#FF9400"
myviolet = "#AD80FF"
myrose = "#CC6675"
myrobin = "#87CCED"
myindigo = "#332187"
mygreen = "#009E73"
mypink = "#FF00CC"
myaqua = "#00C7FF"
myyellow = "#FFC700"

par(mfrow=c(2,2))
#--------------------------------------------------------------
# Panel 1
plot(1:m, reproducibility.rate.poisson.smalln.mle[1,], type='l', lwd=lwd.path, 
     col=alpha(mygrass, alpha=alpha.best),
     xlim=c(1,m+offset.x),ylim=c(0,1), 
     #main=TeX('$Poisson$'),
     main='A. Assumed Model: Poisson',
     cex.main=1.5,
     xlab='number of replications', ylab ='estimated reproducibility rate', cex.lab=1.5)

for(i in 2:k){
  lines(1:m, reproducibility.rate.poisson.smalln.mle[i,],type='l', lwd = lwd.path, col= alpha(mygrass, alpha = alpha.best))
  }
points(m + offset.x, mean.reproducibility.rate.poisson.smalln.mle, pch=16, col=mygrass)

for(i in 1:k){
  lines(1:m, reproducibility.rate.poisson.smalln.bayes[i,],type='l', lwd = lwd.path, col=alpha(myblue, alpha=alpha.best))
  }
points(m + offset.x,mean.reproducibility.rate.poisson.smalln.bayes, pch=16, col=myblue)

for(i in 1:k){
  lines(1:m, reproducibility.rate.poisson.largen.mle[i,],type='l', lwd = lwd.path, col=alpha(myorange, alpha=alpha.best))
  }
points(m + offset.x, mean.reproducibility.rate.poisson.largen.mle, pch=16, col=myorange)

for(i in 1:k){
  lines(1:m, reproducibility.rate.poisson.largen.bayes[i,],type='l', lwd = lwd.path, col=alpha(myviolet, alpha=alpha.best))
  }
points(m + offset.x, mean.reproducibility.rate.poisson.largen.bayes, pch=16, col=myviolet)
#--------------------------------------------------------------
# Panel 2
# Normal
# rgb(0.80, 0.40, 0.46) : salmon
# rgb(0.53, 0.80, 0.93) : robin blue
# rgb(0.20, 0.13, 0.53) : purple
# rgb(0, 0.62, 0.45) : green

plot(1:m, reproducibility.rate.normal.smalln.mle[1,], type='l', lwd=lwd.path, 
     col=alpha(myrose, alpha=alpha.best),
     xlim=c(1,m+offset.x),ylim=c(0,1), 
     #main=TeX('$Poisson$'),
     main='B. Assumed Model: Normal',
     cex.main=1.5,
     xlab='number of replications', ylab ='estimated reproducibility rate', cex.lab=1.5)
for(i in 2:k){
  lines(1:m, reproducibility.rate.normal.smalln.mle[i,],type='l', lwd = lwd.path, col=alpha(myrose, alpha=alpha.best))
}
points(m + offset.x, mean.reproducibility.rate.normal.smalln.mle, pch=16, col=myrose)

for(i in 1:k){
  lines(1:m, reproducibility.rate.normal.smalln.bayes[i,],type='l', lwd = lwd.path, col=alpha(myrobin, alpha=alpha.best))
  }
points(m + offset.x,mean.reproducibility.rate.normal.smalln.bayes, pch=16, col=myrobin)

for(i in 1:k){
  lines(1:m, reproducibility.rate.normal.largen.mle[i,],type='l', lwd = lwd.path, col=alpha(myindigo, alpha=alpha.best))
  }
points(m + offset.x,mean.reproducibility.rate.normal.largen.mle, pch=16, col=myindigo)

for(i in 1:k){
  lines(1:m, reproducibility.rate.normal.largen.bayes[i,],type='l', lwd = lwd.path, col=alpha(mygreen, alpha=alpha.best))
  }
points(m + offset.x,mean.reproducibility.rate.normal.largen.bayes, pch=16, col=mygreen)
#--------------------------------------------------------------
# Panel 3
# mean
# rgb(1, 0, 0.8): pink
# mean low
# rgb(0, 0.78, 1): aqua
# mean high
# rgb(1,0.78,0):yellow
#--------------------------------------------------------------
plot(1:m,  reproducibility.rate.random[1,], type='l', lwd=lwd.path, 
     col=alpha(mypink, alpha=alpha.best),
     xlim=c(1,m+offset.x),ylim=c(0,1), 
     #main=TeX('$Poisson$'),
     main='C. Non-exact replications',
     cex.main=1.5,
     xlab='number of replications', ylab ='estimated reproducibility rate', cex.lab=1.5)

for(i in 2:k){
  lines(1:m, reproducibility.rate.random[i,], type='l', lwd = lwd.path, col=alpha(mypink, alpha=alpha.best))
}
points(m + offset.x,mean.reproducibility.rate.random, pch=16, col=mypink)
# mean low
for(i in 1:k){
  lines(1:m, reproducibility.rate.random.low[i,], type='l', lwd = lwd.path, col=alpha(myaqua, alpha=alpha.best))
}
points(m + offset.x,mean.reproducibility.rate.random.low, pch=16, col=myaqua)
# mean high
for(i in 1:k){
  lines(1:m, reproducibility.rate.random.high[i,], type='l', lwd = lwd.path, col=alpha(myyellow, alpha=alpha.best))
}
points(m + offset.x,mean.reproducibility.rate.random.high, pch=16, col=myyellow)

#--------------------------------------------------------------
# Panel 4
#--------------------------------------------------------------
#
# evaluate at 50 replications (arbitrary)
mm = 50

aug.vec.ind = c(mean.reproducibility.rate.random,
                mean.reproducibility.rate.random.low,
                mean.reproducibility.rate.random.high,
            mean.reproducibility.rate.poisson.smalln.mle,
            mean.reproducibility.rate.poisson.smalln.bayes,
            mean.reproducibility.rate.poisson.largen.mle,
            mean.reproducibility.rate.poisson.largen.bayes,
            mean.reproducibility.rate.normal.smalln.mle,
            mean.reproducibility.rate.normal.smalln.bayes,
            mean.reproducibility.rate.normal.largen.mle,
            mean.reproducibility.rate.normal.largen.bayes)

var.vec = c(var(reproducibility.rate.random[,mm]), 
            var(reproducibility.rate.random.low[,mm]), 
            var(reproducibility.rate.random.high[,mm]),
            var(reproducibility.rate.poisson.smalln.mle[,mm]),
            var(reproducibility.rate.poisson.smalln.bayes[,mm]), 
            var(reproducibility.rate.poisson.largen.mle[,mm]), 
            var(reproducibility.rate.poisson.largen.bayes[,mm]),
            var(reproducibility.rate.normal.smalln.mle[,mm]),
            var(reproducibility.rate.normal.smalln.bayes[,mm]), 
            var(reproducibility.rate.normal.largen.mle[,mm]), 
            var(reproducibility.rate.normal.largen.bayes[,mm]))

# plot
plot(aug.vec.ind, var.vec, type="h", col="gray",
     lwd=3,
     main='D. Estimated variance',
     cex.main=1.5,
     xlab='estimated reproducibility rate', ylab ='estimated variance', cex.lab=1.5)
mycolors = c(mypink, myaqua, myyellow, mygrass, myblue, myorange, myviolet, myrose, myrobin, myindigo, mygreen)
for(i in 1:length(aug.vec.ind)){
  points(aug.vec.ind[i],var.vec[i], col=mycolors[i], pch=16, cex = 2)}
#--------------------------------------------------------------
#--------------------------------------------------------------
dev.off()
#--------------------------------------------------------------

# end all