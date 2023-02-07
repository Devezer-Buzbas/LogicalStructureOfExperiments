k = c(2,8,20,50) # number of possible non-exact replication cases. For Figure 2 of the theory paper it is 8.
      # i.e., Poisson, MLE,n=200, Poisson, MLE,n=30, Poisson, Bayes,n=200,Poisson, Bayes,n=30,
      # Normal, MLE,n=200, Normal, MLE,n=30, Normal, Bayes,n=200, Normal, Bayes,n=30.
lk = length(k)
# phi.vec = list()
var.all = list()

phi.mix = rep(NA, lk)

min.phi = 0.2
max.phi = 0.5

for( i in 1:lk){
  phi.k = runif(k[i],min.phi,max.phi)
# phi.vec[[i]] = phi.k
  phi.mix[i] = mean(phi.k)
  phi.all = c(phi.k, phi.mix[i])
  var.all[[i]] = phi.all*(1-phi.all)
}

m = 1:10 # sequence for number of replications


#par(mfrow=c(1,lk))
par(mfrow=c(2,2))
for( i in 1:lk){
  plot(m, (1/m)*var.all[[i]][k[i]+1], type='l',col='red', lwd=2)
  for(j in 1: k[i]){
  lines(m, (1/m)*var.all[[i]][j], type='l')
    }
}