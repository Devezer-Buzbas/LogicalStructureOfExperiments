# Make Plots
# Four panel 1 x 4 plot of convergence of reproducibility rate 
# Linear regression example
# (1, 1) small sigma, true result 
# (1, 2) small sigma, false result 
# (1, 3) large sigma, true result 
# (1, 4) large sigma, false result 
#--------------------------------------------------------------
library(latex2exp)
#--------------------------------------------------------------
setwd("/home/erkanb/Desktop/lab/Reproducibility/TheoryOfReproducibility/data")

# Global plot parameters
offset.x = 10
alpha.best = 0.1
lwd.path = 1
#--------------------------------------------------------------
# colors:  
# rgb(0.7,1,0.4) : green 
# rgb(1,0.65,0) : yellow 
# rgb(0.68,0.5,1) : violet
# rgb(0.21,0.57,0.90) : blue
# rgb(0.5, 0.5, 0.5) : gray

pdf('/home/erkanb/Desktop/lab/Reproducibility/TheoryOfReproducibility/plots/Figure.ConvergenceOfReproducibilityRate.pdf',
    width=13, height=4)

#--------------------------------------------------------------
# Layout
m <- matrix(c(1,2,3,4,5,5,5,5), nrow = 2,ncol = 4,byrow = TRUE)

layout(mat = m,heights = c(0.8,0.2))
#--------------------------------------------------------------
# From left plots 1 and 2
# small sigma
load("Sigma0.5.Convergence.RData")
# TRUE results l=1
l=1
plot(nrep.vec,phihat.AIC.large[,1,l], type='l', lwd=lwd.path, 
     col=rgb(0.7,1,0.4, alpha=alpha.best),
     xlim=c(1,Nrep+offset.x),ylim=c(0,1), 
     main='Signal:Noise = 1:1, R = TRUE', cex.main=1.5,
     xlab='number of replications', ylab ='estimated reproducibility rate', cex.lab=1.5)
for(j in 2:Nsim){lines(nrep.vec,phihat.AIC.large[,j,l], type='l', lwd = lwd.path, col=rgb(0.7,1,0.4, alpha=alpha.best))}
points(Nrep + offset.x,phihat.AIC.large.mean[Nrep-1,l], pch=8, col=rgb(0.7,1,0.4))

for(j in 1:Nsim){lines(nrep.vec,phihat.BIC.large[,j,l], type='l', lwd = lwd.path, col=rgb(0.68,0.5,1, alpha=alpha.best))}
points(Nrep + offset.x,phihat.BIC.large.mean[Nrep-1,l], pch=8, col=rgb(0.68,0.5,1))

for(j in 1:Nsim){lines(nrep.vec,phihat.AIC.small[,j,l], type='l', lwd = lwd.path, col=rgb(1,0.65,0, alpha=alpha.best))}
points(Nrep + offset.x,phihat.AIC.small.mean[Nrep-1,l], pch=8, col=rgb(1,0.65,0))

for(j in 1:Nsim){lines(nrep.vec,phihat.BIC.small[,j,l], type='l', lwd = lwd.path, col=rgb(0.21,0.57,0.90, alpha=alpha.best))}
points(Nrep + offset.x,phihat.BIC.small.mean[Nrep-1,l], pch=8, col=rgb(0.21,0.57,0.90))

for(j in 1:Nsim){lines(nrep.vec,phihat.mean[,j,l], type='l', lwd = lwd.path, col=rgb(0.5,0.5,0.5, alpha=alpha.best))}
points(Nrep + offset.x,phihat.mean.mean[Nrep-1,l], pch=8, col=rgb(0.5,0.5,0.5))
#-------------------------------------------------------------------------- 

# FALSE results l=2
l=2
plot(nrep.vec,phihat.AIC.large[,1,l], type='l', lwd=lwd.path, 
     col=rgb(0.7,1,0.4, alpha=alpha.best),
     xlim=c(1,Nrep+offset.x),ylim=c(0,1), 
     main='Signal:Noise = 1:1, R = FALSE', cex.main=1.5,
     xlab='number of replications', ylab ='estimated reproducibility rate', cex.lab=1.5)
for(j in 2:Nsim){lines(nrep.vec,phihat.AIC.large[,j,l], type='l', lwd = lwd.path, col=rgb(0.7,1,0.4, alpha=alpha.best))}
points(Nrep + offset.x,phihat.AIC.large.mean[Nrep-1,l], pch=8, col=rgb(0.7,1,0.4))

for(j in 1:Nsim){lines(nrep.vec,phihat.BIC.large[,j,l], type='l', lwd = lwd.path, col=rgb(0.68,0.5,1, alpha=alpha.best))}
points(Nrep + offset.x,phihat.BIC.large.mean[Nrep-1,l], pch=8, col=rgb(0.68,0.5,1))

for(j in 1:Nsim){lines(nrep.vec,phihat.AIC.small[,j,l], type='l', lwd = lwd.path, col=rgb(1,0.65,0, alpha=alpha.best))}
points(Nrep + offset.x,phihat.AIC.small.mean[Nrep-1,l], pch=8, col=rgb(1,0.65,0))

for(j in 1:Nsim){lines(nrep.vec,phihat.BIC.small[,j,l], type='l', lwd = lwd.path, col=rgb(0.21,0.57,0.90, alpha=alpha.best))}
points(Nrep + offset.x,phihat.BIC.small.mean[Nrep-1,l], pch=8, col=rgb(0.21,0.57,0.90))

for(j in 1:Nsim){lines(nrep.vec,phihat.mean[,j,l], type='l', lwd = lwd.path, col=rgb(0.5,0.5,0.5, alpha=alpha.best))}
points(Nrep + offset.x,phihat.mean.mean[Nrep-1,l], pch=8, col=rgb(0.5,0.5,0.5))
#-------------------------------------------------------------------------- 
# From left plots 3 and 4
# large sigma
load("Sigma0.8.Convergence.RData")
# TRUE results l=1
l=1
  plot(nrep.vec,phihat.AIC.large[,1,l], type='l', lwd=lwd.path, 
       col=rgb(0.7,1,0.4, alpha=alpha.best),
       xlim=c(1,Nrep+offset.x),ylim=c(0,1), 
       main='Signal:Noise = 1:4, R = TRUE', cex.main=1.5,
       xlab='number of replications', ylab ='estimated reproducibility rate', cex.lab=1.5)
  for(j in 2:Nsim){lines(nrep.vec,phihat.AIC.large[,j,l], type='l', lwd = lwd.path, col=rgb(0.7,1,0.4, alpha=alpha.best))}
  points(Nrep + offset.x,phihat.AIC.large.mean[Nrep-1,l], pch=8, col=rgb(0.7,1,0.4))
  
  for(j in 1:Nsim){lines(nrep.vec,phihat.BIC.large[,j,l], type='l', lwd = lwd.path, col=rgb(0.68,0.5,1, alpha=alpha.best))}
  points(Nrep + offset.x,phihat.BIC.large.mean[Nrep-1,l], pch=8, col=rgb(0.68,0.5,1))
  
  for(j in 1:Nsim){lines(nrep.vec,phihat.AIC.small[,j,l], type='l', lwd = lwd.path, col=rgb(1,0.65,0, alpha=alpha.best))}
  points(Nrep + offset.x,phihat.AIC.small.mean[Nrep-1,l], pch=8, col=rgb(1,0.65,0))
  
  for(j in 1:Nsim){lines(nrep.vec,phihat.BIC.small[,j,l], type='l', lwd = lwd.path, col=rgb(0.21,0.57,0.90, alpha=alpha.best))}
  points(Nrep + offset.x,phihat.BIC.small.mean[Nrep-1,l], pch=8, col=rgb(0.21,0.57,0.90))
  
  for(j in 1:Nsim){lines(nrep.vec,phihat.mean[,j,l], type='l', lwd = lwd.path, col=rgb(0.5,0.5,0.5, alpha=alpha.best))}
  points(Nrep + offset.x,phihat.mean.mean[Nrep-1,l], pch=8, col=rgb(0.5,0.5,0.5))
#-------------------------------------------------------------------------- 
# FALSE results
  l = 2
  plot(nrep.vec,phihat.AIC.large[,1,l], type='l', lwd=lwd.path, 
       col=rgb(0.7,1,0.4, alpha=alpha.best),
       xlim=c(1,Nrep+offset.x),ylim=c(0,1), 
       main='Signal:Noise = 1:4, R = FALSE', cex.main=1.5,
       xlab='number of replications', ylab ='estimated reproducibility rate', cex.lab=1.5)
  for(j in 2:Nsim){lines(nrep.vec,phihat.AIC.large[,j,l], type='l', lwd = lwd.path, col=rgb(0.7,1,0.4, alpha=alpha.best))}
  points(Nrep + offset.x,phihat.AIC.large.mean[Nrep-1,l], pch=8, col=rgb(0.7,1,0.4))
  
  for(j in 1:Nsim){lines(nrep.vec,phihat.BIC.large[,j,l], type='l', lwd = lwd.path, col=rgb(0.68,0.5,1, alpha=alpha.best))}
  points(Nrep + offset.x,phihat.BIC.large.mean[Nrep-1,l], pch=8, col=rgb(0.68,0.5,1))
  
  for(j in 1:Nsim){lines(nrep.vec,phihat.AIC.small[,j,l], type='l', lwd = lwd.path, col=rgb(1,0.65,0, alpha=alpha.best))}
  points(Nrep + offset.x,phihat.AIC.small.mean[Nrep-1,l], pch=8, col=rgb(1,0.65,0))
  
  for(j in 1:Nsim){lines(nrep.vec,phihat.BIC.small[,j,l], type='l', lwd = lwd.path, col=rgb(0.21,0.57,0.90, alpha=alpha.best))}
  points(Nrep + offset.x,phihat.BIC.small.mean[Nrep-1,l], pch=8, col=rgb(0.21,0.57,0.90))
  
  for(j in 1:Nsim){lines(nrep.vec,phihat.mean[,j,l], type='l', lwd = lwd.path, col=rgb(0.5,0.5,0.5, alpha=alpha.best))}
  points(Nrep + offset.x,phihat.mean.mean[Nrep-1,l], pch=8, col=rgb(0.5,0.5,0.5))
#--------------------------------------------------------------------------
par(mar=c(1,1,1,1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 0,
      legend = c("AIC, n=100", "AIC, n=10", "BIC, n=100","BIC, n=10", 
                 "Non-Exact replications"), 
      col=c(rgb(0.7,1,0.4), rgb(1,0.65,0), rgb(0.68,0.5,1),
            rgb(0.21,0.57,0.90), rgb(0.5,0.5,0.5)), 
            lwd=3, cex=1.3, horiz = TRUE)
dev.off()