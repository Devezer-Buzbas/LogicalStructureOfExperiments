# Make Plots
# 2 x 3 plot of convergence of reproducibility rate for random and nonrandom samples 
# (1, 1) small sample size, true model 
# (1, 2) medium sample size, true model
# (1, 3) large sample size, true model
# (2, 1) small sample size, false model 
# (2, 2) medium sample size, false model
# (2, 3) large sample size, false model
#--------------------------------------------------------------
library(latex2exp)
#--------------------------------------------------------------
setwd("/home/erkan/Desktop/lab/Reproducibility/TheoryOfReproducibility/data")
rm(list=ls())

# Global plot parameters
offset.x = 10
alpha.best = 0.1
lwd.path = 1
#--------------------------------------------------------------
# colors:  rgb(0.8,0.8,0.8):gray, rgb(0.21,0.57,0.90):blue, 
pdf('/home/erkan/Desktop/lab/Reproducibility/TheoryOfReproducibility/plots/Figure.ConvergenceOfReproducibilityRate.NonRandomSample.pdf',
    width=13, height=8)
#--------------------------------------------------------------
# Layout
m <- matrix(c(1,2,3,4,5,6,7,7,7), nrow = 3,ncol = 3,byrow = TRUE)
# 1 2 3
# 4 5 6
# 7 7 7

layout(mat = m,heights = c(0.4,0.4,0.2))
#--------------------------------------------------------------
load("Sigma0.5.Convergence.NonRandom.RData")
# Objects to plot: phihat.BIC[i,j,l,k] (for nonrandom samples)
# and phihat.BIC.rand[i,j,l,k] (for random samples) 
## Estimated reproducibility rate: i x j x l x k matrix
# i: Nrep (length of each simulation i.e., number of replication studies)
# j: Nsim (number of independent simulations, i.e. paths) 
# l: Norg Indicator for the original result (True results l=1, False results l=2) 
# (when the model selected in the original study is true (model 2), 
#   false (model 1), and random (selected from the data))
# k: sample size
# 
# Object to plot: phihat.BIC.rand.truereprate[l,k]
# Summarize estimated true replication rate for random samples (star)
# ---------------------------------------------------------------
l=2 # Original result is the TRUE model selected
for (k in 1:lsample){
plot(nrep.vec,phihat.BIC.rand[,1,l,k], type='l', lwd=lwd.path, 
     col=rgb(0.8, 0.8, 0.8, alpha=alpha.best),
     xlim=c(1,Nrep+offset.x),ylim=c(0,1), 
     main=TeX('$E(Y):\\sigma = 1:1$', paste("n",sampleSize[k])), cex.main=1.5,
     main=paste("mean:signal = 1:1$, n=",sampleSize[k]), cex.main=1.5,
     xlab='number of replications', ylab ='estimated reproducibility rate', cex.lab=1.5)
for(j in 2:Nsim){lines(nrep.vec,phihat.BIC.rand[,j,l,k], type='l', lwd = lwd.path, col=rgb(0.8, 0.8, 0.8, alpha=alpha.best))}
points(Nrep + offset.x,phihat.BIC.rand.truereprate[l,k], pch=8, col=rgb(0.8, 0.8, 0.8))
for(j in 1:Nsim){lines(nrep.vec,phihat.BIC[,j,l,k], type='l', lwd = lwd.path, col=rgb(0.21,0.57,0.90, alpha=alpha.best))}
points(Nrep + offset.x,phihat.BIC.truereprate[l,k], pch=8, col=rgb(0.21,0.57,0.90))
}
# k=2
# plot(nrep.vec,phihat.BIC.rand[,1,l,k], type='l', lwd=lwd.path, 
#      col=rgb(0.8, 0.8, 0.8, alpha=alpha.best),
#      xlim=c(1,Nrep+offset.x),ylim=c(0,1), 
#      main=TeX('$E(Y):\\sigma = 1:1$'), cex.main=1.5,
#      xlab='number of replications', ylab ='estimated reproducibility rate', cex.lab=1.5)
# for(j in 2:Nsim){lines(nrep.vec,phihat.BIC.rand[,j,l,k], type='l', lwd = lwd.path, col=rgb(0.8, 0.8, 0.8, alpha=alpha.best))}
# points(Nrep + offset.x,phihat.BIC.rand.truereprate[l,k], pch=8, col=rgb(0.8, 0.8, 0.8))
# for(j in 1:Nsim){lines(nrep.vec,phihat.BIC[,j,l,k], type='l', lwd = lwd.path, col=rgb(0.21,0.57,0.90, alpha=alpha.best))}
# points(Nrep + offset.x,phihat.BIC.truereprate[l,k], pch=8, col=rgb(0.21,0.57,0.90))
# 
# k=3
# plot(nrep.vec,phihat.BIC.rand[,1,l,k], type='l', lwd=lwd.path, 
#      col=rgb(0.8, 0.8, 0.8, alpha=alpha.best),
#      xlim=c(1,Nrep+offset.x),ylim=c(0,1), 
#      main=TeX('$E(Y):\\sigma = 1:1$'), cex.main=1.5,
#      xlab='number of replications', ylab ='estimated reproducibility rate', cex.lab=1.5)
# for(j in 2:Nsim){lines(nrep.vec,phihat.BIC.rand[,j,l,k], type='l', lwd = lwd.path, col=rgb(0.8, 0.8, 0.8, alpha=alpha.best))}
# points(Nrep + offset.x,phihat.BIC.rand.truereprate[l,k], pch=8, col=rgb(0.8, 0.8, 0.8))
# for(j in 1:Nsim){lines(nrep.vec,phihat.BIC[,j,l,k], type='l', lwd = lwd.path, col=rgb(0.21,0.57,0.90, alpha=alpha.best))}
# points(Nrep + offset.x,phihat.BIC.truereprate[l,k], pch=8, col=rgb(0.21,0.57,0.90))
#-------------------------------------------------------------------------- 
#--------------------------------------------------------------------------
par(mar=c(1,1,1,1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 0,
      legend = c("Random sample", "Non random sample"), 
      col=c(rgb(0.8, 0.8, 0.8), rgb(0.21,0.57,0.90)), lwd=3, cex=1.3, horiz = TRUE)
dev.off()