# Make Plots

#--------------------------------------------------------------
library(latex2exp)
library(ggplot2)
library(corrplot)
#--------------------------------------------------------------
setwd("/home/erkan/Desktop/lab/Reproducibility/TheoryOfReproducibility/data")
#--------------------------------------------------------------
pdf('/home/erkan/Desktop/lab/Reproducibility/TheoryOfReproducibility/plots/Figure.ConvergenceOfReproducibilityRate.Random.Path.pdf',
    width=10, height=10)
#--------------------------------------------------------------
load("ReproducibilityRate.Random.Path.RData")

phihat.AIC <- as.matrix(phihat[,,1])
corrplot(phihat.AIC, method="color", is.corr=FALSE)

var.phihat.AIC <- as.matrix(var.phihat[,,1])
corrplot(var.phihat.AIC, method="color", is.corr=FALSE)

# rows: simulations (Nsim)
# columns: replications (Nrep)
# original data format is "wide"
# i.e., phihat is a matrix

# reshape it to data format "long"
data <- expand.grid(x=sampleSize, y=sigma)
phihat.vec <- c(phihat)
data$Z <- phihat.vec

# filled contour plot
gg <- ggplot(data, aes(x=sampleSize, y=sigma, z=phihat.vec)) +
  geom_tile(aes(fill = phihat.vec)) +
  geom_contour(colour = "white")

gg
#--------------------------------------



dev.off()