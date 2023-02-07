################
##
## @description Theory Of Reproducibility, demonstration of 
## reproducibility rate for a random path: 
# sigma: randomly drawn from a beta distribution with parameters
# alpha=2, beta=2, 
# sampleSize: randomly drawn from discrete uniform (10, 500)
# method: randomly drawn from categorical, equal probability AIC or BIC
# (data for Figure x)  
##
## @param None
##
## @return None
##
## @lastChange 2021-04-25
##
## @changes
##   Fixed ...
##
################
rm(list=ls())
##
 library(caTools)
 library(data.table)
 library(permute)
 library(matrixStats)
 library(MCMCpack)


#############
## PATHS
#############
baseDir <- "/home/erkan/Desktop/lab/Reproducibility/TheoryOfReproducibility"
outputDir <- paste0(baseDir, "/data")
scriptDir <- paste0(baseDir, "/code")

#############
## FUNCTIONS
#############
source(paste0(scriptDir, "/calculateDet.R"))
# source(paste0(scriptDir, "/constants.R"))
source(paste0(scriptDir, "/convertBinary.R"))
source(paste0(scriptDir, "/modelToStr.R"))
source(paste0(scriptDir, "/strToModel.R"))
source(paste0(scriptDir, "/generateModels.R"))
source(paste0(scriptDir, "/compareModels.R"))
source(paste0(scriptDir, "/searchModel.R"))
source(paste0(scriptDir, "/generateY.R"))
source(paste0(scriptDir, "/generateXSet.R"))
source(paste0(scriptDir, "/getBetas.R"))
source(paste0(scriptDir, "/modelExpectation.R"))
#source(paste0(scriptDir, "/modelSimilarByInteraction.R"))
#source(paste0(scriptDir, "/modelSimilarByTerm.R"))
source(paste0(scriptDir, "/GetAICConstant.R"))
#source(paste0(scriptDir, "/GetModelComparison.R"))


###################
## INPUT PARAMETERS
###################
## Number of factors
k <- 3
# Original study parameters
sigma.org <- 0.5 # original study noise
sampleSize.org <- 50 # original study sample size
method.org <-"BIC" # original study method is set to BIC
#--------------------------------------------------------------------
## Error variance
# uniform distributed, parameters:
sigma.min <- 0.1
sigma.max <- 1 
#
# Sample size
sampleSize.min <- 5
sampleSize.max <- 100
#
#--------------------------------------------------------------------
# Generate sigmas for Nsim*Nrep studies
sigma <-seq(sigma.min,sigma.max, length.out=100) 
# Generate sample sizes
sampleSize <- seq(sampleSize.min,sampleSize.max, length.out=100)
# Generate methods
method <-c("AIC","BIC")
#--------------------------------------------------------------------
## Run length for each simulation
Nrep <-30
## Number of simulations (independent paths given models and betas)
Nsample <- length(sampleSize)
Nsigma <- length(sigma)
Nmethod <- length(method)
#grid is <- Nsample*Nsigma* Nmethod
# for each simulation Nrep replications, we are interested 
# in reproducibility rate value at the last iteration, as an estimate
# of the true reproducibility rate that we converge to
# and the variance of iterations (over Nrep) 
#--------------------------------------------------------------------
## Generate all models
modelstemp <- generateModels(k)
models = list(c(),c(),c())
models[[1]] = modelstemp[[2]] # beta1 + beta2
models[[2]] = modelstemp[[3]] # beta1 + beta3
models[[3]] = modelstemp[[4]] # beta1 + beta2+  beta3
L = length(models)

# Same true model for all simulations Nsim and Nrep with same true betas
tModel <- models[[2]] # True Model
# Original result is TRUE for all in these simulations
true.result <- 2
# (as opposed to some other simulations in this project which consider 
# the original result to be TRUE
# as well as FALSE)
#-------------------------------------------------------------------
## Betas
fModel <- as.matrix(models[[1]])
for(model in models){
  if(!is.matrix(model)){
    model <- t(as.matrix(model))
  }
  if(nrow(model) > nrow(fModel)){
    fModel <- model
  }
}
beta <- c(1, rep(1, nrow(fModel)))
betas <- cbind(0, beta[1])
f <- 10^((k-1):0)
for(r in 1:nrow(fModel)){
  index <- sum(fModel[r,] * f)
  betas <- rbind(betas, cbind(as.numeric(index), as.numeric(beta[r+1])))
}

## Correlation between predictor variables 
# (backcompatibility from PLOS ONE. Set to zero because not relevant for this code)
correlation <- 0
#
# -----------------------------------------------------------------------
# Results: Reproducibility rate for each simulation
phihat = var.phihat <- array(NA, dim=c(Nsample, Nsigma, Nmethod))
#------------------------------------------------------------------------
xset.org <- generateXSet(sampleSize.org, k, correlation) # same for all Nsim and Nrep
tModelBetas.org <- getBetas(tModel, betas, sigma.org) # betas are weights, same for all Nsim and Nrep
deterministic.org <- calculateDet(tModel, xset.org, betas, tModelBetas.org) # same for all Nsim and Nrep
#------------------------------------------------------------------------
result.org <- rep(NA, Nrep)
for(j in 1:Nrep){
y.org <- generateY(deterministic.org, sigma.org) 
# independent samples

lm.1.org <- lm(y.org ~ xset.org[,1] + xset.org[,2]) # models[[1]]
lm.2.org <- lm(y.org ~ xset.org[,1] + xset.org[,3]) # models[[2]]
lm.3.org <- lm(y.org ~ xset.org[,1] + xset.org[,2] + xset.org[,3]) # models[[3]]
# fit models

method.stat.org = c(BIC(lm.1.org), BIC(lm.2.org), BIC(lm.3.org)) 
# original study method method.org is defined as BIC in the parameter definitions 
result.org[j] <-  which(method.stat.org == min(method.stat.org))
}
phihat.org = sum(result.org == true.result)/Nrep 

for (i in 1:Nsample){ # for each sample size, the path starts here
  xset <- generateXSet(sampleSize[i], k, correlation)
  for (j in 1:Nsigma){ # for each error variance, the path starts here
    tModelBetas <- getBetas(tModel, betas, sigma[j]) # betas are weights
    ## Generate predictors
    deterministic <- calculateDet(tModel, xset, betas, tModelBetas)
    #----------------------------------------------------------
    for(m in 1:Nmethod){ # for each method, the path starts here
      #----------------------------------------------------------------
result <- rep(NA, Nrep)
#
for(l in 1:Nrep){
    y <- generateY(deterministic, sigma[j])
    # independently generated using original studies parameters
    # for all Nsim and Nrep
    #
    # fit models
    lm.1 <- lm(y ~ xset[,1] + xset[,2]) # models[[1]]
    lm.2 <- lm(y ~ xset[,1] + xset[,3]) # models[[2]]
    lm.3 <- lm(y ~ xset[,1] + xset[,2] + xset[,3]) # models[[3]]
    # 
    # ----------------------------------------------------------------------
    if(method[m]=="AIC"){method.stat = c(AIC(lm.1), AIC(lm.2), AIC(lm.3))}
    if(method[m]=="BIC"){method.stat = c(BIC(lm.1), BIC(lm.2), BIC(lm.3))}
    #----------------------------
    result[l] <- which(method.stat == min(method.stat))
    #
  }
  # ------------------------------------------------------------
  phihat[i,j,m] = sum(result== true.result)/Nrep
  # phihat estimated from quasi-replications, with idealized study in each step 
  # chosen randomly from conditions specified for each step of Nrep
  var.phihat[i,j,m] = var(result)
  }
}
}
setwd("/home/erkan/Desktop/lab/Reproducibility/TheoryOfReproducibility/data")

# HERE WRITE AS A SINGLE FILE
outfile = sprintf("ReproducibilityRate.Random.Path.RData")
save(list = ls(all.names =TRUE), file=outfile)