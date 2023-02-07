################
##
## @description Theory Of Reproducibility, demonstration of 
## convergence to reproducibility rate  and mean reproducibility rate
# strong law of large numbers (i.e., a.s. convergence) (data for Figure 1)  
##
## @param None
##
## @return None
##
## @lastChange 2020-02-06
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
#  library(MCMCpack)


#############
## PATHS
#############
baseDir <- "/home/erkanb/Desktop/lab/Reproducibility/TheoryOfReproducibility"
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
## Error variance
sigma <- c(0.5, 0.8)
## Sample size
sampleSizesmall <- 10
sampleSizelarge <- 100
## Run length for each simulation
Nrep <-300
nrep.vec=c(1:(Nrep-1))
## Number of simulations (independent paths given models and betas)
Nsim <- 100

# Original result
org.result = c("True", "False")
Norg = length(org.result)
#--------------------------------------------------------------------

## Generate all models
modelstemp <- generateModels(k)
models = list(c(),c(),c())
models[[1]] = modelstemp[[2]] # beta1 + beta2
models[[2]] = modelstemp[[3]] # beta1 + beta3
models[[3]] = modelstemp[[4]] # beta1 + beta2+  beta3
L=length(models)

## Correlation between predictor variables 
# (backcompatibility from PLOS ONE. Set to zero because not relevant for this code)
correlation <- 0
#
## Generate predictors
xsetlarge <- generateXSet(sampleSizelarge, k, correlation)
xsetsmall <- generateXSet(sampleSizesmall, k, correlation)
#
## Get AIC Constant (Applies both to AIC and BIC) that depends 
# only on x for all Models
# ICConstantlarge = GetAICConstant(models, xsetlarge)
# ICConstantsmall = GetAICConstant(models, xsetsmall)

tModel=models[[2]] # True Model
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
#------------------------------------------------------------------------
for (s in 1:length(sigma)){ # for each error variance
# Results as vector
result.AIC.large = result.AIC.small = 
result.BIC.large = result.BIC.small = rep(NA,Nrep)

#----------------------------------------------------------------
## Estimated reproducibility rate: i x j x l array
# i: Niter-1 (length of each simulation i.e., number of replication studies)
# j: simulation index. betas are the same throughout, the samples are random only for each run
# i.e., to show individual paths to convergence
# l: Norg Indicator for the original result 
# (when the model selected in the original study is true (model 2), 
# false (model 1))
# ------------------------------------------------------------
phihat.AIC.large = phihat.AIC.small =
phihat.BIC.large = phihat.BIC.small =
phihat.mean = array(NA,c(Nrep-1,Nsim, Norg)) 
# ------------------------------------------------------------
tModelBetas <- getBetas(tModel, betas, sigma[s]) # betas are weights
deterministic.large <- calculateDet(tModel, xsetlarge, betas, tModelBetas)
deterministic.small <- calculateDet(tModel, xsetsmall, betas, tModelBetas)
#----------------------------------------------------------
#----------------------------------------------------------
for(l in 1:Norg){
  # simulation for each original result type (org.result) start here
  # If the original result is a true result (i.e., model 2 selected)
  if(l==1){result.AIC.large[1] = result.AIC.small[1] = result.BIC.large[1] = result.BIC.small[1] = 2}
  # If the original result is a false result (i.e, model 1 selected)
  if(l==2){result.AIC.large[1] = result.AIC.small[1] = result.BIC.large[1] = result.BIC.small[1] = 1}
for(j in 1:Nsim){
  #----------------------------------------------------
  for (i in 2:Nrep){
    # each path starts here
    ylarge <- generateY(deterministic.large, sigma[s])
    ysmall <- generateY(deterministic.small, sigma[s])
    #
    # fit models for large sample size
    lm.large.1 <- lm(ylarge ~ xsetlarge[,1] + xsetlarge[,2]) # models[[1]]
    lm.large.2 <- lm(ylarge ~ xsetlarge[,1] + xsetlarge[,3]) # models[[2]]
    lm.large.3 <- lm(ylarge ~ xsetlarge[,1] + xsetlarge[,2] + xsetlarge[,3]) # models[[3]]
    #
    # fit models for small sample size
    lm.small.1 <- lm(ysmall ~ xsetsmall[,1] + xsetsmall[,2]) # models[[1]]
    lm.small.2 <- lm(ysmall ~ xsetsmall[,1] + xsetsmall[,3]) # models[[2]]
    lm.small.3 <- lm(ysmall ~ xsetsmall[,1] + xsetsmall[,2] + xsetsmall[,3]) # models[[3]]
    # 
    AIC.large = c(AIC(lm.large.1), AIC(lm.large.2), AIC(lm.large.3))
    AIC.small = c(AIC(lm.small.1), AIC(lm.small.2), AIC(lm.small.3))
    BIC.large = c(BIC(lm.large.1), BIC(lm.large.2), BIC(lm.large.3))
    BIC.small = c(BIC(lm.small.1), BIC(lm.small.2), BIC(lm.small.3))
    #----------------------------
    result.AIC.large[i] = which(AIC.large == min(AIC.large))
    result.AIC.small[i] = which(AIC.small == min(AIC.small))
    result.BIC.large[i] = which(BIC.large == min(BIC.large))
    result.BIC.small[i] = which(BIC.small == min(BIC.small))
    # 
  }
  # binary: result reproduced or not with respect to the first result in the sequence
  logical.AIC.large = (result.AIC.large[-1] == result.AIC.large[1]) 
  logical.AIC.small = (result.AIC.small[-1] == result.AIC.small[1])
  logical.BIC.large = (result.BIC.large[-1] == result.BIC.large[1])
  logical.BIC.small = (result.BIC.small[-1] == result.BIC.small[1])
  # Recall
  ## Estimated reproducibility rate: i x l matrix
  # i: Niter (length of each simulation i.e., number of replication studies)
  # l: Norg Indicator for the original result 
  # ( when the model selected in the original study is true (model 2), 
  #   false (model 1), and random (selected from the data))
  # ------------------------------------------------------------
  phihat.AIC.large[,j,l] = cumsum(logical.AIC.large)/nrep.vec
  phihat.AIC.small[,j,l] = cumsum(logical.AIC.small)/nrep.vec
  phihat.BIC.large[,j,l] = cumsum(logical.BIC.large)/nrep.vec 
  phihat.BIC.small[,j,l] = cumsum(logical.BIC.small)/nrep.vec 
  # Mean phihat estimated from quasi-replications, with idealized study in each step 
  # chosen randomly from four type of idealized studies chosen above
  logical.long = c(logical.AIC.large, logical.AIC.small, logical.BIC.large, logical.BIC.small)
  ind = sample(length(logical.long),(Nrep-1),replace = FALSE)
  phihat.mean[,j,l] = cumsum(logical.long[ind])/nrep.vec
  # the first element is the result itself, so no reproducibility rate
  # can be estimated
}
}
# Mean of each idealized experiment (over j)
phihat.AIC.large.mean = apply(phihat.AIC.large, c(1,3), mean)
phihat.AIC.small.mean = apply(phihat.AIC.small, c(1,3), mean)
phihat.BIC.large.mean = apply(phihat.BIC.large, c(1,3), mean)
phihat.BIC.small.mean = apply(phihat.BIC.small, c(1,3), mean)
phihat.mean.mean = apply(phihat.mean, c(1,3), mean)

setwd("/home/erkanb/Desktop/lab/Reproducibility/TheoryOfReproducibility/data")

# HERE WRITE THE SAVE WITH SIGMA INDEX 
outfile = sprintf("Sigma%s.Convergence.RData",sigma[s])
save(list = ls(all.names =TRUE), file=outfile)
}