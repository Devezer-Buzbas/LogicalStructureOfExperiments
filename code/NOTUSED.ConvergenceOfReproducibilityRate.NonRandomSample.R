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
source(paste0(scriptDir, "/constants.R"))
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
source(paste0(scriptDir, "/modelSimilarByInteraction.R"))
source(paste0(scriptDir, "/modelSimilarByTerm.R"))
source(paste0(scriptDir, "/GetAICConstant.R"))
source(paste0(scriptDir, "/GetModelComparison.R"))


###################
## INPUT PARAMETERS
###################
## Number of factors
k <- 3
## Error variance
# sigma <- c(0.5, 0.8)
sigma <- 0.5
## Sample size
sampleSize <- c(10,50,100,250,500,1000)
lsample <- length(sampleSize) 
## Run length for each simulation
# Nrep <-300
Nrep <-100
nrep.vec=c(1:Nrep)
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
models[[3]] = modelstemp[[4]] # beta1 + beta2 + beta3
L=length(models)

## Correlation between predictor variables 
# (backcompatibility from PLOS ONE. Set to zero because not relevant for this code)
correlation = 0
#
## Generate predictors for the largest sample size
#ICConstant <- rep(NA,3)
xset <-list()
for( i in 1:lsample){
  xset[[i]] <- generateXSet(sampleSize[i], k, correlation)
#  ICConstant[i] = GetAICConstant(models, xset[[i]])
  ## Get AIC Constant (Applies both to AIC and BIC) that depends 
  # only on x for all Models
}
#

tModel = models[[2]] # True Model
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
  phihat.BIC <- array(NA,c(Nrep,Nsim, Norg, lsample))
  phihat.BIC.rand <- array(NA,c(Nrep, Nsim, Norg, lsample))
  ## Results: Estimated reproducibility rate: i x j x l x k array
  # i: Nrep (length of each simulation i.e., number of replication studies)
  # j: simulation index. betas are the same throughout, the samples are 
  # random only for each run i.e., to show individual paths to convergence
  # l: Norg Indicator for the original result 
  # k: lsample, the number of sample sizes
  # (when the model selected in the original study is true (model 2), false (model 1))
  # .rand is to estimate the true reproducibility rate when sample is random. See below
  # random eps vs. nonrandom eps
  phihat.BIC.rand.truereprate = phihat.BIC.truereprate = matrix(NA,Norg, lsample)
  # For the last element of phihat.BIC.rand and phihat.BIC in Nrep, averaged over Nsim, it gives the Monte Carlo estimate of the true reproducibility
  # rate to converge under random sample
  # ------------------------------------------------------------
  tModelBetas <- getBetas(tModel, betas, sigma[s]) # betas are weights
  deterministic <- list()
  for(k in 1:lsample){deterministic[[k]] <- calculateDet(tModel, xset[[k]], betas, tModelBetas)}
  #----------------------------------------------------------
  #----------------------------------------------------------
  for(l in 1:Norg){
    # simulation for each original result type (org.result) start here
      if(l==1){result.BIC.original = 2}
      # If the original result is a true result (i.e., model 2 selected)
      if(l==2){result.BIC.original = 1}
      # If the original result is a false result (i.e, model 1 selected)
      #----------------------------------------------------
      for(k in 1:lsample){
        for(j in 1:Nsim){
          result.BIC.rand =result.BIC=logical.BIC.rand=logical.BIC= matrix(NA,lsample,Nrep)
        for (i in 1:Nrep){
        # each path starts here
         eps.rand <- rnorm(sampleSize[k],0,sigma[s]) 
         # random sampling errors
         eps <- c(rnorm(1,0,sigma[s]),rep(NA,sampleSize[k]-1)) 
         for(p in 2:sampleSize[k]){eps[p] <- rnorm(1, eps[p-1], sigma[s])}
         # non random sampling errors
          EY.rand <- deterministic[[k]] + eps.rand
          EY <- deterministic[[k]] + eps
          y.rand <- (EY.rand - mean(EY.rand)) / sd(EY.rand)
          y <- (EY - mean(EY)) / sd(EY)
          # fit models for each sample size
          lm.truemodel.rand <- lm(y.rand ~ xset[[k]][,1]  + xset[[k]][,3]) # models[[2]]
          lm.falsemodel.rand <- lm(y.rand ~ xset[[k]][,1] + xset[[k]][,2]) # models[[1]]
          lm.truemodel <- lm(y ~ xset[[k]][,1]  + xset[[k]][,3]) # models[[2]]
          lm.falsemodel <- lm(y ~ xset[[k]][,1] + xset[[k]][,2]) # models[[1]]
          #----------------------------
          BIC.rand <- c(BIC(lm.truemodel.rand), BIC(lm.falsemodel.rand))
          BIC <- c(BIC(lm.truemodel), BIC(lm.falsemodel))
          #----------------------------
          result.BIC.rand[k,i] = which(BIC.rand == min(BIC.rand))
          result.BIC[k,i] = which(BIC == min(BIC))
          # 
          logical.BIC.rand[k,i] = (result.BIC.rand[k,i] == result.BIC.original)
          logical.BIC[k,i] = (result.BIC[k,i] == result.BIC.original)
          # binary: result reproduced or not with respect to the first result in the sequence
        # ------------------------------------------------------------
        }
        phihat.BIC.rand[,j,l,k] = cumsum(logical.BIC.rand[k,])/nrep.vec
        # the true reproducibility rate to converge
        phihat.BIC[,j,l,k] = cumsum(logical.BIC[k,])/nrep.vec 
        ## Estimated reproducibility rate: i x j x l x k matrix
        # i: Nrep (length of each simulation i.e., number of replication studies)
        # j: Nsim (number of independent simulations, i.e. paths) 
        # l: Norg Indicator for the original result 
        # (when the model selected in the original study is true (model 2), 
        #   false (model 1), and random (selected from the data))
        # k: sample size
      }
  # Summarize estimated true replication rate for random samples
     phihat.BIC.rand.truereprate[l,k] <- mean(phihat.BIC.rand[i,,l,k])
     phihat.BIC.truereprate[l,k] <- mean(phihat.BIC[i,,l,k])
   }
 }
    #--------------------------------------------------------
  setwd("/home/erkan/Desktop/lab/Reproducibility/TheoryOfReproducibility/data")
  # HERE WRITE THE SAVE WITH SIGMA INDEX 
  outfile = sprintf("Sigma%s.Convergence.NonRandom.RData",sigma[s])
  save(list = ls(all.names =TRUE), file=outfile)
}