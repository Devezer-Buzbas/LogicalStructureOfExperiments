################
##
## @description Main reproducibility simulation to run in a cluster
##              (ARRAY functionality)
##
## @param None
##
## @return None
##
## @lastChange 2017-06-17
##
## @changes
##   Parameters output file (06/17/2017)
##   Added a new model X1 + X2 + X1X2 (06/08/2017)
##
################
library(caTools)
library(data.table)
library(permute)
library(matrixStats)
library(MCMCpack, lib="/scratch/nardluis/rpackages")


#############
## COMMAND-LINE
#############
args <- commandArgs(TRUE)
param <- as.integer(args[1])


#############
## PATHS
#############
baseDir <- "/scratch/nardluis"
scriptDir <- paste0(baseDir, "/scripts/reproducibility")
inputDir <- paste0(baseDir, "/data/reproducibility")
outputDir <- paste0(baseDir, "/data/reproducibility")


#############
## FUNCTIONS
#############
source(paste0(scriptDir, "/analysis.R"))
source(paste0(scriptDir, "/calculateDet.R"))
source(paste0(scriptDir, "/calculateDistance.R"))
source(paste0(scriptDir, "/compareModels.R"))
source(paste0(scriptDir, "/constants.R"))
source(paste0(scriptDir, "/convertBinary.R"))
source(paste0(scriptDir, "/generateBetas.R"))
source(paste0(scriptDir, "/generateModels.R"))
source(paste0(scriptDir, "/generateXSet.R"))
source(paste0(scriptDir, "/generateY.R"))
source(paste0(scriptDir, "/getBetas.R"))
source(paste0(scriptDir, "/getPredictors.R"))
source(paste0(scriptDir, "/modelSimilarByInteraction.R"))
source(paste0(scriptDir, "/modelSimilarByTerm.R"))
source(paste0(scriptDir, "/modelToStr.R"))
source(paste0(scriptDir, "/searchModel.R"))
source(paste0(scriptDir, "/seedGenerator.R"))
source(paste0(scriptDir, "/simulator.R"))
source(paste0(scriptDir, "/strToModel.R"))


###################
## INPUT PARAMETERS
###################
# Replications
# Default: 10
replications <- 100

# Length of the simulation
# Default: 10000
timesteps <- 10000

# Number of factors
# Default: 4
k <- 3

# Sample size
# Default: 100
sampleSize <- 100

# Generate all models
models <- generateModels(k)

# Range 0-239
if (k == 3) {
  m <- c(4, 10, 5) ## k = 3
} else if (k == 4) {
  m <- c(5, 15, 7) ## k = 4
}

# True model
tModelIndex <- as.integer(param / 80) + 1
if(tModelIndex == 1){
  tModel <- models[[m[1]]]       # X1 + X2 + X3
} else if(tModelIndex == 2){
  tModel <- models[[m[2]]]      # X1 + X2 + X3 + X1X2 + X1X3 + X2X3 + X1X2X3
} else if(tModelIndex == 3){
  tModel <- models[[m[3]]]        # X1 + X2 + X1X2
}

# Determine the number of predictors of each model
predictors <- getPredictors(models)

# Generate Betas
betaIndex <- as.integer((80 - ((80 * tModelIndex) - param)) / 20) + 1
if(betaIndex == 1){
  beta1 <- 1
  betaO <- 1
} else if(betaIndex == 2){
  beta1 <- 1
  betaO <- 10
} else if(betaIndex == 3){
  beta1 <- 10
  betaO <- 1
} else if(betaIndex == 4){
  beta1 <- 10
  betaO <- 10
}

# Generate weighted betas
betas <- generateBetas(models)

# Sigma (Variance of the error)
sigmaIndex <- as.integer((20 - ((20 * betaIndex) -
            (80 - ((80 * tModelIndex) - param)))) / 10) + 1
if(sigmaIndex == 1){
  sigma <- 0.2
} else if(sigmaIndex == 2){
  sigma <- 0.6
}

# Correlation
correlationIndex <- as.integer((10 - ((10 * sigmaIndex) -
            (20 - ((20 * betaIndex) -
                (80 - ((80 * tModelIndex) - param)))))) / 5) + 1
if(correlationIndex == 1){
  correlation <- 0.2
} else if(correlationIndex == 2){
  correlation <- 0.8
}

# Types of agents
#
# nRay - Number of replicators
# nRob - Number of robusts (term)
# nBob - Number of robusts (interaction)
# nNel - Number of novel
typeIndex <- (5 - ((5 * correlationIndex) -
        (10 - ((10 * sigmaIndex) -
            (20 - ((20 * betaIndex) -
                (80 - ((80 * tModelIndex) - param)))))))) + 1
if(typeIndex == 1){
  nRay <- 17
  nRob <- 1
  nBob <- 1
  nNel <- 1
} else if(typeIndex == 2){
  nRay <- 1
  nRob <- 17
  nBob <- 1
  nNel <- 1
} else if(typeIndex == 3){
  nRay <- 1
  nRob <- 1
  nBob <- 17
  nNel <- 1
} else if(typeIndex == 4){
  nRay <- 1
  nRob <- 1
  nBob <- 1
  nNel <- 17
} else if(typeIndex == 5){
  nRay <- 1
  nRob <- 1
  nBob <- 1
  nNel <- 1
}

# Model comparison
modelCompare <- AIC

# Output filename
outputFile <- paste0("output-", tModelIndex, "-", betaIndex, "-",
    sigmaIndex, "-", correlationIndex, "-", typeIndex, ".csv")

# Parameter filename
paramFile <- paste0("parameter-", tModelIndex, "-", betaIndex, "-",
    sigmaIndex, "-", correlationIndex, "-", typeIndex, ".rds")

# Verbose mode
verbose <- FALSE

# Number of decimal places
ndec <- 4


###################
## SET SEED
###################
seeds <- seedGenerator(replications, paste0(inputDir, "/seeds.csv"))


###################
## SIMULATION
###################
simulator(replications, timesteps, models, k, tModel,
    nRay, nRob, nBob, nNel, betas, sampleSize, correlation, sigma, 
    modelCompare, inputDir, outputDir, outputFile, paramFile,
    verbose, ndec, seeds)
