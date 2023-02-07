################
##
## @description 
# Network plots to show probability of transition
# to a model given another model (Markov chain transition
# probabilities in model-centric approach of Devezer et al. 2018.)
# Model: Markov chain (second order reduced to first order, no learning) 
## @param None
##
## @return None
##
## @lastChange 2018-05-26
##
## @changes
##
################
# begin all
#-------------------------------------------------------------
rm(list=ls())

library("qgraph")
library(permute)
library(matrixStats)
library(MCMCpack)
library(caTools)

#############
## PATHS
#############
# Edit full path to the base directory of CRUSTv1
baseDir <- "/home/erkan/Desktop/lab/Reproducibility/TheoryOfReproducibility/plots" 
scriptDir <- paste0(baseDir, "/Replicator/code") # code directory
outputDir <- paste0(baseDir, "/Replicator/data") # data directory
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
source(paste0(scriptDir, "/modelSimilarByTerm.R"))
source(paste0(scriptDir, "/modelSimilarByInteraction.R"))
source(paste0(scriptDir, "/ReAssign.R"))
###################
## INPUT PARAMETERS
###################
## Number of factors
k <- 3

## Sigma (Error variance)
sigma <- 0.2 #(1:4) error variance to model expectation

## Sample size
sampleSize <- 100

models <- generateModels(k)
## Generate all models
models[[1]]=t(as.matrix(models[[1]]))

## Number of models
L=length(models)
#-----------------------------------
# Scientist model proposal matrices: element ij is the probability of 
# proposing model j if the current global is model i for the first step of the Markov chain

# NEED TO INCORPORATE A DUMMY VARIABLE THAT DETERMINES SOFT VERSUS HARD
# SCIENTISTS (WHETHER ONLY MODELS CONSISTENT WITH EACH SCIENTIST'S
# STRATEGY ARE PROPOSED)

## Tess matrix
tessMatrix <- matrix(data=0.0, L, L)
for(index in 1:L){
  tessMatrix[index,] <- modelSimilarByTerm(models[[index]], models, mode="all")
  nSimModels <- sum(tessMatrix[index,])
  tessMatrix[index,] <- tessMatrix[index,] * (1 / (nSimModels + 1))
  tessMatrix[index, which(tessMatrix[index,] == 0)] <- (1 - sum(tessMatrix[index,])) / (length(tessMatrix[index,]) - nSimModels)
}

## Bo matrix
boMatrix <- matrix(data=0.0, L, L)
for(index in 1:L){
  boMatrix[index,] <- modelSimilarByInteraction(models[[index]], models, mode="all")
  nSimModels <- sum(boMatrix[index,])
  boMatrix[index,] <- boMatrix[index,] * (1 / (nSimModels + 1))
  boMatrix[index, which(boMatrix[index,] == 0)] <- (1 - sum(boMatrix[index,])) / (length(boMatrix[index,]) - nSimModels)
}

## Mave Matrix
maveMatrix <- matrix(data=0.0, L, L )
for(index in 1:L){
  maveMatrix[index,] <- 1 / length(maveMatrix[index,])
} 
# Rey Matrix
# There is no Rey matrix since Rey replicates
# the experiment immediately preeceding her
# incorporated directly in the calculation of transition matrix
# probabilities



# Example usage of library qgraph
# adj <- matrix(c(
#   0,0.8,0.2,
#   0.5,0,0.5,
#   0,0,0),3,3,byrow=TRUE)
# qgraph(adj, layout="circle", theme="gray")

par(mfrow=c(1,3))
qgraph(tessMatrix, layout="circle", theme="gray")
qgraph(boMatrix, layout="circle", theme="gray")
qgraph(maveMatrix, layout="circle", theme="gray")

