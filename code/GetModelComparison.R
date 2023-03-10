GetModelComparison  <- function(xset, sampleSize, tModel, models, Niter, AICConstant){
# Returns ProbMC, whose entry ij is model j beating model i
# calculated based on Niter number of
# independent samples from the true model tModel
# P(S(M_i)>S(M_j)|tModel),  entry ii equals 1 by convention. 
# ----------------------------------------------- 
L=length(models)
ProbMC = matrix(1,L,L)
#pairwise model indices
cModels <- combs(1:L, 2)
#
## Betas
fModel <- t(as.matrix(models[[1]]))
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
#---------------------------
y=matrix(NA,Niter,sampleSize)
for (i in 1:Niter){
  tModelBetas <- getBetas(tModel, betas, sigma) # betas are weights
  deterministic <- calculateDet(tModel, xset, betas, tModelBetas)
  y[i,] <- generateY(deterministic, sigma)
}
#----------------------------
for (i in 1: nrow(cModels)){
  model1 = cModels[i,1]
  model2 = cModels[i,2]
#  pardif=(nrow(models[[model1]])-nrow(models[[model2]]))*log(sampleSize) # BIC 
  pardif=(nrow(models[[model1]])-nrow(models[[model2]]))*2 # AIC 
#  
  # home brew Schwarz Criterion
  # number of parameters for model1 is nrow(models[[model1]])+1
  # where +1 is for the scale parameter, but cancels in pardif
  DeltaAIC=rep(NA,Niter) 
for (j in 1:Niter){
Y=y[j,]
W=t(Y)
DeltaAIC[j]=sampleSize*(log((1/sampleSize)*(W%*%Y-W%*%AICConstant[[model1]]%*%Y))-
                        log((1/sampleSize)*(W%*%Y-W%*%AICConstant[[model2]]%*%Y)))
}
DeltaAIC=pardif+DeltaAIC 
ProbMC[model2,model1]=mean(DeltaAIC<0) # prob. that model 1 beats model 2
ProbMC[model1,model2]=1-ProbMC[model2,model1]
}
return(ProbMC)
}
# end GetModelComparison