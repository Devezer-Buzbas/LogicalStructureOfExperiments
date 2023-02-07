################
##
## @description Find the index of a model inmodels
##
## @param model1 Model 1 in matrix format
## @param model2 models as list each entry a model in matrix format
##
## @return Index of model 1 in models
##
## @lastChange 2018-08-05
##
## @changes
##
################
FindModelIndex <- function(model1, models){
  
  if(!is.matrix(model1)){
    model1 <- t(as.matrix(model1))
  }
  
for(i in 1:length(models)){
  if(identical(models[[i]],model1)){
    break
  }
}
return(i)  
}
