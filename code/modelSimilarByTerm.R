################
##
## @description Generate similar model(s) adding or removing a term
##
## @param model  Model in matrix format
## @param models List of all possible models
## @param mode   "all" returns all possible similar models
##               "random" returns one similar model
##
## @return Similar model(s) adding or removing a term
##
## @lastChange 2017-04-25
##
## @changes
##
################
modelSimilarByTerm <- function(model, models, mode=c("all", "random")){
  
  mode <- match.arg(mode)
  
  if (mode == "all"){
    if(!is.matrix(model)){
      nTerms <- 1
      terms <- which(model == 1)
      model <- t(as.matrix(model))
    } else {
      nTerms <- 0
      terms <- c()
      for(r in 1:nrow(model)){
        if(sum(model[r,]) == 1){
          nTerms <- nTerms + 1
          terms <- c(terms, which(model[r,] == 1))
        }
      }
    }
    
    k <- length(model[1,])
    
    similarModels <- list()
    if(length(setdiff(c(1:k), terms)) > 0){
      values <- setdiff(c(1:k), terms)
      for(index in values){
        newTerm <- rep(0, k)
        newTerm[index] <- 1
        newModel <- rbind(model, newTerm)
        rownames(newModel) <- c()
        similarModels[length(similarModels) + 1] <- list(newModel)
      }
    }
    
    if(length(terms) > 1){
      for(index in terms[terms != 1]){
        newModel <- c()
        for(r in 1:nrow(model)){
          term <- model[r,]
          if(term[index] != 1){
            newModel <- rbind(newModel, term)
          }
        }
        
        rownames(newModel) <- c()
        similarModels[length(similarModels) + 1] <- list(newModel)
      }
    }
    
    result <- array(0, dim=c(1,length(models)))
    for(m in similarModels){
      mIndex <- searchModel(m, models)
      if(mIndex > 0){
        result[mIndex] <- 1
      }
    }
    
    return(result)
  } else if (mode == "random"){
    opAdd <- TRUE
    if(!is.matrix(model)){
      nTerms <- 1
      terms <- which(model == 1)
      
      similarModel <- t(as.matrix(model))
      model <- t(as.matrix(model))
      
      k <- length(model[1,])
    } else {
      nTerms <- 0
      terms <- c()
      for(r in 1:nrow(model)){
        if(sum(model[r,]) == 1){
          nTerms <- nTerms + 1
          terms <- c(terms, which(model[r,] == 1))
        }
      }
      
      k <- length(model[1,])
      
      if((nTerms == k) ||
          ((nTerms > 1) && (runif(1) > 0.5))){
        opAdd <- FALSE
      }
      
      similarModel <- model
    }
    
    if(opAdd){
      values <- setdiff(c(1:k), terms)
      index <- round(runif(1, 1, length(values)))
      newTerm <- rep(0, k)
      newTerm[values[index]] <- 1
      similarModel <- rbind(similarModel, newTerm)
    } else {
      index <- terms[round(runif(1, 1, length(terms)))]
      
      if(index != 1){
        similarModel <- c()
        for(r in 1:nrow(model)){
          term <- model[r,]
          if(term[index] != 1){
            similarModel <- rbind(similarModel, model[r,])
          }
        }
      }
    }
    
    rownames(similarModel) <- c()
    return(similarModel)
  }
}