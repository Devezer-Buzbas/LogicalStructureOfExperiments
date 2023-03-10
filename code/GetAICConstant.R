  GetAICConstant <- function(models, xset){
  xset=as.matrix(xset)
A = list()
for (i in 1:length(models)){
  m=models[[i]]
  X=matrix(NA, nrow=nrow(xset), ncol=nrow(m))
      for (j in 1:nrow(m)){
        if (sum(m[j,])==1){
          X[,j]=xset[,which(m[j,]==1)]}else{
          X[,j] = rowProds(as.matrix(xset[,which(m[j,]==1)]))
        }
      }
  A[[i]] = X%*%solve(t(X)%*%X)%*%t(X) # hat matrix
  }
return(A)
}