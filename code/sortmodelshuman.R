sortmodelshuman <-function(models){
L=length(models)
k=ncol(models[[L]])
r=matrix(c(rep(NA,L*k),c(1:L)),nrow=L, ncol=k+1)
#
# first model has only beta1, 
# and is not a matrix, special case
r[1,1:k]=c(1,rep(0,k-1))

for (i in 2:L){
    for (j in 1:k){
    r[i,j]=sum(rowSums(models[[i]])==j)
    }
}
#
for (i in 1:k){  
r = r[order(r[,i]),]
}
return(r)
}