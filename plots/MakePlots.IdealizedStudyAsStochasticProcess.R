rm(list=ls())


n = 10 # sample size for each idealized study
m = 100 # number of idealized studies used for illustration
y  = seq(from = 0, to = 1, length.out = n+1) # The sample CDF is plotted on n+1 to determine the end points

xtilde = 0.5 
# An arbitrary fixed value that the observable can take, 
# to illustrate the random variable aspect of the idealized study when
# the data value is fixed

ecdf.values = rep(NA,n)
# ecdf values for marginal histogram  at xtilde 

# Using a flexible model for the true model generating the data
# M is a beta distribution (i.e., X ~ Beta(alpha, beta)) 
# where for each idealized study the parameters alpha and beta
# are uniformly chosen on (0,2)

alpha.max = 2 # maximum value that alpha parameter can take
beta.max = 2 # maximum value that beta parameter can take


alpha=runif(m,0, alpha.max) # generate alphas for all idealized studies 
beta=runif(m,0, beta.max) # generate betas for all idealized studies


# generate all samples
x = matrix(NA,m,n)
ECDF = vector(mode="list",length=m)

for( i in 1:m){
x[i,] = sort(rbeta(n, alpha[i], beta[i]))
ECDF[[i]]  = stepfun(x[i,], y, right=T)
# CONTINUE FROM HERE ECDF OR Y VALUES NEED TO BE PULLED FOR EACH ITERATION
ind = findInterval(xtilde,ECDF[i,]) # get the index of x[i,] that is smaller than xtilde
if(ind==0){
  ecdf.value[i]=0}
else{
ecdf.values[i] = x[i,ind]
}# save the value
}


plot(ECDF[[1]], col='blue', main='Idealized study as a stochastic process',
     xlab = 'Values of observable variable',
     ylab = 'Realizations of idealized study',
     do.points = F, xlim = c(0,1), ylim = c(0,1))


for( i in 2: m-1){
  lines(ECDF[[i]], col= 'lightgray', do.points = F)
}


lines(ECDF[[m]], col='red', do.points = F)

abline(v = xtilde,lty=2)
axis(side = 1, at = xtilde, expression(tilde(x)))



