library(qgraph)

# Example usage of library qgraph
 # adj <- matrix(c(
#    0.1,0.7,0.2,
#    0.5,0,0.5,
#    0,0,0),3,3,byrow=TRUE)
# par(mfrow=c(1,3))
 
#  qgraph(adj, layout="circle", theme="gray")

# Experimental conditions from the Toy Example
# Full factorial design
#
# 6 Models (M_A): 
# Binomial, Negative Binomial, Hypergeometric
# Poisson, Exponential, Normal
#
# 2 Data structures (D_s): 
 # Sample sizes n = 30, n =100
 
 # 2 Methods (S_{post}):
 # MLE and Posterior mode (Bayes)
 
 #
 #--------------------------------------------------------------
 # Plot Layout: 8 independent plots
 #--------------------------------------------------------------
 # Number of possible experimental conditions
 m = 6*2*2
model.from = rep(1:m, each = m)
model.to = rep(1:m, m)
#--------------------------------------------
# my colorblind friendly colors
myplum = "#9F0162"
mypluml = "#FF5AAF"
myviolet = "#450270"
myvioletl = "#8400CD"
myopal = "#009175"
myopall = "#00EBC1"
myazure = "#005FCC"
myazurel = "#00C2F9"
mycrims = "#CD022D"
mycrimsl = "#FF6E3A"
mylime = "#004002"
mylimel = "#009503"
#--------------------------------------------
# node symbols and colors
colors = c(myplum, mypluml, myplum, mypluml, myviolet, myvioletl, myviolet, myvioletl, myopal, myopall, 
           myopal, myopall, myazure, myazurel, myazure, myazurel, mycrims, mycrimsl, mycrims, mycrimsl,
           mylime, mylimel, mylime, mylimel)
# EdgeColor = c("#336666")
EdgeColor = c("black")
sh = c(rep("rectangle",2), rep("circle",2))
shapes = rep(sh,12)
#--------------------------------------------------------------
#--------------------------------------------------------------

 # PLOT 1: 
 # Fully open experiment 
 # (3 variables: Open = 3-Open by definition of pi-Open)
 # exact replications
 #
pdf('/home/erkanb/Desktop/lab/Reproducibility/TheoryOfReproducibility/plots/Figure.Opennnes.1.pdf',
    width=5, height=5)
thickness = 10
 temp = thickness*diag(m)
fully.open.width = c(temp)

x <- as.matrix(data.frame(from = model.from , to = model.to, 
                          width = fully.open.width))
qgraph(x, mode = "direct", layout="circle", shape = shapes, colors = colors, edge.color = EdgeColor, labels=FALSE,
       directed = FALSE)
dev.off()
#--------------------------------------------------------------
# PLOT 2: 
# pi = 2 open experiment, M_A and S_{post} open
# D_s randomly chosen with probability 0.5 
pdf('/home/erkanb/Desktop/lab/Reproducibility/TheoryOfReproducibility/plots/Figure.Opennnes.2.pdf',
    width=3, height=3)
model.spost.weight = 0.5*diag(m)*thickness
matrix.block.replace = matrix(0.5*thickness,2,2)
ind = seq(1,m, by=2)
for( i in ind){
model.spost.weight[c(i,i+1),c(i,i+1)] = matrix.block.replace 
}
weights = c(t(model.spost.weight))
x <- as.matrix(data.frame(from = model.from, to = model.to, 
                          width = weights))
qgraph(x, mode = "direct", layout="circle", shape = shapes, colors = colors, edge.color = EdgeColor, labels=FALSE,
       directed = FALSE)
dev.off()
#--------------------------------------------------------------
# PLOT 3: 
# pi = 2 open experiment, M_A and D_s open
# S_post randomly chosen with probability 0.5 
pdf('/home/erkanb/Desktop/lab/Reproducibility/TheoryOfReproducibility/plots/Figure.Opennnes.3.pdf',
    width=3, height=3)
model.ds.weight = 0.5*diag(m)*thickness
matrix.block.replace = matrix(0.5*thickness,2,2)
ind = c(1,2,5,6,9,10,13,14,17,18,21,22)
for( i in ind){
  model.ds.weight[c(i),c(i+2)] = 0.5*thickness
  model.ds.weight[c(i+2),c(i)] = 0.5*thickness 
}
weights = c(t(model.ds.weight))

x <- as.matrix(data.frame(from = model.from, to = model.to, 
                          width = weights))
qgraph(x, mode = "direct", layout="circle", shape = shapes, colors = colors, edge.color = EdgeColor, labels=FALSE,
       directed = FALSE)
dev.off()
#--------------------------------------------------------------
# PLOT 4: 
# pi = 2 open experiment, S_post and D_s open
# M_A randomly chosen with probability 1/6 
pdf('/home/erkanb/Desktop/lab/Reproducibility/TheoryOfReproducibility/plots/Figure.Opennnes.4.pdf',
    width=3, height=3)
spost.ds.weight = matrix(0,m,m)

ind = c(1,5,9,13,17,21)
ind.2 = ind+1
ind.3 = ind+2
ind.4  = ind+3
a = 1/6
spost.ds.weight[ind,ind] = a*thickness 
spost.ds.weight[ind.2,ind.2] = a*thickness 
spost.ds.weight[ind.3,ind.3] = a*thickness 
spost.ds.weight[ind.4,ind.4] = a*thickness 

diag(spost.ds.weight)=a*thickness
weights = c(t(spost.ds.weight))

x <- as.matrix(data.frame(from = model.from, to = model.to, 
                          width = weights))
qgraph(x, mode = "direct", layout="circle", shape = shapes, colors = colors, edge.color = EdgeColor, labels=FALSE,
       directed = FALSE)
dev.off()
#--------------------------------------------------------------
#--------------------------------------------------------------
# PLOT 5: 
# pi = 1 open experiment, M_A open
# S_{post} and D_s are randomly chosen with probability 1/4 
pdf('/home/erkanb/Desktop/lab/Reproducibility/TheoryOfReproducibility/plots/Figure.Opennnes.5.pdf',
    width=3, height=3)
model.weight = matrix(0,m,m)

ind = c(1,2,3,4)
ind.2 = ind+4
ind.3 = ind+8
ind.4  = ind+12
ind.5  = ind+16
ind.6  = ind+20

model.weight[ind,ind] = a*thickness 
model.weight[ind.2,ind.2] = a*thickness 
model.weight[ind.3,ind.3] = a*thickness 
model.weight[ind.4,ind.4] = a*thickness 
model.weight[ind.5,ind.5] = a*thickness 
model.weight[ind.6,ind.6] = a*thickness 

diag(spost.ds.weight)=a*thickness
weights = c(t(spost.ds.weight))

x <- as.matrix(data.frame(from = model.from, to = model.to, 
                          width = weights))
qgraph(x, mode = "direct", layout="circle", shape = shapes, colors = colors, edge.color = EdgeColor, labels=FALSE,
       directed = FALSE)
dev.off()
#--------------------------------------------------------------
# PLOT 6: 
# pi = 1 open experiment, S_{post} open
# M_A and D_s are randomly chosen with probability 1/12
pdf('/home/erkanb/Desktop/lab/Reproducibility/TheoryOfReproducibility/plots/Figure.Opennnes.6.pdf',
    width=3, height=3)
spost.weight = matrix(0,m,m)

ind.mle = c(1,2,5,6,9,10,13,14,17,18,21,22)
ind.bayes = c(3,4,7,8,11,12,15,16,19,20,23,24)
w.ind = 1/length(ind.mle)

spost.weight[ind.mle,ind.mle] = w.ind*thickness 
spost.weight[ind.bayes,ind.bayes] = w.ind*thickness 
diag(spost.weight)=w.ind*thickness
weights = c(t(spost.weight))

x <- as.matrix(data.frame(from = model.from, to = model.to, 
                          width = weights))
qgraph(x, mode = "direct", layout="circle", shape = shapes, colors = colors, edge.color = EdgeColor, labels=FALSE,
       directed = FALSE)
dev.off()
#--------------------------------------------------------------
# PLOT 7: 
# pi = 1 open experiment, D_s open
# M_A and S_{post} are randomly chosen with probability 1/12
pdf('/home/erkanb/Desktop/lab/Reproducibility/TheoryOfReproducibility/plots/Figure.Opennnes.7.pdf',
    width=3, height=3)
ds.weight = matrix(0,m,m)

ind.n30 = c(1,3,5,7,9,11,13,15,17,19,21,23)
ind.n200 = ind.n30+1
w.ind = 1/length(ind.n30)

ds.weight[ind.n30,ind.n30] = w.ind*thickness 
ds.weight[ind.n200,ind.n200] = w.ind*thickness 

diag(ds.weight)=w.ind*thickness
weights = c(t(ds.weight))

x <- as.matrix(data.frame(from = model.from, to = model.to, 
                          width = weights))
qgraph(x, mode = "direct", layout="circle", shape = shapes, colors = colors, edge.color = EdgeColor, labels=FALSE,
       directed = FALSE)
dev.off()
#--------------------------------------------------------------
# PLOT 8: 
# pi = 0 open experiment, 
# Each chosen uniformly randomly with probability 1/24
pdf('/home/erkanb/Desktop/lab/Reproducibility/TheoryOfReproducibility/plots/Figure.Opennnes.8.pdf',
    width=3, height=3)
weights = rep(1/m,m*m)*thickness
x <- as.matrix(data.frame(from = model.from, to = model.to, 
                          width = weights))
qgraph(x, mode = "direct", layout="circle", shape = shapes, colors = colors, edge.color = EdgeColor, labels=FALSE,
       directed = FALSE)
dev.off()
#--------------------------------------------------------------
#--------------------------------------------------------------


# end all