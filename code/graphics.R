################
##
## @description Generate plots
##
## @param None
##
## @return None
##
## @lastChange 2017-01-30
##
## @changes
##   Distribution of global models (01/30/2017)
##   Calculate rate of reproduction (10/10/2016)
##   Distance normalization (08/29/2016)
##   Average distance to the true model (08/27/2016)
##   Show predictors' number (08/27/2016)
##   Number of times the true model is selected (08/27/2016)
##
################
library(data.table)
library(ggplot2)


#############
## PATHS
#############
baseDir <- "/data/Dropbox/Reproducibility"
scriptDir <- paste0(baseDir, "/code")
inputDir <- paste0(baseDir, "/data")
outputDir <- paste0(baseDir, "/data")


#############
## FUNCTIONS
#############
source(paste0(scriptDir, "/constants.R"))


#############
## UPLOAD DATA
#############
# TrueModel-Beta-Sigma-Correlation-Type
inputDir <- "/data/downloads/garbage/reproducibility"
file <- paste0(inputDir, "/output-1-1-1-1-5.csv")
data <- fread(file, sep=";")


###################
## INPUT PARAMETERS
###################
# Replications
replications <- 100

# Length of the simulation
timesteps <- 1000

# Number of factors
k <- 3


#############
## ANALYSIS
#############
# Distance at the end of the simulation run
data$final_global_model_distance[seq(timesteps, timesteps * replications, by=timesteps)]

# Number of times the global model is the true model
data[, list(num=sum(final_global_true_model)), by=replica]

# Distance to the true model when the global model departs from true model
data[which(initial_global_true_model == 1 & final_global_true_model == 0),]$final_global_model_distance

# Rate that RAY replicates the previous result
data[which(strategy == RAY), mean(replicated), by=replica]


#############
## GRAPHICS
#############
## Normalization
smData <- NULL
smDist <- data[,list(maxDist=max(selected_model_distance)), by=replica]
for(r in smDist$replica){
  smData <- rbind(smData,
      data[which(replica == r),
          list(replica,
              distance=selected_model_distance /
                  ifelse(smDist[which(replica == r)]$maxDist == 0,
                      0, smDist[which(replica == r)]$maxDist)),
          by=timestep])
}

igmData <- NULL
igmDist <- data[,list(maxDist=max(initial_global_model_distance)), by=replica]
for(r in igmDist$replica){
  igmData <- rbind(igmData,
      data[which(replica == r),
          list(replica,
              distance=initial_global_model_distance /
                  ifelse(igmDist[which(replica == r)]$maxDist == 0,
                      0, igmDist[which(replica == r)]$maxDist)),
          by=timestep])
}

fgmData <- NULL
fgmDist <- data[,list(maxDist=max(final_global_model_distance)), by=replica]
for(r in fgmDist$replica){
  fgmData <- rbind(fgmData,
      data[which(replica == r),
          list(replica,
              distance=final_global_model_distance /
                  ifelse(fgmDist[which(replica == r)]$maxDist == 0,
                      0, fgmDist[which(replica == r)]$maxDist)),
          by=timestep])
}

## Mean distance of the final global model to the true model over time
aux <- fgmData[,list(mDist=mean(distance)), by=timestep]
ggplot(aux, aes(x=timestep, y=mDist)) +
  geom_line() +
  geom_smooth() +
  xlab("Timestep") + ylab("Distance") +
  ylim(-0.05, 1) +
  theme(axis.title.x = element_text(color='black', size=24, face='bold'),
      axis.title.y = element_text(color='black', size=24, face='bold'),
      axis.text.x = element_text(color='black', size=16, face='bold'),
      axis.text.y = element_text(color='black', size=16, face='bold'),
      axis.line.x = element_line(color='black', size=1, linetype='solid'),
      axis.line.y = element_line(color='black', size=1, linetype='solid'),
      panel.background = element_rect(fill="transparent", color=NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank())

## Distance of a replica of the final global model to the true model over time
REPLICA <- 2
ggplot(fgmData[which(replica == REPLICA),], aes(x=timestep, y=distance)) +
  geom_line() +
  geom_smooth() +
  xlab("Timestep") + ylab("Distance") +
  ylim(-0.05, 1) +
  theme(axis.title.x = element_text(color='black', size=24, face='bold'),
      axis.title.y = element_text(color='black', size=24, face='bold'),
      axis.text.x = element_text(color='black', size=16, face='bold'),
      axis.text.y = element_text(color='black', size=16, face='bold'),
      axis.line.x = element_line(color='black', size=1, linetype='solid'),
      axis.line.y = element_line(color='black', size=1, linetype='solid'),
      panel.background = element_rect(fill="transparent", color=NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank())


## Mean distance of the selected model to the true model over time
aux <- smData[,list(mDist=mean(distance)), by=timestep]
ggplot(aux, aes(x=timestep, y=mDist)) +
  geom_line() +
  geom_smooth() +
  xlab("Timestep") + ylab("Distance") +
  ylim(-0.05, 1) +
  theme(axis.title.x = element_text(color='black', size=24, face='bold'),
      axis.title.y = element_text(color='black', size=24, face='bold'),
      axis.text.x = element_text(color='black', size=16, face='bold'),
      axis.text.y = element_text(color='black', size=16, face='bold'),
      axis.line.x = element_line(color='black', size=1, linetype='solid'),
      axis.line.y = element_line(color='black', size=1, linetype='solid'),
      panel.background = element_rect(fill="transparent", color=NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank())

## Distribution of global models selected over time
ggplot(data, aes(x=timestep, fill=factor(final_global_model))) +
  geom_bar() +
  xlab("Timestep") + ylab("Model") +
  theme(axis.title.x = element_text(color="black", size=24, face="bold"),
      axis.title.y = element_text(color="black", size=24, face="bold"),
      axis.text.x = element_text(color="black", size=16, face="bold"),
      axis.text.y = element_blank(),
      axis.line.x = element_line(color="black", size=1, linetype="solid"),
      axis.line.y = element_line(color="black", size=1, linetype="solid"),
      panel.background = element_rect(fill="transparent", color=NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank())

## Distribution of preditors over time
ggplot(data, aes(x=timestep, fill=factor(predictors))) +
  geom_bar() +
  xlab("Timestep") + ylab("Predictor") +
  theme(axis.title.x = element_text(color="black", size=24, face="bold"),
      axis.title.y = element_text(color="black", size=24, face="bold"),
      axis.text.x = element_text(color="black", size=16, face="bold"),
      axis.text.y = element_blank(),
      axis.line.x = element_line(color="black", size=1, linetype="solid"),
      axis.line.y = element_line(color="black", size=1, linetype="solid"),
      panel.background = element_rect(fill="transparent", color=NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank())


# TrueModel-Beta-Sigma-Correlation-Type
inputDir <- "/data/downloads/garbage/reproducibility"
file <- paste0(inputDir, "/output-1-1-1-1-5.csv")

inputDir <- "/data/Dropbox/Reproducibility/data"
file <- paste0(inputDir, "/output-aic-ALL.csv")

data <- fread(file, sep=";")

## Plot the distribution of global models
p <- list()
for(repli in 1:5){
  p[[repli]] <- ggplot(data[which(replica == repli),], aes(x=final_global_model)) +
    geom_bar(stat = "count")
}
ggsave(paste0("/data/downloads/garbage/plot-replica", repli, ".png"), plot=p)

x <- data[,.N/length(unique(data$replica)), by=final_global_model]
p <- ggplot(data.table(x), aes(x=final_global_model, y=V1)) +
    geom_bar(stat = "identity")
ggsave(paste0("/data/downloads/garbage/plot-average.png"), plot=p)

times <- seq(1000, 5000, 1000)
data[times,]$final_global_model

m <- array(0, dim=c(100, 14))
for(repli in 1:replications){
  d <- data[which(replica == repli), list(n=.N), by=final_global_model]
  for(r in d$final_global_model){
    m[repli, as.numeric(as.character(r))] <- d[which(as.numeric(as.character(final_global_model)) == as.numeric(as.character(r))),]$n
  }
}

for(i in 1:14){
  m[1, i] <- aic[i,10000]
}

for(i in 1:5){
  print(c(i, mean(m[2:6, i]), sd(m[2:6, i])))
}
