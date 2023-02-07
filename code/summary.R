################
##
## @description Generate summary of multiple simulations
##
## @param None
##
## @return None
##
## @lastChange 2017-09-29
##
## @changes
##   Fixed replicatedTGM and replicatedNTGM (06/18/2017)
##   Add Beta Bias (06/08/2017)
##   Fixed changesBeforeFirstTGM (06/08/2017)
##   Fixed replicatedTGM and replicatedNTGM (09/29/2017)
##
################
library(caTools)
library(data.table)
library(ggplot2)


#############
## PATHS
#############
baseDir <- "/data/Dropbox/Reproducibility"
scriptDir <- paste0(baseDir, "/code/v1")
inputDir <- "/data/downloads/garbage/reproducibility/k4"
outputDir <- paste0(baseDir, "/data")


#############
## FUNCTIONS
#############
source(paste0(scriptDir, "/compareModels.R"))
source(paste0(scriptDir, "/constants.R"))
source(paste0(scriptDir, "/convertBinary.R"))
source(paste0(scriptDir, "/generateModels.R"))
source(paste0(scriptDir, "/modelToStr.R"))
source(paste0(scriptDir, "/searchModel.R"))
source(paste0(scriptDir, "/strToModel.R"))


###################
## INPUT PARAMETERS
###################
## Replications
## Default: 100
replications <- 100

## Length of the simulation
## Default: 10000
timesteps <- 10000

k <- 3

## Number of models
numModels <- 14

## Generate all models
models <- generateModels(k)

tModels <- 1:3
model <- c(4, 10, 5)

betas <- 1:4

sigmas <- 1:2

correlations <- 1:2

types <- 1:5

verbose <- TRUE

###################
## SUMMARY
###################
output <- NULL
outputCode <- NULL
frequency <- NULL
frequencyCode <- NULL
for(tModelIndex in tModels){
  for(betaIndex in betas){
    for(sigmaIndex in sigmas){
      for(correlationIndex in correlations){
        for(typeIndex in types){
          
          if(verbose){
            print(paste0(tModelIndex, "-", betaIndex, "-", sigmaIndex,
                    "-", correlationIndex, "-", typeIndex))
          }
          
          ## True model
          tModel <- modelToStr(models[[model[tModelIndex]]])
          
          ## Generate Betas
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
          
          ## Correlation
          if(sigmaIndex == 1){
            sigma <- 0.2
          } else if(sigmaIndex == 2){
            sigma <- 0.6
          }
          
          ## Correlation
          if(correlationIndex == 1){
            correlation <- 0.2
          } else if(correlationIndex == 2){
            correlation <- 0.8
          }
          
          ## Agent types
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
          
          ## Upload data file
          filename <- paste0("output-", tModelIndex, "-", betaIndex, "-",
              sigmaIndex, "-", correlationIndex, "-", typeIndex, ".csv")
          data <- fread(paste0(inputDir, "/", filename), sep=";")
          
          ## Frequency of each model
          freq <- data[which(timestep <= timesteps), list(num=.N),
              by=list(replica, final_global_model)]
          
          for(r in 1:replications){
            aux <- setdiff(1:numModels,
                unique(freq[which(replica == r)]$final_global_model))
            if(length(aux) > 0){
              freq <- rbind(freq, as.data.table(cbind(replica=r,
                          final_global_model=aux, num=rep(0, length(aux)))))
            }
          }
          freq <- freq[order(replica, final_global_model)]
          
          ## Proportion of times True Model was selected as Global Model
          selTrueModel <- data[which((final_global_true_model == 1) &
                      (timestep <= timesteps)), list(nrows=.N / timesteps),
                  by=replica]
          
          aux <- setdiff(1:replications, unique(selTrueModel$replica))
          if(length(aux) > 0){
            selTrueModel <- rbind(selTrueModel, as.data.table(cbind(replica=aux,
                        nrows=rep(0, length(aux)))))
            selTrueModel <- selTrueModel[order(replica)]
          }
          
          ## Proportion of times the Global Model departs from the True Model
          ## when the Global Model was the True Model
          numerator <- data[which((initial_global_true_model == 1) &
                      (final_global_true_model == 0) &
                      (timestep <= timesteps)), list(nrows=.N), by=replica]
          
          aux <- setdiff(1:replications, unique(numerator$replica))
          if(length(aux) > 0){
            numerator <- rbind(numerator, as.data.table(cbind(
                        replica=aux, nrows=rep(0, length(aux)))))
          }
          numerator <- numerator[order(replica)]
          
          denominator <- data[which((initial_global_true_model == 1) &
                      (timestep <= timesteps)), list(nrows=.N), by=replica]
          
          aux <- setdiff(1:replications, unique(denominator$replica))
          if(length(aux) > 0){
            denominator <- rbind(denominator, as.data.table(cbind(
                        replica=aux, nrows=rep(0, length(aux)))))
          }
          denominator <- denominator[order(replica)]
          
          departTrueModel <- as.data.table(cbind(replica=numerator$replica,
              nrows=numerator$nrows / denominator$nrows))
          
          aux <- setdiff(1:replications, unique(departTrueModel$replica))
          if(length(aux) > 0){
            departTrueModel <- rbind(departTrueModel, as.data.table(cbind(
                        replica=aux, nrows=rep(0, length(aux)))))
          }
          departTrueModel <- departTrueModel[order(replica)]
          
          ## Proportion of times that RAY reproduced previous result
          replicated <- data[which((strategy == RAY) &
                      (timestep <= timesteps)),
              list(num=mean(replicated)), by=replica]
          
          aux <- setdiff(1:replications, unique(replicated$replica))
          if(length(aux) > 0){
            replicated <- rbind(replicated, as.data.table(cbind(
                        replica=aux, nrows=rep(0, length(aux)))))
            replicated <- replicated[order(replica)]
          }
          
          ## First time the Global Model is the True Model
          firstTGM <- data[which((initial_global_true_model == 0) &
                      (final_global_true_model == 1) &
                      (timestep <= timesteps)), .SD, by=replica][,
              list(mints=min(timestep)), by=replica]
          
          aux <- setdiff(1:replications, unique(firstTGM$replica))
          if(length(aux) > 0){
            firstTGM <- rbind(firstTGM, as.data.table(cbind(replica=aux,
                        mints=rep(timesteps, length(aux)))))
            firstTGM <- firstTGM[order(replica)]
          }
          
          ## Number of switches
          numSwitches <- data[which(
            (initial_global_true_model != final_global_true_model) &
              (timestep <= timesteps)), list(nrows=.N), by=replica]
          
          aux <- setdiff(1:replications, unique(numSwitches$replica))
          if(length(aux) > 0){
            numSwitches <- rbind(numSwitches, as.data.table(cbind(replica=aux,
                        nrows=rep(0, length(aux)))))
            numSwitches <- numSwitches[order(replica)]
          }
          
          ## Average and Standard Deviation of contiguous time steps with the
          ## Global Model as True Model
          aux <- data[which(timestep <= timesteps),
              list(final_global_true_model), by=replica]
          runLenTGM <- NULL
          for(i in unique(aux$replica)){
            l <- rle(aux[which(replica == i),]$final_global_true_model)
            runLenTGM <- rbind(runLenTGM, cbind(replica=i,
                    avg=mean(l$lengths[l$values == 1]),
                    sd=ifelse(length(l$lengths[l$values == 1]) > 1,
                        sd(l$lengths[l$values == 1]),
                        0)))
          }
          runLenTGM <- as.data.table(runLenTGM)
          
          aux <- setdiff(1:replications, unique(runLenTGM$replica))
          if(length(aux) > 0){
            runLenTGM <- rbind(runLenTGM, as.data.table(cbind(replica=aux,
                        avg=rep(0, length(aux)),
                        sd=rep(0, length(aux)))))
            runLenTGM <- runLenTGM[order(replica)]
          }
          
          ## Number of changes before the True Model is selected as
          ## Global Model for the first time
          changesBeforeFirstTGM <- NULL
          for(r in 1:replications){
            changesBeforeFirstTGM <- rbind(changesBeforeFirstTGM,
                cbind(r, data[
                  which((initial_global_model != final_global_model) &
                          (replica == r) &
                          (timestep <= firstTGM[which(replica == r),]$mints)),
                  .N]))
          }
          changesBeforeFirstTGM <- data.table(changesBeforeFirstTGM)
          names(changesBeforeFirstTGM) <- c("replica", "num")
          
          aux <- setdiff(1:replications, unique(changesBeforeFirstTGM$replica))
          if(length(aux) > 0){
            changesBeforeFirstTGM <- rbind(changesBeforeFirstTGM,
                as.data.table(cbind(replica=aux, num=rep(0, length(aux)))))
          }
          changesBeforeFirstTGM <- changesBeforeFirstTGM[order(replica)]
          
          ## Level of reproducibility when the True Model is the Global Model
          numerator <- data[which((initial_global_true_model == 1) &
                                    (strategy == RAY) & (replicated == 1) &
                                    (timestep <= timesteps)),
              list(nrows=.N), by=replica]
          
          aux <- setdiff(1:replications, unique(numerator$replica))
          if(length(aux) > 0){
            numerator <- rbind(numerator, as.data.table(cbind(
                        replica=aux, nrows=rep(0, length(aux)))))
          }
          numerator <- numerator[order(replica)]
          
          denominator <- data[which((initial_global_true_model == 1) &
                                      (strategy == RAY) &
                      (timestep <= timesteps)), list(nrows=.N), by=replica]
          
          aux <- setdiff(1:replications, unique(denominator$replica))
          if(length(aux) > 0){
            denominator <- rbind(denominator, as.data.table(cbind(
                        replica=aux, nrows=rep(0, length(aux)))))
          }
          denominator <- denominator[order(replica)]
          
          replicatedTGM <- as.data.table(cbind(replica=numerator$replica,
                  nrows=ifelse(denominator$nrows > 0,
                      numerator$nrows / denominator$nrows, 0)))
          
          aux <- setdiff(1:replications, unique(replicatedTGM$replica))
          if(length(aux) > 0){
            replicatedTGM <- rbind(replicatedTGM,
                as.data.table(cbind(replica=aux, nrows=rep(0, length(aux)))))
          }
          replicatedTGM <- replicatedTGM[order(replica)]
          
          ## Level of reproducibility when the True Model is not
          ## the Global Model
          numerator <- data[which((initial_global_true_model == 0) &
                      (strategy == RAY) & (replicated == 1) &
                      (timestep <= timesteps)),
              list(nrows=.N), by=replica]
          
          aux <- setdiff(1:replications, unique(numerator$replica))
          if(length(aux) > 0){
            numerator <- rbind(numerator, as.data.table(cbind(
                        replica=aux, nrows=rep(0, length(aux)))))
          }
          numerator <- numerator[order(replica)]
          
          denominator <- data[which((initial_global_true_model == 0) &
                                      (strategy == RAY) &
                                      (timestep <= timesteps)),
                              list(nrows=.N), by=replica]
          
          aux <- setdiff(1:replications, unique(denominator$replica))
          if(length(aux) > 0){
            denominator <- rbind(denominator, as.data.table(cbind(
                        replica=aux, nrows=rep(0, length(aux)))))
          }
          denominator <- denominator[order(replica)]
          
          replicatedNTGM <- as.data.table(cbind(replica=numerator$replica,
                  nrows=ifelse(denominator$nrows > 0,
                      numerator$nrows / denominator$nrows, 0)))
          
          aux <- setdiff(1:replications, unique(replicatedNTGM$replica))
          if(length(aux) > 0){
            replicatedNTGM <- rbind(replicatedNTGM,
                as.data.table(cbind(replica=aux, nrows=rep(0, length(aux)))))
          }
          replicatedNTGM <- replicatedNTGM[order(replica)]
          
          ## Bias Betas
          biasBeta <- data[which(timestep <= timesteps),
                           list(avg=mean(beta1_true - beta1_estimate)),
                           by=replica]
          
          aux <- setdiff(1:replications, unique(biasBeta$replica))
          if(length(aux) > 0){
            biasBeta <- rbind(biasBeta, as.data.table(cbind(
              replica=aux, nrows=rep(0, length(aux)))))
            biasBeta <- biasBeta[order(replica)]
          }
          
          ## Output
          output <- rbind(output, cbind(1:replications,
                  tModel, beta1, betaO, sigma,
                  correlation, nRay, nRob, nBob, nNel,
                  selTrueModel$nrows, departTrueModel$nrows, replicated$num,
                  firstTGM$mints, numSwitches$nrows, runLenTGM$avg,
                  runLenTGM$sd, changesBeforeFirstTGM$num,
                  replicatedTGM$nrows, replicatedNTGM$nrows, biasBeta$avg))
          
          outputCode <- rbind(outputCode, cbind(1:replications, tModelIndex,
                  betaIndex, sigmaIndex, correlationIndex, typeIndex,
                  selTrueModel$nrows, departTrueModel$nrows, replicated$num,
                  firstTGM$mints, numSwitches$nrows, runLenTGM$avg,
                  runLenTGM$sd, changesBeforeFirstTGM$num,
                  replicatedTGM$nrows, replicatedNTGM$nrows, biasBeta$avg))
          
          for(r in 1:replications){
            frequency <- rbind(frequency, cbind(r, tModel, beta1, betaO,
                    sigma, correlation, nRay, nRob, nBob, nNel,
                    array(freq$final_global_model, dim=c(14, replications))[,r],
                    array(freq$num, dim=c(14, replications))[,r]))
            
            frequencyCode <- rbind(frequencyCode, cbind(r, tModelIndex,
                    betaIndex, sigmaIndex, correlationIndex, typeIndex,
                    array(freq$final_global_model, dim=c(14, replications))[,r],
                    array(freq$num, dim=c(14, replications))[,r]))
          }
        }
      }
    }
  }
}

output <- data.table(output)
names(output) <- c("replica", "tModel", "beta1", "betaO", "sigmaIndex",
    "correlation", "nRay", "nRob", "nBob", "nNel",
    "selTrueModel", "departTrueModel", "replicated",
    "firstTGM", "numSwitches", "avgContTGM", "sdContTGM",
    "changesBeforeFirstTGM", "replicatedTGM", "replicatedNTGM", "biasBeta1")

write.table(output, file=paste0(outputDir, "/summary.csv"),
    append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE, sep=";")

outputCode <- data.table(outputCode)
names(outputCode) <- c("replica", "tModelIndex", "betaIndex", "sigmaIndex",
    "correlationIndex", "typeIndex",
    "selTrueModel", "departTrueModel", "replicated",
    "firstTGM", "numSwitches", "avgContTGM", "sdContTGM",
    "changesBeforeFirstTGM", "replicatedTGM", "replicatedNTGM", "biasBeta1")

write.table(outputCode, file=paste0(outputDir, "/summaryCode.csv"),
    append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE, sep=";")

frequency <- data.table(frequency)
names(frequency) <- c("replica", "tModel", "beta1", "betaO", "sigmaIndex",
    "correlation", "nRay", "nRob", "nBob", "nNel",
    "final_global_model", "frequency")

write.table(frequency, file=paste0(outputDir, "/frequency.csv"),
    append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE, sep=";")

frequencyCode <- data.table(frequencyCode)
names(frequencyCode) <- c("replica", "tModelIndex", "betaIndex", "sigmaIndex",
    "correlationIndex", "typeIndex", "final_global_model", "frequency")

write.table(frequencyCode, file=paste0(outputDir, "/frequencyCode.csv"),
    append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE, sep=";")


###################
## PLOTTING
###################
data <- fread(paste0(outputDir, "/summaryCode.csv"), sep=";")

a <- array(0, dim=c(9,2,4,2,2,5))
for(m in tModels){
  for(b in betas){
    for(e in sigmas){
      for(c in correlations){
        for(t in types){
          a[1, m, b, e, c, t] <- data[which((tModelIndex == m) &
                                              (betaIndex == b) &
                                              (sigmaIndex == e) &
                                              (correlationIndex == c) &
                                              (typeIndex == t)),]$selTrueModel
          
          a[2, m, b, e, c, t] <- data[which((tModelIndex == m) &
                                              (betaIndex == b) &
                                              (sigmaIndex == e) &
                                              (correlationIndex == c) &
                                              (typeIndex == t)),]$departTrueModel
          
          a[3, m, b, e, c, t] <- data[which((tModelIndex == m) &
                                              (betaIndex == b) &
                                              (sigmaIndex == e) &
                                              (correlationIndex == c) &
                                              (typeIndex == t)),]$avgFirstTrueGlobalModel
          
          a[4, m, b, e, c, t] <- data[which((tModelIndex == m) &
                                              (betaIndex == b) &
                                              (sigmaIndex == e) &
                                              (correlationIndex == c) &
                                              (typeIndex == t)),]$replicated
          
          a[5, m, b, e, c, t] <- data[which((tModelIndex == m) &
                                              (betaIndex == b) &
                                              (sigmaIndex == e) &
                                              (correlationIndex == c) &
                                              (typeIndex == t)),]$avgNumSwitches
          
          a[6, m, b, e, c, t] <- data[which((tModelIndex == m) &
                                              (betaIndex == b) &
                                              (sigmaIndex == e) &
                                              (correlationIndex == c) &
                                              (typeIndex == t)),]$avgContTGM
          
          a[7, m, b, e, c, t] <- data[which((tModelIndex == m) &
                                              (betaIndex == b) &
                                              (sigmaIndex == e) &
                                              (correlationIndex == c) &
                                              (typeIndex == t)),]$avgChangesBeforeFirstTGM
          
          a[8, m, b, e, c, t] <- data[which((tModelIndex == m) &
                                              (betaIndex == b) &
                                              (sigmaIndex == e) &
                                              (correlationIndex == c) &
                                              (typeIndex == t)),]$replicatedTGM
          
          a[9, m, b, e, c, t] <- data[which((tModelIndex == m) &
                                              (betaIndex == b) &
                                              (sigmaIndex == e) &
                                              (correlationIndex == c) &
                                              (typeIndex == t)),]$replicatedNTGM
        }
      }
    }
  }
}


g <- list()
for(o in 1:9){
  g[[o]] <- data.table(model=integer(0), t=integer(0), v=numeric(0))
  for(m in tModels){
    
    if(m == 1){
      model <- 4
    } else if(m == 2){
      model = 10
    }
    
    for(t in types){
      v <- mean(a[o, m, 1, , , t], na.rm=TRUE)
      g[[o]] <- rbind(g[[o]], as.data.table(cbind(model,t,v)))
    }
  }
  g[[o]] <- data.table(g[[o]])
}

p <- list()
p[[1]] <- ggplot(g[[1]], aes(x=t, y=v * 100,
            shape=as.factor(model))) +
    geom_point(size=4) +
    xlab("Agent Type") +
    ylab("") +
    #ylim(0, 100) +
    scale_shape_manual("Model", values = c(15, 17),
        labels = c("Simple", "Complex")) +
    scale_x_continuous(breaks = c(1, 2, 3, 4, 5),
        labels = c("Ray", "Rob", "Bob", "Nel", "All")) +
    theme(axis.title.x = element_text(color='black', size=14, face='bold',
            margin=margin(t=0.2, unit = "cm")),
        axis.title.y = element_text(color='black', size=20, face='bold',
            margin=margin(r=0.5, unit = "cm")),
        axis.text.x = element_text(color='black', size=16, face='bold'),
        axis.text.y = element_text(color='black', size=16, face='bold'),
        axis.line.x = element_line(color='black', size=1, linetype='solid'),
        axis.line.y = element_line(color='black', size=1, linetype='solid'),
        panel.background = element_rect(fill="transparent", color=NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"),
        #legend.position = "none",
        legend.title = element_text(color="black", size=14, face="bold"),
        legend.text = element_text(color="black", size=12, face="bold"),
        legend.key = element_rect(fill = "white"))

ggsave("/home/gnardin/trueModel.png", p[[1]], width=11, height=8, units="in")

p[[2]] <- ggplot(g[[2]], aes(x=t, y=v * 100,
            shape=as.factor(model))) +
    geom_point(size=4) +
    xlab("Agent Type") +
    ylab("") +
    #ylim(0, 100) +
    scale_shape_manual("Model", values = c(15, 17),
        labels = c("Simple", "Complex")) +
    scale_x_continuous(breaks = c(1, 2, 3, 4, 5),
        labels = c("Ray", "Rob", "Bob", "Nel", "All")) +
    theme(axis.title.x = element_text(color='black', size=14, face='bold',
            margin=margin(t=0.2, unit = "cm")),
        axis.title.y = element_text(color='black', size=20, face='bold',
            margin=margin(r=0.5, unit = "cm")),
        axis.text.x = element_text(color='black', size=16, face='bold'),
        axis.text.y = element_text(color='black', size=16, face='bold'),
        axis.line.x = element_line(color='black', size=1, linetype='solid'),
        axis.line.y = element_line(color='black', size=1, linetype='solid'),
        panel.background = element_rect(fill="transparent", color=NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"),
        #legend.position = "none",
        legend.title = element_text(color="black", size=14, face="bold"),
        legend.text = element_text(color="black", size=12, face="bold"),
        legend.key = element_rect(fill = "white"))

ggsave("/home/gnardin/departTrueModel.png", p[[2]], width=11, height=8, units="in")

p[[3]] <- ggplot(g[[3]], aes(x=t, y=v,
            shape=as.factor(model))) +
    geom_point(size=4) +
    xlab("Agent Type") +
    ylab("") +
    #ylim(0, 100) +
    scale_shape_manual("Model", values = c(15, 17),
        labels = c("Simple", "Complex")) +
    scale_x_continuous(breaks = c(1, 2, 3, 4, 5),
        labels = c("Ray", "Rob", "Bob", "Nel", "All")) +
    theme(axis.title.x = element_text(color='black', size=14, face='bold',
            margin=margin(t=0.2, unit = "cm")),
        axis.title.y = element_text(color='black', size=20, face='bold',
            margin=margin(r=0.5, unit = "cm")),
        axis.text.x = element_text(color='black', size=16, face='bold'),
        axis.text.y = element_text(color='black', size=16, face='bold'),
        axis.line.x = element_line(color='black', size=1, linetype='solid'),
        axis.line.y = element_line(color='black', size=1, linetype='solid'),
        panel.background = element_rect(fill="transparent", color=NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"),
        #legend.position = "none",
        legend.title = element_text(color="black", size=14, face="bold"),
        legend.text = element_text(color="black", size=12, face="bold"),
        legend.key = element_rect(fill = "white"))

ggsave("/home/gnardin/firstTrueModel.png", p[[3]], width=11, height=8, units="in")

p[[4]] <- ggplot(g[[5]], aes(x=t, y=v,
            shape=as.factor(model))) +
    geom_point(size=4) +
    xlab("Agent Type") +
    ylab("") +
    #ylim(0, 100) +
    scale_shape_manual("Model", values = c(15, 17),
        labels = c("Simple", "Complex")) +
    scale_x_continuous(breaks = c(1, 2, 3, 4, 5),
        labels = c("Ray", "Rob", "Bob", "Nel", "All")) +
    theme(axis.title.x = element_text(color='black', size=14, face='bold',
            margin=margin(t=0.2, unit = "cm")),
        axis.title.y = element_text(color='black', size=20, face='bold',
            margin=margin(r=0.5, unit = "cm")),
        axis.text.x = element_text(color='black', size=16, face='bold'),
        axis.text.y = element_text(color='black', size=16, face='bold'),
        axis.line.x = element_line(color='black', size=1, linetype='solid'),
        axis.line.y = element_line(color='black', size=1, linetype='solid'),
        panel.background = element_rect(fill="transparent", color=NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"),
        #legend.position = "none",
        legend.title = element_text(color="black", size=14, face="bold"),
        legend.text = element_text(color="black", size=12, face="bold"),
        legend.key = element_rect(fill = "white"))

ggsave("/home/gnardin/avgNumberSwitches.png", p[[4]], width=11, height=8, units="in")
