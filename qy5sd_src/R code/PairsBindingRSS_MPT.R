#### Read frequency data, fit with MPT

graphics.off()
rm(list=ls(all=TRUE))

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

require(R2jags)
source("plotPostKO.R")
source("HDIofMCMC.R")
source("plot.credint.R")
#source("Bakeman.R")
source("lineplot.ci3.R")

experiment <- 2
model <- 2
plotoption <- 2 # 0: Plots window; 1: R Graphics device; 2: PDF
traceP <- T  # trace predicted probabilities
traceS <- T  # trace subject parameters

modelname <- paste0("PairsRSS.MPT", model, ".txt")
data <- read.table(paste0("PairsRSS", experiment, "_freq.dat"), header=T)
aggdat <- aggregate(fcorrect ~ setsize+setsize+rsizeList+rsizeNPL, data=data, FUN=mean)
nconds <- dim(aggdat)[1]
nsubj <- length(unique(data$id))
N <- unique(data$N)
Setsizes <- sort(unique(aggdat$setsize))
RSSlist <- sort(unique(aggdat$rsizeList))
RSSnpl <- sort(unique(aggdat$rsizeNPL))

Fdata <- array(NA, dim=c(nconds, nsubj, 3))
SSidx <- NA   # consecutive indices for set sizes (not actual set sizes!)
RSSidx <- matrix(NA, nconds, 2)  # consecutive indices for RSS list / RSS new
GcorrPi <- rep(NA, nconds)
GcorrNoPi <- rep(NA, nconds)
GotherNoPi <- rep(NA, nconds)

for (cond in 1:nconds) {
  subdat <- subset(data, setsize==aggdat[cond, "setsize"] & rsizeList==aggdat[cond, "rsizeList"] & rsizeNPL==aggdat[cond, "rsizeNPL"])
  Fdata[cond,,1] <- subdat[,"fcorrect"]
  Fdata[cond,,2] <- subdat[,"fother"]
  Fdata[cond,,3] <- subdat[,"fnpl"]
  SSidx[cond] <- which(Setsizes == aggdat[cond, "setsize"])
  RSSidx[cond,1] <- which(RSSlist == aggdat[cond, "rsizeList"])
  RSSidx[cond,2] <- which(RSSnpl == aggdat[cond, "rsizeNPL"])
  RSlist <- aggdat[cond, "rsizeList"]
  RSnpl <- aggdat[cond, "rsizeNPL"]
  GcorrPi[cond] <- 1/RSlist             # probability of guessing correct response, given item memory
  GcorrNoPi[cond] <- 1/(RSlist+RSnpl)  # probability of guessing correct respnose, given no item memory
  GotherNoPi[cond] <- (RSlist-1)/(RSlist-1+RSnpl) # p. of guessing "other", given no item memory, and correct was not chosen
}


muparameters <- c("muPb", "muPi")
sgparameters <- c("sgPb", "sgPi")
subjectparameters <- c("Pb", "Pi")
parameters <- c(muparameters, sgparameters, "loglik")
if (model == 2) parameters <- c(parameters, "meanPi", "meanPb", "dPi", "dPb")
if (traceS == T) parameters <- c(parameters, subjectparameters)
if (traceP == T) parameters <- c(parameters, "P")
parnames <- tolower(subjectparameters)

dataList = list(k = Fdata, 
                N = N,
                nsubj = nsubj,
                nconds = nconds,
                nsetsize = length(unique(aggdat$setsize)),
                Setsizes = Setsizes,
                SSidx = SSidx,
                GcorrPi = GcorrPi,
                GcorrNoPi = GcorrNoPi,
                GotherNoPi = GotherNoPi
)
if (model == 2) dataList[["Csetsize"]] <- sort(unique(aggdat$setsize))-5 # centered set sizes

################################# RUN THE CHAINS  ##########################


jags2Model <- jags.parallel(data=dataList, inits=NULL, parameters.to.save=parameters, model.file = modelname,
                            n.chains = 3, n.iter = 10000, n.burnin = 5000,
                            n.thin = 1, n.cluster= 3)
codaSamples <- as.mcmc(jags2Model)
mcmcChain = as.matrix( codaSamples )


# Posterior group mean parameters

x11()
layout(matrix(1:8, 4, 2, byrow=T))
for (s in 1:length(Setsizes)) {
  for (p in 1:length(muparameters)) {
    plotPostKO(mcmcChain[,paste0(muparameters[p], "[", s, "]")], 
               xlab=bquote(paste(mu, "(", .(parnames[p]), "-", .(Setsizes[s]), ")")), breaks=30)
  }
}

if (model == 2) {
  x11()
  layout(matrix(1:4, 2, 2, byrow=T))
  plotPostKO(mcmcChain[,"meanPi"], xlab=bquote(paste(mu, "Pi")), breaks=30)
  plotPostKO(mcmcChain[,"meanPb"], xlab=bquote(paste(mu, "Pb")), breaks=30)
  plotPostKO(mcmcChain[,"dPi"], xlab=bquote(paste(delta, "Pi")), breaks=30)
  plotPostKO(mcmcChain[,"dPb"], xlab=bquote(paste(delta, "Pb")), breaks=30)
}



#### Prepare joint plot of both experiments

if (model == 2) {
  
  parChainFile <- "PairsRSS.MPTparChains.RData"
  chainfilethere <- file.access(parChainFile, mode=0)
  if (chainfilethere==0) {
    load(parChainFile)
  } else {
    MuPi <- array(0, dim=c(length(Setsizes), 2, dim(mcmcChain)[1]))
    MuPb <- array(0, dim=c(length(Setsizes), 2, dim(mcmcChain)[1]))
    meanPi <- array(0, dim=c(length(Setsizes), 2, dim(mcmcChain)[1]))
    meanPb <- array(0, dim=c(length(Setsizes), 2, dim(mcmcChain)[1]))
    DeltaPi <- matrix(0, dim(mcmcChain)[1], 2)
    DeltaPb <- matrix(0, dim(mcmcChain)[1], 2)
  }
  
  CSetsizes <- Setsizes - 5
  
  MuLogitPi <- matrix(0, dim(mcmcChain)[1], c(length(Setsizes)))  # posterior samples of group-level mean of Pi on the logit scale
  MuLogitPb <- matrix(0, dim(mcmcChain)[1], c(length(Setsizes)))
  for (s in 1:length(Setsizes)) {
    MuLogitPi[,s] <- mcmcChain[, paste0("meanPi")] + mcmcChain[, paste0("dPi")] * CSetsizes[s]  
    MuLogitPb[,s] <- mcmcChain[, paste0("meanPb")] + mcmcChain[, paste0("dPb")] * CSetsizes[s]
    # group level Pi/Pb for each set size (on logit scale) as predicted from linear effect of (centered) set size
    MuPi[s, experiment, ] <- 1/(1+exp(-MuLogitPi[,s]))
    MuPb[s, experiment, ] <- 1/(1+exp(-MuLogitPb[,s]))
    # group level Pi/Pb for each set size on probability scale
    PiIdx <- NULL
    for (j in 1:nsubj) PiIdx <- c(PiIdx, which(names(as.data.frame(mcmcChain)) == paste0("Pi[", j, ",", s, "]")))
    PbIdx <- NULL
    for (j in 1:nsubj) PbIdx <- c(PbIdx, which(names(as.data.frame(mcmcChain)) == paste0("Pb[", j, ",", s, "]")))    
    meanPi[s,experiment,] <- rowMeans(mcmcChain[, PiIdx])  # means across subjects of posterior means of Pi for each set size 
    meanPb[s,experiment,] <- rowMeans(mcmcChain[, PbIdx])
  }
  DeltaPi[,experiment] <- mcmcChain[, paste0("dPi")]  # posterior samples of slope of Pi over set size
  DeltaPb[,experiment] <- mcmcChain[, paste0("dPb")]
  print(colMeans(DeltaPi>0))  # proportion of posterior samples > 0
  
  save(MuPi, MuPb, meanPi, meanPb, DeltaPi, DeltaPb, file=parChainFile)
  
}


############### Plot Data and Predictions #######################################3

predFile <- paste0("PairsRSS.postpred", experiment, ".RData")
predfilethere <- file.access(predFile, mode=0)
if (predfilethere==0) {
  load(predFile)
} else {
  meanPropResponse <- list()
  CIupperPropResponse <- list()
  CIlowerPropResponse <- list()
  meanPostPred <- list()
  CIupperPostPred <- list()
  CIlowerPostPred <- list()
}


predP <- array(0,dim=c(nsubj,nconds,3))
predN <- array(0,dim=c(nsubj,nconds,3))
predProb <- array(0,dim=c(nsubj,nconds,3))
datP <- array(0,dim=c(nsubj,nconds,3))
datProb <- array(0,dim=c(nsubj,nconds,3))

for (idx in 1:nsubj) {
  for (c in 1:nconds) {
    for (r in 1:3) {
      condition <- 1
      Pchain <- mcmcChain[, paste("P[", as.character(c), ",", as.character(idx), ",", as.character(r), "]", sep="")]
      predP[idx,c,r] <- mean(Pchain)
      predN[idx,c,r] <- predP[idx,c,r]*N
      datP[idx,c,r] <- Fdata[c,idx,r]/N
    }
  }
}

respCatName <- c("Correct", "Other", "New")
for (respCat in 1:3) {
  
  # Bakeman-McArthur 95% confidence intervals
  
  pdat <- datP[,,respCat]   # pick out the response category (1=correct, 2=other, 3=NPL)  
  rmean = rowMeans(pdat)  # each subject's mean across all conditions
  pcdat <- pdat - rmean + mean(rmean)  
  ppred <- predP[,,respCat]
  rmean = rowMeans(ppred)
  pcpred <- ppred - rmean + mean(rmean)
  
  
  M.data <- array(0,dim=c(4,5,4))  # RSS-npl, RSS-list, setsize
  M.pred <- array(NA,dim=c(4,5,4))  # RSS-npl, RSS-list, setsize
  CI.data.u <- array(NA,dim=c(4,5,4))  # RSS-npl, RSS-list, setsize
  CI.data.l <- array(NA,dim=c(4,5,4))  # RSS-npl, RSS-list, setsize
  CI.pred.u <- array(NA,dim=c(4,5,4))  # RSS-npl, RSS-list, setsize
  CI.pred.l <- array(NA,dim=c(4,5,4))  # RSS-npl, RSS-list, setsize
  
  for (c in 1:nconds) {
    M.data[RSSidx[c,2], RSSidx[c,1], SSidx[c]] <- mean(pcdat[,c])
    std <- sd(pcdat[,c])
    CI.data.u[RSSidx[c,2], RSSidx[c,1], SSidx[c]] <- mean(pcdat[,c]) + 1.96*std/sqrt(nsubj)    
    CI.data.l[RSSidx[c,2], RSSidx[c,1], SSidx[c]] <- mean(pcdat[,c]) - 1.96*std/sqrt(nsubj)   
    M.pred[RSSidx[c,2], RSSidx[c,1], SSidx[c]] <- mean(pcpred[,c])
    std <- sd(pcpred[,c])
    CI.pred.u[RSSidx[c,2], RSSidx[c,1], SSidx[c]] <- mean(pcpred[,c]) + 1.96*std/sqrt(nsubj)    
    CI.pred.l[RSSidx[c,2], RSSidx[c,1], SSidx[c]] <- mean(pcpred[,c]) - 1.96*std/sqrt(nsubj)   
  }
  
  meanPropResponse[[respCat]] <- M.data
  meanPropResponse[[respCat]] <- M.data
  CIupperPropResponse[[respCat]] <- CI.data.u
  CIlowerPropResponse[[respCat]] <- CI.data.l
  meanPostPred[[respCat]] <- M.pred
  CIupperPostPred[[respCat]] <- CI.pred.u
  CIlowerPostPred[[respCat]] <- CI.pred.l
  save(meanPropResponse, CIupperPropResponse, CIlowerPropResponse, meanPostPred, CIupperPostPred, CIlowerPostPred, file=predFile)
  
}

