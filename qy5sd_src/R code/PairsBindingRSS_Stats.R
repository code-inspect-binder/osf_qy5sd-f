# Analysis of Large-Pool Experiment (experiment = 1) and Small-Pool Experiment (experiment = 2) of Oberauer, K. (2018). Working memory capacity limits memory for bindings

setwd(dirname(rstudioapi::getSourceEditorContext()$path))  # sets the directory of location of this script as the current directory

graphics.off()
rm(list=ls(all=TRUE))
library(brms)
library(rstan)
library(metRology)
source("C:/daten/R/Toolbox/ConsecutiveID.R")
source("C:/daten/R/Bayes/plotPostKO.R")

experiment <- 2
scales <- c(0.25, 0.353, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0) 
# scales of the Cauchy prior for sensitivity analysis
# On the choice of scale: Gelman et al. (2008) recommend a scale of 2.5, assuming continuous predictors are standardized
# to an SD = 0.5 (and binary predictors to -0.5, 0.5, which also gives SD=0.5). Here we standardize to SD=1, 
# so Gelman's recommendation would translate to scale = 1.25. 
# Here we use a more informative prior: the Cauchy with scale=sqrt(2)/4 results in a distribution of probabilities 
# (after logistic transformation) spanning the range from 0 to 1, with sd = 0.20 and 
# mean absolute deviation from 0.5 = 0.14 (already a fairly large effect)


### Read Data

dat <- read.table(paste0("PairsRSS", experiment, "_all.dat"), header=F)
names(dat) <- c("id", "session", "block", "trial", "setsize", 
                "rsizeList", "rsizeNPL", "tested", "response", "rcat", "rt")
dat$rsize <- dat$rsizeList + dat$rsizeNPL  # total response set size (NPL = not-presented lures = new words)
dat$correct <- dat$rcat==1
dat$intrusion <- dat$rcat==2
dat$new <- dat$rcat==3

recogdat <- subset(dat, rsize>0)   # data from n-AFC tests
recalldat <- subset(dat, rsize==0)  # data from recall tests

# Write out response frequencies for MPT and MMM
recogdat$rfreq <- 1
recogfreq <- aggregate(rfreq ~ id+setsize+rsizeList+rsizeNPL+rcat, data=recogdat, FUN=length, drop=F)
Rfreq1 <- subset(recogfreq, rcat==1)
Rfreq2 <- subset(recogfreq, rcat==2)
Rfreq3 <- subset(recogfreq, rcat==3)
Rfreq <- cbind(Rfreq1, Rfreq2$rfreq, Rfreq3$rfreq)
names(Rfreq)[names(Rfreq)=="rfreq"] <- "fcorrect"
names(Rfreq)[names(Rfreq)=="Rfreq2$rfreq"] <- "fother"
names(Rfreq)[names(Rfreq)=="Rfreq3$rfreq"] <- "fnpl"
Rfreq$N <- rowSums(Rfreq[,c("fcorrect", "fother", "fnpl")])
Rfreq <- subset(Rfreq, N > 0)  # eliminate the combinations of setsize/rsizeList/rsizeNPL that don't exist (but were kept because drop=F)
write.table(Rfreq, file=paste0(paste0("PairsRSS", experiment, "_freq.dat")), row.names=F)


bffile <- "PairsRSS.BF_Sensitivity.RData"
bffilethere <- file.access(bffile, mode=0)
if (bffilethere==0) {
  load(bffile)
} else {
  BayesFactors <- as.data.frame(matrix(NA, length(scales), 6))
  colnames(BayesFactors) <- c("E1: Set Size", "E1: RS-List", "E1: RS-New", "E2: Set Size", "E2: RS-List", "E2: RS-New")
  rownames(BayesFactors) <- as.character(scales)
}


########## Logistic Regression Model of Recog Data with brms  #############

rstan_options (auto_write=TRUE)
options (mc.cores=parallel::detectCores ()) # Run on multiple cores

# aggregate data, z-standardize all predictor variables for the model
recogfreq <- aggregate(cbind(correct, intrusion, new) ~ id+setsize+rsizeList+rsizeNPL, data=recogdat, FUN=sum)
recogfreq$zsetsize <- (recogfreq$setsize - mean(recogfreq$setsize))/sd(recogfreq$setsize)
recogfreq$zrsizeList <- (recogfreq$rsizeList - mean(recogfreq$rsizeList))/sd(recogfreq$rsizeList)
recogfreq$zrsizeNPL <- (recogfreq$rsizeNPL - mean(recogfreq$rsizeNPL))/sd(recogfreq$rsizeNPL)

for (sc in 1:length(scales))  {
  
  prior <- paste0("cauchy(0, ", as.character(scales[sc]), ")")
  fixefPrior <- c(set_prior(prior, class="b"))
  
  ranefPrior <- set_prior("gamma(1,0.04)", class="sd")
  
  # run the models: full model compared to models removing the main effect of interest
  BRMbinom.Corr.full <- brm(correct ~ zsetsize*zrsizeList*zrsizeNPL + 
                              (1+zsetsize+zrsizeList+zrsizeNPL+zsetsize*zrsizeList+
                                 zsetsize*zrsizeNPL+zrsizeList*zrsizeNPL|id),
                            family=binomial(link="logit"),
                            control=list(adapt_delta=0.8),
                            prior=c(fixefPrior, ranefPrior), chains=3, iter=10000, data=recogfreq, save_all_pars=T)
  
  BRMbinom.Corr.noSS <- brm(correct ~ zsetsize*zrsizeList*zrsizeNPL - zsetsize + 
                              (1+zsetsize+zrsizeList+zrsizeNPL+zsetsize*zrsizeList+
                                 zsetsize*zrsizeNPL+zrsizeList*zrsizeNPL|id),
                            family=binomial(link="logit"),
                            control=list(adapt_delta=0.9),
                            prior=c(fixefPrior, ranefPrior), chains=3, iter=10000, data=recogfreq, save_all_pars=T)
  
  BRMbinom.Corr.noRSL <- brm(correct ~ zsetsize*zrsizeList*zrsizeNPL - zrsizeList + 
                               (1+zsetsize+zrsizeList+zrsizeNPL+zsetsize*zrsizeList+
                                  zsetsize*zrsizeNPL+zrsizeList*zrsizeNPL|id),
                             family=binomial(link="logit"),
                             control=list(adapt_delta=0.8),
                             prior=c(fixefPrior, ranefPrior), chains=3, iter=10000, data=recogfreq, save_all_pars=T)
  
  BRMbinom.Corr.noRSN <- brm(correct ~ zsetsize*zrsizeList*zrsizeNPL - zrsizeNPL + 
                               (1+zsetsize+zrsizeList+zrsizeNPL+zsetsize*zrsizeList+
                                  zsetsize*zrsizeNPL+zrsizeList*zrsizeNPL|id),
                             family=binomial(link="logit"),
                             control=list(adapt_delta=0.8),
                             prior=c(fixefPrior, ranefPrior), chains=3, iter=10000, data=recogfreq, save_all_pars=T)
  
  BFbinom.setsize <- bayes_factor(BRMbinom.Corr.full, BRMbinom.Corr.noSS)
  BFbinom.rsList <- bayes_factor(BRMbinom.Corr.full, BRMbinom.Corr.noRSL)
  BFbinom.rsNew <- bayes_factor(BRMbinom.Corr.full, BRMbinom.Corr.noRSN)
  
  BayesFactors[sc, 3*(experiment-1)+1] <- BFbinom.setsize[[1]]
  BayesFactors[sc, 3*(experiment-1)+2] <- BFbinom.rsList[[1]]
  BayesFactors[sc, 3*(experiment-1)+3] <- BFbinom.rsNew[[1]]
  
  save(BayesFactors, file=bffile)
  rm(BRMbinom.Corr.full,  BRMbinom.Corr.noSS, BRMbinom.Corr.noRSL, BRMbinom.Corr.noRSN)
  
}

# look at range of BFs across Cauchy scales
RangeBF <- matrix(0,6,2)
for (j in 1:6) {
  RangeBF[j,1] <- min(BayesFactors[,j])
  RangeBF[j,2] <- max(BayesFactors[,j])
}


