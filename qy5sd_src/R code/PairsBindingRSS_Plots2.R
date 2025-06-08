# Plots for Pairs-Binding Response-Set-Size (RSS) experiment
# Data for RSS1, subj 9, session 1: 2 separate files; file 1 contains 1 trial from block 6 that was manually deleted (2.9.2018) because otherwise that subject had 1 trial more than all others in one condition

setwd(dirname(rstudioapi::getSourceEditorContext()$path))  # sets the directory of location of this script as the current directory

graphics.off()
rm(list=ls(all=TRUE))
library(Hmisc)
library(gplots)
source("lineplot.ci3.R")
source("plot.credint.R")
source("plotPostKO.R")


plotoption <- 1 # 0: Plots window; 1: R Graphics device; 2: PDF

readdata <- function(experiment) {
  ### Read Data
  dat <- read.table(paste0("PairsRSS", experiment, "_all.dat"), header=F)
  names(dat) <- c("id", "session", "block", "trial", "setsize", 
                  "rsizeList", "rsizeNPL", "tested", "response", "rcat", "rt")
  dat$rsize <- dat$rsizeList + dat$rsizeNPL
  dat$correct <- dat$rcat==1
  dat$intrusion <- dat$rcat==2
  dat$new <- dat$rcat==3
  
  recogdat <- subset(dat, rsize>0)
  recalldat <- subset(dat, rsize==0)
  dat <- list(recogdat,recalldat)
}

cl <- c("red", "blue", "green", "magenta", "orange", "black")
pts <- c(21:25, 21, 0)

if (plotoption == 1) x11(width=7,height=9)
if (plotoption == 2) pdf(paste0("WM.binding.joint.setsize.pdf"),width=7,height=9)
layout(matrix(1:6,3,2,byrow=F))
DV <- c("correct", "intrusion", "new")
Title <- c("Correct", "Other", "New")
Exptitle <- c("Large Pool: ", "Small Pool: ")
Ylim <- list(c(0.5,1), c(0,0.4), c(0, 0.4), c(0, 0.4))

for (experiment in 1:2) {
  
  dat <- readdata(experiment)
  recogdat <- dat[[1]]
  recogagg <- aggregate(cbind(correct, intrusion, new) ~ id+setsize+rsize, data=recogdat, FUN=mean)
  
  Sset <- sort(unique(recogagg$setsize))
  Rset <- sort(unique(recogagg$rsize))
  
  for (dv in DV)
    for (ridx in 1:length(Rset)) {
      r <- Rset[ridx]
      subdat <- subset(recogagg, rsize==r)
      ssize <- sort(unique(subdat$setsize))
      lineplot.ci3(subdat, dv=dv, iv="setsize", id="id", x=1, off=-0.1+0.03*ridx, 
                   xlim=c(0.5,8.5), ylim=Ylim[[which(DV==dv)]], 
                   xlab="Set Size", ylab="P(choice)", pt=pts[ridx], ptcol=cl[ridx], col=cl[ridx], cex=1.3)
      title(paste0(Exptitle[experiment], Title[which(DV==dv)]))
      if (ridx == length(Rset) & dv=="correct") 
        legend(0.5, min(Ylim[[which(DV==dv)]]), c("RSS=2", "RSS=4", "RSS=6", "RSS=8", "RSS=10", "RSS=12"), 
               pch=pts, col=cl, pt.bg=cl, yjust=0)
      if (ridx < length(Rset)) par(new=T)
    }
  
}
if (plotoption == 2) dev.off() 



############ Plot Response Set: List Words vs. New Words


ylim <- c(0.5,1)
Exptitle <- c("Large Pool: ", "Small Pool: ")
if (plotoption == 1) x11(7,7)
if (plotoption == 2) pdf(paste0("WM.binding.joint.RSS.pdf"), width=7, height=7)
layout(matrix(1:4,2,2, byrow=F))

for (experiment in 1:2) {
  
  dat <- readdata(experiment)
  recogdat <- dat[[1]]
  aggdatRecog <- aggregate(correct ~ id+setsize+rsizeList+rsizeNPL, data=recogdat, FUN=mean)
  lgdtxt <- c(expression("RSS"["List"]*"=1"), expression("RSS"["List"]*"=2"), expression("RSS"["List"]*"=4"), 
              expression("RSS"["List"]*"=6"), expression("RSS"["List"]*"=8"))
  
  for (s in 3:length(Sset)) {
    subdat <- subset(aggdatRecog, setsize==Sset[s])
    lineplot.ci3(data=subdat, dv="correct", iv=c("rsizeNPL", "rsizeList"), id="id", 
                 pt=21:25, ptcol=cl, col=cl, cex=1.3, xlim=c(-0.5,4.5), ylim=ylim, xlab=expression("RSS"["New"]), ylab="P(Correct)")
    
    if (experiment == 1 & s == 3) legend(0, min(ylim), lgdtxt, pch=21:(20+length(lgdtxt)), 
                                         col=cl[1:length(lgdtxt)], pt.bg=cl[1:length(lgdtxt)], yjust=0)
    title(paste0(Exptitle[experiment], "Set Size ", Sset[s]))
  }
}

if (plotoption == 2) dev.off()



############# Plot Model Parameters #############################


Setsizes <- sort(unique(recogdat$setsize))
parChainFile <- "PairsRSS.MPTparChains.RData"
chainfilethere <- file.access(parChainFile, mode=0)
if (chainfilethere==0) {
  load(parChainFile)
}

# coordinates for arrows and text for MPT tree diagram
Ax0 <- c(1.0, 1.0, 3.0, 3.0, 5.0, 5.0, 5.0, 5.0, 7.0, 7.0)-1
Ay0 <- c(7.0, 7.0, 4.5, 4.5, 7.0, 7.0, 2.5, 2.5, 1.0, 1.0)+0.3
Ax1 <- c(3.0, 3.0, 5.0, 5.0, 7.0, 7.0, 7.0, 7.0, 9.0, 9.0)-1
Ay1 <- c(9.5, 4.5, 7.0, 2.5, 9.0, 5.8, 4.2, 1.0, 2.5, 0.0)+0.3
#        Pb  1-Pb  Pi  1-Pi  gcor goth gcor 1-gcor goth gnew
Tx <- c(2.1, 2.1, 4.0, 4.0, 5.8, 5.8, 5.8, 5.8, 7.7, 7.7, 
        3.7, 7.7, 7.7, 7.7, 9.7, 9.7)-1
Ty <- c(9.0, 4.8, 6.3, 2.8, 8.8, 5.7, 4.0, 1.2, 2.5, -0.2, 
        9.5, 9.0, 5.8, 4.2, 2.5, 0.0)+0.3
Text <- c("Pb", "1-Pb", "Pi", "1-Pi", expression("1/RSS"["L"]), expression("1-1/RSS"["L"]), 
          "1/RSS", "1-1/RSS", expression("(RSS"["L"]*"-1)/RSS"), expression("1-RSS"["L"]*"/RSS"),
          "Correct", "Correct", "Other", "Correct", "Other", "New")

if (plotoption == 1) x11(7,10)
if (plotoption == 2) pdf(paste0("WM.binding.MPT.pdf"), 7,10)
layout(matrix(c(1,1,2,3,4,5), 3,2, byrow=T))
par(mar=c(1,4,4,2))   #c(bottom, left, top, right), default = c(5,4,4,2)+0.1
plot(0,0, xlim=c(0,10), ylim=c(0,10), type="n", bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
arrows(Ax0, Ay0, Ax1, Ay1, length=0.07)
par(cex=0.8)
text(Tx, Ty, Text)
par(mar=c(4.5,4.5,3,2))   #c(bottom, left, top, right), default = c(5,4,4,2)+0.1
par(cex=0.8)
plot.credint(Setsizes, MuPi, off=0.1, pt=21:22, ptcol=c("black", "red"), col=c("black", "red"),
             ylim=c(0,1), xlab="Set Size", ylab=bquote(paste(mu, "(Pi)")))
plot.credint(Setsizes, MuPb, off=0.1, pt=21:22, ptcol=c("black", "red"), col=c("black", "red"), 
             ylim=c(0,1), xlab="Set Size", ylab=bquote(paste(mu, "(Pb)")))
legend(1,0, c("Large Pool", "Small Pool"), pch=21:22, col=c("black", "red"), pt.bg=c("black", "red"), xjust=0, yjust=0)
par(cex=0.7)
plotPostKO(DeltaPi, xlab=bquote(paste(delta, "(Pi)")), breaks=30)
plotPostKO(DeltaPb, xlab=bquote(paste(delta, "(Pb)")), breaks=30)
if (plotoption == 2) dev.off()


parChainFile <- "PairsRSS.MMMparChains.RData"
chainfilethere <- file.access(parChainFile, mode=0)
if (chainfilethere==0) {
  load(parChainFile)
}

if (plotoption == 1) x11(7,7)
if (plotoption == 2) pdf(paste0("WM.binding.MMM.pdf"), 7,7)
layout(matrix(1:4, 2, 2, byrow=T))
par(mar=c(4.5,4.5,3,2))   #c(bottom, left, top, right), default = c(5,4,4,2)+0.1
par(cex=0.8)
plot.credint(Setsizes, MuA, off=0.1, pt=21:22, ptcol=c("black", "red"), col=c("black", "red"), xlim=c(1, 11), ylim=c(0,1), 
             xlab="Set Size", ylab=bquote(paste(mu, "A | ", mu, "B")), xlabels=as.character(Setsizes))
par(new=T)
plot.credint(max(Setsizes)+2, off=0.1, MuB, pt=21:22, ptcol=c("white", "white"), col=c("black", "red"), 
             xlim=c(1, 11), ylim=c(0,1), xlab="", ylab="", xlabels="All")  
lgdpos <- legend(1,1, c(" Large Pool", " Small Pool"), pch=21:22, col=c("black", "red"), pt.bg=c("black", "red"), xjust=0, yjust=1)
points(x=lgdpos$text$x-0.2, y=lgdpos$text$y, pch=21:22, col=c("black", "red"))
par(new=F)
plot.credint(Setsizes, MuC, off=0.1, pt=21:22, ptcol=c("black", "red"), col=c("black", "red"),
             ylim=c(0,20), xlab="Set Size", ylab=bquote(paste(mu, "C")))
par(cex=0.6)
plotPostKO(DeltaA, xlab=bquote(paste(delta, "A")), breaks=30)
plotPostKO(DeltaC, xlab=bquote(paste(delta, "C")), breaks=30)
if (plotoption == 2) dev.off()


############ Plot K Estimates ######################################3

cl <- c("red", "blue", "green", "magenta", "orange", "black")
#pts <- c(15:20, 0)
Exptitles <- c("Large Pool", "Small Pool")

if (plotoption == 1) x11(10,12)
if (plotoption == 2) pdf(paste0("WM.binding.K.pdf"))
layout(matrix(1:4,2,2, byrow=T))
par(cex=0.75)

for (experiment in 1:2) {
  
  ### Read Data
  dat <- readdata(experiment)
  recogdat <- dat[[1]]
  recalldat <- dat[[2]]
  
  aggdatRecog <- aggregate(correct ~ id+setsize+rsizeList+rsizeNPL, data=recogdat, FUN=mean)
  lgdtxt <- paste0("RSS-List = ", sort(unique(aggdatRecog$rsizeList)))
  
  ########### Computing K #################
  
  aggdatRecog$guessing1 <- 1/(aggdatRecog$rsizeList + aggdatRecog$rsizeNPL)
  aggdatRecog$guessing2 <- 1/aggdatRecog$rsizeList    
  aggdatRecog$Pmem1 <- (aggdatRecog$correct - aggdatRecog$guessing1)/(1 - aggdatRecog$guessing1)
  aggdatRecog$Pmem2 <- (aggdatRecog$correct - aggdatRecog$guessing2)/(1 - aggdatRecog$guessing2)
  aggdatRecog$K1 <- aggdatRecog$Pmem1 * aggdatRecog$setsize
  aggdatRecog$K2 <- aggdatRecog$Pmem2 * aggdatRecog$setsize
  
  
  lineplot.ci3(aggdatRecog, dv="K1", iv=c("setsize", "rsizeList"), id="id", 
               x=1, off=0, xlim=c(0.5,8.5), ylim=c(0,8), main=paste0(Exptitles[experiment], ": Guessing from full RS"),
               xlab="Set Size", ylab="K", pt=21:25, ptcol=cl, col=cl, cex=1.3)
  legend(0.5, 8, lgdtxt, pch=21:(20+length(lgdtxt)), pt.bg=cl[1:length(lgdtxt)], col=cl[1:length(lgdtxt)], yjust=1)
  subdat <- subset(aggdatRecog, !is.nan(K2) & K2>0)
  lineplot.ci3(subdat, dv="K2", iv=c("setsize", "rsizeList"), id="id", 
               x=1, off=0, xlim=c(0.5,8.5), ylim=c(0,8), main=paste0(Exptitles[experiment], ": Guessing from effective RS"),
               xlab="Set Size", ylab="K", pt=21:25, ptcol=cl, col=cl, cex=1.3)
  
  recalldat$guessing2 <- 1/recalldat$setsize
  recalldat$Pmem2 <- (recalldat$correct - recalldat$guessing2)/(1 - recalldat$guessing2) 
  #correction for chance, assuming candidate set = memory list
  recallagg2 <- aggregate(Pmem2 ~ id+setsize, data=recalldat, FUN=mean)
  recallagg2$K2 <- recallagg2$Pmem2 * recallagg2$setsize
  par(new=T)
  lineplot.ci3(recallagg2, dv="K2", iv="setsize", id="id", 
               x=1, off=0, xlim=c(0.5,8.5), ylim=c(0,8), main=paste0(Exptitles[experiment], ": Guessing from effective RS"),
               xlab="Set Size", ylab="K", pt=19, ptcol="white", col="black", cex=1.3)
  legend(0.5, 8, "Recall", pch=19, pt.bg="black", yjust=1)
  
  
}
if (plotoption==2) dev.off()


################ Plot Data and Posterior Predictives ##########################

exptext <- c("Large-Pool", "Small-Pool")
respCatName <- c("Correct", "Other", "New")
titletext = c("Set Size 2", "Set Size 4", "Set Size 6", "Set Size 8")
lgd = c("Data", "Model")
lgd2 = c("RS-New=0", "RS-New=1", "RS-New=2", "RS-New=4")
predcolor <- c("red", "blue")  # for MPT and MMM, respectively
RSSlist <- sort(unique(recogdat$rsizeList))

for (experiment in 1:2) {
  
  predFile <- paste0("PairsRSS.postpred", experiment, ".RData")
  load(predFile)

  for (respCat in 1:3) {
    
    M.data <- meanPropResponse[[respCat]]
    CI.data.u <- CIupperPropResponse[[respCat]]
    CI.data.l <- CIlowerPropResponse[[respCat]]
    if (plotoption == 1) x11(12,9)
    if (plotoption == 2) pdf(paste0("WM.binding.postpred.", experiment, ".", respCatName[respCat] ,".pdf"), width=12, height=9)
    layout(matrix(1:4,2,2, byrow=T))
    par(cex=0.9)
    if (respCat == 1) ylim <- c(0,1.05)
    if (respCat > 1) ylim <- c(0, 0.5)
    #x <- list(1:4, 4.5+(1:4), 9+(1:4), 13.5+(1:4), 18+(1:4))

    
    for (s in 1:4) {
      
      space <- matrix(rep(c(0.5,0,0,0), 5), 4, 5)
      width <- 1
      xlim <- c(0, width*dim(M.data)[1]*dim(M.data)[2]+0.5*dim(M.data)[2])
      barplot2(height=M.data[,,s], width=width, space = space, legend.text=F,
               beside = TRUE, ylim=ylim, xlim=c(0,23), xpd=F,
               density=c(0, 15, 30, 50), angle=45, col=c("white", "black", "black", "black"), 
               prcol="white", border="black",
               names.arg = paste0("RS-L=", RSSlist),
               xlab="", ylab=paste0("Proportion ", respCatName[respCat]),
               plot.ci=TRUE, ci.l=CI.data.l[,,s], ci.u=CI.data.u[,,s], ci.color="black")  
      title(paste0(exptext[experiment], " Experiment, ", titletext[s]))
      if (s == 1) legend(xlim[2],ylim[2]-0.1, legend=lgd2, col=c("white", "black", "black", "black"), 
                         density=c(0, 15, 30, 50), xjust=1)
      
      
      for (MM in 1:2)  {
        
        M.pred <- meanPostPred[[(MM-1)*3+respCat]]
        CI.pred.u <- CIupperPostPred[[(MM-1)*3+respCat]]
        CI.pred.l <- CIlowerPostPred[[(MM-1)*3+respCat]]
        for (r in 1:5) {
          x <- ((r-1)*width*dim(M.data)[1] + 1 + r*0.5):(r*width*dim(M.data)[1] + r*0.5) - 0.5
          par(new=T)
          errbar(x + sign(MM-1.5)*0.2, M.pred[,r,s], yplus=CI.pred.u[,r,s], yminus=CI.pred.l[,r,s], 
                 xlim=xlim, ylim=ylim,add=T,type="p", col=predcolor[MM], errbar.col=predcolor[MM])
        }
        
      } # MM
    } # s
    if (plotoption == 2) dev.off()  
    
  } #respCat
  
} # experiment


