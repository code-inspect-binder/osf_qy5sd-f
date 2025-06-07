# Function for plotting means with 95% Credibility Interval, offset by off on the x-axis
# chains must have the structure array(0, dim=c(UV1, UV2, chainlength)), so the last dimension = chain, 1st and 2nd dimensions = UVs

plot.credint <- function(xaxis, chains, aggfun=mean, off=0, upper=T, lower=T, xdim=1, ylim=c(0,1), 
                         xlim=NULL, xlab=NULL, ylab=NULL, xlabels=NULL, pt=15:25, col=rep("black",11), ptcol=rep("white",11), ...)

{
library(Hmisc)
source("C:\\Daten\\R\\Bayes\\HDIofMCMC.R")
dimensions <- dim(chains)
#chaindim <- which.max(dimensions)  #finds index of longest dimension, which is (most likely) the dimension of samples in the chain
if (is.null(xlim)) xlim=c(xaxis[1]-0.5*(xaxis[2]-xaxis[1]), xaxis[length(xaxis)]+0.5*(xaxis[2]-xaxis[1]))
if (length(dimensions)<3) loop2 = 1 else loop2 = dimensions[2]
if (length(dimensions)<3) {
  ch <- chains
  chains <- array(chains, dim=c(dimensions[1], 1, dimensions[2]))
  chains[,1,] <- ch
}
  

# compute means
m <- matrix(0,dimensions[1], loop2)
for (i in 1:dimensions[1])  {
  for (j in 1:loop2) {
    m[i,j] <- aggfun(chains[i,j,])
  }
}
    
#plot means

if (!is.null(xlabels)) ax = F else ax = T

if (xdim==1)  {
  for (j in 1:loop2)  {
    plot(xaxis + (j-1)*off, m[,j], xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,
         type="b", pch=pt[j], col=col[j], bg=ptcol[j], axes=ax, ...)  
    par(new=T)
  }
}
if (xdim==2)  {
  for (i in 1:dimensions[1])  {
    plot(xaxis + (i-1)*off, m[i,], xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,
         type="b", pch=pt[i], col=col[i], bg=ptcol[i], axes=ax, ...)  
    par(new=T)
  }
}

#compute and plot credible intervals
HDI <- array(0, dim=c(dimensions[1], loop2, 2))
for (i in 1:dimensions[1])  {
  for (j in 1:loop2) {  
    par(new=T)
    if (xdim == 1) {x <- xaxis[i] + (j-1)*off} else {x <- xaxis[j] + (i-1)*off}
    if (!is.na(var(chains[i,j,])))  {
      HDI[i,j,] = HDIofMCMC( chains[i,j,] , 0.95 )
      errbar(x, m[i,j], yplus=HDI[i,j,1], yminus=HDI[i,j,2], errbar.col = col[j], 
             xlim=xlim, ylim=ylim,add=T,type="l", ...)
    }
  }
}
if (!is.null(xlabels)) {
  box(...)
  axis(2, ...)
  axis(1,xaxis,labels=xlabels, ...)
}

return(HDI)

}
  



