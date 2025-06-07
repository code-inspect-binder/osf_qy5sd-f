plotPostKO = function( paramSampleVec , credMass=0.95 , compVal=NULL ,
                       HDItextPlace=0.7 , ROPE=NULL , yaxt=NULL , ylab=NULL ,
                       xlab=NULL , cex.lab=NULL , cex=NULL , cex.axis=NULL, xlim=NULL , xrange=NULL,  main=NULL ,
                       col=NULL , border=NULL , showMean=T, showMode=F, centTendY=0.4, DcentTendY=0.25, showCurve=T , 
                       showHDIText=F, breaks=NULL , ... ) {
  # Override defaults of hist function, if not specified by user:
  # (additional arguments "..." are passed to the hist function)
  if ( is.null(xlab) ) xlab="Parameter"
  if ( is.null(cex.lab) ) cex.lab=1.5
  if ( is.null(cex) ) cex=1.4
  if ( is.null(main) ) main=""
  if ( is.null(yaxt) ) yaxt="n"
  if ( is.null(ylab) ) ylab=""
  if ( is.null(col) ) col="black"
  if ( is.null(border) ) border="black"
  if (is.null(xlim)) {
    if (is.null(xrange)) {
      xlim=range( c( compVal , paramSampleVec ) )
    }
    else  {
      center <- min(paramSampleVec) + (max(paramSampleVec) - min(paramSampleVec))/2
      xlim <- c(center-xrange, center+xrange)
    }  
  }
  # in case xlim is set to 0, and some value has been given to xrange, center plot symmetrically on zero, using maximal extension in case range is set to 0, and given range otherwise
  if (xlim[1]==0 & !is.null(xrange))  {
    if (xrange==0)  {
      maxext <- max(abs(min(paramSampleVec)), abs(max(paramSampleVec)))  #largest extension into positive or negative range
      xlim = c(-maxext, maxext)   #centers plot symmetrically on zero
    }
    else {
      xlim = c(-xrange, xrange)
    }
  }
  
  #source("HDIofMCMC.R") 
  source("C:\\Daten\\R\\Bayes\\KruschkeJAGS\\HDIofMCMC.R")
  
  if (length(dim(paramSampleVec)==2)) {
    nPost <- dim(paramSampleVec)[2] 
  } else {
    nPost = 1
    paramSampleVec <- as.matrix(paramSampleVec) #turn it 2-dimensional so that indexing further down still works!
  }
  colors <- c("black", "red", "blue", "green", "magenta", "orange")
  
  if (is.character(xlab)) parName <- xlab else parName <- "Parameter"
  postSummary = matrix( NA , nrow=nPost , ncol=8 , 
                        dimnames=list( rep(parName, nPost) , 
                                       c("mean","median","mode",
                                         "hdiMass","hdiLow","hdiHigh",
                                         "compVal","pcGTcompVal")))     
  
  maxDensity <- matrix(0,1,nPost)
  cenTendHt <- matrix(0,1,nPost)
  for (p in 1:nPost) {
    histinfo = hist( paramSampleVec[,p] , breaks= dim(paramSampleVec)[1]/100, plot=FALSE )
    maxDensity[p] <- max(histinfo$density)
  }
  yrank <- order(maxDensity)
  for (p in 1:nPost) {
    cenTendHt[p] = (centTendY + DcentTendY*yrank[p]) * max(maxDensity)
  }
  if (showMean==T | showMode==T) ylim <- c( 0, max(c(max(maxDensity), max(cenTendHt))) )
  if (showMean==F & showMode==F) ylim <- c( 0, max(maxDensity) )
  
  for (p in 1:nPost) {
    
    postSummary[p,"mean"] = mean(paramSampleVec[,p])
    postSummary[p,"median"] = median(paramSampleVec[,p])
    mcmcDensity = density(paramSampleVec[,p])
    postSummary[p,"mode"] = mcmcDensity$x[which.max(mcmcDensity$y)]
    
    HDI = HDIofMCMC( paramSampleVec[,p] , credMass )
    postSummary[p,"hdiMass"]=credMass
    postSummary[p,"hdiLow"]=HDI[1]
    postSummary[p,"hdiHigh"]=HDI[2]
    
    # Plot histogram.
    if ( is.null(breaks) ) {
      breaks = length(paramSampleVec[,p])/100
#       breaks = c( seq( from=min(paramSampleVec[,p]) , to=max(paramSampleVec[,p]) ,
#                        by=(HDI[2]-HDI[1])/18 ) , max(paramSampleVec[,p]) )
    }
    if ( !showCurve ) {
      par(xpd=NA)
      histinfo = hist( paramSampleVec[,p] , xlab=xlab , yaxt=yaxt , ylab=ylab ,
                       freq=F , border=border , col=col, density=0.1, 
                       xlim=xlim , main=main , cex=cex , cex.lab=cex.lab ,
                       breaks=breaks , ... )
    }
    #cenTendHt = (centTendY + 0.25*p) * max(histinfo$density)
    cvHt = 0.5*max(histinfo$density)
    if ( showCurve ) {
      par(xpd=NA)
      histinfo = hist( paramSampleVec[,p] , plot=FALSE , breaks)
      densCurve = density( paramSampleVec[,p] , adjust=2 )
      densCurve$x[densCurve$x < xlim[1]] <- NA
      densCurve$x[densCurve$x > xlim[2]] <- NA
      plot( densCurve$x , densCurve$y , type="l" , lwd=1 , col=col , bty="n" ,
            xlim=xlim , ylim=ylim, xlab=xlab , yaxt=yaxt , ylab=ylab ,
            main=main , cex=cex , cex.lab=cex.lab , cex.axis=cex.axis, ... )
    }
    #ROPEtextHt = 0.55*max(histinfo$density)
    # Display mean or mode:
    if ( showMean==T ) {
      meanParam = mean( paramSampleVec[,p] )
      text( meanParam , cenTendHt[p] ,
            bquote(mean==.(signif(meanParam,3))) , adj=c(.5,0) , cex=cex )
    } 
    if (showMean==F & showMode==T) {
      dres = density( paramSampleVec[,p] )
      modeParam = dres$x[which.max(dres$y)]
      text( modeParam , cenTendHt[p] ,
            bquote(mode==.(signif(modeParam,3))) , adj=c(.5,0) , cex=cex )
    }
    # Display the comparison value.
    if ( !is.null( compVal ) ) {
      cvCol = "darkgreen"
      # compute percentage of sample values greater than (gt) and less than (lt) the comparison value
      pcgtCompVal = round( 100 * sum( paramSampleVec[,p] > compVal )
                           / length( paramSampleVec[,p] )  , 1 )
      pcltCompVal = 100 - pcgtCompVal
      lines( c(compVal,compVal) , c(0.96*cvHt,0) ,
             lty="dotted" , lwd=1 , col=cvCol )
      text( compVal , cvHt ,
            bquote( .(pcltCompVal)*"% < " *
                     .(signif(compVal,3)) * " < "*.(pcgtCompVal)*"%" ) ,
            adj=c(pcltCompVal/100,0) , cex=0.8*cex , col=cvCol )
      postSummary[p,"compVal"] = compVal
      postSummary[p,"pcGTcompVal"] = ( sum( paramSampleVec[,p] > compVal ) 
                                       / length( paramSampleVec[,p] ) )
    }
    
    # Display the HDI.
    yy <- 0.03 * max(histinfo$density) * (p-1) #elevation of HDI bar
    lines( HDI , c(yy, yy) , lwd=4, col=colors[p] )
    if (showHDIText==T) {
      text( mean(HDI) , 0 , bquote(.(100*credMass) * "% HDI" ) ,
            adj=c(.5,-1.7) , cex=cex )
      text( HDI[1] , 0 , bquote(.(signif(HDI[1],3))) ,
            adj=c(HDItextPlace,-0.5) , cex=cex )
      text( HDI[2] , 0 , bquote(.(signif(HDI[2],3))) ,
            adj=c(1.0-HDItextPlace,-0.5) , cex=cex )
      par(xpd=F)
    }
    
    if (p < nPost) par(new=T)
    
  }
  
  
  
  #
  return( postSummary )
}
