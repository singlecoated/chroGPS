changewidth <- function(z1,z2,logscale) {
  m <- .5*(start(z1)+end(z1))
  if (logscale) z2 <- exp(z2)
  start(z1) <- m-z2/2; end(z1) <- m+z2/2
  return(z1)
}

setGeneric("adjustPeaks", function(x, adjust, sampleid, logscale=TRUE) standardGeneric("adjustPeaks"))

setMethod("adjustPeaks", signature(x='GRangesList'), function(x, adjust, sampleid, logscale=TRUE) {
  if (length(adjust)!=length(sampleid)) stop('adjust and sampleid must have the same length')
  if (length(adjust)!=length(x)) stop('length(adjust) does not match length(x)')
  allsamples <- unique(sampleid)
  sampleincond <- do.call(cbind,tapply(sampleid,adjust,function(z) allsamples %in% z))
  rownames(sampleincond) <- allsamples
  #Log-scale
  if (logscale) logw <- lapply(as.list(x),function(z) log(width(z))) else logw <- lapply(x,function(z) width(z))
  #Find condition with widest peaks
  meanw <- sapply(logw,mean)
  varw <- sapply(logw,var)
  wcond <- tapply(meanw,adjust,FUN=mean)
  wmax <- names(wcond)[which.max(wcond)]
  #Adjust
  cond2adj <- setdiff(colnames(sampleincond),wmax)
  for (i in cond2adj) {
    sel <- rownames(sampleincond)[(sampleincond[,wmax]==TRUE) & (sampleincond[,i]==TRUE)]
    if (length(sel)>0) { #some common elements
      m1 <- mean(meanw[(sampleid %in% sel) & (adjust==wmax)]); v1 <- mean(varw[(sampleid %in% sel) & (adjust==wmax)])
      m2 <- mean(meanw[(sampleid %in% sel) & (adjust==i)]); v2 <- mean(varw[(sampleid %in% sel) & (adjust==i)])
    } else {  #no common elements
      m1 <- mean(meanw[adjust==wmax]); v1 <- mean(varw[adjust==wmax])
      m2 <- mean(meanw[adjust==i]); v2 <- mean(varw[adjust==i])
    }
    v1 <- ifelse(v1<1e-6,1,v1); v2 <- ifelse(v2<1e-6,1,v2)
    logw[adjust==i] <- lapply(logw[adjust==i],function(z) m1+(z-m2)*sqrt(v1)/sqrt(v2))
    for (j in 1:sum(adjust==i)) x[adjust==i][[j]] <- changewidth(x[adjust==i][[j]],logw[adjust==i][[j]],logscale=logscale)
  }
  return(x)
}
)

### # Adjusting peak width to integrate ChIP-Seq data into ChIP-chip maps
### adjustPeaksInd <- function(x,y) # October 2010
### # adjusts width of peaks in rangedData y to match those of rangedData x
###   {
###     x.mean <- mean(width(x))
###     y.mean <- mean(width(y))
###     x.sd <- sd(width(x))
###     y.sd <- sd(width(y))
###     # New width for y will be
###     y.width <- (((width(y)-y.mean)/y.sd)*x.sd)+x.mean
###     y.shift <- y.width - width(y)
###     start(y) <- start(y) - .5*y.shift
###     end(y) <- end(y) + .5*y.shift
###     return(y)
###   }
###  
### adjustPeaks <- function(x,y) # October 2010
### # adjusts width of peaks in rangedData y to match the numbers of x
###   {
###     x.mean <- mean(unlist(lapply(x,function(z) mean(width(z)))))
###     x.sd <- sd(unlist(lapply(x,function(z) mean(width(z)))))
###     y.mean <- mean(width(y))
###     y.sd <- sd(width(y))
###     # New width for y will be
###     y.width <- (((width(y)-y.mean)/y.sd)*x.sd)+x.mean
###     y.shift <- y.width - width(y)
###     start(y) <- start(y) - .5*y.shift
###     end(y) <- end(y) + .5*y.shift
###     return(y)
###   }
