# Classes
setClass("distGPS", representation(d="matrix",metric='character',type='character'))

## SHOW
setMethod("show","distGPS",function(object) {
  cat("Object of class distGPS with", object@metric, "distances between", nrow(object@d), "objects \n")
}
)

## as.matrix
setMethod("as.matrix",signature(x="distGPS"),function(x) {
  as.matrix(x@d)
}
)

setClass("splitDistGPS", representation(d="list",size="numeric",o="numeric",shuffle="numeric"))
 
## SHOW
setMethod("show","splitDistGPS",function(object) {
cat("Object of class splitDistGPS storing", length(object@d),  "distGPS objects of size", object@size, "with", object@o , "anchor points between them \n")
}
)

setMethod("show","splitDistGPS",function(object) {
cat("Object of class splitDistGPS storing", length(object@d),  "distGPS objects of size", object@size, "with", object@o , "anchor points between them \n")
}
)


# Methods and Functions
setGeneric("distGPS", function(x, metric='tanimoto', weights, uniqueRows=FALSE, genomelength=NULL, mc.cores=1) standardGeneric("distGPS"))

setMethod("distGPS", signature(x='RangedDataList'), function(x, metric='tanimoto', weights, uniqueRows=FALSE, genomelength=NULL, mc.cores=1) {
  # if (metric %in% c('tanimoto','avgdist','realtanimoto')) {
  if (metric %in% c('tanimoto','avgdist')) {
     ans <- rdldist(x=x,metric=metric,mc.cores=mc.cores)
  # } else if (metric %in% c('chisq','chi','chisqDist')) {
  } else if (metric %in% c('chi','chisquare')) {
     ans <- chisqdist(x=x,mc.cores=mc.cores)
     if (metric=='chi') ans <- sqrt(ans)
  } else {
     stop('Invalid metric. Available metrics are tanimoto, avgdist, chi and chisquare')
  }
  new("distGPS",d=ans,metric=metric,type='factors')
}
)

### ## Auxiliary routines
### tanipair <- function(z1,z2) { # Not exactly tanimoto, deprecated
###   n1 <- sum(sum(z1 %in% z2))
###   n2 <- sum(sum(z2 %in% z1))
###   (n1+n2)/(nrow(z1)+nrow(z2))
### }

realtanipair <- function(z1,z2) { # REAL tanimoto pair
  # Replaced %in% with over on April 10th 2013
  #n1 <- sum(sum(z1 %in% z2))
  #n2 <- sum(sum(z2 %in% z1))
  n1 <- sum(sum(z1 %over% z2))
  n2 <- sum(sum(z2 %over% z1))
  (.5*(n1+n2))/(nrow(z1)+nrow(z2) - (.5*(n1+n2)))
}

avgdistpair <- function(z1,z2) {
  # Replaced %in% with over on April 10th 2013
  #n1 <- sum(sum(z1 %in% z2))
  #n2 <- sum(sum(z2 %in% z1))
  n1 <- sum(sum(z1 %over% z2))
  n2 <- sum(sum(z2 %over% z1))
  .5*(n1/nrow(z1) + n2/nrow(z2))
}

realtanipair.old <- function(z1,z2) { # REAL tanimoto pair
  n1 <- sum(sum(z1 %in% z2))
  n2 <- sum(sum(z2 %in% z1))
  (.5*(n1+n2))/(nrow(z1)+nrow(z2) - (.5*(n1+n2)))
}

avgdistpair.old <- function(z1,z2) {
  n1 <- sum(sum(z1 %in% z2))
  n2 <- sum(sum(z2 %in% z1))
  .5*(n1/nrow(z1) + n2/nrow(z2))
}

# This is an alternative function for the chisquare distance since 'chisqdist' gave some problems. In the end the problems came from having a session with GenomicRanges loaded after chroGPS
# Currently deprecated
# chisqDist <- function(z1,z2,genomelength)
#   {
#     cnames <- intersect(names(z1),names(z2)) # matching chromosomes
#     a <- z1[cnames]
#     b <- z2[cnames]
#     matchbases <- sum(unlist(lapply(cnames, function(chr) as.numeric(sum(width(reduce(intersect(ranges(a[chr]),ranges(b[chr])))))))))
#     AnoB <- sum(width(reduce(z1))) - matchbases
#     BnoA <- sum(width(reduce(z2))) - matchbases
#     none <- genomelength - matchbases - AnoB - BnoA
#     return(chisq.test(rbind(c(none,BnoA),c(AnoB,matchbases)))$statistic)
#   }

rdldist <- function(x,metric,genomelength=NULL,mc.cores=1) {
# Compute all pairwise Tanimoto distances between all elements in a RangedDataList object x
# - x: RangedDataList object
# - metric: 'avgdist' or 'tanimoto'
# - genomelength: length of the genome to be used in chisqDist. If missing length is auto calculated to fit everything.
# - mc.cores: number of cores to use in parallel computations (passed on to mclapply)
  # Bypass in case IRanges is older than 1.18
  if (installed.packages()['IRanges','Version']<'1.18.0') { avgdistpair <- avgdistpair.old; realtanipair <- realtanipair.old }
  #if (metric=='avgdist') { distfun <- avgdistpair } else if (metric=='tanimoto') { distfun <- tanipair } else if (metric=='realtanimoto') { distfun <- realtanipair }
  if (metric=='avgdist') { distfun <- avgdistpair } else if (metric=='tanimoto') { distfun <- realtanipair }
  # Branch for alternative chisquare metric deprecated
  # else if (metric=='chisquare')
  #   {
  #     if (genomelength==NULL) genomelength = max(unlist(lapply(x,end))) - min(unlist(lapply(x,start))) + 1
  #     distfun <- chisqDist
  #   }
  #Compute distances
  index <- expand.grid(1:length(x),1:length(x))
  index <- index[index[,1]<index[,2],]
  index <- as.list(data.frame(t(index)))
  if (mc.cores>1) {
    if ('parallel' %in% loadedNamespaces()) {
      d <- parallel::mclapply(index,function(z) distfun(x[[z[1]]], x[[z[2]]]), mc.cores=mc.cores, mc.preschedule=TRUE)
    } else stop('parallel library has not been loaded!')
  } else {
    d <- lapply(index,function(z) distfun(x[[z[1]]], x[[z[2]]]))  
  }
  #Build distance matrix
  ans <- matrix(NA,nrow=length(x),ncol=length(x))
  for (i in 1:length(index)) ans[index[[i]][1],index[[i]][2]] <- ans[index[[i]][2],index[[i]][1]] <- 1-d[[i]]
  diag(ans) <- 0
  rownames(ans) <- colnames(ans) <- names(x)
  ans
}

chisqdist2 <- function(x,seqlen=NULL,mc.cores=1)
  {
    if (is.null(seqlen)) seqlen <- sum(apply(do.call(rbind,lapply(x,function(z) do.call(cbind,lapply(names(z),function(t) max(end(z[t])))))),2,max))
    ans <- do.call(rbind,mclapply(x,function(y) lapply(x, function(z) {
      mnames <- intersect(names(y),names(z)) # Matching names
      only.y <- setdiff(names(y),names(z))
      only.z <- setdiff(names(z),names(y))
      d <- sum(as.data.frame(IRanges::intersect(ranges(y[mnames]),ranges(z[mnames])))$width) # Intersect, just for the common chromosomes
      c <- sum(as.data.frame(IRanges::setdiff(ranges(y[mnames]),ranges(z[mnames])))$width) + ifelse(length(only.y)==0,0,sum(width(reduce(y[only.y])))) # A not in B + chromosomes only in A
      b <- sum(as.data.frame(IRanges::setdiff(ranges(z[mnames]),ranges(y[mnames])))$width) + ifelse(length(only.z)==0,0,sum(width(reduce(z[only.z])))) # B not in A + chromosomes only in B
      a <- seqlen - d - c - b # Rest of the genome
      t <- matrix(c(a,b,c,d),nrow=2,ncol=2,byrow=TRUE)
      return(as.numeric(chisq.test(t)$statistic))
    }),mc.cores=mc.cores,mc.preschedule=FALSE))
  return(ans)
  }


chisqdist <- function(x,mc.cores=1) { # Old chisquare calculation, deprecated
# Compute all pairwise chi-square distances between all elements in a RangedDataList object x
# Note: chi-square dist is between objects (i,j) with frequencies (i.e. coverage) ci and cj is defined as
#         sum (ci-cj)^2/ctot
# where ctot is the total coverage across all samples, i.e. ctot= c1+c2+...+cp, where p is the number of samples
  #Set coverage length to maximum needed by all samples
  warn <- options("warn")$warn
  options(warn= -1)
  if (mc.cores>1) {
    if ('parallel' %in% loadedNamespaces()) { l <- do.call(cbind,parallel::mclapply(x,function(z) sapply(z,function(z) ifelse(is.null(z),0,max(end(z)))),mc.cores=mc.cores,mc.preschedule=FALSE)) }}
  else l <- do.call(cbind,lapply(x,function(z) sapply(z,function(z) ifelse(is.null(z),0,max(end(z))))))
  l[l<=0] <- 0
  lmax <- apply(l,1,max)
  #options(warn=warn)
  #Coverage in each sample
  if (mc.cores>1) {
    if ('parallel' %in% loadedNamespaces()) {
      # covsingle <- parallel::mclapply(x,coverage,width=as.list(lmax), mc.cores=mc.cores, mc.preschedule=FALSE)
      covsingle <- parallel::mclapply(x,function(z) coverage(z,width=as.list(lmax[names(lmax) %in% names(z)])), mc.cores=mc.cores, mc.preschedule=FALSE)
    } else stop('parallel library has not been loaded!')
  } else {
    #covsingle <- lapply(x,coverage,width=as.list(lmax))
    covsingle <- lapply(x,function(z) coverage(z,width=as.list(lmax[names(lmax) %in% names(z)])))
  }
  # Add zero coverage for empty chromosomes to each Rle element
  allnames <- sort(unique(unlist(lapply(x,names)))) # Get all chr names
  covsingle <- lapply(covsingle,function(z) {
    nnames <- setdiff(allnames,names(z))
    Rle2add <- RleList(lapply(nnames,function(t) Rle(0,as.numeric(lmax[t]))))
    names(Rle2add) <- nnames
    if (length(nnames)>0) return(c(z,Rle2add)) else return(z)
  })
  #Coverage sum
  #txt <- paste("covsum=",paste(paste('covsingle[[',1:length(covsingle),']]',sep=''),collapse='+'))
  #eval(parse(text=txt))
  txt <- paste("covsum=sum(",paste(paste('covsingle[[',1:length(covsingle),']]',sep=''),collapse=','),')')
  eval(parse(text=txt))
  #Select bases with coverage sum >0
  sel <- covsum>0
  covsum <- covsum[sel]
  covsingle <- lapply(covsingle,function(z) z[sel])
  #Chi-square distances
  nbases <- sapply(x,function(z) sum(width(z))) #nb bases covered by each sample
  covprop <- mapply('/',covsingle,nbases)
  index <- expand.grid(1:length(x),1:length(x))
  index <- index[index[,1]<index[,2],]
  index <- as.list(data.frame(t(index)))
  dchi <- lapply(index, function(z) sum(sum((covprop[[z[1]]] - covprop[[z[2]]])^2/covsum)))
  options(warn=warn)
  #Build distance matrix
  ans <- matrix(NA,nrow=length(x),ncol=length(x))
  for (i in 1:length(index)) ans[index[[i]][1],index[[i]][2]] <- ans[index[[i]][2],index[[i]][1]] <- dchi[[i]]
  diag(ans) <- 0
  ans
}

setMethod("distGPS", signature(x='data.frame'), function(x, metric='tanimoto', weights, uniqueRows=FALSE, genomelength=NULL, mc.cores=1) {
  distGPS(x=as.matrix(x,rownames.force=TRUE), metric=metric, weights=weights, uniqueRows=uniqueRows, mc.cores=mc.cores)
}
)

setMethod("distGPS", signature(x='matrix'), function(x, metric='tanimoto', weights, uniqueRows=FALSE, genomelength=NULL, mc.cores=1) {
  #Check matrix only has 0's and 1's
  if (any(!(x %in% c(0,1)))) stop('x can only contain 0/1 or FALSE/TRUE')
  #Reduce matrix to unique rows
  if (uniqueRows) {
    xUnique <- uniqueCount(x) 
    x <- as.matrix(xUnique[,c(-1,-ncol(xUnique))])
  }
  #Compute distances
  if (metric=='tanimoto') {
    if (missing(weights)) {
      nboth <- x %*% t(x)
      neach <- matrix(rep(rowSums(x),nrow(x)),nrow=nrow(x),ncol=nrow(x))
    } else {
      x2 <- t(t(x)*sqrt(weights)) # Since they are multiplied by each other, we get 1/weights
      nboth <- x2 %*% t(x2)
      neach <- matrix(rep(rowSums(t(t(x)*weights)),nrow(x)),nrow=nrow(x),ncol=nrow(x))
    }
    nunion <- t(neach) + neach - nboth
    ans <- 1-nboth/nunion
    diag(ans) <- rep(0,ncol(ans))
#  } else if (metric=='avgdist.old') { # This metric didnt return a real symmetric matrix, deprecated and fixed in avgdist below
#    nboth <- x %*% t(x)
#    pincluded <- nboth/matrix(diag(nboth),nrow=nrow(nboth),ncol=ncol(nboth)) #percentage of inclusion (asymmetric)
#    pincluded <- .5*(pincluded[upper.tri(pincluded)]+pincluded[lower.tri(pincluded)])
#    ans <- diag(0,nrow(x))
#    ans[upper.tri(ans)] <- ans[lower.tri(ans)] <- 1-pincluded # not symmetric
  } else if (metric=='avgdist') {
    if (missing(weights)) {
      nboth <- x %*% t(x)
    } else {
      x2 <- t(t(x)*sqrt(weights))
      nboth <- x2 %*% t(x2)      
    }
    pincluded <- nboth/matrix(diag(nboth),nrow=nrow(nboth),ncol=ncol(nboth)) #percentage of inclusion (asymmetric)
    pincluded <- .5*(pincluded[upper.tri(pincluded)]+pincluded[lower.tri(pincluded)])
    ans <- diag(0,nrow(x))
    ans[upper.tri(ans)] <- 1-pincluded
    ans[lower.tri(ans)] <- t(ans)[lower.tri(ans)] # symmetric
    # Correct problem with some distances between different points which should not be zero for isoMDS
    ans[ans==0] <- 0.000001
    diag(ans) <- 0
  } else if (metric %in% c('chisquare','chi')) {
    if (!missing(weights)) warning("Cannot use user-specified weights for chi-square based distances. Ignored the weights")
    coln <- colSums(x)
    x <- x[,coln!=0]; coln <- coln[coln!=0]
    x <- x/rowSums(x)  #turn counts to proportions
    x <- t(t(x)/sqrt(coln))  #weight by column marginals
    ans <- as.matrix(dist(x,method='euclid'))
    if (metric=='chi') ans <- sqrt(ans)
  } else if (metric=='wtanimoto') { # old t.sqrt
    if (missing(weights)) {
      coln <- sqrt(colSums(x)) # Wi = 1/sqrt(ngenes with mark i)
    } else {
      coln <- sqrt(colSums(x)/weights)
    }
    x2 <- t(t(x)/sqrt(coln)) # Since they are multiplicated by each other, sqrt((sqrt(x))^2 becomes sqrt(x)
    nboth <- x2 %*% t(x2)
    neach <- matrix(rep(rowSums(t(t(x)/(coln))),nrow(x)),nrow=nrow(x),ncol=nrow(x))
    nunion <- t(neach) + neach - nboth
    ans <- 1-nboth/nunion
  } else if (metric=='t.dsqrt') { # Tanimoto double weighted, deprecated
    if (missing(weights)) {
      coln <- sqrt(colSums(x)) # Wi = 1/sqrt(ngenes with mark i)
    } else {
      coln <- sqrt(colSums(x)/weights)
    }
    x2.1 <- t(t(x)/sqrt(coln)) # Since they are multiplicated by each other, sqrt((sqrt(x))^2 becomes sqrt(x)
    nboth.1 <- x2.1 %*% t(x2.1)
    x2.2 <- t(t(1-x)/sqrt(nrow(x)-coln))
    nboth.2 <- x2.2 %*% t(x2.2)
    nboth <- (nboth.1 + nboth.2)/2
    neach <- 1*matrix(rep(rowSums(t(t(x)/(coln))),nrow(x)),nrow=nrow(x),ncol=nrow(x))
    nunion <- t(neach) + neach - nboth
    ans <- 1-nboth/nunion
    diag(ans) <- 0
  } else if (metric=='t.lin') { # Tanimoto linear, deprecated
    if (missing(weights)) {
      coln <- sqrt(colSums(x)) # Wi = 1/sqrt(ngenes with mark i)
    } else {
      coln <- sqrt(colSums(x)/weights)
    }
    x2 <- t(t(x)/coln)
    nboth <- x2 %*% t(x2)
    neach <- matrix(rep(rowSums(t(t(x)/(coln^2))),nrow(x)),nrow=nrow(x),ncol=nrow(x))
    nunion <- t(neach) + neach - nboth
    ans <- 1-nboth/nunion
  } else if (metric=='t.perc') { # Wi = 1/%genes with mark i # Tanimoto percentile, deprecated
    coln <- sqrt(100*colSums(x)/max(colSums(x)))
    x2 <- t(t(x)/coln)
    nboth <- x2 %*% t(x2)
    neach <- matrix(rep(rowSums(t(t(x)/(coln^2))),nrow(x)),nrow=nrow(x),ncol=nrow(x))
    nunion <- t(neach) + neach - nboth
    ans <- 1-nboth/nunion
  } else if (metric=='mahalanobis') { # Mahalanobis distance
    pairdiff <- ICSNP::pair.diff(x)
    s=cov(x)
    pairdiff <- pair.diff(x)
    d <- mahalanobis(pairdiff,center=rep(0,ncol(pairdiff)),cov=s,inverted=FALSE)
    ans <- matrix(0,nrow=nrow(x),ncol=nrow(x))
    ans[lower.tri(ans)] <- d   #fill lower triangle
    ans <- ans + t(ans)  #fill upper triangle as well    
  } else if (metric=='euclidean') {
    if (!missing(weights)) warning("Cannot use user-specified weights for euclidean distance. Ignored the weights")
    ans <- as.matrix(dist(x,method='euclidean'))
  } else if (metric=='manhattan') {
    if (!missing(weights)) warning("Cannot use user-specified weights for manhattan distance. Ignored the weights")
    ans <- as.matrix(dist(x,method='manhattan'))
  } else {
    stop('Invalid metric. Available metrics are avgdist, tanimoto, wtanimoto, chi,  chisquare, mahalanobis, euclidean and manhattan')
  }
  if (uniqueRows) {
    rownames(ans) <- colnames(ans) <- xUnique$u    
  } else {
    rownames(ans) <- colnames(ans) <- rownames(x)
  }
  new("distGPS",d=ans,metric=metric,type='genes')
}
)


uniqueCount <- function(x) {
# Internal function to perform row clustering and MDS
  #Find unique rows, count the appearance of each and paste all columns into a single column
  txt <- paste("u <- paste(",paste("x[,",1:ncol(x),"]",collapse=","),", sep=',')",sep="")
  eval(parse(text=txt))
  n <- data.frame(table(u))
  n$u <- as.character(n$u)
  xunique <- unique(x)
  txt <- paste("u <- paste(",paste("xunique[,",1:ncol(xunique),"]",collapse=","),", sep=',')",sep="")
  eval(parse(text=txt))
  xunique <- data.frame(u=u,xunique,stringsAsFactors=FALSE)
  xunique <- merge(xunique, n, by='u')
  return(xunique)
}

# Deprecated
### collapseMarks <- function(x,FUN='OR') {
### # Internal function to perform column merging biological or technical replicates
###   # Currently accepting merging by OR (gene has mark if any of the replicates has mark), AND (gene has mark only if all replicates have mark)
###   id <- as.character(colnames(x))
###   uid <- sort(unique(id))
###   xx <- do.call(cbind,lapply(uid,function(y) {
###     sel <- which(id==y)
###     if (FUN=='OR') if (length(sel)>1) xx <- as.numeric(rowSums(x[,sel])>0) else xx <- as.numeric(x[,sel]>0)
###     else if (FUN=='AND') if (length(sel)>1) xx <- as.numeric(rowSums(x[,sel])==length(sel)) else xx <- as.numeric(x[,sel]>0)
###     xx
###   }))
###   colnames(xx) <- uid
###   xx
### }

# Functions to paralelize distance calculation, remember that splitDistGPS should be an internal object and not visible to the final user...
setGeneric("splitDistGPS", function(x, metric='tanimoto', split=.5, overlap=0.05, reshuffle=TRUE,set.seed=149,mc.cores=1) { standardGeneric("splitDistGPS") } )

setMethod("splitDistGPS", signature(x='data.frame'), function(x, metric, split=.5, overlap=0.05, reshuffle=TRUE,set.seed=149,mc.cores=1) {
  splitDistGPS(x=as.matrix(x,rownames.force=TRUE), metric, split=.5, overlap=0.05, reshuffle=TRUE,set.seed=149,mc.cores=1)
}
)

setMethod("splitDistGPS", signature(x='matrix'), function(x, metric, split=.5, overlap=0.05, reshuffle=TRUE,set.seed=149,mc.cores=1) {
  split <- floor(nrow(x) * split)
  overlap <- floor(nrow(x) * overlap)
  if (reshuffle) {
    set.seed(set.seed)
    sel <- sample(1:nrow(x),nrow(x),replace=FALSE)
    x <- x[sel,]
  }
  width <- split-overlap
  p <- vector('list',ceiling(nrow(x)/width))
  p[[1]] <- 1:split
  for (i in 2:(length(p)-1)) p[[i]] <- p[[i-1]]+width
  p[[length(p)]] <- (max(p[[length(p)-1]])+1-overlap):nrow(x)
  if (mc.cores>1) {
    if ('parallel' %in% loadedNamespaces())
      d <- parallel::mclapply(p,function(y) distGPS(x[y,],metric=metric),mc.cores=mc.cores,mc.preschedule=FALSE)
    else stop('parallel library has not been loaded!')
  } else d <- lapply(p,function(y) distGPS(x[y,],metric=metric))
  new("splitDistGPS",d=d,size=split,o=overlap,shuffle=sel)
})

setMethod("splitDistGPS", signature(x='distGPS'), function(x, split=.5, overlap=0.05, reshuffle=TRUE,set.seed=149,mc.cores=1) {
  # Splits an already existing distGPS object into N distGPS objects across the main distance matrix diagonal
  x <- x@d
  metric <- d@metric
  split <- floor(nrow(x) * split)
  overlap <- floor(nrow(x) * overlap)
  if (reshuffle) {
      set.seed(set.seed)
      sel <- sample(1:nrow(x),nrow(x),replace=FALSE)
      x <- x[sel,sel]
    }
    width <- split-overlap
    p <- vector('list',ceiling(nrow(x)/width))
    p[[1]] <- 1:split
    for (i in 2:(length(p)-1)) p[[i]] <- p[[i-1]]+width
    p[[length(p)]] <- (max(p[[length(p)-1]])+1-overlap):nrow(x)
    if (mc.cores>1) {
      if ('parallel' %in% loadedNamespaces())
        d <- parallel::mclapply(p,function(y) new("distGPS",d=x[y,y],metric=metric,type='genes'),mc.cores=mc.cores,mc.preschedule=FALSE)
      else stop('parallel library has not been loaded!')
    } else d <- lapply(p,function(y) new("distGPS",d=x[y,y],metric=metric,type='genes'))
    new("splitDistGPS",d=d,size=split,o=overlap,shuffle=sel)
  })


