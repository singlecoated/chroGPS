setGeneric("procrustesAdj", function(mds1, d, adjust, sampleid) standardGeneric("procrustesAdj"))

setMethod("procrustesAdj", signature(mds1='mds',d='distGPS'), function(mds1, d, adjust, sampleid) {
  if (length(adjust)!=length(sampleid)) stop('adjust and sampleid must have the same length')
  if (length(adjust)!=nrow(mds1@points)) stop('length(adjust) does not match nrow(mds1@points)')
  # Get contents of distGPS and mds objects
  x <- mds1@points
  d <- d@d
  adjust <- gsub('_','',adjust)
  #If single sample available several times for the same value of the adjustment factor, average them
  g <- paste(sampleid,adjust,sep='_')
  # xmean <- matrix(unlist(by(mds1@points,g,FUN=mean)),byrow=TRUE,ncol=ncol(mds1@points))
  xmean <- matrix(unlist(by(mds1@points,g,FUN=colMeans)),byrow=TRUE,ncol=ncol(mds1@points)) # Use of mean(data.frame) is deprecated
  rownames(xmean) <- unique(g[order(g)])
  adjustmean <- sapply(strsplit(rownames(xmean),split='_'),function(z) z[length(z)])
  allsamples <- unique(sampleid)
  #Set cross-condition distances to NA
  adjid <- unique(adjust)
  for (i in 1:length(adjid)) d[adjust==adjid[i],adjust!=adjid[i]] <- NA; d[adjust!=adjid[i],adjust==adjid[i]] <- NA
  #Repeat adjustment for each pair of conditions
  while (length(unique(adjust))>1) {
    #Find sample x condition matrix
    sampleincond <- do.call(cbind,tapply(sampleid,adjust,function(z) allsamples %in% z))  
    rownames(sampleincond) <- allsamples
    #Find 2 cond with max overlapping samples
    nboth <- t(sampleincond) %*% sampleincond; neach <- diag(nboth); diag(nboth) <- 0
    if (sum(nboth)==0) stop(paste('No common samples found. Cannot use Procrustes.',sep=''))
    #if (sum(nboth)==0) stop(paste('No common samples found between ',sel[1],' and ',sel[2],'. Cannot use Procrustes.',sep=''))
    sel1 <- row(nboth)[which.max(nboth)]; sel2 <- col(nboth)[which.max(nboth)] 
    sel <- c(rownames(nboth)[sel1], rownames(nboth)[sel2])
    if (neach[sel1]<neach[sel2]) sel <- rev(sel)
    #Get data from selected conditions
    insel1 <- adjustmean %in% sel[1]; insel2 <- adjustmean %in% sel[2]
    X1 <- xmean[insel1,]; X2 <- xmean[insel2,]
    n1 <- sub(paste("\\_",sel[1],"$",sep=''),'',rownames(X1))
    n2 <- sub(paste("\\_",sel[2],"$",sep=''),'',rownames(X2))
    idcommon <- n1[n1 %in% n2]
    #Procrustes
    if (length(idcommon)==1) {
      pro <- procrustes(X2[paste(idcommon,sel[2],sep='_'),], X1[paste(idcommon,sel[1],sep='_'),], translation=TRUE, dilation=FALSE, rotation=FALSE)
    } else if (length(idcommon)==2) {
      pro <- procrustes(X2[paste(idcommon,sel[2],sep='_'),], X1[paste(idcommon,sel[1],sep='_'),], translation=TRUE, dilation=TRUE, rotation=FALSE)
    } else {
      pro <- procrustes(X2[paste(idcommon,sel[2],sep='_'),], X1[paste(idcommon,sel[1],sep='_'),], translation=TRUE, dilation=TRUE, rotation=TRUE)
    }
    x[adjust==sel[2],] <- pro$s*x[adjust==sel[2],] %*% pro$R + matrix(pro$tt,nrow=sum(adjust==sel[2]),ncol(x),byrow=TRUE)
    xmean[adjustmean==sel[2],] <- pro$s*xmean[adjustmean==sel[2],] %*% pro$R + matrix(pro$tt,sum(adjustmean==sel[2]),ncol(x),byrow=TRUE)
    adjust[adjust==sel[2]] <- sel[1]
    adjustmean[adjustmean==sel[2]] <- sel[1]
  }
  dapprox <- as.matrix(dist(x,method='euclid'))
  R.square <- cor(d[upper.tri(d)],dapprox[upper.tri(dapprox)],use='complete.obs')
  ans <- new("mds",points=x,R.square=R.square)
  return(ans)
}
)


#Procrustes function modified from MCMCpack
# - Added argument rotation to prevent rotation
# - Adjusted data is not returned, just the Procrustes parameters
procrustes <- function (X, Xstar, translation = FALSE, dilation=FALSE, rotation=FALSE) {
    if (nrow(X) != nrow(Xstar)) {
        cat("X and Xstar do not have same number of rows.\n")
        stop("Check data and call procrustes() again. \n")
    }
    if (ncol(X) != ncol(Xstar)) {
        cat("X and Xstar do not have same number of columns.\n")
        stop("Check data and call procrustes() again. \n")
    }
    n <- nrow(X)
    m <- ncol(X)
    if (translation) { J <- diag(n) - 1/n * matrix(1, n, n) } else { J <- diag(n) }
    C <- t(Xstar) %*% J %*% X
    svd.out <- svd(C)
    if (rotation) { R <- svd.out$v %*% t(svd.out$u) } else { R <- diag(ncol(X)) }
    s <- 1
    if (dilation) {
        mat1 <- t(Xstar) %*% J %*% X %*% R
        mat2 <- t(X) %*% J %*% X
        s.numer <- 0
        s.denom <- 0
        for (i in 1:m) {
            s.numer <- s.numer + mat1[i, i]
            s.denom <- s.denom + mat2[i, i]
        }
        s <- s.numer/s.denom
    }
    tt <- matrix(0, m, 1)
    if (translation) tt <- 1/n * t(Xstar - s * X %*% R) %*% matrix(1, n, 1)
    return(list(R = R, tt = tt, s = s))
}

#Procrustes function copied from MCMCpack package
#     'R', 'tt', and 's' are chosen so that:
# 
#                      s X R + 1 tt' approximately Xstar                 
#     
#     'X.new' is given by:
# 
#                            X.new = s X R + 1 tt'     
#procrustes <- function (X, Xstar, translation = FALSE, dilation = FALSE) {
#    if (nrow(X) != nrow(Xstar)) {
#        cat("X and Xstar do not have same number of rows.\n")
#        stop("Check data and call procrustes() again. \n")
#    }
#    if (ncol(X) != ncol(Xstar)) {
#        cat("X and Xstar do not have same number of columns.\n")
#        stop("Check data and call procrustes() again. \n")
#    }
#    n <- nrow(X)
#    m <- ncol(X)
#    if (translation) { J <- diag(n) - 1/n * matrix(1, n, n) } else { J <- diag(n) }
#    C <- t(Xstar) %*% J %*% X
#    svd.out <- svd(C)
#    R <- svd.out$v %*% t(svd.out$u)
#    s <- 1
#    if (dilation) {
#        mat1 <- t(Xstar) %*% J %*% X %*% R
#        mat2 <- t(X) %*% J %*% X
#        s.numer <- 0
#        s.denom <- 0
#        for (i in 1:m) {
#            s.numer <- s.numer + mat1[i, i]
#            s.denom <- s.denom + mat2[i, i]
#        }
#        s <- s.numer/s.denom
#    }
#    tt <- matrix(0, m, 1)
#    if (translation) tt <- 1/n * t(Xstar - s * X %*% R) %*% matrix(1, n, 1)
#     return(list(R = R, tt = tt, s = s))
#}
