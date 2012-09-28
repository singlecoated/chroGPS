# auxiliary functions

boostMDS <- function(D, Y, rate=.1, maxit=50, tol=0.001, samplesize=1, verbose=TRUE, scale=FALSE, seed=149, plt=FALSE, mc.cores=1) {
# boostMDS, based on HITMDS High-Throughput Dimensional Scaling (HiT-MDS), see http://dig.ipk-gatersleben.de/hitmds/hitmds.html
# Introduction of a quadratic step size grid search 
# Introduced paralellization and estimation of step size by resampling the original points when big MDS are used
#
# Y = hitmds(D, Y, n_dim, rate, maxit, plt)
#
# Embed dissimilarity matrix D into n_dim -dimensional Euclidean vector space.
#
# Arguments:
# D - source dissimilarity matrix
# Y  - initial target point configuration, leave empty {} for auto-init
# rate - try values {0.1 1 10 100} for optimum embedding, start using rate=1
# maxit - maximum number of iterations
# tol - tolerance for the objective function.
#
# Return:
# Y - non-standardized embedded points (apply stdscatter for standardization)
#
# Watch value output after each iteration. The higher, the better.
# 0 means perfect mismatch of embedded data relations with D
# 1 means perfect reconstruction, i.e. most trustful embedding result.
#
# Severals runs are advisable for selecting optimum embedding results.
#
#
# Author:      Marc Strickert
# Institution: Leibniz-Institute of Crop Plant Research, IPK-Gatersleben
# Date:  Wed Feb 25 13:15:25     2009
# Licence:     GPLv2; not for use in critical applications
  if (scale) Y <- scale(Y,center=TRUE,scale=TRUE)
  n_dim <- ncol(Y)
  if(!is.matrix(D))
    stop('D no matrix!')
  n_data <- nrow(D)
  zers <- which(D == 0)
  n_datainvs <- 1. / (n_data * n_data - length(zers))
  if(n_dim == 1 && length(Y) > 0 && !is.matrix(Y))
    Y <- as.matrix(Y)
  if(length(Y) > 0 && (nrow(Y) != n_data))
    stop('Y matrix has wrong size!')
  # initialize via cmdscale if necessary
  if(length(Y) == 0) cmdscale(D,k=n_dim)
  mn_D <-  sum(D) * n_datainvs
  D <- D - mn_D
  D[zers] <- 0
  mo_D <-  sum(D * D)
  pnt_del <- Y
  T <- D
  # Prepare resampling if necessary
  if (samplesize<1)
    {
      set.seed <- seed
      sel <- sample(1:nrow(D),floor(nrow(D)*samplesize),replace=FALSE)
      cat("\nSampling for ",length(sel),"elements")
    } 
  targetf <- function(r,resample=FALSE) {
    if (resample) {
      D <- D[sel,sel]
      zers <- which(D == 0)
      Y <- Y[sel,]
      n_data <- nrow(D)
      n_datainvs <- 1. / (n_data * n_data - length(zers))
      mo_D <-  sum(D * D)      
      pnt_del <- pnt_del[sel,]
    }
    Ynew <- Y + r * pnt_del / sqrt(abs(pnt_del)+.001)
    T <- as.matrix(dist(Ynew))
    T[zers] <- 0
    mn_T <- sum(T) * n_datainvs
    T <- T - mn_T
    T[zers] <- 0
    mi_T <- sum(T * D)
    mo_T <- sum(T * T)
    f <- 2 / (abs(mi_T) + abs(mo_T))
    mi_T <- mi_T * f
    mo_T <- mo_T * f
    #sqrt(mi_T * mi_T / (mo_D * mo_T * f))
    (mi_T * mi_T / (mo_D * mo_T * f)) # return r square
  }
  i <- 1; fup <- tol+.1
  while(i<=maxit && fup>tol) {
    T <- as.matrix(dist(Y))
    #T <- as.matrix(dist(Y, diag=T, upper=T))
    T[zers] <- 0
    mn_T <- sum(T) * n_datainvs
    T <- T - mn_T
    T[zers] <- 0 # unknowns: zero force
    mi_T <- sum(T * D)
    mo_T <- sum(T * T)
    f <- 2 / (abs(mi_T) + abs(mo_T))
    mi_T <- mi_T * f
    mo_T <- mo_T * f     
    tmpT <- T * mi_T - D * mo_T
    T <- T + (0.1 + mn_T)
    tmpT <- tmpT / T 
    # calc point i update strength 
    for(j in 1:n_dim) {
      tmp <- Y[,rep(j,n_data)]
      tmp <- t(tmp) - tmp
      tmp = tmp * tmpT
      #pnt_del[,j] = apply(tmp,1,sum)
      pnt_del[,j] = rowSums(tmp)
    }
    if (i==1) {
      #cat("\nResampling is",(samplesize<1))
      r2 <- targetf(0,resample=(samplesize<1))
      if (verbose) cat('   Correl','   Step size\n',r2,'\n',sep='')
    }    
    # Grid step size search
    rateseq <- seq(.1*rate,2*rate,length=5)
    #cat("\nResampling is",(samplesize<1))
    if (mc.cores>1) {
      if ('multicore' %in% loadedNamespaces()) r2new <- unlist(multicore::mclapply(rateseq,targetf,resample=(samplesize<1),mc.cores=mc.cores,mc.preschedule=FALSE))
      else stop('multicore library has not been loaded!')
    } else r2new <- unlist(lapply(rateseq,targetf,resample=(samplesize<1)))
    #Consider quadratic step size
    fit <- lm(r2new ~ rateseq + I(rateseq^2))
    ropt <- -coef(fit)[2]/(2*coef(fit)[3])
    rateseq <- c(rateseq, ropt)
    #cat("\nResampling is",(samplesize<1))
    r2new <- c(r2new,targetf(ropt,resample=(samplesize<1)))
    #Update solution, objective function and optimal rate
    fup <- r2new[which.max(r2new)] - r2
    if (fup > 0) {
      r2 <- r2new[which.max(r2new)]
      rate <- rateseq[which.max(r2new)]
      Y <- Y + rate * pnt_del / sqrt(abs(pnt_del)+.001)
      #Y <- Y + rate * i * .25 * (1+ i %% 2) / maxit * pnt_del / sqrt(abs(pnt_del)+.001)
    }
    if (verbose) cat(r2,rate,'\n')
    if (plt) plot(Y)
    i <- i+1
  } # while
  return(Y)
  #return(new("mds",points=Y,R.square=r2)) # Since resampling is allowed, this r2 would be the one for the sample and not for the whole MDS
} # function

setGeneric("mds", function(d,m=NULL,k=2,type='classic',add=FALSE,cor.method='pearson',splitMDS=FALSE,split=0.05,overlap=0.0025,reshuffle=TRUE,set.seed=149,mc.cores=1,...) standardGeneric("mds"))

setMethod("mds", signature=c(d="distGPS",m="missing"), function(d,m=NULL,k=2,type='classic',add=FALSE,cor.method='pearson',splitMDS=FALSE,split=0.05,overlap=0.0025,reshuffle=TRUE,set.seed=149,mc.cores=1,...) {
  if (splitMDS)
    {
      dsplit <- splitDistGPS(d,split=split,overlap=overlap,reshuffle=reshuffle,set.seed=set.seed,mc.cores=mc.cores)
      #ans <- mds(dsplit,k=k,type=type,add=add,splitMDS=FALSE,mc.cores=mc.cores)
      ans <- splitMDS(dsplit,k=k,type=type,plt=FALSE,mc.cores=mc.cores)
      ans <- ans@points
      d <- d@d
      d <- d[rownames(ans),rownames(ans)] # Reindex internally with rownames from m@points for correct calculation of r2
    }
  else {
    d <- d@d
    if (type=='classic') {
      ans <- cmdscale(as.dist(d), k=k, add=add)
      if (add) ans <- ans$points
    } else if (type=='isoMDS') {
      ans <- isoMDS(d, k=k, maxit=50, trace=FALSE)$points
    } else if (type=='boostMDS') {
      cat('\nPerforming boostMDS without initial configuration. Classic MDS will be used\n')
      ans <- boostMDS(d, Y=cmdscale(as.dist(d), k=k, add=add),...)
    } else {
      stop('Available MDS type for a distance matrix are: classic, isoMDS, boostMDS')
    }
  }
  dapprox <- dist(ans, method='euclidean')
  dapprox <- as.matrix(dapprox)
  R.square <- ifelse(nrow(d)==2, 1, cor(d[upper.tri(d)], dapprox[upper.tri(dapprox)], method=cor.method)^2)
  new("mds", points=ans, R.square=R.square)
}
)

setMethod("mds", signature=c(d="distGPS",m="mds"), function(d,m,type='boostMDS',...) {
# boostMDS
  if (type != c('boostMDS')) stop('Only boostMDS type is available if an input MDS configuration is given')
  d@d <- d@d[rownames(m@points),rownames(m@points)] # Reindex internally with rownames from m@points for correct calculation of r2
  ans <- boostMDS(d@d, Y=m@points, ...)
  dapprox <- dist(ans, method='euclidean')
  dapprox <- as.matrix(dapprox)
  R.square <- ifelse(nrow(d@d)==2, 1, cor(d@d[upper.tri(d@d)], dapprox[upper.tri(dapprox)], method=cor.method)^2)
  new("mds", points=ans, R.square=R.square)
}
)

setMethod("mds", signature=c(d="splitDistGPS",m="missing"), function(d,k=2,type='classic',plt=FALSE,mc.cores=1,...) {    
  ans <- splitMDS(d,k=k,type=type,plt=plt,mc.cores=mc.cores)
  ans
}
)

# Functions to paralelize MDS, adapted on april 2012 for returning a mds class object, R.square is the average r2 correlation of all independent MDSs
splitMDS <- function(d,k=2,type='classic',plt=FALSE,mc.cores=1)
{
  if (mc.cores>1) {
    if ('multicore' %in% loadedNamespaces()) m <- multicore::mclapply(d@d,mds,k=k,type=type,splitMDS=FALSE,mc.cores=mc.cores,mc.preschedule=FALSE)
    else stop('multicore library has not been loaded!')
  } else m <- lapply(d@d,mds,k=k,type=type,splitMDS=FALSE)
  # Ahora hay que juntar, es mejor hacerlo secuencialmente o por parejas y luego eliminar los puntos comunes ?
  # Version secuencial, funciona bien
  mpr <- m[[1]]@points # Init procrustes mds
  for (i in 1:(length(m)-1))
    {
      pr <- mergeMDS(tail(mpr,d@o),m[[i+1]]@points[1:d@o,],m[[i+1]]@points[(d@o+1):nrow(m[[i+1]]@points),],scale=FALSE,symmetric=FALSE)
      mpr <- rbind(mpr,pr$Y2rot)
      if (plt) plot(mpr)
    }
  new("mds",points=mpr,R.square=mean(unlist(lapply(m,function(x) x@R.square))))
  #new("mds",points=mpr,R.square=NA)
}

# Procrustes function to match anchor points in different MDSs (Consider using the ProcrustesAdj)
mergeMDS <- function (X, Y, Y2, scale = TRUE, symmetric = FALSE, scores = "sites", ...)
# Procrustes modified to rotate additinal matrix. X and Y are the matching elements, Y2 has also the Y elements no matching X
{
    if (ncol(X) < ncol(Y)) {
        warning("X has fewer axes than Y: X adjusted to comform Y\n")
        addcols <- ncol(Y) - ncol(X)
        for (i in 1:addcols) X <- cbind(X, 0)
    }
    ctrace <- function(MAT) sum(diag(crossprod(MAT)))
    c <- 1
    if (symmetric) {
        X <- scale(X, scale = FALSE)
        Y <- scale(Y, scale = FALSE)
        X <- X/sqrt(ctrace(X))
        Y <- Y/sqrt(ctrace(Y))
        #Y2 <- Y2/sqrt(ctrace(Y2))
    }
    xmean <- apply(X, 2, mean)
    ymean <- apply(Y, 2, mean)
    #y2mean <- apply(Y2, 2, mean)
    if (!symmetric) {
        X <- scale(X, scale = FALSE)
        Y <- scale(Y, scale = FALSE)
        #Y2 <- scale(Y2, scale = FALSE)
    }
    XY <- crossprod(X, Y)
    sol <- svd(XY)
    A <- sol$v %*% t(sol$u)
    if (scale) {
        c <- sum(sol$d)/ctrace(Y)
    }
    Yrot <- c * Y %*% A
    Y2rot <- c * Y2 %*% A
    b <- xmean - t(A %*% ymean)
    R2 <- ctrace(X) + c * c * ctrace(Y) - 2 * c * sum(sol$d)
    reslt <- list(Yrot = Yrot, Y2rot = Y2rot, X = X, ss = R2, rotation = A,
        translation = b, scale = c, symmetric = symmetric, call = match.call())
    reslt$svd <- sol
    class(reslt) <- "procrustes"
    return(reslt)
}

# Jan's function to add quantitative or whatever variable to an existing MDS, currently only for 2D MDS
addVar <- function(mds1,z,plot=TRUE,label='z',pos=3,...) {
  if (ncol(mds1@points)!=2) stop('Currently only available for 2D MDS')
  x <- mds1@points[,1]
  y <- mds1@points[,2]
  #zs <- scale(z,scale=FALSE)
  zs <- z
  lm.out <- lm(zs~-1+x+y)
  b <- lm.out$coefficients
  yhat <- fitted(lm.out)
  fa <- sqrt(var(yhat))
  bstar <- fa*b
  gof <- summary(lm.out)$r.squared
  ans <- list(b=b,bstar=bstar,gof=gof)
  if (plot)
    {
      arrows(0,0,ans$bstar[1],ans$bstar[2],...)
      text(ans$bstar[1],ans$bstar[2],label,pos=pos,...)
    }
  ans
}
