setClass("mds", representation(points="matrix", R.square="numeric"))

## SHOW
setMethod("show","mds",function(object) {
cat("Object of class MDS approximating distances between", nrow(object@points), "objects \n")
cat("R-squared=",round(object@R.square,4), "\n")
}
)

plotContour <- function(x,probContour=.5,labels='',xlim,ylim,labcex=1,col='red',lwd=2,lty=1,...)
          {
            prob <- x$dens/sum(x$dens)
            probo <- prob[order(prob)]
            cutoff <- probo[match(TRUE,cumsum(probo) > 1-probContour)]
            contour(x=x$grid1,y=x$grid2,z=matrix(x$dens,nrow=length(x$grid1),ncol=length(x$grid2)),levels=cutoff*sum(x$dens),labels=labels,xlim=xlim,ylim=ylim,
              labcex=labcex,col=col,lwd=lwd,axes=FALSE,lty=lty,add=TRUE,plot.axes=FALSE)
          }

contour2dDP <- function(x, ngrid, grid=NULL, probContour=0.5,xlim,ylim,labels='',labcex=.01,col=colors()[393],lwd=4,lty=1,contour.type='single',contour.fill=FALSE,...) {
  #Contour of 2-dimensional distribution based on Dirichlet Process density estimation
  if (missing(xlim)) xlim <- range(x[,1])
  if (missing(ylim)) ylim <- range(x[,2])
  #DP density estimate
  s2 <- var(x); m2 <- colMeans(x); psiinv2 <- solve(s2)
  prior <- list(a0=1,b0=1/5,nu1=4,nu2=4,s2=s2,m2=m2,psiinv2=psiinv2,tau1=0.01,tau2=0.01)
  mcmc <- list(nburn=500,nsave=1000,nskip=1,ndisplay=1001)
  fit1 <- DPdensity(y=x,ngrid=ngrid,grid=grid,prior=prior,mcmc=mcmc,state=state,status=TRUE,na.action=na.omit,method='neal')
  #Select set of points with probContour probability
  prob <- fit1$dens/sum(fit1$dens)
  probo <- prob[order(prob)]
  cutoff <- probo[match(TRUE,cumsum(probo) > 1-probContour)]
  # If single or multiple contour, plot it
  if (contour.type=='single') {
    if (contour.fill==FALSE) {
      contour(x=fit1$grid1,y=fit1$grid2,z=matrix(fit1$dens,nrow=length(fit1$grid1),ncol=length(fit1$grid2)),levels=cutoff*sum(fit1$dens),labels=labels,xlim=xlim,ylim=ylim,
              labcex=labcex,col=col,lwd=lwd,axes=FALSE,lty=lty,add=TRUE,plot.axes=FALSE)
    }
    else filled.contour(x=fit1$grid1,y=fit1$grid2,z=matrix(fit1$dens,nrow=length(fit1$grid1),ncol=length(fit1$grid2)),levels=cutoff*sum(fit1$dens),labels=labels,xlim=xlim,ylim=ylim,
                        labcex=labcex,col=col,lwd=lwd,axes=FALSE,lty=lty,plot.axes=FALSE)
    return(fit1)
  }
  else if (contour.type=='multiple') {
    if (contour.fill==FALSE) {
      contour(fit1$grid1,fit1$grid2,z=matrix(fit1$dens,nrow=length(fit1$grid1),ncol=length(fit1$grid2)),xlim=xlim,ylim=ylim,col=col,add=TRUE,plot.axes=FALSE)
    }
    else filled.contour(fit1$grid1,fit1$grid2,z=matrix(fit1$dens,nrow=length(fit1$grid1),ncol=length(fit1$grid2)),xlim=xlim,ylim=ylim,col=col,plot.axes=FALSE)
    return(fit1)
  }
  else if (contour.type=='none') { return(fit1) }
  else warning('Invalid contour type')
}

contour3dDP <- function(x, col='red', probContour=0.5, ngrid=30, contour.type='none') {
  #Multivariate density fit
  xdf <- data.frame(x)
  fit <- ssden(~X1*X2*X3,data=xdf)
  #Estimate density on a grid
  xseq1 <- seq(min(x[,1]),max(x[,1]),length=ngrid)
  xseq2 <- seq(min(x[,2]),max(x[,2]),length=ngrid)
  xseq3 <- seq(min(x[,3]),max(x[,3]),length=ngrid)
  xgrid <- expand.grid(xseq1,xseq2,xseq3)
  names(xgrid) <- names(xdf)
  dens <- dssden(fit, x=xgrid)
  #Select set of points with probContour probability
  probCutoff <- function(cutoff) sum(prob[prob>cutoff])
  prob <- dens/sum(dens)
  cutseq <- seq(min(prob),max(prob),length=100)
  probcutseq <- sapply(as.list(cutseq),probCutoff)
  cutoff <- cutseq[which.min(abs(probcutseq-probContour))]
  xsel <- xgrid[prob>cutoff,]
  #Find convex hull of selected points  
  contourx <- convhulln(xsel, options = "Tv")
  contourx <- list(highDensPoints=xsel,convhull=contourx)
  if (type=='contours') for (i in 1:nrow(contourx$convhull)) rgl.triangles(contourx$highDensPoints[contourx$convhull[i,],],col=col,alpha=.5)
  return(list(dens=dens,prob=prob,contour=contourx))
}

### EXAMPLE
### ### 3D
### x <- rmvnorm(1000,rep(1,2,3))
### contourx <- contour3d(x, probContour=.5)
###  
### plot3d(x)
### points3d(contourx$highDensPoints,col=2)
### for (i in 1:nrow(contourx$convhull)) rgl.triangles(contourx$highDensPoints[contourx$convhull[i,],],col=2,alpha=.5)

setMethod("plot", signature(x="mds"), function(x,drawlabels=FALSE,labels,plantar=FALSE, # Type of plot
                                      #contour=FALSE,contour.only=FALSE,contour.type='single',contour.fill=FALSE,clus,k,contour.sel,contour.col=rainbow(length(unique(contour.sel))),contour.prob=0.75, # Contour info
                                      point.cex=1,text.cex=1,text.pos=NULL,point.col='blue',text.col='black',point.pch=20,type.3d='p',radius=point.cex/(.05*nrow(x@points)),app='fill',alpha=0.5, # Graphics info
                                      scalecol=FALSE,scale,palette=terrain.colors(100,alpha=alpha),...) { # Scale color info
  # Limits for the plot window
  xlim <- 1.2*range(x@points)
  # Labels if needed
  if (drawlabels & missing(labels)) labels <- rownames(x@points)
  # Overriding colors if scale color is chosen (for expression, mark density, etc)
  if (scalecol)
    {
      if (missing(scale)) stop('You must provide a numeric scale value for each point to use for coloring')
      scale <- as.integer((scale*length(palette))/max(scale,na.rm=TRUE))
      point.col <- palette[scale]
    }
  # One-dimensional MDS (points along one axis)
  if (ncol(x@points)==1) {
    if (plantar) warning('Plantar views only available for 3D maps')
    plot(x=x@points[,1], y=rep(0,nrow(x@points)), xlim=xlim, ylim=xlim, xlab='', ylab='', col=point.col, cex=point.cex, pch=point.pch,...)
    if (drawlabels) text(x@points, rep(0,2), labels,  col=text.col, cex=text.cex, pos=text.pos,...)
  }
  # Two-dimensional MDS
  else if (ncol(x@points)==2) {
    if (plantar) warning('Plantar views only available for 3D maps')
    plot(x@points, xlim=xlim,ylim=xlim,xlab='',ylab='', col=point.col, cex=point.cex, pch=point.pch,...)
    if (drawlabels) text(x@points[,1],x@points[,2],labels, col=text.col, cex=text.cex, pos=text.pos,...)
  }
  # Three-dimensional MDS (3D with RGL and multiple 2D with plantar views of 2-axis combinations)
  else if (ncol(x@points)==3) {
    if (plantar) {
      x.xy <- new('mds',points=x@points[,1:2],R.square=x@R.square)
      x.yz <- new('mds',points=x@points[,2:3],R.square=x@R.square)
      x.xz <- new('mds',points=x@points[,-2],R.square=x@R.square)
      plot(x.xy,drawlabels=drawlabels,labels=labels,plantar=FALSE,point.cex=point.cex,text.cex=text.cex,text.pos=text.pos,point.col=point.col,text.col=text.col,point.pch=point.pch,type.3d=type.3d,app=app,alpha=alpha,...)
      plot(x.yz,drawlabels=drawlabels,labels=labels,plantar=FALSE,point.cex=point.cex,text.cex=text.cex,text.pos=text.pos,point.col=point.col,text.col=text.col,point.pch=point.pch,type.3d=type.3d,app=app,alpha=alpha,...)
      plot(x.xz,drawlabels=drawlabels,labels=labels,plantar=FALSE,point.cex=point.cex,text.cex=text.cex,text.pos=text.pos,point.col=point.col,text.col=text.col,point.pch=point.pch,type.3d=type.3d,app=app,alpha=alpha,...)
    }
    else {
      if ('rgl' %in% loadedNamespaces()) {
        plot3d(x@points, type='p', radius=radius, xlim=xlim,ylim=xlim,zlim=xlim,xlab='',ylab='',zlab='', col=point.col, cex=point.cex, alpha=0.5, ...) # point.cex/(2*nrow(x@points))
        if (drawlabels) text3d(x@points[,1],x@points[,2],x@points[,3],texts=labels, col=text.col, cex=text.cex, pos=text.pos,...) # Draw labels before spheres to make them visible
        if (type.3d=='s') spheres3d(x@points, radius=radius, xlim=xlim,ylim=xlim,zlim=xlim,xlab='',ylab='',zlab='', col=point.col, cex=point.cex, alpha=0.5, ...) # point.cex/(2*nrow(x@points))
        #if (type.3d=='spheres') spheres3d(x@points, xlim=xlim,ylim=xlim,zlim=xlim,xlab='',ylab='',zlab='', col=point.col, radius=point.cex/(2*nrow(x@points)),alpha=alpha,front='fill',back='fill',...)
      } else stop('rgl library has not been loaded!')
    }
  }
}
)

