# Adopting the S3 hclust class into an S4 one
setOldClass("hclust")

# Classes, methods and functions for clustering and clustering reproducibility assessment of distGPS objects

# Classes
setClass("clusGPS", representation(h="hclust",clus="list",adjusted="logical"))

## SHOW
setMethod("show","clusGPS",function(object) {
cat("Object of class clusGPS with clustering for",length(object@h$height)+1,"elements\n")
}
)

pden.adjust <- function(clus,mc.cores=1)
  {
    if (!is(clus,"clusGPS")) stop('Object must be of "clusGPS" class')
    if (clus@adjusted) stop('Posterior density for clustering is already adjusted')
    if (mc.cores>1) {
      if ('multicore' %in% loadedNamespaces()) {
        clus@clus <- multicore::mclapply(clus@clus,function(x) { # For each k
          id <- x$id
          priorprob <- as.numeric(prop.table(table(x$id)))
          pointprob <- x$pointprob[,-1:-2] # Get just prob values, leave coordinates aside
          pointprob <- t(t(pointprob) * priorprob) # A = T( T(A) * PriorProb)
          pointprob <- pointprob/rowSums(pointprob) # A = A / rowSums(A)
          pointprob <- cbind(x$pointprob[,1:2],pointprob) # Put coordinates back
          # PP of accurate class. for each element/cluster in partition for k=x
          propclass <- multicore::mclapply(1:(ncol(pointprob)-2),function(y) { pointprob[id==y,paste('C',y,sep='')] / rowSums(pointprob[id==y,-1:-2]) },mc.cores=mc.cores,mc.preschedule=FALSE)
          return(list(id=x$id,pden=x$pden,propclass=propclass,pointprob=pointprob)) # return clus list
        },mc.cores=mc.cores,mc.preschedule=FALSE)
      }
      else stop('multicore library has not been loaded!')
    } else {
      clus@clus <- lapply(clus@clus,function(x) { # For each k
        id <- x$id
        priorprob <- as.numeric(prop.table(table(x$id)))
        pointprob <- x$pointprob[,-1:-2] # Get just prob values, leave coordinates aside
        pointprob <- t(t(pointprob) * priorprob) # A = T( T(A) * PriorProb)
        pointprob <- pointprob/rowSums(pointprob) # A = A / rowSums(A)
        pointprob <- cbind(x$pointprob[,1:2],pointprob) # Put coordinates back
        propclass <- lapply(1:(ncol(pointprob)-2),function(y) { pointprob[id==y,paste('C',y,sep='')] / rowSums(pointprob[id==y,-1:-2]) }) # PP of accurate class. for each element/cluster in partition for k=x
        return(list(id=x$id,pden=x$pden,propclass=propclass,pointprob=pointprob)) # return clus list
      })
    }
    clus@adjusted <- TRUE
    return(clus)
  }

setGeneric("clusGPS", function(d,m,h,grid,ngrid=1000,densgrid=TRUE,type='agnes',method='average',samplesize=1,p.adjust=TRUE,k=2:5,mc.cores=1,set.seed=149,verbose=TRUE,...) standardGeneric("clusGPS"))

setMethod("clusGPS", signature=c(d='distGPS',m='mds'), function(d,m,h,grid,ngrid=1000,densgrid=TRUE,type='agnes',method='average',samplesize=1,p.adjust=TRUE,k=2:5,mc.cores=1,set.seed=149,verbose=TRUE,...) {
  if (ncol(m@points)!=2) stop('Currently only supported for 2 dimensional maps')
  m@points <- m@points[rownames(d@d),] # Reindex internally to revert possible splitMDS
  #d@d <- d@d[rownames(m@points),rownames(m@points)] # Reindex internally with rownames from m@points for correct calculation of r2
  # Resampling if wanted
  if (samplesize<1)
    {
      set.seed(set.seed)
      sel <- sample(1:as.integer(nrow(d@d)*samplesize),replace=FALSE)
      d@d <- d@d[sel,sel]
      m@points <- m@points[sel,]
    }
  # Clustering
  if (type=='agnes') {
    if (missing(h)) {
      if (verbose) cat('\nPerforming Agglomerative Nesting\n')
      h <- as.hclust(cluster::agnes(d@d,method=method,...))
    }
    # Since we must evaluate the density of all clusters over a common grid, we first obtain the grid for all the map points unless a specific grid is given
    if (missing(grid))
      {
        if (verbose) cat('\nPrecalculating Grid\n')
        if (densgrid) # Adequate grid step points to data density to account for higher resolution where needed
          {
            gridx <- quantile(m@points[,1],probs=seq(0,1,length.out=as.integer(sqrt(ngrid))))
            gridy <- quantile(m@points[,2],probs=seq(0,1,length.out=as.integer(sqrt(ngrid))))
            grid <- cbind(gridx,gridy)
          }
        else # Do not adequate grid step points to data density automatically
          {
            left <- rep(0,2) # Grid precalculation taken from DPdensity source code
            right <- rep(0,2)
            left[1] <- min(m@points[,1])-0.5*sqrt(var(m@points[,1]))
            right[1] <- max(m@points[,1])+0.5*sqrt(var(m@points[,1]))      
            left[2] <- min(m@points[,2])-0.5*sqrt(var(m@points[,2]))
            right[2] <- max(m@points[,2])+0.5*sqrt(var(m@points[,2]))
            ngrid <- as.integer(sqrt(ngrid))
            grid1 <- seq(left[1],right[1],length=ngrid)
            grid2 <- seq(left[2],right[2],length=ngrid)
            grid <- cbind(grid1,grid2)
          }
      }
    names(k) <- as.character(k)
    clus <- lapply(k,function(x,m,p.adjust) {
      id <- cutree(h,k=x) # ID of cluster to which each element is assigned
      if (mc.cores>1) {
        if ('multicore' %in% loadedNamespaces()) {
          pden <- multicore::mclapply(1:x,function(y) {
            mpoints <- m@points[as.numeric(id)==y,]
            if (verbose) cat("\nCalculating posterior density of mis-classification for cluster",y,"of",x,"\n")
            contour2dDP(mpoints,grid=grid,xlim=range(mpoints),ylim=range(mpoints),contour.type='none',...)
          },mc.cores=mc.cores,mc.preschedule=FALSE)
        }
        else stop('multicore library has not been loaded!')
      } else {
        pden <- lapply(1:x,function(y) {
          mpoints <- m@points[as.numeric(id)==y,]
          if (verbose) cat("\nCalculating posterior density of mis-classification for cluster",y,"of",x,"\n")
          contour2dDP(mpoints,grid=grid,xlim=range(mpoints),ylim=range(mpoints),contour.type='none',...)
        })
      }
      names(pden) <- as.character(1:x) #DPdensity object for each cluster 1:k in the partition for k=x
      probs <- lapply(pden,function(x) x$dens) # Posterior density for each element in clusters 1:k in the partition for k=x
      #normpoints <- ceiling(normCoords(m@points,c(1,floor(sqrt(ngrid))))+.1) # Assigning each point in the MDS to its corresponding grid cell (coordinate transformation)
      #normpoints <- ceiling(normCoords(m@points,c(1,nrow(pden[[1]]$dens)))+.1) # Assigning each point in the MDS to its corresponding grid cell (coordinate transformation, ONLY VALID FOR REGULAR GRIDS)
      # New assignation of mds points to gridpoints, now valid for regular and irregular grids
      normx <- unlist(lapply(m@points[,1],function(x) sum(x>grid[,1])+1))
      normy <- unlist(lapply(m@points[,2],function(x) sum(x>grid[,2])+1))
      normpoints <- cbind(normx,normy)
      pointprob <- t(apply(normpoints,1,function(x) unlist(lapply(probs,function(y) y[x[1],x[2]]))))
      pointprob <- cbind(normpoints,pointprob)
      colnames(pointprob) <- c('X','Y',paste('C',1:x,sep=''))
      # if (p.adjust) pointprob[,-1:-2] <- esteps(pointprob[,-1:-2])$pp # Posterior probability table for # Old adjustment for prior probabilities esteps
      propclass <- lapply(1:x,function(y) { pointprob[id==y,paste('C',y,sep='')] / rowSums(pointprob[id==y,-1:-2]) } ) # Posterior probability of accurate classification for each element/cluster in partition for k=x
      return(list(id=id,pden=pden,propclass=propclass,pointprob=pointprob)) # return clus list
    },m,p.adjust)
    ans <- new("clusGPS",h=h,clus=clus,adjusted=FALSE)
    # PP adjustment
    if (p.adjust) ans <- pden.adjust(ans,mc.cores=mc.cores)
    return(ans)
  }
  else stop('Currently supporting agnes clustering only')
})

setMethod("clusGPS", signature=c(d='distGPS',m='missing'), function(d,type='agnes',method='average',samplesize=1,k=2:5,...) {
  # Resampling if wanted
  if (samplesize<1)
    {
      sel <- sample(1:as.integer(nrow(d@d)*sample))
      d@d <- d@d[sel,sel]
    }
  if (type=='agnes') {
    h <- as.hclust(cluster::agnes(d@d,method=method,...))
    clus <- lapply(k,function(x,m,p.adjust) {
      clus <- cutree(h,k=x)
    })
    new("clusGPS",h=h,clus=list(clus=clus,pden=NA,propclass=NA,pointprob=NA),adjusted=FALSE)
  }
})

setMethod("plot", signature=c(x="clusGPS"), function(x,type='stats',cut=.7,cut.col='red',cut.lwd=4,cut.lty=2,k,probContour=.7,contour.col=NULL,labels=as.character(1:k),labcex=1,...) {
  # If type=='stats' prints stats of correct classification for each individual cluster
  # If type=='stats' prints stats of avg correct classification for clustering
  # If type=='contours' draws contour information over existing map
  # If type=='density' draws full contour lines and information about the DPdensity object
  if (type=='stats') {
    i <- as.character(k)
    #for (i in 1:length(x@clus)) {
    #plot(unlist(lapply(x@clus[[i]]$propclass,mean)),type='o',main=sprintf('Posterior density for %s clusters',names(x@clus[i])),xlab='Cluster ID',ylab='avg agreement score',...)
    plot(unlist(lapply(x@clus[[i]]$propclass,mean)),type='o',...)
    abline(cut,0,col=cut.col,lty=cut.lty,lwd=cut.lwd)
    #}
    }
  else if (type=='avgstat')
    {
      ans <- lapply(x@clus,function(x) lapply(x$propclass,mean))
      plot(c(1,as.numeric(unlist(lapply(ans,function(x) mean(unlist(x)))))),type='o',...)
      #plot(c(1,as.numeric(unlist(lapply(ans,function(x) mean(unlist(x)))))),type='o',main='Avg Posterior density / k',xlab='Number of clusters (k)',ylab='avg agreement score',...)
      abline(cut,0,col=cut.col,lty=cut.lty,lwd=cut.lwd)
    }
  else if (type=='contours')
    {
      if (!(as.character(k) %in% names(x@clus))) stop('No clustering information for this number of clusters')
      pden <- x@clus[[as.character(k)]]$pden
      if (is.null(contour.col)) contour.col <- rainbow(length(pden))
      for (i in 1:length(pden)) {
        den <- pden[[i]]
        prob <- den$dens/sum(den$dens)
        probo <- prob[order(prob)]
        cutoff <- probo[match(TRUE,cumsum(probo) > 1-probContour)]
        contour(x=den$grid1,y=den$grid2,z=matrix(den$dens,nrow=length(den$grid1),ncol=length(den$grid2)),levels=cutoff*sum(den$dens),col=contour.col[i],add=TRUE,axes=FALSE,labels=labels,labcex=labcex,...)
      }
    }
  else if (type=='density')
    {
      if (!(as.character(k) %in% names(x@clus))) stop('No clustering information for this number of clusters')
      pden <- x@clus[[as.character(k)]]$pden
      if (is.null(contour.col)) contour.col <- rainbow(length(pden))
      plot(0,col=NA,...)
      for (i in 1:length(pden)) { den <- pden[[i]]; contour(x=den$grid1,y=den$grid2,z=matrix(den$dens,nrow=length(den$grid1),ncol=length(den$grid2)),col=contour.col[i],add=TRUE,axes=FALSE,...) }
      for (i in 1:length(pden)) { den <- pden[[i]]; plot(den,...) }      
    }
  else stop("Invalid type of contour plot. Valid types are 'avgstat', 'stats', 'contours', 'density' ")
})

normCoords <- function(coords,newrange)
# First make origin in 0,0
# Then divide by max to turn them to 0,1
# Finally expand or contract to the new range
{
  coords <- apply(coords,2,function(x) x-min(x))
  coords <- apply(coords,2,function(x) x/max(x))
  newrange <- newrange[2] - newrange[1]
  coords <- coords * newrange
}

enrichmentClusGPS <- function(x,clus,k,index=NULL,p.adjust=TRUE,method='BH',plt=FALSE,signif=0.05,mc.cores=1,...)
  # Computes enrichment or depletion of marks in a given cluster
  # X is table with epigenes / factors
  {
    if (!is.null(index)) x <- x[index,] # Reindex if necessary (if we have used reshuffling for split MDS)
    id <- clus@clus[[as.character(k)]]$id # Cluster / gene assignment
    if (mc.cores>1) {
      if ('multicore' %in% loadedNamespaces()) {
        # Real proportion of factor abundance in each cluster  
        ans1 <- do.call(rbind,multicore::mclapply(colnames(x),function(myfactor) {
          as.numeric(unlist(lapply(1:k,function(mycut) prop.table(table(x[id==mycut,myfactor]))['1']))) # Proportion of genes with mark myfactor in cluster mycut
        },mc.cores=mc.cores,mc.preschedule=FALSE))
        rownames(ans1) <- colnames(x)
        # Pvalue for observed vs expected proportion
        ans2 <- do.call(rbind,multicore::mclapply(colnames(x),function(myfactor) {
          as.numeric(unlist(lapply(1:k,function(mycut) fisher.test(rbind(table(x[id==mycut,myfactor]),table(x[id!=mycut,myfactor])))$p.value))) # Proportion of genes with mark myfactor not in in cluster mycut
        },mc.cores=mc.cores,mc.preschedule=FALSE))  
        rownames(ans2) <- colnames(x)
      }
      else stop('multicore library has not been loaded!')
    } else {
      # Real proportion of factor abundance in each cluster  
      ans1 <- do.call(rbind,lapply(colnames(x),function(myfactor) {
        as.numeric(unlist(lapply(1:k,function(mycut) prop.table(table(x[id==mycut,myfactor]))['1']))) # Proportion of genes with mark myfactor in cluster mycut
      }))
      rownames(ans1) <- colnames(x)
      # Pvalue for observed vs expected proportion
      ans2 <- do.call(rbind,lapply(colnames(x),function(myfactor) {
        as.numeric(unlist(lapply(1:k,function(mycut) fisher.test(rbind(table(x[id==mycut,myfactor]),table(x[id!=mycut,myfactor])))$p.value))) # Proportion of genes with mark myfactor not in in cluster mycut
      }))  
      rownames(ans2) <- colnames(x)
    }
    if (p.adjust) ans2 <- apply(ans2,2,p.adjust,method=method)    
    fc <- (ans1 / colSums(x/nrow(x))) - 1
    par(las=3)
    if (plt) {
      for (i in 1:k)
        {
          alpha <- rep(0,ncol(x))
          cols <- rep('grey',ncol(x))
          cols[ans2[,i]<signif] <- 'yellow'
          cols[ans2[,i]<1e-05] <- 'orange'
          cols[ans2[,i]<1e-10] <- 'red'
          cols[ans2[,i]<1e-20] <- 'purple'
          par(mar=c(5,10,5,5),las=2)
          barplot2(fc[,i],col=cols,xlim=c(min(fc[,i],na.rm=TRUE),max(fc[,i],na.rm=TRUE)),...)
          #legend('topleft',pch=18,col=c('purple','red','orange','yellow'),title=sprintf('Cluster/Factor distribution pvalue for %s Fisher Test',ifelse(p.adjust,paste(method,"adjusted"),'unadjusted')),
          #       legend=c('<1e-20','<1e-10','<1e-05',as.character(signif)))
        }
    }
    list(props=ans1,p.value=ans2,fc=fc)
  }
