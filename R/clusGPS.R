# Adopting the S3 hclust class into an S4 one
setOldClass("hclust")

# Classes, methods and functions for clustering and clustering reproducibility assessment of distGPS objects

# Classes
setClass("clusGPS", representation(h="hclust",clus="list",adjusted="logical"))

# Some accessor and helper functions
setMethod("show","clusGPS",function(object) {
cat("Object of class clusGPS with clustering for",length(object@h$height)+1,"elements.\n")
cat(length(object@clus),"clustering configuration(s) with name(s)",names(object@clus),"\n")
}
)

setGeneric("clusNames",function(clus) standardGeneric("clusNames")); setMethod("clusNames","clusGPS",function(clus) { names(clus@clus) })

setGeneric("tabClusters",function(clus,name) standardGeneric("tabClusters"))
setMethod("tabClusters","clusGPS",function(clus,name) {if (!as.character(name) %in% names(clus@clus)) stop('No cluster entry with this name')
  else table(clus@clus[[as.character(name)]]$id)
})

setGeneric("clusterID",function(clus,name) standardGeneric("clusterID"))
setMethod("clusterID","clusGPS",function(clus,name) {
  if (!as.character(name) %in% names(clus@clus)) stop('No cluster entry with this name')
  else as.numeric(clus@clus[[as.character(name)]]$id)
})

# Auxiliary functions
pden.adjust <- function(clus,mc.cores=1)
  {
    if (!is(clus,"clusGPS")) stop('Object must be of "clusGPS" class')
    if (clus@adjusted) stop('Posterior density for clustering is already adjusted')
    if (mc.cores>1) {
      if ('parallel' %in% loadedNamespaces()) {
        clus@clus <- parallel::mclapply(clus@clus,function(x) { # For each k
          id <- x$id
          priorprob <- as.numeric(prop.table(table(x$id)))
          pointprob <- x$pointprob[,-1:-2] # Get just prob values, leave coordinates aside
          pointprob <- t(t(pointprob) * priorprob) # A = T( T(A) * PriorProb)
          pointprob <- pointprob/rowSums(pointprob) # A = A / rowSums(A)
          pointprob <- cbind(x$pointprob[,1:2],pointprob) # Put coordinates back
          # PP of accurate class. for each element/cluster in partition for k=x
          propclass <- parallel::mclapply(1:(ncol(pointprob)-2),function(y) { pointprob[id==y,paste('C',y,sep='')] / rowSums(pointprob[id==y,-1:-2]) },mc.cores=mc.cores,mc.preschedule=FALSE)
          return(list(id=x$id,pden=x$pden,propclass=propclass,pointprob=pointprob)) # return clus list
        },mc.cores=mc.cores,mc.preschedule=FALSE)
      }
      else stop('parallel library has not been loaded!')
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

preCalcGrid <- function(m,densgrid,ngrid)
  {
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
    return(grid)
  }

# Skeleton for premergeClusters methods for the future. For Now we will use it pedestrian way
### setGeneric("premergeClusters", function(clus,m,n,method='manual',plot=FALSE,mc.cores=mc.cores,...) standardGeneric("premergeClustersGPS"))
###  
### setMethod("premergeClusters", signature=c(clus='clusGPS',m='mds'), function(clus,m,n,method='manual',plot=FALSE,mc.cores=mc.cores,...) {
###   if ('parallel' %in% loadedNamespaces() & mc.cores>1)
###     clus@clus <- parallel::mclapply(clus@clus,premergeClusters,m=m,n=n,method=method,plot=plot,mc.cores=mc.cores)
###   else clus@clus <- lapply(clus@clus,premergeClusters,m=m,n=n,method=method,plot=plot)
###  
### })
###  
###  
### setMethod("premergeClusters", signature=c(clus='list',m='mds'), function(clus,m,n,method='manual',plot=FALSE,...) 

###premergeClusters <- function(m,id,n,method='manual',recalcDist=TRUE,retCentroids=FALSE,plt=FALSE,verbose=FALSE)
#### Clean the noisy data (ie tiny clusters generated during the bottom-up hierarchical clustering process, and merge them into bigger ones until we can calculate their density)
###  {
###    cat(sprintf('\nPre-merging non-clustered points in nodules of size %d...\n',n))
###    if (method=='manual') {
###      sqdist <- function(x,y) apply(y,1,function(r) sqrt((x[1]-r[1])^2 + (x[2]-r[2])^2))
###        #centroids <- do.call(rbind,by(m@points,id,mean)) # Compute mean of coordinates for every cluster
###        centroids <- do.call(rbind,by(m@points,id,colMeans)) # Compute mean of coordinates for every cluster
###        dist.centroids <- as.matrix(dist(centroids,method='euclidean')) # First call
###        diag(dist.centroids) <- Inf
###        conf <- 1
###        while ((min(table(id))<=n)) # While there are still small clusters
###          {
###            if (plt==TRUE) { plot(m,point.col='white'); points(centroids,col=1:nrow(centroids),cex=2,lwd=2,pch=19) }
###            id.small <- names(table(id)[table(id)<n]) # ID of small clusters        
###            # Look for minimum distance (out of the diagonal) and reassign clusters, smaller to bigger
###            if (length(id.small)>1)
###              {
###                inds <- which(dist.centroids[id.small,]==min(dist.centroids[id.small,]),arr.ind=TRUE)[1,] # In case of tie, get first one
###                oldclus <- id.small[inds[1]]
###                newclus <- names(table(id))[inds[2]]
###              }
###            else {
###              oldclus <- id.small
###              newclus <- names(table(id))[which.min(dist.centroids[id.small,])]
###            }
###            id[id==as.numeric(oldclus)] <- as.numeric(newclus)
###            # recompute centroid for changed cluster
###            centroids[as.character(newclus),] <- colMeans(m@points[id == as.numeric(newclus),]) # Just necessary to compute mean of coordinates for the changed cluster
###            # recompute distances for changed cluster...
###            #centroids <- do.call(rbind,by(m@points,id,mean)) # Compute mean of coordinates for every cluster
###            centroids <- centroids[rownames(centroids) %in% as.character(id),] # Remove centroids from old cluster
###            if (!recalcDist) { dist.centroids <- dist.centroids[rownames(centroids),rownames(centroids)] # Remove distance row and column for old cluster
###              dist.centroids[newclus,] <- dist.centroids[,newclus] <- sqdist(centroids[newclus,],centroids)
###              diag(dist.centroids) <- Inf }
###            else { dist.centroids <- as.matrix(dist(centroids,method='euclidean')) # First call
###                   diag(dist.centroids) <- Inf }         
###            if (verbose) cat(sprintf('\nConfiguration %d: Assigned cluster %s to %s',conf,oldclus,newclus))
###            conf <- conf+1
###          }
###        #ifelse(retCentroids,return(list(centroids=centroids,id=id)),return(id))
###      if (retCentroids==TRUE) return(list(centroids=centroids,id=id)) else return(id)
###      }
###    else if (method=='auto') {
###      cat('\nNot yet implemented...')
###    }
###  }

premergeClusters <- function(m,id,n,method='manual',recalcDist=TRUE,retCentroids=FALSE,plt=FALSE,verbose=TRUE)
# Clean the noisy data (ie tiny clusters generated during the bottom-up hierarchical clustering process, and merge them into bigger ones until we can calculate their density)
  {
    cat(sprintf('\nPre-merging non-clustered points in nodules of size %d...\n',n))
    if (method=='manual') {
      sqdist <- function(x,y) apply(y,1,function(r) sqrt((x[1]-r[1])^2 + (x[2]-r[2])^2))
        #centroids <- do.call(rbind,by(m@points,id,mean)) # Compute mean of coordinates for every cluster
        centroids <- do.call(rbind,by(m@points,id,colMeans)) # Compute mean of coordinates for every cluster
        dist.centroids <- as.matrix(dist(centroids,method='euclidean')) # First call
        diag(dist.centroids) <- Inf
        conf <- 1
        while (min(table(id))<n) # While there are still small clusters
          {
            if (plt==TRUE) { plot(m,point.col='white'); points(centroids,col=1:nrow(centroids),cex=2,lwd=2,pch=19) }
            id.small <- names(table(id)[table(id)<n]) # ID of small clusters        
            # Look for minimum distance (out of the diagonal) and reassign clusters, smaller to bigger
            if (length(id.small)>1)
              {
                inds <- which(dist.centroids[id.small,]==min(dist.centroids[id.small,]),arr.ind=TRUE)[1,] # In case of tie, get first one
                oldclus <- id.small[inds[1]]
                newclus <- names(table(id))[inds[2]]
                #cat(sprintf('\nConfiguration %d: Assigning cluster %s to %s',conf,oldclus,newclus))
              }
            else {
              oldclus <- id.small
              newclus <- names(table(id))[which.min(dist.centroids[id.small,])]
              #newclus <- names(table(id))[which(dist.centroids[id.small,])==min(dist.centroids[id.small,])]
              #cat(sprintf('\nConfiguration %d: Assigning cluster %s to %s',conf,oldclus,newclus))
            }
            #cat('\nReassigning cluster')
            id[id==as.numeric(oldclus)] <- as.numeric(newclus)
            # recompute centroid for changed cluster
            #cat('\nRecompute centroids')
            centroids[as.character(newclus),] <- colMeans(m@points[id == as.numeric(newclus),]) # Just necessary to compute mean of coordinates for the changed cluster
            # recompute distances for changed cluster...
            #centroids <- do.call(rbind,by(m@points,id,mean)) # Compute mean of coordinates for every cluster
            centroids <- centroids[rownames(centroids) %in% as.character(id),] # Remove centroids from old cluster
            if (!recalcDist) { dist.centroids <- dist.centroids[rownames(centroids),rownames(centroids)] # Remove distance row and column for old cluster
                               #cat('\nRecalcDist1')
              dist.centroids[newclus,] <- dist.centroids[,newclus] <- sqdist(centroids[newclus,],centroids)
              diag(dist.centroids) <- Inf }
            else { dist.centroids <- as.matrix(dist(centroids,method='euclidean')) # First call
                   #cat('\nRecalcDist2')
                   diag(dist.centroids) <- Inf }
            if (verbose) cat(sprintf('\nConfiguration %d: Assigned cluster %s to %s',conf,oldclus,newclus))
            conf <- conf+1
            #print(table(id)[table(id)>=70])
          }
        #ifelse(retCentroids,return(list(centroids=centroids,id=id)),return(id))
      if (retCentroids==TRUE) return(list(centroids=centroids,id=id)) else return(id)
      }
    else if (method=='auto') {
      cat('\nNot yet implemented...')
    }
  }

setGeneric("clusGPS",
           function(d,m,h,sel=NULL,id=NULL,grid,ngrid=1000,densgrid=FALSE,preMerge=TRUE,type='hclust',method='average',samplesize=1,p.adjust=TRUE,k,mc.cores=1,set.seed=149,verbose=TRUE,minpoints=70,...) standardGeneric("clusGPS"))

setMethod("clusGPS", signature=c(d='distGPS',m='mds'),
          function(d,m,h,sel=NULL,id=NULL,grid,ngrid=1000,densgrid=FALSE,preMerge=TRUE,type='hclust',method='average',samplesize=1,p.adjust=TRUE,k,mc.cores=1,set.seed=149,verbose=TRUE,minpoints=70,...) {
  # MDS density control, currently only possible for 2D maps
  if (ncol(m@points)!=2) stop('Currently only supported for 2 dimensional maps')
  # Iteration control, either multiple cluster cuts (k) or multiple cluster size (minpoints), you can't have it both ways
  if (length(k)>1 & length(minpoints)>1) stop('Multiple clustering is possible either at cluster number or cluster size, not both.')
  m@points <- m@points[rownames(d@d),] # Reindex internally to revert possible splitMDS
  # Resampling if wanted, consider deprecating...
  if (samplesize<1)
    {
      set.seed(set.seed)
      sel <- sample(1:as.integer(nrow(d@d)*samplesize),replace=FALSE)
      d@d <- d@d[sel,sel]
      m@points <- m@points[sel,]
    }
  # Clustering, consider deprecating internal clustering and performing it outside of this function...
  if (type=='hclust') {
    if (missing(h)) {
      if (verbose) cat(sprintf('\nPerforming Hierarchical Clustering, method=%s\n',method))
      h <- as.hclust(hclust(as.dist(d@d),method=method,...))
    }
    # If we want to work within a certain cluster only to subdivide it, subset accordingly
    if (!is.null(sel)) {
      cat(sprintf('\nSubselecting points from cluster'))    
      #sel <- cutree(h,k=sel$k)==sel$clus # Subselect elements from the specified cluster at the specified cut
      m@points <- m@points[sel,]; d@d <- d@d[sel,sel]
    }    
    # Since we must evaluate the density of all clusters over a common grid, we first obtain the grid for all the map points unless a specific grid is given
    if (missing(grid)) {
        if (verbose) cat('\nPrecalculating Grid\n')
        grid <- preCalcGrid(m,densgrid,ngrid)
      }
    names(k) <- as.character(k)
    # Establishing iteration item, either k or minpoints
    if (length(k)>1 & length(minpoints)==1) iter <- k else if (length(k)==1 & length(minpoints)>1) iter <- minpoints else iter <- k
    names(iter) <- as.character(iter)
    clus <- lapply(iter,function(x,m,p.adjust) {
      id <- cutree(h,k=x) # ID of cluster to which each element is assigned
      if (!is.null(sel)) id <- id[sel] # If we want to work just with a subcluster, subselect only those points
      if (preMerge) id <- premergeClusters(m,id,minpoints,recalcDist=TRUE,retCentroids=FALSE,plt=FALSE,verbose=FALSE)    
      if (length(k)>1 & length(minpoints)==1) idx <- as.numeric(names(table(id))[table(id)>=minpoints]) # We iterate over k
      else if (length(k)==1 & length(minpoints)>1) idx <- as.numeric(names(table(id))[table(id)>=x]) # We iterate over minpoints
      else idx <- as.numeric(names(table(id)[table(id)>=minpoints])) # By default we iterate over k
      if (mc.cores>1) {
        if ('parallel' %in% loadedNamespaces()) {
          pden <- parallel::mclapply(idx,function(y) {
            mpoints <- m@points[as.numeric(id)==y,]
            if (verbose)
              if (length(k)>1) cat("Calculating posterior density of mis-classification for cluster:",y)
              else if (length(minpoints)>1) cat("Calculating posterior density of mis-classification for cluster:",y)
              else cat("Calculating posterior density of mis-classification for cluster:",y)
            contour2dDP(mpoints,grid=grid,xlim=range(mpoints),ylim=range(mpoints),contour.type='none',...)
          },mc.cores=mc.cores,mc.preschedule=FALSE)
        }
        else stop('parallel library has not been loaded!')
      } else {
        pden <- lapply(idx,function(y) {
          mpoints <- m@points[as.numeric(id)==y,]
          if (verbose)
            if (length(k)>1) cat("Calculating posterior density of mis-classification for cluster:",y)
            else if (length(minpoints)>1) cat("Calculating posterior density of mis-classification for cluster:",y)
            else cat("\nCalculating posterior density of mis-classification for cluster:",y)
          contour2dDP(mpoints,grid=grid,xlim=range(mpoints),ylim=range(mpoints),contour.type='none',...)
        })
      }
      names(pden) <- as.character(idx) #DPdensity object for each cluster in the partition for k=x where ngenes >= minpoints
      probs <- lapply(pden,function(x) x$dens) # Posterior density for each element in clusters 1:k in the partition for k=x
      # New assignation of mds points to gridpoints, now valid for regular and irregular grids
      normx <- unlist(lapply(m@points[,1],function(x) sum(x>grid[,1])+1))
      normy <- unlist(lapply(m@points[,2],function(x) sum(x>grid[,2])+1))
      normpoints <- cbind(normx,normy)
      pointprob <- t(apply(normpoints,1,function(x) unlist(lapply(probs,function(y) y[x[1],x[2]]))))      
      # On the fly p.adjust, takes less time and memory
      if (p.adjust)
        {
          cat('\nAdjusting posterior probabilities...\n')
          priorprob <- as.numeric(prop.table(table(idx)))
          pointprob <- t(t(pointprob) * priorprob) # A = T( T(A) * PriorProb)
          pointprob <- pointprob/rowSums(pointprob) # A = A / rowSums(A)
        }
      pointprob <- cbind(normpoints,pointprob)
      colnames(pointprob) <- c('X','Y',paste('C',idx,sep=''))
      propclass <- lapply(idx,function(y) { pointprob[id==y,paste('C',y,sep='')] / rowSums(pointprob[id==y,-1:-2]) } ) # Posterior probability of accurate classification for each element/cluster in partition for k=x
      names(propclass) <- names(pden)
      return(list(id=id,pden=pden,propclass=propclass,pointprob=pointprob)) # return clus list
    },m,p.adjust)
    ans <- new("clusGPS",h=h,clus=clus,adjusted=p.adjust)
    # PP adjustment
    return(ans)
  }
  else stop('Currently supporting hclust only, supply an external clustering object (h) if you want to use other clustering functions')
})

### # DEPRECATED
### setMethod("clusGPS", signature=c(d='distGPS',m='missing'), function(d,h,id=NULL,type='hclust',method='average',samplesize=1,k,...) {
###   # Resampling if wanted
###   if (samplesize<1)
###     {
###       if (!missing(h)) stop('Cannot resample if an initial clustering is already provided')
###       sel <- sample(1:as.integer(nrow(d@d)*sample))
###       d@d <- d@d[sel,sel]
###     }
###   if (type=='hclust') {
###     if (missing(h)) h <- as.hclust(hclust(as.dist(d@d),method=method,...))
###     clus <- lapply(k,function(x,m,p.adjust) {
###       id <- cutree(h,k=x)
###       list(id=id)
###     })
###     names(clus) <- as.character(k)
###     #new("clusGPS",h=h,clus=list(clus=clus,pden=NA,propclass=NA,pointprob=NA),adjusted=FALSE)
###     new("clusGPS",h=h,clus=clus,adjusted=FALSE)
###   }
###   else stop('Currently supporting hclust only, supply an external clustering object (h) if you want to use other clustering functions')
### })

setMethod("plot", signature=c(x="clusGPS"), function(x,m,type='stats',cut=.7,cut.col='red',cut.lwd=4,cut.lty=2,k,value=FALSE,probContour=.7,contour.col=NULL,drawlabels=FALSE,labels='',labcex=1,...) {
  # If type=='stats' prints stats of correct classification for each individual cluster
  # If type=='avgstat' prints stats of avg correct classification for clustering
  # If type=='contours' draws contour information over existing map
  # If type=='density' draws full contour lines and information about the DPdensity object
  if (type=='stats') {
    i <- as.character(k)
    #for (i in 1:length(x@clus)) {
    #plot(unlist(lapply(x@clus[[i]]$propclass,mean)),type='o',main=sprintf('Posterior density for %s clusters',names(x@clus[i])),xlab='Cluster ID',ylab='avg agreement score',...)
    plot(unlist(lapply(x@clus[[i]]$propclass,mean)),type='o',...)
    wmean <- weighted.mean(unlist(lapply(x@clus[[i]]$propclass,mean)),unlist(lapply(x@clus[[i]]$propclass,length)))
    #print(wmean)
    abline(cut,0,col=cut.col,lty=cut.lty,lwd=cut.lwd)
    abline(wmean,0,col='black',lty=2,lwd=2)
    #}
    }
  else if (type=='avgstat')
    {
      #ans <- lapply(x@clus,function(x) lapply(x$propclass,mean))
      #plot(c(1,as.numeric(unlist(lapply(ans,function(x) mean(unlist(x)))))),type='o',...)
      #plot(c(1,as.numeric(unlist(lapply(ans,function(x) mean(unlist(x)))))),type='o',main='Avg Posterior density / k',xlab='Number of clusters (k)',ylab='avg agreement score',...)
      wmean <- lapply(x@clus, function(y) weighted.mean(unlist(lapply(y$propclass,mean)),unlist(lapply(y$propclass,length))))
      #plot(c(1,as.numeric(unlist(wmean))),type='o',...)
      plot(as.numeric(unlist(wmean)),type='o',...)
      abline(cut,0,col=cut.col,lty=cut.lty,lwd=cut.lwd)
    }
  else if (type=='contours')
    {
      if (!(as.character(k) %in% names(x@clus))) stop('No clustering information for this number of clusters')
      pden <- x@clus[[as.character(k)]]$pden
      if (length(labels)==1) if (labels=='') labels <- names(pden)      
      if (is.null(contour.col)) contour.col <- rainbow(length(pden))
      for (i in 1:length(pden)) {
        den <- pden[[i]]
        prob <- den$dens/sum(den$dens)
        probo <- prob[order(prob)]
        cutoff <- probo[match(TRUE,cumsum(probo) > 1-probContour)]
        contour(x=den$grid1,y=den$grid2,z=matrix(den$dens,nrow=length(den$grid1),ncol=length(den$grid2)),levels=cutoff*sum(den$dens),col=contour.col[i],add=TRUE,axes=FALSE,drawlabels=drawlabels,labels=labels[i],labcex=labcex,...)
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

### profileClusters <- function(x,uniqueCount=TRUE,clus,i,minpoints,merged=FALSE,log2=TRUE,plt=FALSE)
###   # Computes enrichment or depletion of marks in a given cluster
###   # X is table with epigenes / factors, h is a valid clustering for that matrix, k is a cluster cut configuration
###   # n is used to consider only clusters with more than or n points, 0 for all points
###   # Univ is used to compute average factor distribution, if TRUE all genes are considered, if FALSE only genes in clusters where ngenes>=n are considered
###   # Result is given in the shape of a heatmap
###   {
###     if (missing(minpoints) & !merged) stop('For non-merged clusterings the minimum cluster size considered in the clustering is needed')
###     if (uniqueCount) { x <- uniqueCount(x); x <- x[,-c(1,ncol(x))] }
###     id <- clus@clus[[as.character(i)]]$id # If this come from merged clusters, cluster 0 is the one for elements from clusters below the minpoint number specified when creating
###     if (merged) idx <- as.numeric(names(table(id))) else idx <- as.numeric(names(table(id)[table(id)>=minpoints])) # For non-merged clusters, ignore cluster ID for clusters below minpoints
###     # Compute observed proportion of factor distribution in full dataset
###     tprop <- colSums(x)/nrow(x); names(tprop) <- colnames(x)
###     # Compute same proportion for each cluster and return ratios against tprop
###     cprop <- do.call(rbind,lapply(idx, function(y) {
###       xx <- x[id==y,]
###       ans <- colSums(xx)/nrow(xx)
###       ans <- ans / tprop
###       if (log2)
###         {
###           ans <- log2(ans)
###           ans[ans==-Inf] <- min(ans[is.finite(ans)]) - .1 # To remove -Inf
###           ans[ans==Inf] <- max(ans[is.finite(ans)]) + .1 # To remove +Inf
###         }     
###       ans
###     }))
###     if (merged) rownames(cprop) <- c('None',colnames(clus@clus[[as.character(i)]]$pointprob)[-1:-2]) else rownames(cprop) <- colnames(clus@clus[[as.character(i)]]$pointprob)[-1:-2]
###     cprop
###   }

profileClusters <- function(x,uniqueCount=TRUE,weights,clus,i,minpoints,merged=FALSE,log2=TRUE,plt=FALSE)
  # Computes enrichment or depletion of marks in a given cluster
  # X is table with epigenes / factors, h is a valid clustering for that matrix, k is a cluster cut configuration
  # n is used to consider only clusters with more than or n points, 0 for all points
  # Univ is used to compute average factor distribution, if TRUE all genes are considered, if FALSE only genes in clusters where ngenes>=n are considered
  # Result is given in the shape of a heatmap
  {
    if (missing(minpoints) & !merged) stop('For non-merged clusterings the minimum cluster size considered in the clustering is needed')
    if (uniqueCount) { x <- uniqueCount(x); x <- x[,-c(1,ncol(x))] }
    id <- clus@clus[[as.character(i)]]$id # If this come from merged clusters, cluster 0 is the one for elements from clusters below the minpoint number specified when creating
    if (merged) idx <- as.numeric(names(table(id))) else idx <- as.numeric(names(table(id)[table(id)>=minpoints])) # For non-merged clusters, ignore cluster ID for clusters below minpoints
    # Compute observed proportion of factor distribution in full dataset
    tprop <- colSums(x)/nrow(x); names(tprop) <- colnames(x)
    # Compute same proportion for each cluster and return log2 ratios against tprop    
    cprop <- do.call(rbind,lapply(idx, function(y) {
      xx <- x[id==y,]
      ans <- colSums(xx)/nrow(xx)
      if (log2)
        {
          #ans <- ifelse(ans>=tprop,ans/tprop,-1*(1/(ans/tprop)))
          ans <- log2(ans/tprop)
          ans[ans==-Inf] <- min(ans[is.finite(ans)]) - .1 # To remove -Inf
          ans[ans==Inf] <- max(ans[is.finite(ans)]) + .1 # To remove +Inf
          #ans <- sign(ans) * log2(abs(ans))
        }
      else ans <- ans-tprop
      ans
    }))
    if (merged) rownames(cprop) <- c('None',colnames(clus@clus[[as.character(i)]]$pointprob)[-1:-2]) else rownames(cprop) <- colnames(clus@clus[[as.character(i)]]$pointprob)[-1:-2]
    if (!missing(weights))
      {
        colnames(cprop) <- names(weights)
        cprop <- do.call(cbind,lapply(sort(unique(colnames(cprop))),function(i) rowMeans(as.data.frame(cprop[,colnames(cprop) %in% i]))))
        colnames(cprop) <- sort(unique(names(weights)))
      }
    cprop
  }

clusOverlap <- function(cluspoints,clusgenes,i,j,method='unweighted',overlap=TRUE)
# Compute cluster overlap, that is the correct classification rate just for clusters i and j
  {
  if (method=='unweighted') {
    sel.i <- clusgenes[,i]==1 # Points from cluster i
    sel.j <- clusgenes[,j]==1 # Points from cluster j
    pclus1 <- mean(cluspoints[sel.i,i] / rowSums(cluspoints[sel.i,c(i,j)]))
    pclus2 <- mean(cluspoints[sel.j,j] / rowSums(cluspoints[sel.j,c(i,j)]))
    clus.ov <- (pclus1 + pclus2) / 2
    # This is the CCR for two clusters, the higher the CCR, the better separated the clusters are
    if (overlap) clus.ov <- 2*(1-clus.ov) # Transform CCR in two cluster overlap, 0 no overlap at all, 1 total overlap
  }
  else stop('Currently only unweighted overlap method is supported')
  return(clus.ov)
}

setGeneric("mergeClusters", function (clus, clus.method = "unweighted", cpt.method = "mean", logscale = TRUE, brake = rep(1,length(clus@clus)), plt = TRUE, mc.cores = 1) standardGeneric("mergeClusters"))

setMethod("mergeClusters", signature=c(clus="clusGPS"), function (clus, clus.method = "unweighted", cpt.method = "mean", logscale = TRUE, brake = rep(1,length(clus@clus)), plt = TRUE, mc.cores = 1) {
  #names(brake) <- names(clus@clus)
  cnames <- names(clus@clus) # We iterate using indexes and not names, thus, save them
  clus@clus <- lapply(1:length(clus@clus),function(i) mergeClusters(clus@clus[[i]],clus.method=clus.method,cpt.method=cpt.method,logscale=logscale,brake=brake[i],plt=plt,mc.cores=mc.cores))
  names(clus@clus) <- cnames # Restore cluster names
  clus
})

setMethod("mergeClusters", signature=c(clus="list"), function (clus, clus.method = "unweighted", cpt.method = "mean", logscale = TRUE, brake = 1, plt = TRUE, mc.cores = 1) {
# Init of mergeClusters function code
 doMerge <- TRUE
  x <- clus$pointprob
  xx <- cluspoints <- cbind(x[, c("X", "Y")], x[, paste("C", names(sort(table(clus$id),decreasing = TRUE)), sep = "")[1:(ncol(x) - 2)]])
  #xx <- cluspoints <- cbind(x[, c("X", "Y")], x[, paste("C", order(table(clus$id),decreasing = TRUE), sep = "")[1:(ncol(x) - 2)]])
  #xx <- cluspoints <- x
  clusgenes <- cbind(xx[, 1:2], do.call(cbind, lapply(as.numeric(gsub("C","", colnames(xx)[-1:-2])), function(x) as.numeric(x == clus$id))))
  colnames(clusgenes) <- colnames(cluspoints)
  if (doMerge) { pden <- clus$pden[gsub("C", "", colnames(xx)[-1:-2])]; names(pden) <- paste("C", names(pden), sep = ""); mclus <- vector("list", ncol(cluspoints)); z <- 1 }
  maxOverlap <- vector("numeric", 0)
  clusovs <- vector("list", 0)
  while (ncol(cluspoints) > 3) {
    if ("parallel" %in% loadedNamespaces()) 
      clusov <- do.call(cbind, parallel::mclapply(colnames(cluspoints[,-1:-2]), function(x) unlist(lapply(colnames(cluspoints[,-1:-2]), function(y) clusOverlap(cluspoints, clusgenes, x, y, method = clus.method))), mc.cores = mc.cores, mc.preschedule=FALSE))
    else clusov <- do.call(cbind, lapply(colnames(cluspoints[,-1:-2]), function(x) unlist(lapply(colnames(cluspoints[,-1:-2]), function(y) clusOverlap(cluspoints, clusgenes, x, y, method = clus.method)))))
    if ((nrow(clusov) == ncol(clusov)) & nrow(clusov) >= 2) diag(clusov) <- -1
    maxOverlap <- c(maxOverlap, max(clusov))
    clusovs <- c(clusovs, clusov)
    i <- which(clusov == max(clusov), arr.ind = TRUE)[1,1] + 2
    j <- which(clusov == max(clusov), arr.ind = TRUE)[1,2] + 2
    sel.i <- clusgenes[, i] == 1
    sel.j <- clusgenes[, j] == 1
    # Update cluspoints
    cluspoints[, i] <- (cluspoints[, i] * (sum(sel.i)/(sum(sel.i) + sum(sel.j)))) + (cluspoints[, j] * (sum(sel.j)/(sum(sel.i) +  sum(sel.j))))
    colnames(cluspoints)[i] <- colnames(clusgenes)[i] <- paste(colnames(cluspoints)[i], colnames(cluspoints)[j], sep = ".")
    cluspoints <- cluspoints[, -c(j)]
    # Update clusgenes
    clusgenes[, i] <- as.numeric(clusgenes[, i] | clusgenes[,j])
    clusgenes <- clusgenes[, -c(j)]
    # Update pden if needed (only if doMerge==TRUE)
    if (doMerge) {
      pden[[i - 2]]$dens <- (pden[[i - 2]]$dens * (sum(sel.i)/(sum(sel.i) + sum(sel.j)))) + (pden[[j - 2]]$dens * (sum(sel.j)/(sum(sel.i) + sum(sel.j))))
      names(pden)[i - 2] <- colnames(cluspoints)[i]
      pden <- pden[-c(j - 2)]
      clus$pden <- pden
    }
    # Assign to current clus
    clus$pointprob <- cluspoints
    clus$id <- colSums(1:ncol(as.data.frame(cluspoints[, -1:-2])) * t(as.data.frame(clusgenes[,-1:-2])))
    clus$propclass <- lapply(colnames(as.data.frame(clusgenes)[-1:-2]), function(x) cluspoints[clusgenes[, x] == 1, x]/rowSums(as.data.frame(cluspoints[clusgenes[,x] == 1, -1:-2])))
    names(clus$pden) <- colnames(cluspoints)[-1:-2] # Fix bug in cluster names for density information
    if (doMerge) { mclus[[z]] <- clus; z <- z + 1 }
  }
  # Computing changepoint to identify 'cliff' in the maximum overlap between all possible pairwise merging
  if (cpt.method == "mean") cpt <- cpt.mean(log(maxOverlap),Q=1) else if (changepoint::cpt.method == "var") cpt <- changepoint::cpt.var(log(maxOverlap),Q=1)
  mycpt <- maxOverlap[[cpt@cpts[1] - brake]]
  mycpt.auto <- maxOverlap[[cpt@cpts[1] - 1]]
  if (plt==TRUE) {
    cptdata <- data.set(cpt)
    cptpoints <- cpts(cpt)
    changepoint::plot(cpt,pch=19,cex=1.5,ylim=c(min(cptdata),1),cpt.col=NA,xlab='Merging step',ylab=ifelse(logscale,'maxOverlap (log)','maxOverlap')) # Draw points but not cpt line
    text(x=1:length(maxOverlap),y=cptdata,round(maxOverlap,2),pos=3)
    text(x=.5*length(maxOverlap),y=0.5,'# of clusters by overlap threshold')
    axis(3,at=1:length(maxOverlap),labels=length(maxOverlap):1,padj=5,tick=FALSE)
    points(x = cptpoints[1], cex = 1.5, y = log(maxOverlap)[cptpoints[1]], col = "red", pch = 19);
    abline(v=cptpoints[1],lty=2,col='red')
  }
  return(mclus[[cpt@cpts[1] - brake]])
})
