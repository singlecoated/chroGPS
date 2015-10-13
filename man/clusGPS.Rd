\name{clusGPS}
\alias{clusGPS}
\alias{clusGPS-methods}
\alias{clusGPS,distGPS,mds-method}
\alias{clusNames}
\alias{tabClusters}
\alias{clusterID}
\alias{clusNames,clusGPS-method}
\alias{tabClusters,clusGPS-method}
\alias{clusterID,clusGPS-method}
\alias{contour2dDP}
\alias{plotContour}
\title{
Computation of cluster density estimates for cluster contour representation and correct-classification rates (cluster robustness). A pre-computed clustering of elements used in the map has to be given as an input, which is useful to explore results using different clustering
algorithms and methodologies (top-down, bottom-up, etc).
}
\description{
After performing a pre-merging step so that all clusters have a minimum size, semiparametric bayesian
density is estimated using a Dirichlet process mixture of normals. This
is used both to compute bayesian mis-classification posterior
probabilities (correct classification rates) and to estimate probability contours which can be
visualized on the MDS map.

The functions \code{contour2dDP} and \code{plotContour} functions can be used to compute bayesian density estimates for a given set of elements (points) 
from a pre-generated 2D MDS object. These functions are used internally by clusGPS to draw cluster contours but are also useful to visualize other type of
contours over the map (ie genes from a given Gene Ontology term, having
a specific epigenetic mark of interest, etc).

The S4 accessors \code{clusNames,tabClusters} and \code{clusterID}
retrieve information stored within a \code{clusGPS} object.
}
\section{Methods}{
\describe{
\item{signature(d='distGPS',m='mds')}{ Hierarchical clustering is
  performed for the elements whose pairwise distances are given in
  \code{d}. For each cluster partition given in \code{k}, cluster
  identity for each element is returned, and semiparametric
  bayesian density estimation is computed using the point density
  information from \code{m}.
}
\item{plot}{\code{signature(m = "clusGPS")}: S4 plot method for
  \code{clusGPS} objects. 
}
\item{clusNames}{\code{signature(m = "clusGPS")}: Retrieves names of the
  clustering configurations stored in \code{clusGPS} objects, one for
  each distance threshold indicated in \code{k}, that get automatically named accordingly.
}
\item{tabClusters}{\code{signature(m = "clusGPS")}: Returns a table with
  the number of elements in each of the clusters found for an existing
  clustering configuration with name \code{name} within the \code{clusGPS} object. 
}
\item{clusterID}{\code{signature(m = "clusGPS")}: Returns a vector of
  cluster assignments for all the elements in an existing
  clustering configuration \code{name} within the \code{clusGPS} object.
  }}}
\usage{
clusGPS(d, m, h, sel=NULL, id=NULL, grid, ngrid=1000, densgrid=FALSE, preMerge=TRUE, type = "hclust", method =
"average", samplesize = 1, p.adjust = TRUE, k, mc.cores = 1,
set.seed = 149, verbose=TRUE, minpoints=70,...)
contour2dDP(x, ngrid, grid = NULL, probContour = 0.5, xlim, ylim, 
    labels = "", labcex = 0.01, col = colors()[393], lwd = 4, 
    lty = 1, contour.type = "single", contour.fill = FALSE,
minpoints=100, ...)
clusNames(clus)
tabClusters(clus,name)
clusterID(clus,name)
}
\arguments{
  \item{d}{
    Object of class \code{distGPS} with the pairwise observed dissimilarities
    between elements.
  }
  \item{m}{
    (Optional). Object of class \code{mds} with a MDS object generated from the
    distances in \code{d}. Only MDS type "boostMDS" is available. The mds
    function performs an optimization of the approximated distances in
    \code{m} in order to improve R-square correlation between them and the
    observed dissimilarities en \code{d}, maximizing goodness of fit.
}
\item{h}{
  (Optional). Object of class \code{hclust} with a pre-calculated clustering for the
  elements in \code{d}.
}
\item{sel}{
  (Optional). Logical vector indicating which elements from \code{d} will be used for performing hierarchical clustering with average linkage. This is useful if we want to focus on a given set of points only (i.e. those from a big cluster which we want to study in more detail).
}
\item{id}{
  (Optional). Label of the cluster which we want to further subdivide, ignoring points from all other clusters. Deprecated, use parameter \code{sel} specified above.
}
\item{grid}{
    Matrix of dimension ngrid*nvar giving the diagonal points of the
    grid where the density estimate is evaluated. The default value is NULL: grid dimensions are
    chosen according to the range of the data, and granularity is
    automatically determined according to data density, in order to
    provide a more accurate estimation in high density areas, where more
    resolution is needed.
  }
  \item{ngrid}{
    Number of grid points where the density estimate is evaluated. This
    argument is ignored if a grid is specified. The default value is
    1000. Higher values are recommended if data presents very high
    density areas.
  }
  \item{densgrid}{
    Set to true to generate grid points from the quantile distribution
    of the data using the grid size defined by ngrid. This is useful if
    the data presents areas of very different density, ranging from very
    sparse to extremely dense areas, optimizing grid
    granularity where is necessary, therefore improving resolution of density
    estimation and reducing computation time. 
  }
  \item{preMerge}{
    If TRUE will perform a first pre-merging step so that any cluster smaller than \code{minpoints} gets merged with its closest cluster based on their centroid distances. This is performed until no clusters < minpoints exist.
  }
  \item{type}{
    Type of clustering to be performed. Currently only "hclust"
    (Agglomerative Nesting) is supported, but any other clustering type
    can be used by providing a pre-calculated object \code{h}. This variable is to become deprecated, since clusGPS will only work with a precomputed clustering.
  }
  \item{method}{
    Clustering method. See \code{\link{hclust}} for details. This variable is to become deprecated, since clusGPS will only work with a precomputed clustering.
  }
  \item{samplesize}{
    Proportion of elements to sample for computing clustering and
    density estimation. This is useful to generate density contours from
    a subset of the data, speeding up computation.
  }
  \item{p.adjust}{
    Set to TRUE to adjust the bayesian posterior probabilities of mis-classification.
  }
  \item{k}{
    Integer vector indicating the number of clusters on which density
    estimation will be computed for mis-classification or contour calculation.
  }
  \item{mc.cores}{
    Number of cores to be used for parallel computation with the
    \code{parallel} package.
  }
  \item{set.seed}{
    If samplesize<1, random seed to be used to perform random sampling
    of the data.
  }
  \item{verbose}{
    Set to TRUE to output clustering process information.
  }
  \item{minpoints}{
    If preMerge is FALSE, then the algorithm will ignore clusters with fewer than \code{minpoints} elements. This is
    useful if the clustering method used tends to generate many very
    small clusters of limited use and difficult interpretation and for which density estimates may not be correctly computed. The default method is to preMerge clusters since this ensures density estimation is available for all clusters and helps interpreting the map, since no elements are ignored.
  }
  \item{x}{
    Numeric matrix indicating coordinates of the points for which a
    probability contour is calculated in contour2dDP.
  }
  \item{probContour}{
    Numeric matrix indicating coordinates of the points for which a
    probability contour is calculated in contour2dDP.
  }
  \item{contour.type}{
    For contour2dDP, type of contour, either 'single' (surrounding the points within the
    given probContour probability) or 'multiple' to generate 
    terrain-like density contour lines.
  }  
  \item{contour.fill}{
    Deprecated.
  }
  \item{xlim,ylim,labels,labcex,col,lwd,lty}{
    Graphical parameters given to contour2dDP.
  }
  \item{clus}{
    A valid \code{clusGPS} object from which we want to extract information.
  }
  \item{name}{
    Character indicating a valid name within a \code{clusGPS} object,
    from which we want to extract information.
  }
  \item{\dots}{
    Additional parameters.
  }
}
\value{
The function \code{clusGPS} returns an object of class
\code{clusGPS}. See help for \code{clusGPS-methods} for
details. \code{contour2dDP} returns a \code{DPdensity} object with
density contour information which can be plotted as 2D contours with our
\code{plotContour} function, as well as with the \code{plot} function
from the \code{DPpackage} package.
}
\author{
Oscar Reina
}
\examples{
# Not run
# data(s2)
# # Computing distances
# d <- distGPS(s2.tab,metric='tanimoto',uniqueRows=TRUE)
# # Creating MDS object
# mds1 <- mds(d,type='isoMDS')
# mds1
# plot(mds1)
# Precomputing clustering
# h <- hclust(as.dist(d@d),method='average')
# # Calculating densities (contours and probabilities), takes a while
# clus <- clusGPS(d,mds1,preMerge=TRUE,k=max(cutree(h,h=0.5)))
# # clus contains information for contours and probabilities
# plot(clus,type='contours',k=125,lwd=3,probContour=.75)
# plot(clus,type='stats',k=125,ylim=c(0,1))
# plot(clus,type='avgstat')
# plot(clus,type='density',k=3,ask=TRUE,xlim=range(mds1@points),ylim=range(mds1@points))
}
\keyword{ cluster }
