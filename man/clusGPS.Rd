\name{clusGPS}
\alias{clusGPS}
\alias{clusGPS-methods}
\alias{clusGPS,distGPS,mds-method}
\alias{clusGPS,distGPS,missing-method}
\alias{contour2dDP}
\alias{contour3dDP}
\alias{plotContour}
\title{
Hierarchical clustering (Agglomerative Nesting) of elements based on their input pairwise distances.
}
\description{
Elements are clustered using Agglomerative Nesting of their pairwise
dissimilarities (distances). Additionally, semiparametric bayesian
density is estimated using a Dirichlet process mixture of normals. This
is used both to compute bayesian mis-classification posterior
probabilities and to estimate probability contours for each cluster in any given cut of
the clustering tree.

contour2dDP, contour3dDP and plotContour functions can be used to compute and plot bayesian density probability contours for a given set of elements (points)
from a pre-generated MDS object. These functions are used internally by clusGPS to draw cluster contours but are also useful to visualize other type of
contours over the map (ie genes from a given Gene Ontology term, having a specific epigenetic mark of interest, etc).
}
\section{Methods}{
\describe{
\item{signature(d='distGPS',m='mds')}{ Agglomerative nesting is
  performed for the elements whose pairwise distances are given in
  \code{d}. For each cluster partition given in \code{k}, cluster
  identity for each element is returned, and semiparametric
  bayesian density estimation is computed using the point density
  information from \code{m}.
}
\item{signature(d='distGPS',m='missing')}{ Agglomerative nesting is
  performed for the elements whose pairwise distances are given in
  \code{d}. For each cluster partition given in \code{k}, cluster
  identity for each element is returned.}
\item{plot}{\code{signature(m = "clusGPS")}: S4 plot method for
  \code{mds} objects. }
}}
\usage{
clusGPS(d, m, h, grid, ngrid=1000, densgrid=TRUE, type = "agnes", method =
"average", samplesize = 1, p.adjust = TRUE, k = 2:5, mc.cores = 1,
set.seed = 149, verbose=TRUE, ...)
contour2dDP(x, ngrid, grid = NULL, probContour = 0.5, xlim, ylim, 
    labels = "", labcex = 0.01, col = colors()[393], lwd = 4, 
    lty = 1, contour.type = "single", contour.fill = FALSE, ...)
contour3dDP(x, col = "red", probContour = 0.5, ngrid = 30, contour.type = "none") 
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
    \code{m} in order to improve r-square correlation between them and the
  observed dissimilarities en \code{d}, maximizing goodness of fit.
}
\item{h}{
  (Optional). Object of class \code{hclust} with a pre-calculated clustering for the
  elements in \code{d}.
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
  \item{type}{
    Type of clustering to be performed. Currently only "agnes"
    (Agglomerative Nesting) is supported, but any other clustering type
    can be used by providing a pre-calculated object \code{h}.
  }
  \item{method}{
    Clustering method. See \code{\link{agnes}} for details.
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
    \code{multicore} package.
  }
  \item{set.seed}{
    If samplesize<1, random seed to be used to perform random sampling
    of the data.
  }
  \item{verbose}{
    Set to TRUE to output clustering process information.
  }
  \item{x}{
    Numeric matrix indicating coordinates of the points for which a
    probability contour is calculated in contour2dDP and contour3dDP.
  }
  \item{probContour}{
    Numeric matrix indicating coordinates of the points for which a
    probability contour is calculated in contour2dDP and contour3dDP.
  }
  \item{contour.type}{
    For contour2dDP, type of contour, either 'single' (surrounding the points within the
    given probContour probability) or 'multiple' to generate 
    terrain-like density contour lines. For contour3dDP currently only
    'contours' (to drawn contours as tridimensional polygons) and 'none'
    are available.
  }  
  \item{contour.fill}{
    Deprecated.
  }
  \item{xlim,ylim,labels,labcex,col,lwd,lty}{
    Graphical parameters given to contour2dDP.
  }
  \item{\dots}{
    Additional parameters.
  }
}
\value{
Object of class \code{clusGPS}. See help for \code{clusGPS-methods} for details.  
}
\author{
Oscar Reina
}
\examples{
# Not run
# data(geneSample)
# # Computing distances
# d <- distGPS(geneSample,metric='tanimoto',uniqueRows=TRUE)
# # Creating MDS
# mds1 <- mds(d)
# mds1
# plot(mds1)
# # Calculating densities (contours and probabilities), takes a while
# c=clusGPS(d,mds1,k=2:3,mc.cores=1)
# # c contains information for contours and probabilities
# plot(c,type='contours',k=3,lwd=3,probContour=.75)
# plot(c,type='stats',k=3)
# plot(c,type='avgstat',k=3)
# plot(c,type='density',k=3,ask=TRUE,xlim=range(mds1@points),ylim=range(mds1@points))
}
\keyword{ cluster }
