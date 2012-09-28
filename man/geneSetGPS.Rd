\name{geneSetGPS}
\alias{geneSetGPS}
\alias{geneSetGPS-methods}
\alias{geneSetGPS,data.frame,mds,character-method}
\alias{geneSetGPS,matrix,mds,character-method}
\title{
  Highlight point (gene) position over a Multi-dimensional Scaling plot.
}
\description{
Given a list of genes of interest, the function highlights their position over a Multi-dimensional Scaling plot.
}
\usage{
geneSetGPS(x, m, genes, uniqueCount = TRUE, ...)
}
\arguments{
  \item{x}{
    Matrix or data frame of observations x variables (tipically genes x
    epigenetic factors), with gene identifiers as rownames.
  }
  \item{m}{
    Object of class \code{mds} with a valid Multidimensional Scaling
    representation for the elements in \code{x}.
  }
  \item{genes}{
    Character vector containing gene identifiers, matching those on
    rownames(x).
  }
  \item{uniqueCount}{
    Set to FALSE if the MDS has been generated directly from the data in
    \code{x}, otherwise set to TRUE to match gene identifiers with their
    unique pattern of observed variables.
  }
  \item{\dots}{
    Additional parameters given to the generic function \code{plot}.
}
}
\value{
  Matrix with coordinates on the given input MDS object for the genes selected.
}
\author{
Oscar Reina
}
\examples{
# Not run, see further examples in chroGPS-manual.pdf document
# data(geneSample)
# # Non-unique rows
# d <- distGPS(geneSample,metric='tanimoto',uniqueRows=FALSE)
# mds1 <- mds(d)
# pts <- geneSetGPS(geneSample,mds1,c('FBgn0051048','FBgn0010406','FBgn0034642'),uniqueCount=FALSE)@points
# plot(mds1)
# points(pts,col='red',cex=3)
# # Unique rows
# d <- distGPS(geneSample,metric='tanimoto',uniqueRows=TRUE)
# mds1 <- mds(d)
# plot(mds1)
# pts <- geneSetGPS(geneSample,mds1,c('FBgn0051048','FBgn0010406','FBgn0034642'),uniqueCount=TRUE)@points
# points(pts,col='red',cex=3)
}
\keyword{ graphics }
