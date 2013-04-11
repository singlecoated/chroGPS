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
# Not run
# data(s2)
# d <- distGPS(s2.tab,metric='tanimoto',uniqueRows=TRUE)
# mds1 <- mds(d)
# set.seed(149)
# sampleGenes <- rownames(s2.tab)[sample(1:nrow(s2.tab),10,rep=FALSE)]
# pts <- geneSetGPS(s2.tab,mds1,genes=sampleGenes,uniqueCount=TRUE)
# plot(mds1)
# points(getPoints(pts),col='red',cex=3)
}
\keyword{ graphics }
