\name{domainDist}
\alias{domainDist}
\title{Overview of intra and inter-domain distances.}
\description{
  Given a distance of pairwise distances or dissimilarities between
  elements, return intra and inter-group sets of distances based on a
  given group definition. This is useful to get an insight on domain
  robustness for functional related genes or factors.
}
\usage{
domainDist(d, gps='factors', domain, type='intra', col='white', avg=FALSE,
plot=TRUE, ...)
}
\arguments{
  \item{d}{
    Distance/Dissimilarities matrix, usually the slot \code{d} on
    a \code{distGPS} object, but any distance matrix can be given as input.
  }
  \item{gps}{
    'factors' for a chroGPS-factors distance matrix, 'genes' for a
    chroGPS-genes one.
  }
  \item{domain}{
    Character vector with group identity for each element \code{d}. It can be a functional domain classification (i.e.
    'Activation', 'Repression', etc), given for each factor on a
    chroGPS-factors map or for each gene in a chroGPS-genes
    map. However, any classification of interest can be used (pathways,
    gene ontology, etc.)
  }
  \item{type}{
    Intradomain ('intra') or Interdomain ('inter') distance overview.
  }
  \item{col}{
    Character vector with colors to be passed to \code{plot}.
  }
  \item{avg}{
    TRUE to return also the average inter or intra domain distance.
  }
  \item{plot}{
    TRUE to generate inter/intra domain boxplots.
  }
  \item{\dots}{
    Additional parameters given to the generic function \code{plot}.
  }
}
\value{
  List of inter or intra domain distances.
}
\examples{
# Not run
# data(s2)
# d <- distGPS(s2,metric='avgdist',mc.cores=1)
# d.intra <- domainDist(as.matrix(d),domain=s2names$Color,type='intra',plot=TRUE)
# d.inter <- domainDist(as.matrix(d),domain=s2names$Color,type='inter',plot=TRUE)
}
\keyword{multivariate,cluster}
