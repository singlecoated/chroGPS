\name{mds}
\alias{mds}
\alias{mds-methods}
\alias{mds,distGPS,missing-method}
\alias{mds,splitDistGPS,missing-method}
\alias{plot,mds,ANY-method}
\title{
  Metric and non-metric Multidimensional Scaling for a distGPS object. 
}
\description{
  Generation of Multidimensional Scaling objects for the dissimilarities
  between elements given as an input in a \code{distGPS} object. Metric
  and non-metric algorithms are available, as well as an optimization
  algorithm for improving r-square correlation between observed and
  approximated distances. The MDS calculation for a given distance
  matrix can be splitted into smaller individual tasks and run in
  parallel, greatly improving CPU time and system memory usage.
}
\section{Methods}{
  \describe{
    \item{mds}{\code{signature(d = "distGPS", m = "missing")}: Creates a
      \code{mds} object with points in a k-dimensional space approximating
      the pairwise distances in \code{d}. }
    \item{mds}{\code{signature(d = "distGPS", m = "mds")}: For the
      observed dissimilarities in \code{d} and a valid spatial
      representation of them in \code{m}, the function returns a
      \code{mds} object with an optimized representation of \code{d} in
      terms of R-square. See help for \code{boostMDS} for details.}
    \item{plot}{\code{signature(m = "mds")}: S4 plot method for
      \code{mds} objects. }
}}
\usage{
mds(d, m = NULL, k = 2, type = "classic", add = FALSE, cor.method = "pearson", splitMDS = FALSE, split = 0.05, overlap = 0.0025, reshuffle = TRUE, set.seed = 149, mc.cores = 1, ...)
}
\arguments{
  \item{d}{
    Object of class \code{distGPS} with the pairwise observed dissimilarities
    between elements, a distance matrix.
  }
  \item{m}{
    (Optional). Object of class \code{mds} with a MDS object generated from the
    distances in \code{d}. Only MDS type "boostMDS" is available. The mds
  function performs an optimization of the approximated distances in
  \code{m} in order to improve r-square correlation between them and the
  observed dissimilarities en \code{d}, maximizing goodness of fit.
}
  \item{k}{
    Dimensionality of the reconstructed space, typically set to 2
    or 3.
  }
  \item{type}{
    Set to \code{"classic"} to perform classical MDS (uses
    function \code{cmdscale} from package \code{stats}). Set to
    \code{"isoMDS"} to use Kruskal's non-metric MDS (uses function
    \code{isoMDS} from package \code{MASS}) Set to "boostMDS" to perform
    r-square optimization of a pre-computed input MDS for that distance
    matrix.
  }
  \item{add}{
    Logical indicating if an additive constant c* should be
    computed, and added to the non-diagonal dissimilarities such
    that all n-1 eigenvalues are non-negative in \code{cmdscale}.
  }
  \item{cor.method}{
    A character string indicating which correlation
    coefficient (or covariance) is to be computed.  One of "pearson"
    (default), "kendall", or "spearman", can be abbreviated.
  }
  \item{splitMDS}{
    Set to \code{TRUE} to perform computation of the MDS in parallel
    (see parameters below).
  }
  \item{split}{
    Proportion of elements to include in each (but last) distance submatrix.
  }
  \item{overlap}{
    Proportion of elements to be used as common anchor points between
    two adjacent distance submatrixes. These points will be used as
    spatial references to stitch each two adjacent MDS objects by Procrustes.
  }
  \item{reshuffle}{
    Set to TRUE to perform random resampling of the input distance
    matrix before splitting it for parallel computation. This is often
    necessary to sufficiently capture the inherent variability of the
    data in each distance submatrix so that the stitching process can
    work properly, as the original data may present an arbitrary sorting
    of some kind. If a previous resampling of the data has been
    performed, this is not necessary.
}
\item{set.seed}{
  Random seed to perform the resampling.
}
\item{mc.cores}{
  Number of cores to be passed to the \code{mclapply} function from
  the \code{multicore} package, used to perform the parallel MDS computations.
}
\item{\dots}{
  Additional parameters passed to \code{cmdscale}, \code{isoMDS} or
  \code{boostMDS}, see each individual help file for details.
}
}
\value{
  The function returns a \code{mds} object. See help ("mds-Class") for details.
}
\author{
Oscar Reina
}
\references{
See functions cmdscale, isoMDS from package \code{MASS}
}
\examples{
x <- rbind(c(rep(0,15),rep(1,5)),c(rep(0,15),rep(1,5)),c(rep(0,19),1),c(rep(1,5),rep(0,15)))
rownames(x) <- letters[1:4]
d <- distGPS(x,metric='tanimoto',uniqueRows=TRUE)
mds1 <- mds(d)
mds1
plot(mds1)
}
\keyword{ graphs }
\keyword{ mds }
