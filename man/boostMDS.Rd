\name{boostMDS}
\alias{boostMDS}
\title{Improve goodness-of-fit of a given MDS solution in terms of R-square.}
\description{
  Given a distance matrix and a valid MDS representation for it, improve
  the R-square correlation between observed and approximated distances
  until converged is reached for a given threshold.
}
\usage{
boostMDS(D, Y, rate = 0.01, maxit = 50, tol = 0.001, samplesize,
verbose = TRUE, scale = FALSE, seed = 149, plt = FALSE, mc.cores = 1)
}
\arguments{
  \item{D}{
    Distance matrix.
  }
  \item{Y}{
    Matrix with points from a valid MDS solution for the distances in D.
  }
  \item{rate}{
    Grid step rate, start with 0.1 which usually is a good compromise,
    try also 0.01, 1, 10.
  }
  \item{maxit}{
    Maximum number of iterations.
  }
  \item{tol}{
    Tolerace for R-square convergence.
  }
  \item{samplesize}{
    When there are over 100 points to represent, the gradiend descent step size is determined
    using a fraction \code{samplesize} of the original data points.
    By default 0.01 with a minimum of 100 points, which typically gives
    very stable results. Setting large \code{samplesize} can
    significantly increase the computational cost.
  }
  \item{verbose}{
    Give details of the gains in R-square and step size.
  }
  \item{scale}{
    Whether to scale the MDS coordinates in the output MDS.
  }
  \item{seed}{
    A random seed to be used in the resampling process if samplesize < 1.
  }
  \item{plt}{
    Whether to plot the intermediate solutions or not.
  }
  \item{mc.cores}{
    Number of cores to use in parallelized grid step size search.
  }
}
\value{
  The function returns a matrix with the coordinates of a valid MDS
  solution for distance matrix D where the R-square correlation has been
  improved. However, have in mind that an MDS solution with better
  R-square does not necessarily mean the solution is easier to
  interpret. As with any MDS approach, a balance must be found between
  pure 'technical' goodness-of-fit and usefulness of the delivered
  solution in terms of answering the original hypothesis.
}
\references{
  boostMDS is based on hitMDS (High-Throughput Multidimensional Scaling,
  see see http://dig.ipk-gatersleben.de/hitmds/hitmds.html for details)
}
\examples{
# Not run, see also chroGPS-manual.pdf file for examples
#data(geneSample)
#d = distGPS(geneSample,uniqueRows=TRUE)
#m = mds(d,type='isoMDS')
#m
#plot(m)
#m = boostMDS(d@d,m@points)
#plot(m)
}
\keyword{multivariate,cluster}
