\name{distGPS-class}
\Rdversion{1.1}
\docType{class}
\alias{distGPS-class}
\alias{show,distGPS-method}
\alias{splitDistGPS-class}

\title{Class \code{"distGPS"}}
\description{
  Pairwise distances between elements. Function \code{distGPS} creates
  objects of this class.
  \code{splitDistGPS} in an private class used internally for parallel
  Multidimensional Scaling.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("distGPS", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{d}:}{Object of class \code{"matrix"} with pairwise
      dissimilarities (distances) between elements. }
    \item{\code{metric}:}{Object of class \code{"character"} indicating
      the metric type used for calculating distances. See function \code{distGPS}. }
    \item{\code{type}:}{Object of class \code{"character"}, deprecated. }
  }
}
\author{
  Oscar Reina
}
\examples{
showClass("distGPS")
}
\keyword{classes}
