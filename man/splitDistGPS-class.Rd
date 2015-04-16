\name{splitDistGPS-class}
\Rdversion{1.1}
\docType{class}
\alias{splidDistGPS-class}
\alias{show,splitDistGPS-method}
\alias{splitDistGPS,distGPS-method}

\title{Class \code{"splitDistGPS"}}
\description{
  Set of pairwise distances between elements. This is an internal class 
  to be used with the parallel version of \code{mds}, and
  should not be used on its own.
}
\section{Objects from the Class}{
Objects from this class are used internally for parallel
Multidimensional Scaling. See \code{mds} for details.
}
\section{Slots}{
  \describe{
    \item{\code{d}:}{ List of \code{distGPS} objects. }
    \item{\code{size}:}{Object of class \code{"numeric"} indicating
      the size of the individual \code{distGPS} objects in the list. See function \code{mds}. }
    \item{\code{o}:}{Object of class \code{"numeric"} with the overlap
      (anchor points) between adjacent \code{distGPS} objects. See function \code{mds}. }
    \item{\code{shuffle}:}{Object of class \code{"numeric"}, deprecated. }
  }
}
\author{
  Oscar Reina
}
\examples{
showClass("splitDistGPS")
}
\keyword{classes}
