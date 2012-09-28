\name{addVar}
\alias{addVar}
\title{Plot vector of a quantitative variable over a MDS map.}
\description{
  Given a quantitative variable as a numeric vector with one element for each point on a MDS map,
  calculate and plot the weight vector corresponding to that variable.
}
\usage{
addVar(mds1, z, plot = TRUE, label = "z", pos = 3, ...) 
}
\arguments{
  \item{mds1}{
    An object of class \code{mds} with the MDS object.
  }
  \item{z}{
    Numeric vector with the quantitative variable, one element for each
    point.
  }
  \item{plot}{
    Set to TRUE to calculate and draw the resulting vector on the MDS.
  }
  \item{label}{
    Something to be printed on the tip of the vector arrow, usually the
    name of the given variable.
  }
  \item{pos}{
    Graphical position where the label is drawn respect to the vector
    arrow tip.
  }
  \item{\dots}{
    Additional parameters given to the generic function \code{plot}.
  }
}
\value{
  A named list with the vector components.
}
\examples{
# Not run
# See chroGPS-manual.pdf for examples.
}
\keyword{multivariate,cluster}
