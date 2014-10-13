\name{clusGPS-class}
\Rdversion{1.1}
\docType{class}
\alias{clusGPS-class}
\alias{clusGPS-method}
\alias{hclust-class}
\alias{plot,clusGPS-method}
\alias{plot,clusGPS,ANY-method}
\alias{show,clusGPS-method}

\title{Class \code{"clusGPS"}}
\description{
Agglomerative Nesting for a distGPS object. Contains probability
contours and bayesian posterior probability of mis-classification for the clusters evaluated.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("clusGPS", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{h}:}{Object of class \code{"hclust"} with Agglomerative
      Nesting or user-provided cluster object. }
    \item{\code{clus}:}{Object of class \code{"list"} with probability
      contour and bayesian posterior probability of mis-classification
      information for the clusters evaluated. }
    \item{\code{adjusted}:}{Object of class \code{"logical"} indicating
      if bayesian posterior probabilities of mis-classification are
      adjusted for multiple testing. }
  }
}
\details{
  \describe{
  Parameters for the S4 plot method for \code{mds} objects.
  \item{\code{x}:}{
    Object of class \code{"mds"} with a 2D or 3D Multidimensional
    Scaling to be plotted.
  }
  \item{\code{drawlabels}:}{
    TRUE to use rownames of the MDS points as text labels.
  }
  \item{\code{labels}:}{
    Alternative character vector giving the text labels for the MDS
    points.  
  }
  \item{\code{plantar}:}{
    If a 3D MDS is used, set plantar to TRUE to plot projected views of
    the MDS using XY, YZ and XZ axis decomposition.
  }
  \item{\code{point.cex}:}{
    Size of the points / spheres for the MDS plot.
  }
  \item{\code{text.cex}:}{
    Size of text labels for the MDS points.
  }
  \item{\code{text.pos}:}{
    Alignment position of the text labels respective to its points (1,2,3,4).
  }
  \item{\code{point.col}:}{
    Color for the MDS points / spheres.
  }
  \item{\code{text.col}:}{
    Color for the MDS text labels.
  }
  \item{\code{point.pch}:}{
    PCH type for the MDS points.
  }
  \item{\code{type.3d}:}{
    Use 'p' for points, 's' for spheres.
  }
  \item{\code{radius}:}{
    Radius for the spheres on a 3D MDS plot. Automatically generated
    from point.cex and the number of points in the MDS.
  }
  \item{\code{app}:}{
    Appearance of the 3D spheres on a 3D MDS plot, can be 'fill',
    'lines', 'grid'.
  }
  \item{\code{alpha}:}{
    Number between 0 and 1 with the level of transparency to be used on
    spheres on a 3D MDS.
  }
  \item{\code{scalecol}:}{
    Set to TRUE to use a color scale for points, for instance to color
    points (genes) according to their expression level on a
    chroGPS-genes MDS plot.
  }
  \item{\code{scale}:}{
    Scale to use to generate scale colors (for instance normalized gene
    expression for each element (gene) on chroGPS-genes MDS).
  }
  \item{\code{palette}:}{
    Color palette to be used for scale colors.
  }
}}
\author{
  Oscar Reina
}
\examples{
showClass("clusGPS")
}
\keyword{classes}
