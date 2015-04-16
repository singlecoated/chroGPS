\name{gps2xgmml}
\alias{gps2xgmml}
\alias{gps2xgmml,mds-method}
\alias{gps2xgmml,mds,ANY-method}
\title{Export an 'mds' object to Cytoscape .xgmml format }
\description{
  gps2xgmml creates a .xgmml file for visualizing MDS results
  in Cytoscape. Two-dimensional MDS maps can be visualized in Cytoscape
  as usual. For three-dimensional maps Cytoscape's 3D Renderer
  (http://wiki.cytoscape.org/Cytoscape_3/3D_Renderer) is required.
}
\usage{
gps2xgmml(x, fname='out.xgmml', names.arg, fontSize=4, col=col2hex('steelblue'), cex)
}
\arguments{
  \item{x}{Object of class \code{mds}}
  \item{fname}{Name of output file}
  \item{names.arg}{Names for each point. If missing, they're taken from \code{x}. }
  \item{fontSize}{Font size}
  \item{col}{Fill colour(s) for the plotting symbols. Should be given in
  hexadecimal, e.g. as returned by function \code{col2hex} from
  \code{gplots}. Tips: col2hex('steelblue') looks nice in 2D/3D plots,
  col2hex('steelblue') looks nice in 2D plots and a bit faded on 3D plots.}
  \item{cex}{Expansion factor for plotting symbols. By default,
    \code{cex=12} for 2D plots and \code{cex=100} for 3D plots.}
}
\details{
  The .xgmml file contains the map co-ordinates in 2 or 3 dimensions,
  depending on the number of dimensions stored in the input \code{mds}
  object.
  To visualize properly a file with 3D co-ordinates, you need to install
  Cytoscape's 3D Renderer
  (http://wiki.cytoscape.org/Cytoscape_3/3D_Renderer) and start
  Cytoscape following the instructions provided therein.
  
  An .xgmml file with 3D co-ordinates can still be visualized in regular
  Cytoscape but the z-axis will be ignored.
}
\value{
 Generates an .xgmml file that can be opened in Cytoscape (File ->
 Import -> Network).
}
\examples{
#See help(mds) for an example
}
\keyword{manip}
