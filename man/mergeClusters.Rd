\name{mergeClusters}
\alias{mergeClusters}
\alias{mergeClusters-methods}
\alias{mergeClusters,clusGPS-method}
\alias{mergeClusters,list-method}
\title{
Unsupervised cluster merging based on their observed overlap with
automatic changepoint detection.
}
\description{
The function uses contour density estimation as computed by the
\code{clusGPS} function to merge significantly overlapping clusters in
an unsupervised manner. In
each step, clusters with highest overlap are merged, their individual
density estimates are updated in a computational feasible manner, and
the process continues until the maximum overlap between any given pair
of clusters drops swiftly, as detected by the \code{cpt.mean} function
in the \code{changepoint} package.
}
\usage{
mergeClusters(clus, clus.method = "unweighted", cpt.method = "mean", logscale = TRUE, brake = rep(1, length(clus@clus)), plt = TRUE, mc.cores = 1)
}
\arguments{
  \item{clus}{
A \code{clusGPS} object from which we want to merge clusters with a
significant overlap when visualized on a chroGPS MDS map. This is quite
useful when a clustering method (i.e. hierarchical clustering with
average linkage) tends to return a high number of overlapping clusters.
}
  \item{clus.method}{
Currently only 'unweighted' method is supported, that is, cluster
overlap is computed based on spatial location of contours, but the
computed overlaps are not weighted for cluster size.
}
  \item{cpt.method}{
Use 'mean' for using \code{cpt.mean} function in \code{changepoint}
package for computing overlap changepoint. Use 'var' for
\code{cpt.var}. See specific function help for details.
}
  \item{logscale}{
Defaults to TRUE. Whether to use decimal or log scale values for computing overlap changepoint.
}
  \item{brake}{
(Optional). By default, the function returns the clusters from the
optimal merging step as detected by the changepoint functions
(brake=1). By using smaller values (0, -1, -2, ...) or bigger ones (2,
3, 4, ...) the algorithm can be forced to return the result from any
previous or later merging step respectively.
}
  \item{plt}{
Set to TRUE to visualize maximum cluster overlap for each merging step
and changepoint detection (optimal merging step).
}
\item{mc.cores}{
Numbers of cores to use in parallel computation.
}
}
\value{
A \code{clusGPS} object where significantly overlapping clusters
are merged, highly improving visualization, cluster robustness and
further study of the epigenetic configuration of the chroGPS map.
}
\references{
Changepoint package from Killick et al, 2012.
}
\author{
Oscar Reina.
}
\seealso{
See documentation for package \code{changepoint}, \code{clusGPS} for
epigenetic cluster generation.
}
\examples{
# Not run
# data(s2)
# # Computing distances
# d <- distGPS(s2.tab,metric='tanimoto',uniqueRows=TRUE)
# # Creating MDS object
# mds1 <- mds(d,type='isoMDS')
# mds1
# plot(mds1)
# Precomputing clustering
# h <- hclust(as.dist(d@d),method='average')
# # Calculating densities (contours and probabilities), takes a while
# clus <- clusGPS(d,mds1,preMerge=TRUE,k=300) # Generating a high number of clusters
# clus <- mergeClusters(clus)
}
\keyword{ ~clustering }
\keyword{ ~changepoint }
