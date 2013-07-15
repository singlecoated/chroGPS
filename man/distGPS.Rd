\name{distGPS}
\alias{distGPS}
\alias{distGPS-methods}
\alias{distGPS,RangedDataList-method}
\alias{distGPS,data.frame-method}
\alias{distGPS,matrix-method}
\alias{splitDistGPS,data.frame-method}
\alias{splitDistGPS,matrix-method}
\alias{uniqueCount}
\title{Compute matrix with pairwise distances between objects. Several
  GPS metrics are available.}
\description{
The function computes pairwise distances between invididuals
(e.g. samples or genes) according to a user-specified metric.
Several metrics are available. The precise definition of each metric
depends on the class of the first argument (see details section).
}
\section{Methods}{
\describe{
  
  distGPS:
  
  \item{signature(x='RangedDataList')}{ Each element in \code{x} is
  assumed to indicate the binding sites for a different sample,
  e.g. epigenetic factor. Typically \code{space(x)} indicates the
  chromosome, \code{start(x)} the start position and \code{end(x)} the
  end position (in bp). Strand information is ignored.  }

\item{signature(x='matrix')}{ Rows in \code{x} contain individuals for
  which we want to compute distances. Columns in \code{x} contain the
  variables, and should only contain either 0's and 1's or \code{FALSE}
  and \code{TRUE}.}

splitDistGPS:

This is a set of internal classes and functions to be used in the
parallel computation of Multidimensional Scaling.

uniqueCount:

This function collapses a chroGPS-genes matrix or data frame so that
elements with the same combination of variables are aggregated into a
single entry. Elements become then identified by their unique pattern
and a frequency count is also returned.

}}
\usage{
distGPS(x, metric='tanimoto', uniqueRows=FALSE, genomelength=NULL, mc.cores=1)
}
\arguments{
  \item{x}{Object for which we want to compute distances}
  \item{metric}{Desired distance metric. Valid options for
    chroGPS-factors map are 'tanimoto',
    'avgdist', 'chisquare' and 'chi' (see details). For chroGPS-genes
    maps, metrics 'wtanimoto', 'euclidean' and 'manhattan' are also
    available. }
  \item{uniqueRows}{If set to \code{TRUE} and \code{x} is a
    \code{matrix} or \code{data.frame}, duplicated rows are removed
    prior to distance calculation. This can save substantial computing
    time and memory. Notice however that the dimension of the distance
    matrix is equal to the number of unique rows in \code{x}, instead of
    \code{nrow..
      (x)}.}
  \item{genomelength}{For 'chi' and 'chisquare' metrics, numeric value
    indicating the length of the genome. If not given the function
    uses the minimum length necessary to fit the total length of the result.}
  \item{mc.cores}{If \code{mc.cores}>1 and \code{parallel} package is
    loaded, computations are performed in parallel with \code{mc.cores}
    processors when possible.}
}
\value{
  Object of class \code{distGPS}, with matrix of pairwise
  dissimilarities (distances) between objects.
}
\details{
  For \code{RangedDataList} objects, distances are defined as follows.
  
  Let \code{a1} and \code{a2} be two \code{RangedData} objects.
  Define as \code{n1} the number of \code{a1} intervals overlapping with
  some interval in \code{a2}. Define \code{n2} analogously.
  The Tanimoto distance between \code{a1} and \code{a2} is defined as 
  \code{(n1+n2)/(nrow(z1)+nrow(z2))}.
  The average distance between \code{a1} and \code{a2} is defined as
  \code{.5*(n1/nrow(z1) + n2/nrow(z2))}.
  The wtanimoto distance in chroGPS-genes weights each epigenetic factor
  (table columns) according to its frequency (table rows).
  The chi-square distance is defined as the usual chi-square distance on
  a binary matrix \code{B} which is automatically computed by
  \code{distGPS}. 
  The binary matrix \code{B} is the
  matrix with \code{length(x)} rows and number of columns equal to the
  genome length, where \code{B[i,j]==1} indicates that element \code{i}
  has a binding site at base pair \code{j}.
  The chi distance is simply defined as the square root of the
  chi-square distance.
  Finally, euclidean and manhattan metrics have the same definition than
  in the base R function \code{dist}.

  When choosing a metric one should consider the effect of outliers,
  i.e. samples with large distance to all other samples.
  Tanimoto and Average Distance take values between 0 and 1, and
  therefore outlying distances have a limited effect.
  Chi-square and Chi distances are not limited between 0 and 1,
  i.e. some distances may be much larger than others. The Chi metric is
  slightly more robust to outliers than the Chi-square metric.
  
  For \code{matrix} or \code{data.frame} objects, \code{x} must be a
  matrix with 0's and 1's (or \code{FALSE} and \code{TRUE}).
  The usual definitions
  are used for Tanimoto (which is equivalent to Jaccard's index),
  Chi-square and Chi.
  Average overlap between rows \code{i} and \code{j} is simply the
  average between the proportion of elements in \code{i} also in
  \code{j} and the proportion of elements in \code{j} also in \code{i}.

}
\seealso{ \code{\link{mds}} to create MDS-oriented objects, \code{\link{procrustesAdj}} for
  Procrustes adjustment. }
\examples{
x <- rbind(c(rep(0,15),rep(1,5)),c(rep(0,15),rep(1,5)),c(rep(0,19),1),c(rep(1,5),rep(0,15)))
rownames(x) <- letters[1:4]
d <- distGPS(x,metric='tanimoto')
du <- distGPS(x,metric='tanimoto',uniqueRows=TRUE)
mds1 <- mds(d)
mds1
plot(mds1)
d <- distGPS(x,metric='chisquare')
mds1 <- mds(d)
mds1
plot(mds1)
}
\keyword{multivariate}
\keyword{clustering}
