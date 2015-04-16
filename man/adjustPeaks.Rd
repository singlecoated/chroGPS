\name{adjustPeaks}
\alias{adjustPeaks}
\alias{adjustPeaks-methods}
\alias{adjustPeaks,RangedDataList-method}
\title{Adjust peak width so that samples obtained under different
  conditions become comparable.}
\description{
  Peaks obtained under different conditions (e.g. chip-chip, chip-seq,
  mnase-seq) are typically not comparable in terms of their width.
  \code{adjustPeaks} modifies the mean and SD of the peak width distribution
  for each condition, so that they become equivalent to the condition
  with widest peaks. See details.
}
\section{Methods}{
\describe{

\item{signature(x='RangedDataList')}{ Each element in \code{x} contains
  the binding sites for a different sample. The start, end and
  chrosomome of each binding sites should be accessed via  \code{start},
  \code{end} and \code{space}. }

}}
\usage{
adjustPeaks(x, adjust, sampleid, logscale = TRUE) 
}
\arguments{
  \item{x}{\code{RangedDataList} indicating the binding sites for each sample/experiment.}
  \item{adjust}{Vector indicating the adjustment factor, i.e. the
    condition under which each sample has been obtained.}
  \item{sampleid}{Vector containing the sample
    identifier. \code{sampleid} should take the same value for samples
    obtained under different conditions, as this is used to detect
    the samples to be used for Procrustes adjustment.}
  \item{logscale}{If set to \code{TRUE} the mean and SD are matched for
    log width, otherwise the original widths are used. Working in log
    scale can help reduce the effect of outliers (e.g. an usually long
    binding site).}
}
\value{
  \code{RangedDataList} object with adjusted widths.
}
\details{
  In a sense, the peak calling resolution is decreased so that they
  become comparable to the less precise technology
  (notice that there is no reliable way to increase the precision given
  by a low-resolution technology).
  
}
\seealso{ \code{\link{procrustesAdj}} for an alternative, more general, adjustment
  based on Procrustes. \code{\link{distGPS}} for computing distances,
  \code{\link{mds}} to create MDS-oriented objects.
}
\examples{
#See examples in help(procrustesAdj)
}
\keyword{multivariate,cluster}
