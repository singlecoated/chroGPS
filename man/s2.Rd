\name{s2}
\alias{s2}
\alias{s2names}
\alias{s2Seq}
\alias{s2SeqNames}
\alias{s2.tab}
\alias{s2.wt}
\alias{toydists}
\alias{d}
\alias{d2}
\alias{d3}
\docType{data}
\title{
  Sample binding site and related data from S2 and BG3 cell lines in Drosophila melanogaster.
}
\description{
  chroGPS example dataset including ChIP-CHIP (modEncode) and 
  ChIP-Seq (NCBI GEO GSE19325) data for Drosophila melanogaster S2 and BG3 cell
  lines as well as S2 wildtype gene expression values coming from Affymetrix
  Drosophila2 arrays. The object \code{toydists} stores precomputed
  \code{distGPS} objects (called \code{d, d2, d3}) for the epigenetic factors used in the dynamic
  vignette that comes with the package.
}
\usage{data(s2)}
\source{
  \code{http://www.modencode.org}
  \code{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19325}
}
\references{
  \code{http://www.modencode.org}
  \code{http://www.ncbi.nlm.nih.gov/geo/}
}
\examples{
data(s2)
class(s2)
s2
s2names$Factor
data(s2Seq)
s2Seq
# See vignette examples for several uses of these datasets.
}
\keyword{datasets}
\keyword{modEncode}
