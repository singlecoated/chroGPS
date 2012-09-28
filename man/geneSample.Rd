\name{geneSample}
\alias{geneSample}
\alias{S2Controls}
\docType{data}
\title{
  Sample binding site data from Drosophila melanogaster S2-DSRC cell
  lines (ChIP-chip), coming from modEncode Release 21.
}
\description{
  Sample binding site data from Drosophila melanogaster S2-DSRC cell
  lines (ChIP-chip), coming from modEncode Release 21. Binding sites are
  annotated with the Drosophila melanogaster genome version on Biomart
  ENSEMBL 66 (February 2012) using the \code{annotatePeakInBatch}
  function from the \code{ChIPpeakAnno} package. A gap of 1000bp is
  considered for feature overlap, filtering the resulting features to
  keep only stricts overlaps within those 1000bp. 
}
\usage{data(geneSample)}
\source{
  \code{http://www.modencode.org}
}
\references{
  \code{http://www.modencode.org}
}
\examples{
data(geneSample)
## maybe str(geneSample) ; plot(geneSample) ...
}
\keyword{datasets}
