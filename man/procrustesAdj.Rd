\name{procrustesAdj}
\alias{procrustesAdj}
\alias{procrustesAdj-methods}
\alias{procrustesAdj,mds,distGPS-method}
\title{Use Procrustes to adjust an MDS map containing samples obtained under different
  conditions, e.g. technology or genetic backgrounds.}
\description{
The function adjusts a previous \code{mds} to take into account that
samples were obtained under different conditions, e.g. technological or
genetic.
Pairwise adjustments are performed by identifying samples present in
both conditions and using Procrustes.
When there are more than two conditions, sequential pairwise adjustments
are applied (in the order that maximizes the number of common samples in
each pairwise adjustment).
}
\section{Methods}{
\describe{

\item{signature(x='mds')}{ \code{x} is a \code{mds} object with
  the results of an MDS analysis. }


}}
\usage{
procrustesAdj(mds1, d, adjust, sampleid)
}
\arguments{
  \item{mds1}{Object of class \code{mds} with a Multi-dimensional scaling analysis on a
    distance matrix, typically obtained by a previous call to \code{mds}.}
  \item{d}{Object of class \code{distGPS} with the matrix used to create the Multidimensional Scaling
    object usually through a call to \code{mds}.}
  \item{adjust}{Vector indicating the adjustment factor, i.e. the
    condition under which each sample has been obtained.}
  \item{sampleid}{Vector containing the sample
    identifier. \code{sampleid} should take the same value for samples
    obtained under different conditions, as this is used to detect
    the samples to be used for Procrustes adjustment.}
}
\value{
  Adjusted \code{mds} object. Have in mind that only original distances
  between samples obtained under the same condition should be conserved,
  as the adjusted distances manipulated by Procrustes no longer
  correlate with the distances between their points in the adjusted MDS.
}
\details{
  We implement the Procrustes adjustment as follows.
  First we identify common samples, i.e. those obtained both under conditions A and B.
  Second, we use Procrustes to estimate the shift, scale and rotation that best
  matches the position of the samples in B to those in A.
  If only 1 sample was obtained under both conditions, only the shift
  is estimated.
  Last, we apply the estimated shift, scale and rotation to all B
  samples.
  That is, the Procruses parameters are estimated using common samples
  only, which are then applied to all samples to perform the adjustment.

  Notice that the R square of the adjusted \code{mds} is typically
  improved after Procrustes adjustment, since distances between samples
  obtained under different conditions are set to \code{NA} and therefore
  MDS needs to approximate distances between less points.
  
  When several replicates are available for a given \code{sampleid}
  under the same condition (\code{adjust}), the average position of all
  replicates is used.
}
\seealso{ \code{\link{distGPS}} for computing distances,
  \code{\link{mds}} to create MDS-oriented objects.
}
\examples{
st1 <- runif(100,1,1000); st2 <- runif(100,500,1500)  #Peak starts
st3 <- runif(100,1000,2000); st4 <- runif(100,1500,2000)
#cond1: more precise technology
cond1 <- RangedDataList(s1=RangedData(IRanges(st1,st1+100)),s2=RangedData(IRanges(st2,st2+100)),s3=RangedData(IRanges(st3,st3+100)))
#cond2: less precise
cond2 <- RangedDataList(s1=RangedData(IRanges(st1-200,st1+300)),s2=RangedData(IRanges(st2-200,st2+300)),s5=RangedData(IRanges(st4-200,st4+300)))
x <- c(cond1,cond2)
d <- distGPS(x,metric='tanimoto')  #compute distances
mds1 <- mds(d)  #MDS
#Adjust via Procrustes
mds2 <- procrustesAdj(mds1,d,adjust=rep(c('seq','chip'),each=3),sampleid=names(x))
plot(mds1)
plot(mds2)
#Adjust via peak width
xadj <- adjustPeaks(x,adjust=rep(c('seq','chip'),each=3),sampleid=names(x))
dadj <- distGPS(xadj)
mds3 <- mds(dadj)
plot(mds3)
}
\keyword{multivariate,cluster}
