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
# Unadjusted map
data(s2)
data(s2Seq)
data(toydists) # precomputed distances
# d2 <- distGPS(c(reduce(s2),reduce(s2Seq)),metric='avgdist') # not run
mds2 <- mds(d2,k=2,type='isoMDS')
cols <- c(as.character(s2names$Color),as.character(s2SeqNames$Color))
sampleid <-
  c(as.character(s2names$Factor),as.character(s2SeqNames$Factor))
pchs <- rep(c(20,17),c(length(s2),length(s2Seq)))
point.cex <- rep(c(8,5),c(length(s2),length(s2Seq)))
par(mar=c(2,2,2,2))
plot(mds2,drawlabels=TRUE,point.pch=pchs,point.cex=point.cex,text.cex=.7,
point.col=cols,text.col='black',labels=sampleid,font=2)
#legend('topleft',legend=sprintf('R2=%.3f - %stress=%.3f',getR2(mds2),getStress(mds2)),bty='n',cex=1)
legend('topright',legend=c('ChIP-Chip','ChIP-Seq'),pch=c(20,17),pt.cex=c(1.5,1))

# Procrustes Adjusted map
adjust <- rep(c('chip','seq'),c(length(s2),length(s2Seq)))
sampleid <-
c(as.character(s2names$Factor),as.character(s2SeqNames$Factor))
mds3 <- procrustesAdj(mds2,d2,adjust=adjust,sampleid=sampleid)
par(mar=c(0,0,0,0),xaxt='n',yaxt='n')
plot(mds3,drawlabels=TRUE,point.pch=pchs,point.cex=point.cex,text.cex=.7,
point.col=cols,text.col='black',labels=sampleid,font=2)
#legend('topleft',legend=sprintf('R2=%.3f - %stress=%.3f',getR2(mds3),getStress(mds3)),bty='n',cex=1)
legend('topright',legend=c('ChIP-Chip','ChIP-Seq'),pch=c(20,17),pt.cex=c(1.5,1))

# Peak Width Adjusted map
s2.pAdj <-
adjustPeaks(c(reduce(s2),reduce(s2Seq)),adjust=adjust,sampleid=sampleid,logscale=TRUE)
# d3 <- distGPS(s2.pAdj,metric='avgdist')
mds4 <- mds(d3,k=2,type='isoMDS')
par(mar=c(0,0,0,0),xaxt='n',yaxt='n')
plot(mds4,drawlabels=TRUE,point.pch=pchs,point.cex=point.cex,text.cex=.7,
point.col=cols,text.col='black',labels=sampleid,font=2)
#legend('topleft',legend=sprintf('R2=%.3f - %s=%.3f',getR2(mds4),getStress(mds4)),bty='n',cex=1)
legend('topright',legend=c('ChIP-Chip','ChIP-Seq'),pch=c(20,17),pt.cex=c(1.5,1))

}
\keyword{multivariate,cluster}
