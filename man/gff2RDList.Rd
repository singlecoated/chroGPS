\name{gff2RDList}
\alias{gff2RDList}
\title{Retrieve binding site information from
  GFF3 files. }
\description{
An auxiliary function to retrieve binding site information from
  GFF3 format files (for instance those downloaded from modEncode, see
  function \code{getURL}).}
\usage{
gff2RDList(filenames,listnames,dir,quote=NULL,chrprefix='') 
}
\arguments{
  \item{filenames}{
    GFF3 filenames to read.
  }
  \item{listnames}{
    Names for each read filename, will be used as names of the returned
    \code{RangedDataList}. If not given, filenames will be used as listnames.
  }
  \item{dir}{
    Directory where the GFF3 files are located.
  }
  \item{quote}{
    Quote character used in the GFF3 files around fields, if any.
  }
  \item{chrprefix}{
    Prefix to be placed before the chromosome names if desired, for
    instance 'chr'.
  }
}
\value{
  A list with Enriched and Depleted binding sites, each one is an object of class \code{RangedDataList} with the \code{RangedData}
  objects containing the respective enriched or depleted binding sites from each GFF3 file.
}
\examples{
# Not run
#getURL('http://intermine.modencode.org/release-30/features.do?type=submission&action=export&format=gff3&submission=modENCODE_2984&feature=BindingSite&UCSC','test')
#test <- gff2RDList('test.gff3',dir=getwd())
#test
#test$Enriched[[1]]
#test$Depleted[[1]]
}
\keyword{multivariate,cluster}
