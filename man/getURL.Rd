\name{getURL}
\alias{getURL}
\title{Retrieve file from URL.}
\description{
  A function that can be used to retrieve any file of interest from the
  internet, in our case, modEncode binding site information GFF files
  into the working directory. See also help for function \code{gff2RDList}.
}
\usage{
getURL(urls, filenames, extension='.gff3', method='internal') 
}
\arguments{
  \item{urls}{
    Character vector with one or more target URLs to download.
  }
  \item{filenames}{
    Character vector with the filename for each URL target.
  }
  \item{extension}{
    If desired, an extension to append to filenames.
  }
  \item{method}{
    Either 'internal' to use the system's default or 'wget' if it is installed.
  }
}
\value{
  Message indicating the path to downloaded file(s).
}
\examples{
# Not run
#getURL('http://www.google.com/index.html','index','.html')
}
\keyword{multivariate,cluster}
