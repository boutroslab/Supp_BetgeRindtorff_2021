\name{matList}

\alias{matList}

\docType{data}

\title{List of gene expression matrices for CRIS classification}

\description{A list containing example gene expression data
  for CRIS classification using the kTSP classifier.
}

\usage{
data(matList)
}

\format{
  The \code{matList$Training} matrix contains gene epressions
  for 80 genes and 416 samples.
  The \code{matList$Testing} matrix contains gene epressions
  for 72 genes and 208 samples.
  Row names for these two matrices are valid gene symbols
  corresponding to those used in the kTSP classifier.
}


\seealso{
  \code{\link{phenoList}}
}

\author{
  Luigi Marchionni \email{marchion@jhu.edu}
}

\examples{

### Load gene expression datasets
data(matList)

### Show the class and names of the ``matList'' object
class(matList)
names(matList)

### Show the dimentions of the ``matList'' matrices
sapply(matList, dim)

}

\keyword{datasets}
