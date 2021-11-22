\name{phenoList}

\alias{phenoList}

\docType{data}

\title{List of factors containing CRIS predictions from the NTP classifier}

\description{A list containing CRIS predictions based on NTP for the
  example gene expression dataset available in the \code{\link{matList}}
  object.
}

\usage{
data(phenoList)
}

\format{
  The factor \code{phenoList$Training} matrix contains NTP predictions
  for 416 samples.
  The factor \code{phenoList$Testing} matrix contains NTP predictions
  for 208 samples.
  Factor levels for the CRIS classes are encoded as follow:
  \code{CRISA}, \code{CRISB}, \code{CRISC}, \code{CRISD}, and \code{CRISE}.
}


\seealso{
  \code{\link{matList}}
}

\author{
  Luigi Marchionni \email{marchion@jhu.edu}
}

\examples{

### Load gene expression datasets
data(phenoList)

### Show the class and names of the ``phenoList'' object
class(phenoList)
names(phenoList)

### Show the length of the ``phenoList'' factors
sapply(phenoList, length)

### Summary for the ``phenoList'' factors
sapply(phenoList, summary)

}

\keyword{datasets}
