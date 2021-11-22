\name{predictCRISclassKTSP}

\alias{predictCRISclassKTSP}

\title{Function to classify samples using a kTSP predictor for CRIS classification.}

\description{
  \code{predictCRISclassKTSP} classifies new test samples
  using kTSP predictor for CRIS classification.
}

\usage{
  predictCRISclassKTSP(mat)
}

\arguments{
  \item{mat}{the test data: a numeric matrix in which the rows
    represent the genes and the columns the samples.
    The rownames of this matrix must be valid gene symbols corresponding
    to the genes used in the kTSP classifier.}
}


\author{
  Luigi Marchionni \email{marchion@jhu.edu}
}


\examples{

##################################################
### Load gene expression data for the training set
data(matList)
data(phenoList)

### Data structure
sapply(matList, class)
sapply(matList, dim)

### Valid gene expression matrix with all CRIS genes
newMat <- matList$Training

### To make predictions on 1 matrix
newPreds <- predictCRISclassKTSP(newMat)

### Counts classifications
summary(newPreds)

###  Load predictions based on NPT
data(matList)

### NPT classification
refClass <- phenoList$Training

### Agreement

}

\keyword{kTSP, classification, prediction}
