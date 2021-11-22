\name{cris_classifier}
\alias{cris_classifier}
\title{cris_classifier}
\description{
the function perfoms cris classifier on CRC data.
Data should be provided as .gct or .txt file, with first column corresponding to gene symbol
}
\usage{
cris_classifier(input.exp.filename, output.name, nresmpl, rnd.seed)
}
\arguments{
  \item{input.exp.filename}{gene expression data table in linear format. Data should be provided .gct or .txt file, with first column corresponding to gene symbol}
  \item{output.name}{ define the file path and name of the output}
  \item{nresmpl}{ number of permutations for montecarlo simulation }
  \item{rnd.seed}{ random sed for simulation }
}
\details{
The function performs cris classifier on CRC data. Data should be provided as linear values in .gct or .txt file format, with first column corresponding to gene symbol
}
\examples{

demo <- list.files(pattern="txt.gz$", system.file("data",package="CRISclassifier"), full.names=TRUE)

cris_classifier(input.exp.filename = demo, output.name="cris", nresmpl=1)
}

\value{
\item{input.exp.filename}{the expression set to classify}
\item{output.name}{the the output file name}
\item{nresmpl}{number of permutations for the montecarlo simulation}
\item{rnd.seed}{randome seed}
}
\author{
Claudio Isella
}
\keyword{CRISclassifier}

