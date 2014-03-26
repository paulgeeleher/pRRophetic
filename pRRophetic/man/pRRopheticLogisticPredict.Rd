\name{pRRopheticLogisticPredict}
\alias{pRRopheticLogisticPredict}
\title{Take an expression matrix and if duplicate genes are measured, summarize them by their means}
\usage{
pRRopheticLogisticPredict(testMatrix, drug, tissueType = "all",
  doCV = FALSE, batchCorrect = "eb", minNumSamples = 10, selection = -1,
  printOutput = TRUE, numGenesSelected = 1000, numSens = 15,
  numRes = 55)
}
\arguments{
  \item{testMatrix}{a gene expression matrix with gene
  names as row ids and sample names as column ids.}
}
\value{
a gene expression matrix that does not contain duplicate
gene ids
}
\description{
This function accepts two expression matrices, with gene
ids as rownames() and sample ids as colnames(). It will
find all duplicate genes and summarize their expression by
their mean.
}
\keyword{by}
\keyword{duplicate}
\keyword{gene}
\keyword{ids}
\keyword{mean.}
\keyword{summarize}
\keyword{their}

