\name{doVariableSelection}
\alias{doVariableSelection}
\title{a function to do variable selection on expression matrices.}
\usage{
doVariableSelection(exprMat, removeLowVaryingGenes)
}
\arguments{
  \item{a}{matrix of expression levels, rows contain genes
  and columns contain samples.}

  \item{removeLowVaryingGenes}{the proportion of low
  varying genes to be removed.}
}
\value{
a vector of row ids to keep
}
\description{
This funtino will I.e. remove genes with low variation It
returns a vector of row ids to keep. Note, rownames() must
be specified.
}

