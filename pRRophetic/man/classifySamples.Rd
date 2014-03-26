\name{classifySamples}
\alias{classifySamples}
\title{Dichotimize a training expression set and fit a logistic ridge regression model which is applied to the test expression matirx.}
\usage{
classifySamples(trainingExprData, trainingPtype, testExprData,
  batchCorrect = "eb", minNumSamples = 10, selection = -1,
  printOutput = TRUE, numGenesSelected = 1000, numSens = 15,
  numRes = 55)
}
\arguments{
  \item{trainingExprData}{- Gene expression matrix for
  samples for which we the phenotype is already known.}

  \item{trainingPtype}{The known phenotype, a vector in the
  same order as the columns of "trainingExprData" or with
  the same names as colnames of "trainingExprData".}

  \item{testExprData}{Gene expression matrix for samples on
  which we wish to predict a phenotype. Gene names as rows,
  samples names as columns.}

  \item{batchCorrect}{The type of batch correction to be
  used. Options are "eb", "none", .....}

  \item{removeLowVaryingGenes}{The proportion of genes with
  lowest variation to be removed from the model, by default
  20 precent. Set to 0 to use all genes.}

  \item{variableSelectionMethod}{To improve performance,
  fit the logistic regression dataset with only this number
  of genes. By default 2000}

  \item{numGenesSelected}{Specifies how genes are selected
  for "variableSelectionMethod". Options are "tTests",
  "pearson" and "spearman".}

  \item{minNumTestSamples}{The minimum number of test
  samples, print an error if the number of columns of
  "testExprData" is below this threshold. A large number of
  test samples may be necessary to correct for batch
  effects.}
}
\description{
Dichotimize a training expression set and fit a logistic
ridge regression model which is applied to the test
expression matirx. This function will return a set of
probabilities.
}

