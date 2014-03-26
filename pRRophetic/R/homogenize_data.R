#' Take two expression matrices and return homogenized versions of the matrices.
#' 
#' This function accepts two expression matrices, with gene ids as rownames() and 
#' sample ids as colnames(). It will deal with duplicate gene ids.
#' subset and order the matrices correctly.
#' and perform homogenize the data using whatever method is specified (by default Combat from the sva library).
#'
#' @param testExprData Gene expression matrix for samples on which we wish to predict a phenotype. Gene names as rows, samples names as columns.
#' @param trainingExprData Gene expression matrix for samples for which we the phenotype is already known.
#' @param batchCorrect The type of batch correction to be used. Options are "eb" for Combat, "none", or "qn" for quantile normalization.
#' @param selection parameter can be used to specify how duplicates are handled, by default value -1 means ask the user. 1 means summarize duplictes by their mean and 2 means to disguard all duplicate genes.
#' @param Set to FALSE to supress output
#'
#' @return a list containing two entries $train and $test, which are the homogenized input matrices.
#'
#' @import sva
#' @import preprocessCore
#'
#' @keywords homogenize expression data
#'
#' @export
homogenizeData <- function(testExprMat, trainExprMat, batchCorrect="eb", selection=-1, printOutput=TRUE)
{
  # Check batchCorrect paramter
  if(!(batchCorrect %in% c("eb", "qn", "none")))
    stop("\"batchCorrect\" must be one of \"eb\", \"qn\" or \"none\"")

  # check if both row and column names have been specified
  if(is.null(rownames(trainExprMat)) || is.null(rownames(testExprMat)))
  {
    stop("ERROR: Gene identifiers must be specified as \"rownames()\" on both training and test expression matrices. Both matices must have the same type of gene identifiers.")
  }
  
  # check that some of the row names overlap between both datasets (print an error if none overlap.
  if(sum(rownames(trainExprMat) %in% rownames(testExprMat)) == 0)
  {
    stop("ERROR: The rownames() of the supplied expression matrices do not match. Note that these are case-sensitive.")
  } else {
    if(printOutput) cat(paste("\n", sum(rownames(trainExprMat) %in% rownames(testExprMat)), " gene identifiers overlap between the supplied expression matrices... \n", paste=""));
  }

  # if there are duplicate gene names, give the option of removing them or summarizing them by their mean.
  if((sum(duplicated(rownames(trainExprMat))) > 0) || sum(sum(duplicated(rownames(testExprMat))) > 0))
  {
    if(selection == -1) #print the following if we're asking the user how to proceed.
    {  
      cat("\nExpression matrix contain duplicated gene identifiers (i.e. duplicate rownames()), how would you like to proceed:")
      cat("\n1. Summarize duplicated gene ids by their mean value (acceptable in most cases)")
      cat("\n2. Disguard all duplicated genes (recommended if unsure)")
      cat("\n3. Abort (if you want to deal with duplicate genes ids manually)\n")
    }
    
    while(is.na(selection) | selection <= 0 | selection > 3 )
    {
      selection <- readline("Selection: ")
      selection <- ifelse(grepl("[^1-3.]", selection), -1 , as.numeric(selection))
    }
    
    cat("\n")
    
    if(selection == 1) # summarize duplicates by their mean
    {
      if((sum(duplicated(rownames(trainExprMat))) > 0))
      {
	trainExprMat <- summarizeGenesByMean(trainExprMat)
      }
      if((sum(duplicated(rownames(testExprMat))) > 0))
      {
	testExprMat <- summarizeGenesByMean(testExprMat)
      }
    }
    else if(selection == 2) # disguard all duplicated genes
    {
      if((sum(duplicated(rownames(trainExprMat))) > 0))
      {
	keepGenes <- names(which(table(rownames(trainExprMat)) == 1))
	trainExprMat <- trainExprMat[keepGenes, ]
      }

      if((sum(duplicated(rownames(testExprMat))) > 0))
      {
	keepGenes <- names(which(table(rownames(testExprMat)) == 1))
	testExprMat <- testExprMat[keepGenes, ]
      }      
    } else {
      stop("Exectution Aborted!")
    }
    
  }
  
  # subset and order gene ids on the expression matrices
  commonGenesIds <- rownames(trainExprMat)[rownames(trainExprMat) %in% rownames(testExprMat)]
  trainExprMat <- trainExprMat[commonGenesIds, ]
  testExprMat <- testExprMat[commonGenesIds, ]
  
  # subset and order the two expresison matrices
  if(batchCorrect == "eb")
  {
    # subset to common genes andbatch correct using ComBat
    dataMat <- cbind(trainExprMat, testExprMat)
    mod <- data.frame("(Intercept)"=rep(1, ncol(dataMat)))
    rownames(mod) <- colnames(dataMat)
    whichbatch <- as.factor(c(rep("train", ncol(trainExprMat)), rep("test", ncol(testExprMat))))
    combatout <- ComBat(dataMat, whichbatch, mod=mod)
    return(list(train=combatout[, whichbatch=="train"], test=combatout[, whichbatch=="test"], selection=selection))
  }
  else if(batchCorrect == "qn")
  {
    # library("preprocessCore")
    dataMat <- cbind(trainExprMat, testExprMat)
    dataMatNorm <- normalize.quantiles(dataMat)
    whichbatch <- as.factor(c(rep("train", ncol(trainExprMat)), rep("test", ncol(testExprMat))))
    return(list(train=dataMatNorm[, whichbatch=="train"], test=dataMatNorm[, whichbatch=="test"], selection=selection))
  } else {
    return(list(train=trainExprMat, test=testExprMat, selection=selection))
  }
}