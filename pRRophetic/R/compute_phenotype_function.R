## Functions to predict a phenotype from microarray expression data of different platforms.
##

#' Calculates phenotype from microarray data.
#'
#' This function uses ridge regression to calculate a phenotype from an gene expression,
#' given a gene expression matrix where the phenotype is already known. The function 
#' also integrates the two datasets using a user-defined procedure, power transforms
#' the known phenotype and provides several other options for flexible and powerful prediction
#' from a gene expression matrix.
#' 
#'
#' @param trainingExprData The training data. A matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "testExprData"
#' @param trainingPtype The known phenotype for "trainingExprData". A numeric vector which MUST be the same length as the number of columns of "trainingExprData".
#' @param testExprData The test data where the phenotype will be estimted. It is a matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "trainingExprData".
#' @param batchCorrect How should training and test data matrices be homogenized. Choices are "eb" (default) for ComBat, "qn" for quantiles normalization or "none" for no homogenization.
#' @param powerTransformPhenotype Should the phenotype be power transformed before we fit the regression model? Default to TRUE, set to FALSE if the phenotype is already known to be highly normal.
#' @param removeLowVaryingGenes What proportion of low varying genes should be removed? 20 percent be default
#' @param minNumSamples How many training and test samples are requried. Print an error if below this threshold
#' @param selection How should duplicate gene ids be handled. Default is -1 which asks the user. 1 to summarize by their or 2 to disguard all duplicates.
#' @param printOutput Set to FALSE to supress output
#'
#' @return A vector of the estimated phenotype, in the same order as the columns of "testExprData".
#'
#' @import sva
#' @import ridge
#' @import car
#'
#' @keywords predict, phenotype
#' @export
calcPhenotype <- function(trainingExprData, trainingPtype, testExprData, batchCorrect="eb", powerTransformPhenotype=TRUE, removeLowVaryingGenes=.2, minNumSamples=10, selection=-1, printOutput=TRUE)
{
  
  # check if the supplied data are of the correct classes
  if(class(testExprData)[1] != "matrix") stop("ERROR: \"testExprData\" must be a matrix.");
  if(class(trainingExprData)[1] != "matrix") stop("ERROR: \"trainingExprData\" must be a matrix.");
  if(class(trainingPtype)[1] != "numeric") stop("ERROR: \"trainingPtype\" must be a numeric vector.");
  if(ncol(trainingExprData) != length(trainingPtype)) stop("The training phenotype must be of the same length as the number of columns of the training expressin matrix.");
  
  # check if an adequate number of training and test samples have been supplied.
  if((ncol(trainingExprData) < minNumSamples) || (ncol(testExprData) < minNumSamples))
  {
    stop(paste("There are less than", minNumSamples, "samples in your test or training set. It is strongly recommended that you use larger numbers of samples in order to (a) correct for batch effects and (b) fit a reliable model. To supress this message, change the \"minNumSamples\" parameter to this function."))
  }

  # Get the homogenized data
  homData <- homogenizeData(testExprData, trainingExprData, batchCorrect=batchCorrect, selection=selection, printOutput=printOutput)
  
  # Do variable selection if specified. By default we remove 20% of least varying genes.
  # Otherwise, keep all genes.
  if(removeLowVaryingGenes > 0 && removeLowVaryingGenes < 1)
  {
    keepRows <- doVariableSelection(cbind(homData$test, homData$train), removeLowVaryingGenes=removeLowVaryingGenes)
  }
  else 
    keepRows <- seq(1:nrow(homData$train))
  
  
  # PowerTranform phenotype if specified.
  offset = 0
  if(powerTransformPhenotype)
  {
    if(min(trainingPtype) < 0) # all numbers must be postive for a powerTranform to work, so make them positive.
    {
      offset <- -min(trainingPtype) + 1
      trainingPtype <- trainingPtype + offset
    }
    
    transForm <- powerTransform(trainingPtype)[[6]]
    trainingPtype <- trainingPtype^transForm
  }
  
  # create the Ridge Regression model on our training data
  if(printOutput) cat("\nFitting Ridge Regression model... ");
  trainFrame <- data.frame(Resp=trainingPtype, t(homData$train[keepRows, ]))
  rrModel <- linearRidge(Resp ~ ., data=trainFrame)
  if(printOutput) cat("Done\n\nCalculating predicted phenotype...");
  
  # calculate the relative contribution of each gene to the prediction
  # i might report these, I don't know if there's any point.
  totBeta <- sum(abs(coef(rrModel)))
  eachBeta <- abs(coef(rrModel))
  eachContribution <- eachBeta/totBeta
  
  # predict the new phenotype for the test data.
  # if there is a single test sample, there may be problems in predicting using the predict() function for the linearRidge package
  # This "if" statement provides a workaround
  if(class(homData$test)[1] == "numeric")
  {
    n <- names(homData$test)
    homData$test <- matrix(homData$test, ncol=1)
    rownames(homData$test) <- n
    testFrame <- data.frame(t(homData$test[keepRows, ]))
    preds <- predict(rrModel, newdata=rbind(testFrame, testFrame))[1]
  }
  else #predict for more than one test sample
  {
    testFrame <- data.frame(t(homData$test[keepRows, ]))
    preds <- predict(rrModel, newdata=testFrame)
  }
  
  # if the response variable was transformed, untransform it now...
  if(powerTransformPhenotype)
  {
    preds <- preds^(1/transForm)
    preds <- preds - offset
  }
  if(printOutput) cat("Done\n\n");
  
  return(preds)
}



#' Cross validation on training dataset
#'
#' This function does cross validation on a training set to estimate prediction accuracy on a training set.
#' If the actual test set is provided, the two datasets can be subsetted and homogenized before the 
#' cross validation analysis is preformed. This may improve the estimate of prediction accuracy.
#'
#' @param trainingExprData The training data. A matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "testExprData"
#' @param trainingPtype The known phenotype for "trainingExprData". A numeric vector which MUST be the same length as the number of columns of "trainingExprData".
#' @param testExprData The test data where the phenotype will be estimted. It is a matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "trainingExprData".
#' @param cvFold Specify the "fold" requried for cross validation. "-1" will do leave one out cross validation (LOOCV)
#' @param powerTransformPhenotype Should the phenotype be power transformed before we fit the regression model? Default to TRUE, set to FALSE if the phenotype is already known to be highly normal.
#' @param batchCorrect How should training and test data matrices be homogenized. Choices are "eb" (default) for ComBat, "qn" for quantiles normalization or "none" for no homogenization.
#' @param removeLowVaryingGenes What proportion of low varying genes should be removed? 20 precent be default
#' @param minNumSamples How many training and test samples are requried. Print an error if below this threshold
#' @param selection How should duplicate gene ids be handled. Default is -1 which asks the user. 1 to summarize by their or 2 to disguard all duplicates.
#' @param printOutput Set to FALSE to supress output
#'
#' @return An object of class "pRRopheticCv", which is a list with two members, "cvPtype" and "realPtype", which correspond to the cross valiation predicted phenotype and the  user provided measured phenotype respectively.
#'
#' @import sva
#' @import ridge
#' @import car
#'
#' @keywords predict phenotype
#' @export
predictionAccuracyByCv <- function(trainingExprData, trainingPtype, testExprData=-1, cvFold=-1, powerTransformPhenotype=TRUE, batchCorrect="eb", removeLowVaryingGenes=.2, minNumSamples=10, selection=1)
{

  # check if an adequate number of training and test samples have been supplied.
  if((ncol(trainingExprData) < minNumSamples))
  {
    stop(paste("There are less than", minNumSamples, "samples in your test or training set. It is strongly recommended that you use larger numbers of samples in order to (a) correct for intrinsic difference in your training and test sets and (b) fit a reliable model. To supress this message, change the \"minNumSamples\" parameter to this function."))
  }

  # check if a test matrix was supplied, if not, don't homogenize the data (but store it in a list called homData anyway for convenience)
  if(is.null(testExprData))
  {
    homData <- list()
    homData$selection <- selection
    homData$train <- trainingExprData
  } else if(!is.null(testExprData)) {
    homData <- homogenizeData(testExprData, trainingExprData, batchCorrect=batchCorrect, selection=selection) # homogenize the data.
  }
  
  nTrain <- ncol(trainingExprData)
  predPtype <- numeric() # a numeric vector to hold the predicted phenotypes for the CV subgroups

  # Perform either N fold cross validation or LOOCV, depending on the "cvFold" variable.
  if(cvFold == -1) # if we are doing leave-one-out cross-validation (LOOCV).
  {
    for(i in 1:nTrain)
    {
      testMatTemp <- matrix(homData$train[,i], ncol=1)
      rownames(testMatTemp) <- rownames(homData$train)
      #predPtype[i] <- calcPhenotype(testMatTemp, trainCvSet[,-i], trainingPtype[-i], batchCorrect="none", minNumSamples=0, selection=homData$selection, removeLowVaryingGenes=removeLowVaryingGenes, powerTransformPhenotype=powerTransformPhenotype)
      predPtype[i] <- calcPhenotype(homData$train[,-i], trainingPtype[-i], testMatTemp, batchCorrect="none", minNumSamples=0, selection=homData$selection, removeLowVaryingGenes=removeLowVaryingGenes, powerTransformPhenotype=powerTransformPhenotype, printOutput=FALSE)
      
      # print an update for each 20% of the this, this is slow, so we should give some updates....
      if(i %% as.integer(nTrain/5) == 0)
      cat(paste(i, "of" , nTrain, "iterations complete. \n"))
    }
  }
  else if(cvFold > 1) # if we are doing N-fold cross validation
  {
    randTestSamplesIndex <- sample(1:nTrain) # create a random order for samples
  
    # create a vector which indicates which samples are in which group... and split into list of groups
    sampleGroup <- rep(cvFold, nTrain)
    groupSize <- as.integer(nTrain / cvFold)
    for(j in 1:(cvFold-1)) { sampleGroup[(((j-1)*groupSize)+1):(j*groupSize)] <- rep(j, groupSize) }
    cvGroupIndexList <- split(randTestSamplesIndex, sampleGroup)
    
    # predict on each of the groups....
    for(j in 1:cvFold)
    {
      
      # create the ranomdly chosen "training" and "test" sets for cross validation
      testCvSet <- homData$train[, cvGroupIndexList[[j]]]
      trainCvSet <- homData$train[, unlist(cvGroupIndexList[-j])]
      trainPtypeCv <- trainingPtype[unlist(cvGroupIndexList[-j])]
      
      predPtype <- c(predPtype, calcPhenotype(trainCvSet, trainPtypeCv, testCvSet, batchCorrect="none", minNumSamples=0, selection=homData$selection, removeLowVaryingGenes=removeLowVaryingGenes, powerTransformPhenotype=powerTransformPhenotype, printOutput=FALSE))
      
      cat(paste("\n", j, " of ", cvFold, " iterations complete.", sep=""))
    }
    
    # re-order the predicted phenotypes correctly, as they were ranomized when this started.
    predPtype <- predPtype[order(randTestSamplesIndex)]
    
  } else {
    stop("Unrecognised value of \"cvFold\"")
  }
  
  finalData <- list(cvPtype=predPtype, realPtype=trainingPtype)
  class(finalData) <- "pRRopheticCv"
  
  return(finalData)
}


#' R^2 from "pRRopheticCv" object.
#'
#' Given an object of class "pRRopheticCv", i.e. the output of cross validation, calculate 
#' the R^2 value for the prediction (an estimate of prediction accuracy).
#'
#' @param cvOutData an object of class "pRRopheticCv", i.e. the outpout of the "predictionAccuracyByCv()" function
#'
#' @return a numeric vector containing the R squared value from the cross validation.
#'
#' @import car
#'
#' @keywords r-squared
#' @export
estimateRsqr.pRRopheticCv <- function(cvOutData, powerTranform=TRUE)
{
  # calculate the R^2
  return(summary(lm(cvOutData$realPtype~cvOutData$cvPtype))$r.squared)
}


#' Confidence intervals from "pRRopheticCv" object.
#'
#' Given an object of class "pRRopheticCv", i.e. the output of cross validation, calculate 
#' an average confidence interval for the predictions.
#'
#' @param cvOutData an object of class "pRRopheticCv", i.e. the outpout of the "predictionAccuracyByCv()" function
#' @param conf the confidence interval required, by default 95 precent confidence interval.
#'
#' @return a numeric vector containing the average upper and lower confidence interval.
#'
#'
#' @keywords confidence interval
#' @export
estimateCI.pRRopheticCv <- function(cvOutData, conf=.95)
{

  # report 95% confidence intervals and 50% confidence intervals (estimated from the untransformed data)
  allDeviationsRaw <- cvOutData$cvPtype - cvOutData$realPtype
  allDeviations <- abs(allDeviationsRaw)
  inCi <- which(allDeviations < quantile(allDeviations, conf))
  ci <- c(min(allDeviationsRaw[inCi]), max(allDeviationsRaw[inCi]))
  return(ci)
}


#' Mean prediction error from "pRRopheticCv" object.
#'
#' Given an object of class "pRRopheticCv", estiamte the mean prediction error,
#' i.e. the mean difference between the predicted and measured phenotype.
#'
#' @param cvOutData an object of class "pRRopheticCv", i.e. the outpout of the "predictionAccuracyByCv()" function
#'
#' @return a numeric vector containing the mean prediction error from the cross validation.
#'
#' @keywords prediction error
#' @export
estimateMeanPredictionError.pRRopheticCv <- function(cvOutData)
{
  allDeviationsRaw <- abs(cvOutData$cvPtype - cvOutData$realPtype)
  return(mean(allDeviationsRaw))
}

#' Median prediction error from "pRRopheticCv" object.
#'
#' Given an object of class "pRRopheticCv", estiamte the median prediction error,
#' i.e. the median difference between the predicted and measured phenotype.
#'
#' @param cvOutData an object of class "pRRopheticCv", i.e. the outpout of the "predictionAccuracyByCv()" function
#'
#' @return a numeric vector containing the median prediction error from the cross validation.
#'
#' @keywords prediction error
#' @export
estimateMedianPredictionError.pRRopheticCv <- function(cvOutData)
{
  allDeviationsRaw <- abs(cvOutData$cvPtype - cvOutData$realPtype)
  return(mean(median))
}


#' Summary of "pRRopheticCv" object.
#'
#' Given an object of class "pRRopheticCv", print various metrics that 
#' summarize the performance of the cross validataion analysis
#'
#' @param cvOutData an object of class "pRRopheticCv", i.e. the outpout of the "predictionAccuracyByCv()" function
#'
#' @keywords summary
#' @export
summary.pRRopheticCv <- function(cvOutData)
{
  cat("\nSummary of cross-validation results:\n\n")
  
  corOut <- cor.test(cvOutData[[1]], cvOutData[[2]])
  cat(paste("Pearsons correlation:", round(corOut$estimate, digits=2), ", P = ", corOut$p.value, "\n"))
  cat(paste("R-squared value: ", round(estimateRsqr.pRRopheticCv(cvOutData), digits=2), "\n", sep=""))
  cis <- estimateCI.pRRopheticCv(cvOutData)
  cat(paste("Estimated 95% confidence intervals: ", round(cis[1], digits=2), ", ", round(cis[2], digits=2), "\n", sep=""))
  cat(paste("Mean prediction error: ", round(estimateMeanPredictionError.pRRopheticCv(cvOutData), digits=2), "\n\n", sep=""))
}

#' Plot "pRRopheticCv" object.
#'
#' Given an object of class "pRRopheticCv", plot the cross validation
#' predicted values against the measured values. Also plots a regression
#' line.
#'
#' @param cvOutData an object of class "pRRopheticCv", i.e. the outpout of the "predictionAccuracyByCv()" function
#'
#' @keywords plot
#' @export
plot.pRRopheticCv <- function(cvOutData)
{
  coefs <- coef(lm(cvOutData$realPtype~cvOutData$cvPtype))
  plot(cvOutData$cvPtype, cvOutData$realPtype, main="Measured phenotype Vs. C.V. predicted phenotype", xlab="Predicted Phenotype", ylab="Measured Phenotype")
  abline(coefs[1], coefs[2], col="red")
}

