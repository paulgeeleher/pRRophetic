#' a function to do variable selection on expression matrices.
#'
#' This funtino will I.e. remove genes with low variation
#' It returns a vector of row ids to keep. Note, rownames() must be specified.
#' 
#' @param a matrix of expression levels, rows contain genes and columns contain samples.
#' @param removeLowVaryingGenes the proportion of low varying genes to be removed.
#' 
#' @return a vector of row ids to keep
doVariableSelection <- function(exprMat, removeLowVaryingGenes)
{
  vars <- apply(exprMat, 1, var)
  return(order(vars, decreasing=TRUE)[seq(1:as.integer(nrow(exprMat)*(1-removeLowVaryingGenes)))])
}










