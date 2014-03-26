#' Take an expression matrix and if duplicate genes are measured, summarize them by their means
#' 
#' This function accepts two expression matrices, with gene ids as rownames() and 
#' sample ids as colnames(). It will find all duplicate genes and summarize their
#' expression by their mean.
#'
#' @param exprMat a gene expression matrix with gene names as row ids and sample names as column ids.
#'
#' @return a gene expression matrix that does not contain duplicate gene ids
#'
#' @keywords summarize duplicate gene ids by their mean.
#'
#' @export
summarizeGenesByMean <- function(exprMat)
{
  geneIds <- rownames(exprMat)
  t <- table(geneIds) # how many times is each gene name duplicated
  allNumDups <- unique(t)
  allNumDups <- allNumDups[-which(allNumDups == 1)]

  # create a *new* gene expression matrix with everything in the correct order....
  # start by just adding stuff that isn't duplicated
  exprMatUnique <- exprMat[which(geneIds %in% names(t[t == 1])), ]
  gnamesUnique <- geneIds[which(geneIds %in% names(t[t == 1]))]

  # add all the duplicated genes to the bottom of "exprMatUniqueHuman", summarizing as you go
  for(numDups in allNumDups) 
  {
    geneList <- names(which(t == numDups))
    
    for(i in 1:length(geneList))
    {
      exprMatUnique <- rbind(exprMatUnique, colMeans(exprMat[which(geneIds == geneList[i]), ]))
      gnamesUnique <- c(gnamesUnique, geneList[i])
      # print(i)
    }
  }
  
  if(class(exprMatUnique) == "numeric")
  {
    exprMatUnique <- matrix(exprMatUnique, ncol=1)
  }
  
  rownames(exprMatUnique) <- gnamesUnique
  return(exprMatUnique)
}
