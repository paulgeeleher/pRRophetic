## This file contains functions for prediction and classification from the CGP cell line data....

#' Given a gene expression matrix, predict drug senstivity for a drug in CGP
#' 
#' Given a gene expression matrix, predict drug senstivity for a drug in CGP.
#'
#' @param testMatrix a gene expression matrix with gene names as row ids and sample names as column ids.
#' @param drug the name of the drug for which you would like to predict sensitivity, one of A.443654, A.770041, ABT.263, ABT.888, AG.014699, AICAR, AKT.inhibitor.VIII, AMG.706, AP.24534, AS601245, ATRA, AUY922, Axitinib, AZ628, AZD.0530, AZD.2281, AZD6244, AZD6482, AZD7762, AZD8055, BAY.61.3606, Bexarotene, BI.2536, BIBW2992, Bicalutamide, BI.D1870, BIRB.0796, Bleomycin, BMS.509744, BMS.536924, BMS.708163, BMS.754807, Bortezomib, Bosutinib, Bryostatin.1, BX.795, Camptothecin, CCT007093, CCT018159, CEP.701, CGP.082996, CGP.60474, CHIR.99021, CI.1040, Cisplatin, CMK, Cyclopamine, Cytarabine, Dasatinib, DMOG, Docetaxel, Doxorubicin, EHT.1864, Elesclomol, Embelin, Epothilone.B, Erlotinib, Etoposide, FH535, FTI.277, GDC.0449, GDC0941, Gefitinib, Gemcitabine, GNF.2, GSK269962A, GSK.650394, GW.441756, GW843682X, Imatinib, IPA.3, JNJ.26854165, JNK.9L, JNK.Inhibitor.VIII, JW.7.52.1, KIN001.135, KU.55933, Lapatinib, Lenalidomide, LFM.A13, Metformin, Methotrexate, MG.132, Midostaurin, Mitomycin.C, MK.2206, MS.275, Nilotinib, NSC.87877, NU.7441, Nutlin.3a, NVP.BEZ235, NVP.TAE684, Obatoclax.Mesylate, OSI.906, PAC.1, Paclitaxel, Parthenolide, Pazopanib, PD.0325901, PD.0332991, PD.173074, PF.02341066, PF.4708671, PF.562271, PHA.665752, PLX4720, Pyrimethamine, QS11, Rapamycin, RDEA119, RO.3306, Roscovitine, Salubrinal, SB.216763, SB590885, Shikonin, SL.0101.1, Sorafenib, S.Trityl.L.cysteine, Sunitinib, Temsirolimus, Thapsigargin, Tipifarnib, TW.37, Vinblastine, Vinorelbine, Vorinostat, VX.680, VX.702, WH.4.023, WO2009093972, WZ.1.84, X17.AAG, X681640, XMD8.85, Z.LLNle.CHO, ZM.447439.
#' @param tissueType specify if you would like to traing the models on only a subset of the CGP cell lines (based on the tissue type from which the cell lines originated). This be one any of "all" (for everything, default option), "allSolidTumors" (everything except for blood), "blood", "breast", "CNS", "GI tract" ,"lung", "skin", "upper aerodigestive"
#' @param doCV
#'
#' @return a gene expression matrix that does not contain duplicate gene ids
#'
#' @keywords summarize duplicate gene ids by their mean.
#'
#' @export
pRRopheticPredict <- function(testMatrix, drug, tissueType="all", doCV=FALSE, batchCorrect="eb", powerTransformPhenotype=TRUE, removeLowVaryingGenes=.2, minNumSamples=10, selection=-1, printOutput=TRUE)
{
  cgpTrainData <- getCGPinfo(drug, tissueType) # get the IC50 and expression data for this drug/tissueType
  
  predictedPtype <- calcPhenotype(cgpTrainData$trainDataOrd, cgpTrainData$ic50sOrd, testMatrix, batchCorrect=batchCorrect, powerTransformPhenotype=powerTransformPhenotype, removeLowVaryingGenes=removeLowVaryingGenes, minNumSamples=minNumSamples, selection=selection, printOutput=printOutput)

  return(predictedPtype)
  
}


#' This function uses X fold cross validation on the TrainingSet to estimate the accuracy of the 
#' phenotype prediction fold: How many fold cross-validation to use.
#'
#' This function does cross validation on a training set to estimate prediction accuracy on a training set.
#' If the actual test set is provided, the two datasets can be subsetted and homogenized before the 
#' cross validation analysis is preformed. This may improve the estimate of prediction accuracy.
#'
#' @param testExprData The test data where the phenotype will be estimted. It is a matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "trainingExprData".
#' @param drug the name of the drug for which you would like to predict sensitivity, one of A.443654, A.770041, ABT.263, ABT.888, AG.014699, AICAR, AKT.inhibitor.VIII, AMG.706, AP.24534, AS601245, ATRA, AUY922, Axitinib, AZ628, AZD.0530, AZD.2281, AZD6244, AZD6482, AZD7762, AZD8055, BAY.61.3606, Bexarotene, BI.2536, BIBW2992, Bicalutamide, BI.D1870, BIRB.0796, Bleomycin, BMS.509744, BMS.536924, BMS.708163, BMS.754807, Bortezomib, Bosutinib, Bryostatin.1, BX.795, Camptothecin, CCT007093, CCT018159, CEP.701, CGP.082996, CGP.60474, CHIR.99021, CI.1040, Cisplatin, CMK, Cyclopamine, Cytarabine, Dasatinib, DMOG, Docetaxel, Doxorubicin, EHT.1864, Elesclomol, Embelin, Epothilone.B, Erlotinib, Etoposide, FH535, FTI.277, GDC.0449, GDC0941, Gefitinib, Gemcitabine, GNF.2, GSK269962A, GSK.650394, GW.441756, GW843682X, Imatinib, IPA.3, JNJ.26854165, JNK.9L, JNK.Inhibitor.VIII, JW.7.52.1, KIN001.135, KU.55933, Lapatinib, Lenalidomide, LFM.A13, Metformin, Methotrexate, MG.132, Midostaurin, Mitomycin.C, MK.2206, MS.275, Nilotinib, NSC.87877, NU.7441, Nutlin.3a, NVP.BEZ235, NVP.TAE684, Obatoclax.Mesylate, OSI.906, PAC.1, Paclitaxel, Parthenolide, Pazopanib, PD.0325901, PD.0332991, PD.173074, PF.02341066, PF.4708671, PF.562271, PHA.665752, PLX4720, Pyrimethamine, QS11, Rapamycin, RDEA119, RO.3306, Roscovitine, Salubrinal, SB.216763, SB590885, Shikonin, SL.0101.1, Sorafenib, S.Trityl.L.cysteine, Sunitinib, Temsirolimus, Thapsigargin, Tipifarnib, TW.37, Vinblastine, Vinorelbine, Vorinostat, VX.680, VX.702, WH.4.023, WO2009093972, WZ.1.84, X17.AAG, X681640, XMD8.85, Z.LLNle.CHO, ZM.447439.
#' @param cvFold Specify the "fold" requried for cross validation. "-1" will do leave one out cross validation (LOOCV)
#' @param powerTransformPhenotype Should the phenotype be power transformed before we fit the regression model? Default to TRUE, set to FALSE if the phenotype is already known to be highly normal.
#' @param batchCorrect How should training and test data matrices be homogenized. Choices are "eb" (default) for ComBat, "qn" for quantiles normalization or "none" for no homogenization.
#' @param removeLowVaryingGenes What proportion of low varying genes should be removed? 20 percent by default.
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
#' @keywords predict, phenotype
#' @export
pRRopheticCV <- function(drug, tissueType="all", testExprData=NULL, cvFold=-1, powerTransformPhenotype=TRUE, batchCorrect="eb", removeLowVaryingGenes=.2, minNumSamples=10, selection=1)
{
  cgpTrainData <- getCGPinfo(drug, tissueType) # get the IC50 and expression data for this drug/tissueType

  # I may need to alter this function so it can either take the test data or not take the test data....
  cvOut <- predictionAccuracyByCv(cgpTrainData$trainDataOrd, cgpTrainData$ic50sOrd, testExprData=testExprData, cvFold=cvFold, powerTransformPhenotype=powerTransformPhenotype, batchCorrect=batchCorrect, removeLowVaryingGenes=removeLowVaryingGenes, minNumSamples=minNumSamples, selection=selection)
  return(cvOut)
}


#' Given a drug and tissue type, return CGP expression and drug sensitivity data.
#'
#' Given a drug and tissue type, return CGP expression and drug sensitivity data.
#' 
#' @param drug the name of the drug for which you would like to predict sensitivity, one of A.443654, A.770041, ABT.263, ABT.888, AG.014699, AICAR, AKT.inhibitor.VIII, AMG.706, AP.24534, AS601245, ATRA, AUY922, Axitinib, AZ628, AZD.0530, AZD.2281, AZD6244, AZD6482, AZD7762, AZD8055, BAY.61.3606, Bexarotene, BI.2536, BIBW2992, Bicalutamide, BI.D1870, BIRB.0796, Bleomycin, BMS.509744, BMS.536924, BMS.708163, BMS.754807, Bortezomib, Bosutinib, Bryostatin.1, BX.795, Camptothecin, CCT007093, CCT018159, CEP.701, CGP.082996, CGP.60474, CHIR.99021, CI.1040, Cisplatin, CMK, Cyclopamine, Cytarabine, Dasatinib, DMOG, Docetaxel, Doxorubicin, EHT.1864, Elesclomol, Embelin, Epothilone.B, Erlotinib, Etoposide, FH535, FTI.277, GDC.0449, GDC0941, Gefitinib, Gemcitabine, GNF.2, GSK269962A, GSK.650394, GW.441756, GW843682X, Imatinib, IPA.3, JNJ.26854165, JNK.9L, JNK.Inhibitor.VIII, JW.7.52.1, KIN001.135, KU.55933, Lapatinib, Lenalidomide, LFM.A13, Metformin, Methotrexate, MG.132, Midostaurin, Mitomycin.C, MK.2206, MS.275, Nilotinib, NSC.87877, NU.7441, Nutlin.3a, NVP.BEZ235, NVP.TAE684, Obatoclax.Mesylate, OSI.906, PAC.1, Paclitaxel, Parthenolide, Pazopanib, PD.0325901, PD.0332991, PD.173074, PF.02341066, PF.4708671, PF.562271, PHA.665752, PLX4720, Pyrimethamine, QS11, Rapamycin, RDEA119, RO.3306, Roscovitine, Salubrinal, SB.216763, SB590885, Shikonin, SL.0101.1, Sorafenib, S.Trityl.L.cysteine, Sunitinib, Temsirolimus, Thapsigargin, Tipifarnib, TW.37, Vinblastine, Vinorelbine, Vorinostat, VX.680, VX.702, WH.4.023, WO2009093972, WZ.1.84, X17.AAG, X681640, XMD8.85, Z.LLNle.CHO, ZM.447439.
#' @param tissueType specify if you would like to traing the models on only a subset of the CGP cell lines (based on the tissue type from which the cell lines originated). This be one any of "all" (for everything, default option), "allSolidTumors" (everything except for blood), "blood", "breast", "CNS", "GI tract" ,"lung", "skin", "upper aerodigestive"#' 
#' 
#' @return a list with two entries, trainDataOrd the ordered expression data and ic50sOrd the drug sensitivity data. 
getCGPinfo <-  function(drug, tissueType="all")
{

    # was a valid tissue type specified; tissue types represeted by > 40 cell lines
  if(!tissueType %in% c("all", "allSolidTumors", "blood", "breast", "CNS", "GI tract" ,"lung", "skin", "upper aerodigestive")) stop("ERROR: the tissue type specified must be one of \"all\", \"allSolidTumors\", \"blood\", \"breast\", \"CNS\", \"GI tract\", \"lung\", \"skin\" or \"upper aerodigestive\"");
  
  # was a valid drug specified
  possibleDrugs <- c("A.443654", "A.770041", "ABT.263", "ABT.888", "AG.014699", "AICAR", "AKT.inhibitor.VIII", "AMG.706", "AP.24534", "AS601245", "ATRA", "AUY922", "Axitinib", "AZ628", "AZD.0530", "AZD.2281", "AZD6244", "AZD6482", "AZD7762", "AZD8055", "BAY.61.3606", "Bexarotene", "BI.2536", "BIBW2992", "Bicalutamide", "BI.D1870", "BIRB.0796", "Bleomycin", "BMS.509744", "BMS.536924", "BMS.708163", "BMS.754807", "Bortezomib", "Bosutinib", "Bryostatin.1", "BX.795", "Camptothecin", "CCT007093", "CCT018159", "CEP.701", "CGP.082996", "CGP.60474", "CHIR.99021", "CI.1040", "Cisplatin", "CMK", "Cyclopamine", "Cytarabine", "Dasatinib", "DMOG", "Docetaxel", "Doxorubicin", "EHT.1864", "Elesclomol", "Embelin", "Epothilone.B", "Erlotinib", "Etoposide", "FH535", "FTI.277", "GDC.0449", "GDC0941", "Gefitinib", "Gemcitabine", "GNF.2", "GSK269962A", "GSK.650394", "GW.441756", "GW843682X", "Imatinib", "IPA.3", "JNJ.26854165", "JNK.9L", "JNK.Inhibitor.VIII", "JW.7.52.1", "KIN001.135", "KU.55933", "Lapatinib", "Lenalidomide", "LFM.A13", "Metformin", "Methotrexate", "MG.132", "Midostaurin", "Mitomycin.C", "MK.2206", "MS.275", "Nilotinib", "NSC.87877", "NU.7441", "Nutlin.3a", "NVP.BEZ235", "NVP.TAE684", "Obatoclax.Mesylate", "OSI.906", "PAC.1", "Paclitaxel", "Parthenolide", "Pazopanib", "PD.0325901", "PD.0332991", "PD.173074", "PF.02341066", "PF.4708671", "PF.562271", "PHA.665752", "PLX4720", "Pyrimethamine", "QS11", "Rapamycin", "RDEA119", "RO.3306", "Roscovitine", "Salubrinal", "SB.216763", "SB590885", "Shikonin", "SL.0101.1", "Sorafenib", "S.Trityl.L.cysteine", "Sunitinib", "Temsirolimus", "Thapsigargin", "Tipifarnib", "TW.37", "Vinblastine", "Vinorelbine", "Vorinostat", "VX.680", "VX.702", "WH.4.023", "WO2009093972", "WZ.1.84", "X17.AAG", "X681640", "XMD8.85", "Z.LLNle.CHO", "ZM.447439")
  if(!drug %in% possibleDrugs) stop(paste("ERROR: the drug specified was not found. Note drug names are case sensitive. Please select from: ", (paste(possibleDrugs, collapse=", "))));
  
  ######################################################################
  ########## i need to create a valid data file velow and show
  ######################################################################
  
  # load("/home/pgeeleher/postdoc_stuff/HDAC_project/Scripts/r_package_files/data/drugAndPhenoCgp.RData") # contains drugToCellLineDataCgp (to map cell lines to .CEL file names), gdsc_brainarray_syms (the gene expression data), drugSensitivityDataCgp (the drug sensitivity data)
  data(drugAndPhenoCgp) # contains drugToCellLineDataCgp (to map cell lines to .CEL file names), gdsc_brainarray_syms (the gene expression data), drugSensitivityDataCgp (the drug sensitivity data)
  
  colIc50Name <- paste(drug, "_IC_50", sep="")
  ic50s <- as.numeric(drugSensitivityDataCgp[, colIc50Name])
  names(ic50s) <- drugSensitivityDataCgp[ ,"Cell.Line"]
  whichNas <- which(is.na(ic50s))
  ic50s <- ic50s[-whichNas]
  tissue <- drugSensitivityDataCgp[ ,"Cancer.Type"]
  names(tissue) <- drugSensitivityDataCgp[ ,"Cell.Line"]
  tissue <- tissue[-whichNas]

  # if a tissue type has been specified, use only tissues of that type.
  if(tissueType != "all")
  {
    if(tissueType == "allSolidTumors")
    {
      tissueType <- ic50s <- ic50s[!(tissue %in% "blood")]
    }
    else
    {
      ic50s <- ic50s[tissue %in% tissueType]
    }
  }

  # map the drug sensitivity and expression data
  pDataUnique <- drugToCellLineDataCgp[drugToCellLineDataCgp$Source.Name %in% 
  names(which(table(drugToCellLineDataCgp$Source.Name) == 1)), ]
  rownames(pDataUnique) <- pDataUnique$Source.Name
  commonCellLines <- rownames(pDataUnique)[rownames(pDataUnique) %in% names(ic50s)]
  pDataUniqueOrd <- pDataUnique[commonCellLines, ]
  ic50sOrd <- ic50s[commonCellLines]
  trainDataOrd <- gdsc_brainarray_syms[, pDataUniqueOrd$"Array.Data.File"]
  
  return(list(ic50sOrd=ic50sOrd, trainDataOrd=trainDataOrd))

}

#' Take an expression matrix and if duplicate genes are measured, summarize them by their means
#' 
#' This function accepts two expression matrices, with gene ids as rownames() and 
#' sample ids as colnames(). It will find all duplicate genes and summarize their
#' expression by their mean.
#'
#' @param testMatrix a gene expression matrix with gene names as row ids and sample names as column ids.
#'
#' @return a gene expression matrix that does not contain duplicate gene ids
#'
#' @keywords summarize duplicate gene ids by their mean.
#'
#' @export
pRRopheticLogisticPredict <- function(testMatrix, drug, tissueType="all", doCV=FALSE, batchCorrect="eb", minNumSamples=10, selection=-1, printOutput=TRUE, numGenesSelected=1000, numSens=15, numRes=55)
{
  cgpTrainData <- getCGPinfo(drug, tissueType) # get the IC50 and expression data for this drug/tissueType

  predictedPtype <- classifySamples(cgpTrainData$trainDataOrd, cgpTrainData$ic50sOrd, testMatrix, batchCorrect=batchCorrect,minNumSamples=minNumSamples, selection=selection, printOutput=printOutput, numGenesSelected=numGenesSelected, numSens=numSens, numRes=numRes)

  return(predictedPtype[,1])
}


#' Check the distribution of the drug response (IC50) data using a QQ-plot.
#' 
#' Visualize the distribution of the transformed IC50 data for a drug of interest using a QQ plot. If the distribution of the IC50 values deviates wildly from a normal distribtion, it is likely not suitalbe for prediction using a linear model (like linear ridge regression). This drug may be more suitable to constructing a model using a logistic or other type of model.
#
#' @param drug the name of the drug for which you would like to predict sensitivity, one of A.443654, A.770041, ABT.263, ABT.888, AG.014699, AICAR, AKT.inhibitor.VIII, AMG.706, AP.24534, AS601245, ATRA, AUY922, Axitinib, AZ628, AZD.0530, AZD.2281, AZD6244, AZD6482, AZD7762, AZD8055, BAY.61.3606, Bexarotene, BI.2536, BIBW2992, Bicalutamide, BI.D1870, BIRB.0796, Bleomycin, BMS.509744, BMS.536924, BMS.708163, BMS.754807, Bortezomib, Bosutinib, Bryostatin.1, BX.795, Camptothecin, CCT007093, CCT018159, CEP.701, CGP.082996, CGP.60474, CHIR.99021, CI.1040, Cisplatin, CMK, Cyclopamine, Cytarabine, Dasatinib, DMOG, Docetaxel, Doxorubicin, EHT.1864, Elesclomol, Embelin, Epothilone.B, Erlotinib, Etoposide, FH535, FTI.277, GDC.0449, GDC0941, Gefitinib, Gemcitabine, GNF.2, GSK269962A, GSK.650394, GW.441756, GW843682X, Imatinib, IPA.3, JNJ.26854165, JNK.9L, JNK.Inhibitor.VIII, JW.7.52.1, KIN001.135, KU.55933, Lapatinib, Lenalidomide, LFM.A13, Metformin, Methotrexate, MG.132, Midostaurin, Mitomycin.C, MK.2206, MS.275, Nilotinib, NSC.87877, NU.7441, Nutlin.3a, NVP.BEZ235, NVP.TAE684, Obatoclax.Mesylate, OSI.906, PAC.1, Paclitaxel, Parthenolide, Pazopanib, PD.0325901, PD.0332991, PD.173074, PF.02341066, PF.4708671, PF.562271, PHA.665752, PLX4720, Pyrimethamine, QS11, Rapamycin, RDEA119, RO.3306, Roscovitine, Salubrinal, SB.216763, SB590885, Shikonin, SL.0101.1, Sorafenib, S.Trityl.L.cysteine, Sunitinib, Temsirolimus, Thapsigargin, Tipifarnib, TW.37, Vinblastine, Vinorelbine, Vorinostat, VX.680, VX.702, WH.4.023, WO2009093972, WZ.1.84, X17.AAG, X681640, XMD8.85, Z.LLNle.CHO, ZM.447439.
#
#' @import car
#
#' @export
pRRopheticQQplot <- function(drug)
{
  possibleDrugs <- c("A.443654", "A.770041", "ABT.263", "ABT.888", "AG.014699", "AICAR", "AKT.inhibitor.VIII", "AMG.706", "AP.24534", "AS601245", "ATRA", "AUY922", "Axitinib", "AZ628", "AZD.0530", "AZD.2281", "AZD6244", "AZD6482", "AZD7762", "AZD8055", "BAY.61.3606", "Bexarotene", "BI.2536", "BIBW2992", "Bicalutamide", "BI.D1870", "BIRB.0796", "Bleomycin", "BMS.509744", "BMS.536924", "BMS.708163", "BMS.754807", "Bortezomib", "Bosutinib", "Bryostatin.1", "BX.795", "Camptothecin", "CCT007093", "CCT018159", "CEP.701", "CGP.082996", "CGP.60474", "CHIR.99021", "CI.1040", "Cisplatin", "CMK", "Cyclopamine", "Cytarabine", "Dasatinib", "DMOG", "Docetaxel", "Doxorubicin", "EHT.1864", "Elesclomol", "Embelin", "Epothilone.B", "Erlotinib", "Etoposide", "FH535", "FTI.277", "GDC.0449", "GDC0941", "Gefitinib", "Gemcitabine", "GNF.2", "GSK269962A", "GSK.650394", "GW.441756", "GW843682X", "Imatinib", "IPA.3", "JNJ.26854165", "JNK.9L", "JNK.Inhibitor.VIII", "JW.7.52.1", "KIN001.135", "KU.55933", "Lapatinib", "Lenalidomide", "LFM.A13", "Metformin", "Methotrexate", "MG.132", "Midostaurin", "Mitomycin.C", "MK.2206", "MS.275", "Nilotinib", "NSC.87877", "NU.7441", "Nutlin.3a", "NVP.BEZ235", "NVP.TAE684", "Obatoclax.Mesylate", "OSI.906", "PAC.1", "Paclitaxel", "Parthenolide", "Pazopanib", "PD.0325901", "PD.0332991", "PD.173074", "PF.02341066", "PF.4708671", "PF.562271", "PHA.665752", "PLX4720", "Pyrimethamine", "QS11", "Rapamycin", "RDEA119", "RO.3306", "Roscovitine", "Salubrinal", "SB.216763", "SB590885", "Shikonin", "SL.0101.1", "Sorafenib", "S.Trityl.L.cysteine", "Sunitinib", "Temsirolimus", "Thapsigargin", "Tipifarnib", "TW.37", "Vinblastine", "Vinorelbine", "Vorinostat", "VX.680", "VX.702", "WH.4.023", "WO2009093972", "WZ.1.84", "X17.AAG", "X681640", "XMD8.85", "Z.LLNle.CHO", "ZM.447439")
  if(!drug %in% possibleDrugs) stop(paste("ERROR: the drug specified was not found. Note drug names are case sensitive. Please select from: ", (paste(possibleDrugs, collapse=", "))));
  
  # load("/home/pgeeleher/postdoc_stuff/HDAC_project/Scripts/r_package_files/data/drugAndPhenoCgp.RData") # contains drugToCellLineDataCgp (to map cell lines to .CEL file names), gdsc_brainarray_syms (the gene expression data), drugSensitivityDataCgp (the drug sensitivity data)
  data(drugAndPhenoCgp) # contains drugToCellLineDataCgp (to map cell lines to .CEL file names), gdsc_brainarray_syms (the gene expression data), drugSensitivityDataCgp (the drug sensitivity data)
  
  colIc50Name <- paste(drug, "_IC_50", sep="")
  ic50s <- as.numeric(drugSensitivityDataCgp[, colIc50Name])
  names(ic50s) <- drugSensitivityDataCgp[ ,"Cell.Line"]
  whichNas <- which(is.na(ic50s))
  ic50s <- ic50s[-whichNas]
  
  offset = 0
  if(min(ic50s) < 0) # all numbers must be postive for a powerTranform to work, so make them positive.
  {
    offset <- -min(ic50s) + 1
    ic50s <- ic50s + offset
  }
    
  transForm <- powerTransform(ic50s)[[6]]
  ic50s <- ic50s^transForm
  
  qqnorm(ic50s, main=paste("QQplot on power-transformed IC50 values for", drug))
  qqline(ic50s, col="red")
}




