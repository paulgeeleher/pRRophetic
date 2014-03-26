setwd("/home/pgeeleher/postdoc_stuff/HDAC_project/Scripts/r_package_files/vignettes/")

## load the library
library(pRRophetic)
sapply(paste("/home/pgeeleher/postdoc_stuff/HDAC_project/Scripts/r_package_files/R/", dir("../R/"), sep=""), source)
library(car)
library(sva)
library(ridge)
library(ggplot2)
library(genefilter)

#################################################################
## PREDICTING CLINICAL OUTCOME FROM THE CGP CELL LINES DATA
#################################################################

# Load the bortezomib expression data. This is stored in a matrix with Gene Symbols as row names. Note that in order for "pRRopheticPredict" to work, the matrix MUST be annotated with official gene symbols.
data("bortezomibData.RData") # load("/home/pgeeleher/postdoc_stuff/HDAC_project/Scripts/r_package_files/data/bortezomibData.RData") # exprDataBortezomib, bortIndex, studyResponse, studyIndex,

# Assess the normality of the data. For many of the drug in CGP, the data deviate wildly from a normal distribtion and should not be fit as the response variable in a linear model (i.e. pRRopheticPredict is based on linear ridge regression).
pRRopheticQQplot("Bortezomib")

cvOut <- pRRopheticCV("Bortezomib", cvFold=5, testExprData=exprDataBortezomib)
summary(cvOut)
plot(cvOut)

# Based on the qqplot and the result of the Shapiro test, it is likely acceptable to use these data for prediction of bortezomib sensitivity. Predict bortezomib sensitivity using all cell lines, then only cell lines from haematological cancers and then only cell lines from derived from solid tumors
predictedPtype <- pRRopheticPredict(exprDataBortezomib, "Bortezomib", selection=1)
predictedPtype_blood <- pRRopheticPredict(exprDataBortezomib, "Bortezomib", "blood", selection=1)
predictedPtype_solid <- pRRopheticPredict(exprDataBortezomib, "Bortezomib", "allSolidTumors", selection=1)

# compare these three types of models.
t.test(predictedPtype[((studyResponse == "PGx_Responder = NR") & bortIndex)], predictedPtype[((studyResponse == "PGx_Responder = R") & bortIndex)], alternative="greater")
t.test(predictedPtype_blood[((studyResponse == "PGx_Responder = NR") & bortIndex)], predictedPtype_blood[((studyResponse == "PGx_Responder = R") & bortIndex)], alternative="greater")
t.test(predictedPtype_solid[((studyResponse == "PGx_Responder = NR") & bortIndex)], predictedPtype_solid[((studyResponse == "PGx_Responder = R") & bortIndex)], alternative="greater")

# make a boxplot of the results of the blood-only model...
df <- data.frame(Pred=c(predictedPtype_blood[((studyResponse == "PGx_Responder = NR") & bortIndex)], predictedPtype_blood[((studyResponse == "PGx_Responder = R") & bortIndex)]), Response=factor(c(rep("NR", sum((studyResponse == "PGx_Responder = NR") & bortIndex)), rep("R", sum((studyResponse == "PGx_Responder = R") & bortIndex)))))
ggplot(data=df, aes(y=Pred, x=Response, fill=Response)) + geom_boxplot(alpha=.3) + theme_bw() + ylab("Predicted Bortezomib Sensitivity")

#################################################################
## Predict on CCLE data? Using both linear and logistic models ##
#################################################################
# lets predict for PD0332991 using pRRopheticPredict
data(ccleData) # load(file="/home/pgeeleher/postdoc_stuff/HDAC_project/Scripts/r_package_files/data/ccleData.RData") # sensDataCcle, exprMatCcle, 


cvOut_pd <- pRRopheticCV("PD.0325901", cvFold=10, testExprData=exprMatCcle)
summary(cvOut_pd)
plot(cvOut_pd)

# run the prediction for PD0325901 on the CCLE data
predictedPtype_ccle <- pRRopheticPredict(exprMatCcle, "PD.0325901", selection=1)

# get the ActArea for the CCLE cell lines for which we have just predicted IC50
cellLinesWithCcleIc50s <- names(predictedPtype_ccle)[names(predictedPtype_ccle) %in% sensDataCcle$CCLE.Cell.Line.Name]
predCcleOrd <- predictedPtype_ccle[names(predictedPtype_ccle)]
ccleActArea_pd <- -sensDataCcle$"ActArea"[sensDataCcle$Compound == "PD-0325901"]
names(ccleActArea_pd) <- sensDataCcle$"CCLE.Cell.Line.Name"[sensDataCcle$Compound == "PD-0325901"]
ccleActAreaord <- ccleActArea_pd[cellLinesWithCcleIc50s]

# compare prediction to measured IC50, it is actually higher than the correlation achieved for remeasuring the drug sensitivity...
cor.test(predictedPtype_ccle[cellLinesWithCcleIc50s], ccleActAreaord, method="spearman")

plot(predictedPtype_ccle[cellLinesWithCcleIc50s], ccleActAreaord, xlab="", ylab="")


## ERLOTINIB
pRRopheticQQplot("Erlotinib")

# lets predict for erlotinib using pRRopheticLogisiticPredict. This function will return the log odds of sensitivity
predictedPtype_ccle_erlotinib <- pRRopheticLogisticPredict(exprMatCcle, "Erlotinib", selection=1)

# get the ActArea for the CCLE cell lines for which we have just predicted IC50
cellLinesWithCcleIc50s <- names(predictedPtype_ccle_erlotinib)[names(predictedPtype_ccle_erlotinib) %in% sensDataCcle$CCLE.Cell.Line.Name]
predCcleOrd <- predictedPtype_ccle_erlotinib[names(predictedPtype_ccle_erlotinib)]
ccleActArea_pd <- sensDataCcle$"ActArea"[sensDataCcle$Compound == "Erlotinib"]
names(ccleActArea_pd) <- sensDataCcle$"CCLE.Cell.Line.Name"[sensDataCcle$Compound == "Erlotinib"]
ccleActAreaord <- ccleActArea_pd[cellLinesWithCcleIc50s]

# There are a very large number of cell lines resistant to Erlotinib (within the drug screening window), so a correlation is not an appropriate measure of concordance. So lets do a t-test between some of the most sensitve and resistant cell lines to assess whether signal is being captured by the predictions.
resistant <- names(sort(ccleActAreaord))[1:55] # select 55 highly resistant cell lines.
sensitive <- names(sort(ccleActAreaord, decreasing=TRUE))[1:15] # select the 15 most sensitive cell lines.
t.test(predictedPtype_ccle_erlotinib[resistant], predictedPtype_ccle_erlotinib[sensitive])

# despite the fact that IC50 values are not correlated for this drug between these studies, the most sensitive/resistant samples are separated highly significantly with this logistic models.
## MAKE THIS A GGPLOT BOXPLOT....
boxplot(list(Resistant=predictedPtype_ccle_erlotinib[resistant], Sensitive=predictedPtype_ccle_erlotinib[sensitive]), pch=20, vertical=TRUE, method="jitter", ylab="Log-odds of sensitivity")


#################################################################
###### PREDICTION USING NON-CGP DATA AS A TRAINING SET ##########
#################################################################

## Include, as an example, prediction from the bortezomib clinical data. i.e. try to predict CR, PR, MR, NC, PD from CR, PR, MR, NC, PD 
## This serves as both an example of prediction directly from clinical data and 
## of using a dataset other than the CGP from which to predict.

# prepare the training data and test expression data.
trainExpr <- exprDataBortezomib[, (detailedResponse %in% c(1,2,3,4,5)) & studyIndex %in% c("studyCode = 25", "studyCode = 40")]
trainPtype <- detailedResponse[(detailedResponse %in% c(1,2,3,4,5)) & studyIndex %in% c("studyCode = 25", "studyCode = 40")]
testExpr <- exprDataBortezomib[, (detailedResponse %in% c(1,2,3,4,5)) & studyIndex %in% c("studyCode = 39")]

ptypeOut <- calcPhenotype(trainExpr, trainPtype, testExpr, selection=1)

testPtype <- detailedResponse[(detailedResponse %in% c(1,2,3,4,5)) & studyIndex %in% c("studyCode = 39")]
cor.test(testPtype, ptypeOut, alternative="greater")

# This t-test allows us to compare results directly to the cell line model. The cell line model outperforms the clinical model.
t.test(ptypeOut[testPtype %in% c(3,4,5)], ptypeOut[testPtype %in% c(1,2)], alternative="greater")










