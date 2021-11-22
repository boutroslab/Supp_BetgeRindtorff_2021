### R code from vignette source 'cris.Rnw'

###################################################
### code chunk number 1: Library
###################################################
library(CRISclassifier)


###################################################
### code chunk number 2: cris.Rnw:115-119
###################################################
demo <- list.files(pattern="txt.gz$", system.file("data",package="CRISclassifier"), full.names=TRUE)
print(demo)
cris_classifier(input.exp.filename = demo, output.name="cris", nresmpl=1)



###################################################
### code chunk number 3: cris.Rnw:141-142 (eval = FALSE)
###################################################
## require(predictCRIS)


###################################################
### code chunk number 4: loadData
###################################################
data(matList)
data(phenoList)


###################################################
### code chunk number 5: cris.Rnw:161-166
###################################################
sapply(matList, class)
sapply(matList, dim)
sapply(phenoList, class)
sapply(phenoList, length)
sapply(phenoList, summary)


###################################################
### code chunk number 6: kTSPclassify
###################################################
### Valid gene expression matrix with all CRIS genes
newMat <- matList$Training
### To make predictions on 1 matrix
newPreds <- predictCRISclassKTSP(newMat)
### Counts classifications
summary(newPreds)
### NPT classification
refClass <- phenoList$Training


###################################################
### code chunk number 7: kTSPmultiClassify
###################################################
### For all matrices
newPredsList <- lapply(matList, predictCRISclassKTSP)
### Count classifications
lapply(newPredsList, summary)


###################################################
### code chunk number 8: sessioInfo
###################################################
toLatex(sessionInfo())


