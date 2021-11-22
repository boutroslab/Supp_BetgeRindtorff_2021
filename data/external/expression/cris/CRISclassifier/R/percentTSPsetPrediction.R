percentTSPsetPrediction <- function(predMat) {
   classLabels <- paste("CRIS", LETTERS[1:5], sep="")
   colnames(predMat) <- gsub("^k[1234]\\.", "", colnames(predMat))
   classMax <- summary(factor(unlist(strsplit(colnames(predMat), split="_")),
                       levels=classLabels))
   out <- t(apply(predMat, 1, function(x, y) {
      out <- summary(factor(x, levels=y))}, 
      y=classLabels))
   colnames(out) <- classLabels
   out <- sweep(out, 2, classMax, "/")
}

