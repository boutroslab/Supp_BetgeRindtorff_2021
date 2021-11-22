countTSPsetPrediction <- function(predMat) {
   classLabels <- paste("CRIS", LETTERS[1:5], sep="")
   out <- t(apply(predMat, 1, function(x, y) {
      out <- summary(factor(x, levels=y))}, 
      y=classLabels))
   colnames(out) <- classLabels
   out
}

