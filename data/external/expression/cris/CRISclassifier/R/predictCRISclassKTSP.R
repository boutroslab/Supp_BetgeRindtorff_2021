predictCRISclassKTSP <- function(mat) {
   ## Check inout
   if ( !is.matrix(mat) ){
      stop("The input mast be a numeric matrix.\n", 
	   "Row names for the input matrix must be valid gene symbols.\n")
   }
   ## Get individual predictions for each pair
   tspSetClassAllPreds <- getTSPsetPrediction(mat)
   ## Get count for each class
   tspSetClassPercent <- percentTSPsetPrediction(tspSetClassAllPreds)
   ## Classify according to max rule with TIES
   tspSetClassPredsFinal <- tspSetClassifyTiesDub(tspSetClassPercent, tspSetClassAllPreds)
   ## Make into factor: Handling TIES
   classLabels <- paste("CRIS", LETTERS[1:5], sep="")
   tspSetClassPredsFinal <- factor(tspSetClassPredsFinal, levels=classLabels)
   out <- list()
   out[[1]] <- tspSetClassPercent
   out[[2]] <- tspSetClassPredsFinal
   names(out) <- c("tspSetClassPercent","tspSetClassPredsFinal")
   return(out)
}


