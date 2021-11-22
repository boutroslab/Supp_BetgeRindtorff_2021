tspSetClassifyTiesDub <- function(countMat, preds) {
   ## Classify according to max rule (NA is not possible)
   countList <- lapply(apply(countMat, 1, list), unlist)
   colnames(preds) <- gsub("^k[1234]\\.", "", colnames(preds))
   predsList <- lapply(apply(preds, 1, list), unlist)
   out <- mapply(oneTspSetClassifyTiesDub, countList, predsList,
		 SIMPLIFY=TRUE)
}

