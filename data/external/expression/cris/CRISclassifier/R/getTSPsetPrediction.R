getTSPsetPrediction <- function(mat) {
   tsp <- makeCRIStspSet()
   gns <- as.vector(rownames(mat))
   tsp <- tsp[ sapply(tsp, function(x, y) all(x %in% y), y=gns ) ]
   availableTsp <- unique(gsub("^k[1234]\\.", "", names(tsp)))
   if ( length(availableTsp) < 10 ) {
      stop("You need at least one TSP for CRIS class comparison.\n",
           "Row names for the input matrix must be valid gene symbols.\n")
   }
   sapply(tsp, getOneTSPprediction, mat)
}


