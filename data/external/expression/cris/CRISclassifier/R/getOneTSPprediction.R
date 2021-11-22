getOneTSPprediction <- function(tsp, mat) {
   out <- factor(mat[ tsp[1] , ] > mat[ tsp[2] , ], levels=c(TRUE, FALSE))
   levels(out) <- names(tsp)
   out
}

