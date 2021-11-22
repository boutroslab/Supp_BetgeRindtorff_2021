oneTspSetClassifyTiesDub <- function(x, y) {
   out <- x[ which(x == max(x))]
   if ( length(out) == 1 ) ans <- names(which.max(out))
   else if ( length(out) == 2 ) {
      nm <- names(out == max(out))
      nm <- combn(nm, 2, paste, collapse="_")
      ans <- unique(y[names(y) %in% nm])
      if ( length(ans) > 1 ) ans <- NA
      } else ans <- NA
      ans
}

