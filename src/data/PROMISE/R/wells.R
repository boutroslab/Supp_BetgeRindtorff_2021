
#' @title Returns a data.frame with row names, col names and well names for the specified mult-well plate
#' 
#' @description This function returns a data.frame with well names, row names, and col names for the specified
#' multi-well plate. Currently supported are 96-well plates, 384 well plates, and 1536 well plates
#' 
#' @param nrWells The design of the multi-well plate. Either 96, 384 (default), or 1536
#' 
#' @return A data.frame with nrWells rows with well names, row names and col names.
#' 
#' @author Bernd Fischer
#' 
#' @examples print(wells)
#' @export
wells <- function(nrWells=384) {
    if (nrWells == 384) {
        wells = data.frame(names = sprintf("%s%0.2d",
                                           rep(LETTERS[1:16],each=24),
                                           rep(1:24,times=16)),
                           rows = rep(LETTERS[1:16],each=24),
                           cols = sprintf("%0.2d",
                                          rep(1:24,times=16)),
                           stringsAsFactors = FALSE)
        attr(wells,"nrwells") = nrWells
        attr(wells,"nrows") = 16
        attr(wells,"ncols") = 24
        attr(wells,"rows") = LETTERS[1:16]
        attr(wells,"cols") = sprintf("%0.2d",1:24)
    } else {
        if (nrWells == 96) {
            wells = data.frame(names = sprintf("%s%0.2d",
                                               rep(LETTERS[1:8],each=12),
                                               rep(1:12,times=8)),
                               rows = rep(LETTERS[1:8],each=12),
                               cols = sprintf("%0.2d",
                                              rep(1:12,times=8)),
                               stringsAsFactors = FALSE)
            attr(wells,"nrwells") = nrWells
            attr(wells,"nrows") = 8
            attr(wells,"ncols") = 12
            attr(wells,"rows") = LETTERS[1:8]
            attr(wells,"cols") = sprintf("%0.2d",1:12)
        } else {
            if (nrWells == 1536) {
                wells = data.frame(names = sprintf("%s%s%0.2d",
                                                   rep(rep(c("A","B"),times=c(26,6)),each=48),
                                                   rep(LETTERS[c(1:26,1:6)],each=48),
                                                   rep(1:48,times=32)),
                                   rows = sprintf("%s%s",
                                                  rep(rep(c("A","B"),times=c(26,6)),each=48),
                                                  rep(LETTERS[c(1:26,1:6)],each=48)),
                                   cols = sprintf("%0.2d",
                                                  rep(1:48,times=32)),
                                   stringsAsFactors = FALSE)
                attr(wells,"nrwells") = nrWells
                attr(wells,"nrows") = 32
                attr(wells,"ncols") = 48
                attr(wells,"rows") = LETTERS[1:32]
                attr(wells,"cols") = sprintf("%0.2d",1:48)
            } else {
                stop("platedesign with ",nrWells," wells unknown")
            }
        }
    }
    wells
}

