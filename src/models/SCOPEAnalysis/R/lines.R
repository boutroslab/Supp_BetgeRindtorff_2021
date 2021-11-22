
#' @title Get Replicate IDs for Plates
#'
#' @description Get hard-coded replicate IDs for all the plates.
#'
#' @usage get_replicate_ids_for_plates()
#'
#' @return A named vector containing the replicate IDs
#'
#' @author Jan Sauer
#'
#' @examples print(get_replicate_ids_for_plates)
#' @export
get_replicate_ids_for_plates = function(plates = NULL) {
  rep_ids = c(
    "D004T01P003L02"= 1, "D004T01P004L03"= 1, "D004T01P005L08"= 1,
    "D004T01P007L02"= 2, "D004T01P008L03"= 2, "D004T01P009L08"= 2,
    "D007T01P003L02"= 1, "D007T01P004L03"= 1, "D007T01P005L08"= 1,
    "D007T01P007L02"= 2, "D007T01P008L03"= 2, "D007T01P009L08"= 2,
    "D010T01P017L02"= 1, "D010T01P018L03"= 1, "D010T01P019L08"= 1,
    "D010T01P021L02"= 2, "D010T01P022L03"= 2, "D010T01P023L08"= 2,
    "D013T01P001L02"= 1, "D013T01P002L03"= 1, "D013T01P003L08"= 1,
    "D013T01P905L02"= 2, "D013T01P906L03"= 2, "D013T01P907L08"= 2,
    "D018T01P901L02"= 1, "D018T01P902L03"= 1, "D018T01P003L08"= 1,
    "D018T01P905L02"= 2, "D018T01P906L03"= 2, "D018T01P907L08"= 2,
    "D019T01P001L02"= 1, "D019T01P002L03"= 1, "D019T01P003L08"= 1,
    "D019T01P005L02"= 2, "D019T01P006L03"= 2, "D019T01P007L08"= 2,
    "D020T01P001L02"= 1, "D020T01P002L03"= 1, "D020T01P003L08"= 1,
    "D020T01P905L02"= 2, "D020T01P906L03"= 2, "D020T01P907L08"= 2,
    "D020T02P009L02"= 1, "D020T02P010L03"= 1, "D020T02P011L08"= 1,
    "D020T02P013L02"= 2, "D020T02P014L03"= 2, "D020T02P015L08"= 2,
    "D021T01P901L02"= 1, "D021T01P902L03"= 1, "D021T01P003L08"= 1,
    "D021T01P905L02"= 2, "D021T01P906L03"= 2, "D021T01P907L08"= 2,
    "D022T01P001L02"= 1, "D022T01P002L03"= 1, "D022T01P003L08"= 1,
    "D022T01P005L02"= 2, "D022T01P906L03"= 2, "D022T01P907L08"= 2,
    "D027T01P001L02"= 1, "D027T01P002L03"= 1, "D027T01P003L08"= 1,
    "D027T01P905L02"= 2, "D027T01P906L03"= 2, "D027T01P907L08"= 2,
    "D030T01P001L02"= 1, "D030T01P002L03"= 1, "D030T01P003L08"= 1,
    "D030T01P905L02"= 2, "D030T01P906L03"= 2, "D030T01P907L08"= 2,
    "D046T01P001L02"= 1, "D046T01P002L03"= 1, "D046T01P003L08"= 1,
    "D046T01P005L02"= 2, "D046T01P006L03"= 2, "D046T01P007L08"= 2,
    # "D052T01P001L08"= 1, "D052T01P003L08"= 2,
    "D054T01P004L08"= 1, "D054T01P006L08"= 2,
    "D055T01P007L02"= 1, "D055T01P008L03"= 1, "D055T01P009L08"= 1,
    "D055T01P011L02"= 2, "D055T01P012L03"= 2, "D055T01P013L08"= 2)

  if(is.null(plates)) {
    return(rep_ids)
  } else {
    return(rep_ids[plates])
  }
}

#' @title Get All Lines
#'
#' @description Get the names of all lines
#'
#' @usage get_lines()
#'
#' @return A vector of line IDs
#'
#' @author Jan Sauer
#'
#' @examples print(get_lines)
#' @export
get_lines = function() {
  return(unique(substr(names(get_replicate_ids_for_plates()), 1, 7)))
}

#' @title Get All Plates
#'
#' @description Get the names of all plates
#'
#' @usage get_plates()
#'
#' @return A vector of plate IDs
#'
#' @author Jan Sauer
#'
#' @examples print(get_plates)
#' @export
get_plates = function() {
  return(unique(names(get_replicate_ids_for_plates())))
}
