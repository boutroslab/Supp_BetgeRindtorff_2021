#' @title Load layout library
#' 
#' @description Loads the layout library
#' 
#' @param library The library ID (e.g. L02)
#' @param configdir The directory of the configuration files
#' @param cols A custom vector of columns to return. Column names not found 
#' in the library file are silently ignored.
#' 
#' @return A library dataframe
#' 
#' @author Jan Sauer
#' 
#' @examples print(loadLibrary)
#' @export
loadLibrary <- function(
  library, configdir, 
  cols = c("Well_ID_384", "Product.Name", "concentration")) {
  
  source(file.path(configdir, "watchdogConfig.R"))
  
  # Load layout
  layout = read.csv2(file = file.path(layoutdir, sprintf("%s.txt", library)))
  if("Well_ID_384" %in% colnames(layout)) row.names(layout) = layout$Well_ID_384
  
  # Keep only columns also in the layout file
  cols = cols[cols %in% colnames(layout)]
  
  layout = layout[,cols]
  return(layout)
}

#' @title Load layout library fron Plate
#' 
#' @description Loads the layout library corresponding to the plate
#' 
#' @param platename The plate ID
#' @param configdir The directory of the configuration files
#' @param cols A custom vector of columns to return. Column names not found 
#' in the library file are silently ignored.
#' 
#' @return A library dataframe
#' 
#' @author Jan Sauer
#' 
#' @examples print(loadLibraryFromPlate)
#' @export
loadLibraryFromPlate <- function(
  platename, configdir, 
  cols = c("Well_ID_384", "Product.Name", "concentration")) {
  
  source(file.path(configdir, "watchdogConfig.R"))
  
  # Load layout
  regex = regexpr(pattern = "L[0-9]{2}", text = platename)
  layout_id = substr(
    x = platename, start = regex, stop = regex + attr(regex, "match.length"))
  layout = read.csv2(file = file.path(layoutdir, sprintf("%s.txt", layout_id)))
  if("Well_ID_384" %in% colnames(layout)) row.names(layout) = layout$Well_ID_384
  
  # Keep only columns also in the layout file
  cols = cols[cols %in% colnames(layout)]
  
  layout = layout[,cols]
  return(layout)
}