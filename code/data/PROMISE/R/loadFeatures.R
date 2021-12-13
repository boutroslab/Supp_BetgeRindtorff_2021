#' @title Load features for a well
#' 
#' @description Retrieve all features for a well
#' 
#' @param platename The plate ID
#' @param row The plate row
#' @param col The plate column
#' @param configdir The directory of the configuration files
#' 
#' @return A list of feature dataframes. The entries are as follows: 
#'   "organoid_features" correspond to features for individual organoids.
#'   "clump_features" corespond to features of entire contiguous foreground 
#'   patches without any (i.e. groups of connected organoids)
#'   "texture_features" correspond to the shape-independent haralick features 
#'   of the foreground without any object separation.
#' 
#' @author Jan Sauer
#' 
#' @examples print(loadFeaturesForWell)
#' @export
loadFeaturesForWell <- function(platename, row, col, configdir) {
  source(file.path(configdir, "watchdogConfig.R"))
  
  # The feature file could be in the main folder or in the 'wells' subfolder
  # if it's already been compressed
  feature_fn = featuresHdf5Filename(
    filedir = file.path(featuresdir, platename), platename = platename, 
    row = row, col = col)
  if(!file.exists(feature_fn)) {
    feature_fn = featuresHdf5Filename(
      filedir = file.path(featuresdir, platename, "wells"), 
      platename = platename, row = row, col = col)
  }
  if(!file.exists(feature_fn)) {
    warning(sprintf(
      "Feature file for '%s / %s / %s' doesn't exist",
      platename, row, col))
    return(NULL)
  }
  
  fo_key = ifelse(
    test = "features" %in% h5ls(feature_fn)$name, 
    yes = "features", no = "features_organoids")
  fon_key = ifelse(
    test = "feature_names" %in% h5ls(feature_fn)$name, 
    yes = "feature_names", no = "feature_names_organoids")
  
  features_organoids = data.frame(h5read(file = feature_fn, name = fo_key))
  colnames(features_organoids) = h5read(file = feature_fn, name = fon_key)
  
  return(features_organoids)
}