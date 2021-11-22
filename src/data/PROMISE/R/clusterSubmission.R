
#' @title A helper script to submit an entire plate to the queue
#' 
#' @description A helper script to submit an entire plate to the queue
#' 
#' @param workflow The workflow function to call (e.g. PROMISEconversion, PROMISEfeatures)
#' @param plateIndir The incoming directory for the plate. Can be a dummy value for workflows that don't need it.
#' @param platename The plate that the well is on
#' @param configdir The directory of the configuration files
#' 
#' @return A boolean value indicating if the features were correctly calculated and saved
#' 
#' @author Jan Sauer
#' 
#' @examples print(submitPlateToQueue)
#' @export
submitPlateToQueue <- function(workflow, plateIndir, platename, configdir) {
  source(file.path(configdir, "watchdogConfig.R"))
  all_wells = wells()
  for(i in seq_len(nrow(all_wells))) {
    submitClusterQueue(
      plateIndir = plateIndir, platename = platename, 
      row = all_wells[i, "rows"], col = all_wells[i, "cols"], 
      configdir = configdir, custom_function = workflow)
  }
}