#' @title Dummy function to test the watchdog
#' 
#' @description This function is a dummy to test the watchdog. It writes an empty text file per well.
#' 
#' @param plateIndir The base directory on the server, which contains a folder for each plate
#' @param platename The ID (and folder name) of the plate
#' @param row The row of the well
#' @param col The column of the well
#' @param configdir The configuration file, which defines the number of z-stacks, channels, and fields
#' 
#' @return A boolean value indicating if the projection was successful
#' 
#' @author Bernd Fischer
#' 
#' @examples print(dummyfct)
#' @export
dummyfct <- function(plateIndir, platename, row, col, configdir) {
    # TODO: Define from config file
  source(file.path(configdir, "watchdogConfig.R"))

  filenames = imageFilenames(plateIndir, platename, row, col, configdir)
  
  # Create hdf5 file (and folders if necessary)
  dir.create(file.path(hdf5dir, platename), showWarnings = FALSE, recursive = TRUE)
  h5FileName = hdf5Filename(filedir = file.path(hdf5dir, platename), platename = platename, row = row, col = col)
  h5FileName = gsub(".h5",".txt",h5FileName)
  file.create(h5FileName)

  return(TRUE)
}

