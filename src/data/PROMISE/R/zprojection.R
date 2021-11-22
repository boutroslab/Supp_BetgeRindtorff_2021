#' @title Project all images in a well
#' 
#' @description This function projects all images in a well and stores the results in an hdf5 file
#' 
#' @param plateIndir The base directory on the server, which contains a folder for each plate
#' @param platename The ID (and folder name) of the plate
#' @param row The row of the well
#' @param col The column of the well
#' @param configdir The configuration file, which defines the number of z-stacks, channels, and fields
#' 
#' @return A boolean value indicating if the projection was successful
#' 
#' @author Jan Sauer
#' 
#' @examples print(projection)
#' @export
projection <- function(plateIndir, platename, row, col, configdir) {
  source(file.path(configdir, "watchdogConfig.R"))

  # Load input images
  inputFile = hdf5Filename(filedir = file.path(hdf5dir, platename), platename = platename, row = row, col = col)
  metadata = h5read(file = inputFile, name = "metadata")
  metadata[metadata[,1] == "PackageVersion",2] = as.character(packageVersion("PROMISE"))

  h5f = H5Fopen(inputFile)
  h5d = H5Dopen(h5f,name = "images")
  h5s = H5Dget_space(h5d)
  Dims = H5Sget_simple_extent_dims(h5s)
  H5Sclose(h5s)
  H5Dclose(h5d)
  H5Fclose(h5f)

  if(Dims$rank != 5) {
      warning("The package expects data with 5 dimensions: [x, y, z-stacks, channels, fields]")
      return(FALSE)
  }

  xDim = Dims$size[1]
  yDim = Dims$size[2]
  nStacks = Dims$size[3]
  nChannels = Dims$size[4]
  nFields = Dims$size[5]
  
  # Create hdf5 file (and directory if necessary)
  dir.create(file.path(hdf5projection, platename), showWarnings = FALSE, recursive = TRUE)
  h5FileName = projectionHdf5Filename(filedir = file.path(hdf5projection, platename), platename = platename, row = row, col = col, configdir = configdir)
  file_made = h5createFile(h5FileName)
  if(!file_made) {H5close(); return(FALSE)}

  dataset_made = h5createDataset(file = h5FileName, dataset = "images", dims = c(xDim, yDim, nChannels, nFields), 
                                 H5type = "H5T_NATIVE_UINT16", chunk = c(hdf5_chunksize, hdf5_chunksize, 1, 1), level = hdf5_compression)
  if(!dataset_made) {H5close(); return(FALSE)}
  h5write(obj = metadata, file = h5FileName, name = "metadata", level = hdf5_compression)

  t = Sys.time()  
  t1 = t
  for(c in seq_len(nChannels)) {
    for(f in seq_len(nFields)) {
        cat("c = ",c," f = ",f)
        
        # TODO: Find proper parameters
        img_data = h5read(file = inputFile, name = "images", index=list(NULL,NULL,NULL,c,f))[,,,1,1]
        output_array = contrastProjection(imageStack = img_data, w_x = 35, w_y = 35, smoothing = 5, 
                                               brushShape = "disc", interpolation = 3, fix.gaussian.blur = TRUE, 
                                               blur.size = 100)
        h5write(obj = output_array, file = h5FileName, name = "images", index=list(NULL,NULL,c,f))
        t2 = Sys.time()
        cat(" time=",format(t2 - t1),"\n")
        t1 = t2
    }
  }
  t2 = Sys.time()
  cat("total time=",format(t2 - t),"\n")
  
  # Create hdf5 file (and directory if necessary)
  # dir.create(file.path(hdf5projection, platename), showWarnings = FALSE, recursive = TRUE)
  # h5FileName = projectionHdf5Filename(filedir = file.path(hdf5projection, platename), platename = platename, row = row, col = col)
  # file_made = h5createFile(h5FileName)
  # if(!file_made) {H5close(); return(FALSE)}
  # 
  # # Create dataset
  # dataset_made = h5createDataset(file = h5FileName, dataset = "images", dims = c(2048, 2048, nChannels, nFields), 
  #                                H5type = "H5T_NATIVE_UINT16", chunk = c(hdf5_chunksize, hdf5_chunksize, 1, 1), level = hdf5_compression)
  # if(!dataset_made) {H5close(); return(FALSE)}
  # 
  # metadata_made = h5createDataset(file = h5FileName, dataset = "metadata", dims = dim(metadata), storage.mode = "character", size = 50, level = hdf5_compression)
  # if(!metadata_made) {H5close(); return(FALSE)}
  # 
  # # Write data
  # h5write(obj = output_array, file = h5FileName, name = "images", index=list(NULL, NULL, NULL, NULL))
  # h5write(obj = metadata, file = h5FileName, name = "metadata")
  
  H5close()
  return(TRUE)
}
