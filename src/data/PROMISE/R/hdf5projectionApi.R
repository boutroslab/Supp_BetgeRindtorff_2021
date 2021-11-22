#' @title Load or save a projection from hdf5 to an image
#' 
#' @description Loads the hdf5projection file for a given well and returns it 
#' as an array and optionally saves it to a TIF file. The fields are saved as 
#' a single image. Note that there is NO GUARANTEE that the TIFF images have 
#' the true intensities. Analysis should be performed exclusively on the hdf5 
#' files. 
#' 
#' The function can optionally preprocess the intensities. The raw images are 
#' most likely too dark to see the the naked eye. To fix this, the brightest
#' 0.1% of all pixels can be cropped for each channel independently. This 
#' makes the picture more easily viewable with the naked eye. This 
#' preprocessing only applies to the saved file. The function returns the raw 
#' integer values in the hdf5 file.
#' 
#' @param platename The plate that the well is on
#' @param row The row of the well
#' @param col The column of the well
#' @param configdir The directory of the configuration files
#' @param outfile The file to write to. Leave as NULL to not save
#' @param preprocess A boolean indicating if the intensities should be 
#' preprocessed.
#' 
#' @return A 4D array with the dimensions (X, Y, channel, field)
#' 
#' @author Jan Sauer
#' 
#' @examples print(loadProjection)
#' @export
loadProjection <- function(
  platename, row, col, configdir, outfile=NULL, preprocess=TRUE) {
  source(file.path(configdir, "watchdogConfig.R"))
  
  proj_fn = projectionHdf5Filename(
    filedir = file.path(hdf5projection, platename), 
    platename = platename, row = row, col = col)
  
  img = h5read(proj_fn, "images")
  if(!is.null(outfile)) {
    if(!endsWith(outfile, "tif") & !endsWith(outfile, "tiff")) {
      warning(paste0(
        "R is a bit special about how it saves images, so this function only ",
        "writes TIFF images. Please make sure that the output filename ends ",
        "with *.tif or *.tiff"))
      return(img)
    }
    dir.create(path = dirname(outfile), showWarnings = FALSE, recursive = TRUE)
    xdim = dim(img)[1]
    ydim = dim(img)[2]
    num_channels = dim(img)[3]
    num_fields = dim(img)[4]
    output_img = NULL
    if(num_fields == 4) {
      output_img = Image(data = 0, dim = c(xdim * 2, ydim * 2, num_channels), colormode = "color")
      output_img[1:xdim, 1:ydim,] = img[,,,as.numeric(fields_layout[1])]
      output_img[(xdim+1):(2*xdim), 1:ydim,] = img[,,,as.numeric(fields_layout[2])]
      output_img[1:xdim, (ydim+1):(2*ydim),] = img[,,,as.numeric(fields_layout[3])]
      output_img[(xdim+1):(2*xdim), (ydim+1):(2*ydim),] = img[,,,as.numeric(fields_layout[4])]
    } else {
      warning(sprintf("This function cannot handle layouts with %s fields", num_fields))
      return(img)
    }
    
    # Preprocess
    if(preprocess) {
      for(ch in seq_len(num_channels)) {
        q = quantile(x=output_img[,,ch], probs=c(0.1, 0.999))
        output_img@.Data[,,ch] = output_img@.Data[,,ch] - q[1]
        output_img@.Data[,,ch][output_img@.Data[,,ch] < 0] = 0
        output_img@.Data[,,ch][output_img@.Data[,,ch] > q[2]] = q[2]
      }
      output_img = normalize(output_img)
    } else {
      output_img = output_img / 65535
    }
    
    writeImage(x = output_img, files = outfile, type = "tiff")
  }
  return(img)
}

#' @title Load or save a segmentation mask from hdf5 to an image
#' 
#' @description Loads the segmentation mask for a given well and returns it 
#' as an array and optionally saves it to a TIF file. The fields are saved as 
#' a single image.
#' 
#' @param platename The plate that the well is on
#' @param row The row of the well
#' @param col The column of the well
#' @param configdir The directory of the configuration files
#' @param outfile The file to write to. Leave as NULL to not save
#' 
#' @return A 3D array with the dimensions (X, Y, field)
#' 
#' @author Jan Sauer
#' 
#' @examples print(loadSegmentationMask)
#' @export
loadSegmentationMask <- function(
  platename, row, col, configdir, outfile=NULL) {
  source(file.path(configdir, "watchdogConfig.R"))
  
  segmask_fn = segmentationHdf5Filename(
    filedir = file.path(segmentationdir, platename), 
    platename = platename, row = row, col = col)
  
  img = h5read(segmask_fn, "mask")
  if(!is.null(outfile)) {
    if(!endsWith(outfile, "tif") & !endsWith(outfile, "tiff")) {
      warning(paste0(
        "R is a bit special about how it saves images, so this function only ",
        "writes TIFF images. Please make sure that the output filename ends ",
        "with *.tif or *.tiff"))
      return(img)
    }
    dir.create(path = dirname(outfile), showWarnings = FALSE, recursive = TRUE)
    xdim = dim(img)[1]
    ydim = dim(img)[2]
    num_fields = dim(img)[3]
    output_img = NULL
    if(num_fields == 4) {
      output_img = Image(data = 0, dim = c(xdim * 2, ydim * 2))
      output_img[1:xdim, 1:ydim] = img[,,as.numeric(fields_layout[1])]
      output_img[(xdim+1):(2*xdim), 1:ydim] = img[,,as.numeric(fields_layout[2])]
      output_img[1:xdim, (ydim+1):(2*ydim)] = img[,,as.numeric(fields_layout[3])]
      output_img[(xdim+1):(2*xdim), (ydim+1):(2*ydim)] = img[,,as.numeric(fields_layout[4])]
    } else {
      warning(sprintf("This function cannot handle layouts with %s fields", num_fields))
      return(img)
    }
    
    writeImage(x = output_img, files = outfile, type = "tiff")
  }
  return(img)
}