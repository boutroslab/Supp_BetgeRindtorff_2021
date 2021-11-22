#' @title Thumbnail images of z-projecion
#' 
#' @description This function generates thumbnail images of the z-rpojection
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
#' @examples print(thumbnailImages)
#' @export
thumbnailImages <- function(plateIndir, platename, row, col, configdir) {
  source(file.path(configdir, "watchdogConfig.R"))

  # Load input images
  inputFile = projectionHdf5Filename(filedir = file.path(hdf5projection, platename), platename = platename, row = row, col = col, configdir = configdir)
  metadata = h5read(file = inputFile, name = "metadata")
  Img = h5read(inputFile, name="images")

  if(length(dim(Img)) != 4) {
      warning("The package expects data with 4 dimensions: [x, y, channels, fields]")
      return(FALSE)
  }

  # Quantile-truncate and normalize each channel separately
  for(c in seq_len(dim(Img)[3])) {
    Img[,,c,] = pmax(Img[,,c,] - median(Img[,,c,]), 0)
    Img[,,c,] = pmin(Img[,,c,], quantile(Img[,,c,], 0.999))
    Img[,,c,] = normalize(Img[,,c,])
  }
  
  Img2 = resize(Img, w = 2048 / 8, h = 2048 / 8)
  # If Img2 has only two channels then insert a blank green channel
  if(dim(Img2)[3] == 2) {
    Img2 = abind(
      Img2[,,1,,drop=FALSE], 
      array(data = 0, dim = dim(Img2[,,1,,drop=FALSE])),
      Img2[,,2,,drop=FALSE], along = 3)
  }
  Img3 = rgbImage(tile(Img2[,,1,as.integer(fields_layout)], nx = sqrt(length(fields))),
                  tile(Img2[,,2,as.integer(fields_layout)], nx = sqrt(length(fields))),
                  tile(Img2[,,3,as.integer(fields_layout)], nx = sqrt(length(fields))))

  Img2 = resize(Img, w = 2048 / 64, h = 2048 / 64)
  if(dim(Img2)[3] == 2) {
    Img2 = abind(
      Img2[,,1,,drop=FALSE], 
      array(data = 0, dim = dim(Img2[,,1,,drop=FALSE])),
      Img2[,,2,,drop=FALSE], along = 3)
  }
  Img4 = rgbImage(tile(Img2[,,1,as.integer(fields_layout)], nx = sqrt(length(fields))),
                  tile(Img2[,,2,as.integer(fields_layout)], nx = sqrt(length(fields))),
                  tile(Img2[,,3,as.integer(fields_layout)], nx = sqrt(length(fields))))

  if (!file.exists(file.path(htmldir, platename))) {
      dir.create(file.path(htmldir, platename), showWarnings = FALSE, recursive = TRUE)
  }
  if (row == "A" & col == "01") {
      makeIndexPage(platename, configdir)
  }
  fn1 = thumbnailFilename(filedir = file.path(htmldir, platename), platename = platename, row = row, col = col, configdir = configdir, level=1)
  writeImage(Img3, fn1, quality = 99)
  fn2 = thumbnailFilename(filedir = file.path(htmldir, platename), platename = platename, row = row, col = col, configdir = configdir, level=2)
  writeImage(Img4, fn2, quality = 99)

  H5close()
  return(TRUE)
}


#' @title Thumbnail images of DNN Segmentation Masks
#' 
#' @description This function generates thumbnail images of the segmentation masks
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
#' @examples print(thumbnailSegmentation)
#' @export
thumbnailSegmentation <- function(plateIndir, platename, row, col, configdir) {
  source(file.path(configdir, "watchdogConfig.R"))
  
  # Load segmentation masks
  inputFile = segmentationHdf5Filename(filedir = file.path(segmentationdir, platename), platename = platename, row = row, col = col)
  metadata = h5read(file = inputFile, name = "metadata")
  Img = h5read(inputFile, name="mask")
  
  if(length(dim(Img)) != 3) {
    warning("The package expects data with 3 dimensions: [x, y, fields]")
    return(FALSE)
  }
  
  # Normalize
  Img = Img / 2
  Img2 = resize(Img, w = 2048 / 8, h = 2048 / 8)
  Img3 = tile(Img2[,,as.integer(fields_layout)], nx = sqrt(length(fields)))
  
  Img2 = resize(Img, w = 2048 / 64, h = 2048 / 64)
  Img4 = tile(Img2[,,as.integer(fields_layout)], nx = sqrt(length(fields)))
  
  if (!file.exists(file.path(segmentationhtmldir, platename))) {
    dir.create(file.path(segmentationhtmldir, platename), showWarnings = FALSE, recursive = TRUE)
  }
  if (row == "A" & col == "01") {
    makeIndexPage_segmentation(platename, configdir)
  }
  fn1 = thumbnailFilename(filedir = file.path(segmentationhtmldir, platename), platename = platename, row = row, col = col, level=1)
  writeImage(Img3, fn1, quality = 99)
  fn2 = thumbnailFilename(filedir = file.path(segmentationhtmldir, platename), platename = platename, row = row, col = col, level=2)
  writeImage(Img4, fn2, quality = 99)
  
  H5close()
  return(TRUE)
}


#' @title Writes an index html-page for the thumbnail images of z-projecion
#' 
#' @description This function writes an index.html page for the thumbnail images of the z-projection
#' 
#' @param platename The ID (and folder name) of the plate
#' @param configdir The configuration file, which defines the number of z-stacks, channels, and fields
#' 
#' @return A boolean value indicating if the index file was written successfully
#' 
#' @author Bernd Fischer
#' 
#' @examples print(makeIndexPage)
#' @export
makeIndexPage <- function(platename, configdir) {
    source(file.path(configdir, "watchdogConfig.R"))
    
    row = "A"
    col="01"
    
    # Load input images
    inputFile = projectionHdf5Filename(filedir = file.path(hdf5projection, platename), platename = platename, row = row, col = col, configdir = configdir)
    metadata = h5read(file = inputFile, name = "metadata")

    W = wells(nrWells)
    fn1 = thumbnailFilename(filedir = file.path(htmldir, platename), platename = platename, row = W$row, col = W$col, configdir = configdir, level=1, addPath=FALSE)
    fn2 = thumbnailFilename(filedir = file.path(htmldir, platename), platename = platename, row = W$row, col = W$col, configdir = configdir, level=2, addPath=FALSE)
    M = matrix(fn2, nrow=attr(W,"nrows"), byrow = TRUE)
    L = matrix(fn1, nrow=attr(W,"nrows"), byrow = TRUE)
    page = openPage(file.path(htmldir, platename, "index.html"), link.css="hwriter.css")
    hwrite(metadata[metadata[,1] == "Plate",2], heading=1, page=page)
    M = hwriteImage(M, link=L, table=FALSE)
    row.names(M) = attr(W, "rows")
    colnames(M) = attr(W, "cols")
    hwrite(M, link=L, page=page, br=TRUE)
    hwrite(metadata, page=page)
    closePage(page, splash=FALSE)
    file.copy(system.file("images","hwriter.css", package="hwriter"),file.path(htmldir, platename, "hwriter.css"))
    return(TRUE)   
}


#' @title Writes an index html-page for the thumbnail images of the DNN segmentation
#' 
#' @description This function writes an index.html page for the thumbnail images of the z-projection
#' 
#' @param platename The ID (and folder name) of the plate
#' @param configdir The configuration file, which defines the number of z-stacks, channels, and fields
#' 
#' @return A boolean value indicating if the index file was written successfully
#' 
#' @author Jan Sauer
#' 
#' @examples print(makeIndexPage_segmentation)
#' @export
makeIndexPage_segmentation <- function(platename, configdir) {
  source(file.path(configdir, "watchdogConfig.R"))
  
  row = "A"
  col = "01"
  
  # Load input images
  inputFile = segmentationHdf5Filename(filedir = file.path(segmentationdir, platename), platename = platename, row = row, col = col)
  metadata = h5read(file = inputFile, name = "metadata")
  
  W = wells(nrWells)
  fn1 = thumbnailFilename(filedir = file.path(segmentationhtmldir, platename), platename = platename, row = W$row, col = W$col, level=1, addPath=FALSE)
  fn2 = thumbnailFilename(filedir = file.path(segmentationhtmldir, platename), platename = platename, row = W$row, col = W$col, level=2, addPath=FALSE)
  M = matrix(fn2, nrow=attr(W,"nrows"), byrow = TRUE)
  L = matrix(fn1, nrow=attr(W,"nrows"), byrow = TRUE)
  page = openPage(file.path(segmentationhtmldir, platename, "index.html"), link.css="hwriter.css")
  hwrite(metadata[metadata[,1] == "Plate",2], heading=1, page=page)
  M = hwriteImage(M, link=L, table=FALSE)
  row.names(M) = attr(W, "rows")
  colnames(M) = attr(W, "cols")
  hwrite(M, link=L, page=page, br=TRUE)
  hwrite(metadata, page=page)
  closePage(page, splash=FALSE)
  file.copy(system.file("images","hwriter.css", package="hwriter"),file.path(segmentationhtmldir, platename, "hwriter.css"))
  return(TRUE)   
}


#' @title Writes a sorted index html-page for the thumbnail images of z-projecion
#' 
#' @description This function writes an index.html page for the thumbnail images of the z-projection, which is sorted according to drug and concentration
#' 
#' @param platename The ID (and folder name) of the plate
#' @param configdir The configuration file, which defines the number of z-stacks, channels, and fields
#' 
#' @return A boolean value indicating if the index file was written successfully
#' 
#' @author Jan Sauer
#' 
#' @examples print(makeSortedIndexPage)
#' @export
makeSortedIndexPage <- function(platename, configdir) {
  source(file.path(configdir, "watchdogConfig.R"))
  
  # Load metadata
  row = "A"
  col= "01"
  inputFile = projectionHdf5Filename(filedir = file.path(hdf5projection, platename), platename = platename, row = row, col = col, configdir = configdir)
  metadata = h5read(file = inputFile, name = "metadata")
  
  # Load well names
  W = wells(nrWells)
  fn1 = thumbnailFilename(filedir = file.path(htmldir, platename), platename = platename, row = W$row, col = W$col, configdir = configdir, level=1, addPath=FALSE)
  fn2 = thumbnailFilename(filedir = file.path(htmldir, platename), platename = platename, row = W$row, col = W$col, configdir = configdir, level=2, addPath=FALSE)
  
  # Load layout and append the full well filename as a column
  layoutID = metadata[metadata[,1] == "Layout", 2]
  # layout = read.table(file.path(layoutdir, sprintf("%s.txt", layoutID)), header = TRUE, stringsAsFactors = FALSE)
  layout = read.csv2(file.path(layoutdir, sprintf("%s.txt", layoutID)), header = TRUE, stringsAsFactors = FALSE)
  layout$filename_large = fn1[match(layout$Well_ID_384, W$names)]
  layout$filename_small = fn2[match(layout$Well_ID_384, W$names)]
  
  # Extract the negative and positive controls first
  treatments_NegCTRL = sort(unique(layout[startsWith(layout$substance, "DMSO"), "Product.Name"]))
  M_neg = lapply(treatments_NegCTRL, function(x) layout[layout$Product.Name == x, "filename_small"])
  M_neg = do.call(what = rbind, args = M_neg)
  rownames(M_neg) = treatments_NegCTRL
  L_neg = lapply(treatments_NegCTRL, function(x) layout[layout$Product.Name == x, "filename_large"])
  L_neg = do.call(what = rbind, args = L_neg)
  rownames(L_neg) = treatments_NegCTRL
  I_neg = lapply(treatments_NegCTRL, function(x) layout[layout$Product.Name == x, "Well_ID_384"])
  I_neg = do.call(what = rbind, args = I_neg)
  rownames(I_neg) = treatments_NegCTRL
  
  treatments_PosCTRL = sort(unique(layout[startsWith(layout$substance, "Stauro"), "Product.Name"]))
  M_pos = lapply(treatments_PosCTRL, function(x) layout[layout$Product.Name == x, "filename_small"])
  M_pos = do.call(what = rbind, args = M_pos)
  rownames(M_pos) = treatments_PosCTRL
  L_pos = lapply(treatments_PosCTRL, function(x) layout[layout$Product.Name == x, "filename_large"])
  L_pos = do.call(what = rbind, args = L_pos)
  rownames(L_pos) = treatments_PosCTRL
  I_pos = lapply(treatments_PosCTRL, function(x) layout[layout$Product.Name == x, "Well_ID_384"])
  I_pos = do.call(what = rbind, args = I_pos)
  rownames(I_pos) = treatments_PosCTRL
  
  # Extract the well names for each drug and sort them by concentration (-9, -3, 1, 3, 9)
  treatments = sort(unique(layout[startsWith(layout$substance, "Drug"), "Product.Name"]))
  M = list()
  L = list()
  I = list()
  for(treatment in treatments) {
    sublayout = layout[layout$Product.Name == treatment,]
    sublayout = sublayout[order(sublayout$concentration, decreasing = FALSE),]
    M[[treatment]] = sublayout$filename_small
    names(M[[treatment]]) = sublayout$concentration
    L[[treatment]] = sublayout$filename_large
    names(L[[treatment]]) = sublayout$concentration
    I[[treatment]] = sublayout$Well_ID_384
    names(I[[treatment]]) = sublayout$concentration
  }
  
  if(length(unique(lapply(M, names))) != 1 | length(unique(lapply(L, names))) != 1 | length(unique(lapply(I, names))) != 1) {
    warning("Layout file must have the same number of concentration replicates for each non-control drug")
    return(NULL)
  }
  
  M = do.call(rbind, M)
  L = do.call(rbind, L)
  I = do.call(rbind, I)
  
  # Write html page
  page = openPage(file.path(htmldir, platename, "index_sorted.html"), link.css="hwriter.css")
  hwrite(metadata[metadata[,1] == "Plate",2], heading=1, page=page)
  M_neg = hwriteImage(M_neg, link=L_neg, table=FALSE)
  row.names(M_neg) = row.names(L_neg)
  colnames(M_neg) = colnames(L_neg)
  for(i in seq_along(M_neg)) { M_neg[i] = paste0(M_neg[i], "<br>", I_neg[i]) }
  M_pos = hwriteImage(M_pos, link=L_pos, table=FALSE)
  row.names(M_pos) = row.names(L_pos)
  colnames(M_pos) = colnames(L_pos)
  for(i in seq_along(M_pos)) { M_pos[i] = paste0(M_pos[i], "<br>", I_pos[i]) }
  M = hwriteImage(M, link=L, table=FALSE)
  row.names(M) = row.names(L)
  colnames(M) = colnames(L)
  for(i in seq_along(M)) { M[i] = paste0(M[i], "<br>", I[i]) }
  hwrite("Negative Controls", heading=2, page=page)
  hwrite(M_neg, page=page, br=TRUE, style="text-align:center;")
  hwrite("Positive Controls", heading=2, page=page)
  hwrite(M_pos, page=page, br=TRUE, style="text-align:center;")
  hwrite("Treatments", heading=2, page=page)
  hwrite(M, page=page, br=TRUE, style="text-align:center;")
  hwrite(metadata, page=page)
  closePage(page, splash=FALSE)
  file.copy(system.file("images","hwriter.css", package="hwriter"),file.path(htmldir, platename, "hwriter.css"))
  
  return(TRUE)   
}
