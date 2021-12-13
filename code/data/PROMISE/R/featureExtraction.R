# Note: most of the code in this file is outdated. It's in here because it 
# might be useful again for one reason or another. Be cautious with any function
# not explictly marked with '@export', they most likely don't even work anymore.

#' @title Extract organoid features
#' 
#' @description Calculates the features of a projection. It calculates three 
#' feature sets: 
#'     * A primitive segmentation of organoid clumps
#'     * An attempt to separate connected organoids via watershedding
#'     * The texture features of only the foreground / background
#' 
#' @param plateIndir This parameter does nothing, but makes it easier to integrate the function into the workflow
#' @param platename The plate that the well is on
#' @param row The row of the well
#' @param col The column of the well
#' @param configdir The directory of the configuration files
#' 
#' @return A boolean value indicating if the features were correctly calculated and saved
#' 
#' @author Jan Sauer
#' 
#' @examples print(extractOrganoidFeatures)
#' @export
extractOrganoidFeatures <- function(plateIndir, platename, row, col, configdir) {
  # This is a necessary manual import as there seems to be a bug with EBImage:
  # computeFeatures.moment isn't loaded into the namespace. Or something. Who 
  # knows. R is a terrible language.
  library(EBImage)
  source(file.path(configdir, "watchdogConfig.R"))
  # Load segmap
  segmap_fn = segmentationHdf5Filename(
    filedir = file.path(segmentationdir, platename), 
    platename = platename, row = row, col = col)
  segmap_mask = h5read(segmap_fn, "mask")
  
  # Load projection
  proj_fn = projectionHdf5Filename(
    filedir = file.path(hdf5projection, platename), 
    platename = platename, row = row, col = col)
  proj = h5read(proj_fn, "images")
  
  # Calculate features for every field
  num_fields = dim(proj)[4]
  features = vector(mode = "list", length = num_fields)
  # features_clumps = vector(mode = "list", length = num_fields)
  # features_noseg = vector(mode = "list", length = num_fields)
  for(fld in seq_len(num_fields)) {
    img = Image(data = proj[,,,fld], colormode = "color")
    # The edge detection doesn't really work for human organoids, so just
    # set "edges" and "foreground" to 1
    mask = matrix(
      as.integer(segmap_mask[,,fld] > 0), 
      nrow = nrow(segmap_mask), 
      ncol = ncol(segmap_mask))
    
    # Small erosion to make the watershedding of touching objects better
    mask_erode = erode(mask, kern = makeBrush(11, "disc"))

    # Find distance to background and perform primitive watershedding. The 
    # tolerance of 15 was chosen empirically and is bound to change as the 
    # segmentation gets better
    dmap = distmap(mask_erode, "euclidean")
    wshed = watershed(dmap, tolerance = 15)
    
    # Undo the erosion by voronoi propagating the watershed labels into the 
    # original mask
    mask_labeled = propagate(x = mask, seeds = wshed, mask = mask)
    
    # This is a hacky solution to avoid errors due to missing segmentation.
    # WARNING: If anybody changes the number of features to compute (haralick scales), 
    # then the size of the array must be modified.
    if(sum(mask_labeled != 0) == 0) {
      features[[fld]] = matrix(0, nrow = 0, ncol = 1572)
      # features_clumps[[fld]] = matrix(0, nrow = 0, ncol = 1572)
      # features_noseg[[fld]] = matrix(0, nrow = 0, ncol = 118)
      next
    }
    
    features[[fld]] = computeFeatures(
      x = mask_labeled, ref = normalize(img), 
      haralick.scales = c(1, 2, 4, 8, 16, 32, 64, 128, 256))
    features[[fld]] = cbind(features[[fld]], "FIELD" = fld)
    
    # features_clumps[[fld]] = computeFeatures(
    #   x = bwlabel(mask_labeled > 0), ref = normalize(img), 
    #   haralick.scales = c(1, 2, 4, 8, 16, 32, 64, 128, 256))
    # features_clumps[[fld]] = cbind(features_clumps[[fld]], "FIELD" = fld)
    # 
    # features_noseg[[fld]] = computeFeatures.haralick(
    #   x = mask, ref = normalize(img),
    #   haralick.scales = c(1, 2, 4, 8, 16, 32, 64, 128, 256))
    # features_noseg[[fld]] = cbind(features_noseg[[fld]], "FIELD" = fld)
  }
  
  features = do.call(what = rbind, args = features)
  # features_clumps = do.call(what = rbind, args = features_clumps)
  # features_noseg = do.call(what = rbind, args = features_noseg)
  feature_names = colnames(features)
  # feature_names_clumps = colnames(features_clumps)
  # feature_names_noseg = colnames(features_noseg)
  
  # Write features to file for well
  feature_dir = file.path(featuresdir, platename, "wells")
  dir.create(feature_dir, recursive = TRUE, showWarnings = FALSE)
  h5FileName = featuresHdf5Filename(
    filedir = feature_dir, platename = platename, 
    row = row, col = col, configdir = configdir)
  file_made = h5createFile(h5FileName)
  if(!file_made) {H5close(); return(FALSE)}
  
  dataset_made = h5createDataset(
    file = h5FileName, dataset = "features", dims = dim(features), 
    storage.mode = "double", level = hdf5_compression, showWarnings = FALSE)
  if(!dataset_made) {H5close(); return(FALSE)}
  # dataset_made = h5createDataset(
  #   file = h5FileName, dataset = "features_clumps", dims = dim(features_clumps), 
  #   storage.mode = "double", level = hdf5_compression, showWarnings = FALSE)
  # if(!dataset_made) {H5close(); return(FALSE)}
  # dataset_made = h5createDataset(
  #   file = h5FileName, dataset = "features_noseg", dims = dim(features_noseg), 
  #   storage.mode = "double", level = hdf5_compression, showWarnings = FALSE)
  # if(!dataset_made) {H5close(); return(FALSE)}
  h5write(obj = feature_names, file = h5FileName, name = "feature_names", level = hdf5_compression)
  # h5write(obj = feature_names_clumps, file = h5FileName, name = "feature_names_clumps", level = hdf5_compression)
  # h5write(obj = feature_names_noseg, file = h5FileName, name = "feature_names_noseg", level = hdf5_compression)
  h5write(obj = features, file = h5FileName, name = "features")
  # h5write(obj = features_clumps, file = h5FileName, name = "features_clumps")
  # h5write(obj = features_noseg, file = h5FileName, name = "features_noseg")
  H5close()
  return(TRUE)
}

#' @title Combine well features
#' 
#' @description Combines the features for each individual well of a plate. 
#' This will only combine features if the features for every well have been 
#' calculated. The output is written to a file.
#' 
#' @param platename The plate that the well is on
#' @param configdir The directory of the configuration files
#' 
#' @return A boolean value indicating if the features were correctly combined 
#' and written to disk
#' 
#' @author Jan Sauer
#' 
#' @examples print(combineWellFeatures)
#' @export
combineWellFeatures <- function(platename, configdir) {
  warning("Feature combination should be performed with the FeatureAnalysis 
          python package.")
  invisible(NULL)
  
  source(file.path(configdir, "watchdogConfig.R"))
  
  # Make sure all files are present
  feature_files = sort(list.files(file.path(featuresdir, platename)))
  all_wells = wells(nrWells = nrWells)
  expected_feature_files = character(nrWells)
  for(i in seq_len(nrow(all_wells))) {
    r = all_wells[i, "rows"]
    c = all_wells[i, "cols"]
    expected_feature_files[i] = featuresHdf5Filename(
      filedir = file.path(featuresdir, platename), platename = platename, 
      row = r, col = c, configdir = configdir)
  }
  if(!identical(feature_files, basename(expected_feature_files))) {
    warning(sprintf("Not all feature files present for '%s'", platename))
    return(FALSE)
  }
  
  features_organoids = setNames(
    object = vector(mode = "list", length = nrWells), 
    nm = all_wells$names)
  feature_names_organoids = setNames(
    object = vector(mode = "list", length = nrWells), 
    nm = all_wells$names)
  features_noseg = setNames(
    object = vector(mode = "list", length = nrWells), 
    nm = all_wells$names)
  feature_names_noseg = setNames(
    object = vector(mode = "list", length = nrWells), 
    nm = all_wells$names)
  features_clumps = setNames(
    object = vector(mode = "list", length = nrWells), 
    nm = all_wells$names)
  feature_names_clumps = setNames(
    object = vector(mode = "list", length = nrWells), 
    nm = all_wells$names)
  well_names_organoids = c()
  well_names_clumps = c()
  well_names_noseg = c()
  for(i in seq_len(nrow(all_wells))) {
    r = all_wells[i, "rows"]
    c = all_wells[i, "cols"]
    n = all_wells[i, "names"]
    full_path = featuresHdf5Filename(
      filedir = file.path(featuresdir, platename), platename = platename, 
      row = r, col = c, configdir = configdir)
    keys = h5ls(full_path)$name
    
    features_organoids[[n]] = data.frame(h5read(file = full_path, name = "features"))
    features_noseg[[n]] = data.frame(h5read(file = full_path, name = "features_noseg"))
    features_clumps[[n]] = data.frame(h5read(file = full_path, name = "features_clumps"))
    
    if("feature_names" %in% keys) {
      feature_names_organoids[[n]] = h5read(file = full_path, name = "feature_names")
      colnames(features_organoids[[n]]) = feature_names_organoids[[n]]}
    if("feature_names_noseg" %in% keys) {
      feature_names_noseg[[n]] = h5read(file = full_path, name = "feature_names_noseg")
      colnames(features_noseg[[n]]) = feature_names_noseg[[n]]}
    if("feature_names_clumps" %in% keys) {
      feature_names_clumps[[n]] = h5read(file = full_path, name = "feature_names_clumps")
      colnames(features_clumps[[n]]) = feature_names_clumps[[n]]}
    
    well_names_organoids = c(well_names_organoids, rep_len(n, nrow(features_organoids[[n]])))
    well_names_clumps = c(well_names_clumps, rep_len(n, nrow(features_clumps[[n]])))
    well_names_noseg = c(well_names_noseg, rep_len(n, nrow(features_noseg[[n]])))
    H5close()
  }
  
  # Test that all feature names are identical
  unique_fno = unique(feature_names_organoids)
  unique_fno = unique_fno[!sapply(unique_fno, is.null)]
  if(length(unique_fno) != 1) {
    warning(paste0(
      "Feature names for '", platename, 
      "' (organoids) are not identical across all wells"))
    return(FALSE)
  }
  unique_fnc = unique(feature_names_clumps)
  unique_fnc = unique_fnc[!sapply(unique_fnc, is.null)]
  if(length(unique_fnc) != 1) {
    warning(paste0(
      "Feature names for '", platename, 
      "' (clumps) are not identical across all wells"))
    return(FALSE)
  }
  unique_fnns = unique(feature_names_noseg)
  unique_fnns = unique_fnns[!sapply(unique_fnns, is.null)]
  if(length(unique_fnns) != 1) {
    warning(paste0(
      "Feature names for '", platename, 
      "' (foreground haralick) are not identical across all wells"))
    return(FALSE)
  }
  
  # Apply column names to features because empty feature sets have no column names
  for(well in names(features_organoids)) {
    colnames(features_organoids[[well]]) = as.character(unique_fno[[1]])
  }
  for(well in names(features_clumps)) {
    colnames(features_clumps[[well]]) = as.character(unique_fnc[[1]])
  }
  for(well in names(features_noseg)) {
    colnames(features_noseg[[well]]) = as.character(unique_fnns[[1]])
  }
  
  # Combine features (rbindlist is MUCH faster)
  features_noseg = as.data.frame(data.table::rbindlist(
    features_noseg, use.names = TRUE, fill = FALSE))
  features_organoids = as.data.frame(data.table::rbindlist(
    features_organoids, use.names = TRUE, fill = FALSE))
  features_clumps = as.data.frame(data.table::rbindlist(
    features_clumps, use.names = TRUE, fill = FALSE))
  
  # Sanity check for well names
  if(length(well_names_organoids) != nrow(features_organoids)) {
    warning(paste0(
      "The number of wells and the size of the feature matrix for '", 
      platename, "' don't match (organoids)"))
    return(FALSE)
  }
  if(length(well_names_clumps) != nrow(features_clumps)) {
    warning(paste0(
      "The number of wells and the size of the feature matrix for '", 
      platename, "' don't match (clumps)"))
    return(FALSE)
  }
  if(length(well_names_noseg) != nrow(features_noseg)) {
    warning(paste0(
      "The number of wells and the size of the feature matrix for '", 
      platename, "' don't match (foreground haralick)"))
    return(FALSE)
  }
  
  # Move all existing well features to a subfolder
  dir.create(path = file.path(
    featuresdir, platename, "wells"), recursive = FALSE, showWarnings = FALSE)
  for(fn in feature_files) {
    from = file.path(featuresdir, platename, fn)
    to = file.path(featuresdir, platename, "wells", fn)
    file.rename(from = from, to = to)
  }
  
  # Save plate features
  h5FileName = file.path(
    featuresdir, platename, sprintf("%s_features.h5", platename))
  file_made = h5createFile(h5FileName)
  if(!file_made) {H5close(); return(FALSE)}
  dataset_made = h5createDataset(
    file = h5FileName, dataset = "features", dims = dim(features_organoids), 
    storage.mode = "double", chunk = c(1, ncol(features_organoids)), 
    level = hdf5_compression)
  if(!dataset_made) {H5close(); return(FALSE)}
  dataset_made = h5createDataset(
    file = h5FileName, dataset = "features_clumps", dims = dim(features_clumps), 
    storage.mode = "double", chunk = c(1, ncol(features_clumps)), 
    level = hdf5_compression)
  if(!dataset_made) {H5close(); return(FALSE)}
  dataset_made = h5createDataset(
    file = h5FileName, dataset = "features_noseg", dims = dim(features_noseg), 
    storage.mode = "double", chunk = c(1, ncol(features_noseg)), 
    level = hdf5_compression)
  if(!dataset_made) {H5close(); return(FALSE)}
  # Feature names
  h5write(obj = colnames(features_organoids), file = h5FileName, 
          name = "feature_names", level = hdf5_compression)
  h5write(obj = colnames(features_clumps), file = h5FileName, 
          name = "feature_names_clumps", level = hdf5_compression)
  h5write(obj = colnames(features_noseg), file = h5FileName, 
          name = "feature_names_noseg", level = hdf5_compression)
  
  # Well names
  h5write(obj = well_names_organoids, file = h5FileName, 
          name = "well_names", level = hdf5_compression)
  h5write(obj = well_names_clumps, file = h5FileName, 
          name = "well_names_clumps", level = hdf5_compression)
  h5write(obj = well_names_noseg, file = h5FileName, 
          name = "well_names_noseg", level = hdf5_compression)
  
  # Features
  h5write(
    obj = as.matrix(features_organoids), 
    file = h5FileName, name = "features")
  h5write(
    obj = as.matrix(features_clumps), 
    file = h5FileName, name = "features_clumps")
  h5write(
    obj = as.matrix(features_noseg), 
    file = h5FileName, name = "features_noseg")
  H5close()
  return(TRUE)
}

