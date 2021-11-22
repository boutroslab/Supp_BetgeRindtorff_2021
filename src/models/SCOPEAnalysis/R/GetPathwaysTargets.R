
#' @title Get Mode of Action / Pathway
#'
#' @description Loads the annotation data and retrieves the pathway on which 
#'     drugs act.
#'
#' @usage get_mode_of_action()
#'
#' @return A list of pathways
#'
#' @author Jan Sauer
#'
#' @examples print(get_mode_of_action)
#' @export
get_mode_of_action = function(drugs) {
  data("Layouts", package = "SCOPEAnalysis")
  # This package deals exclusively with HUMAN DATA
  Layouts = Layouts[Layouts$Library_ID %in% c(2, 3, 8), ]
  Layouts = Layouts[, c("Product.Name", "Pathway")]
  Layouts$Product.Name = as.character(Layouts$Product.Name)
  Layouts$Pathway = as.character(Layouts$Pathway)
  Layouts = unique(Layouts)
  
  # Some entries might be duplicated. For those that are, remove any without pathway
  # information
  dup_drugs = names(which(table(Layouts$Product.Name) == 2))
  remove_entries = which(
    is.na(Layouts$Pathway) &
      Layouts$Product.Name %in% dup_drugs)
  if(length(remove_entries) > 0) Layouts = Layouts[-remove_entries, ]
  rownames(Layouts) = Layouts$Product.Name
  return(Layouts[drugs, "Pathway"])
}

#' @title Get Drug Targets
#'
#' @description Loads the annotation data and retrieves the targets on which 
#'     drugs act. There is some preprocessing involved here as a drug can 
#'     have multiple targets. In this case, a comma-separated list is returned.
#'
#' @usage get_targets()
#'
#' @return A list of drug targets
#'
#' @author Jan Sauer
#'
#' @examples print(get_targets)
#' @export
get_targets = function(drugs) {
  data("Layouts", package = "SCOPEAnalysis")
  # This package deals exclusively with HUMAN DATA
  Layouts = Layouts[Layouts$Library_ID %in% c(2, 3, 8), ]
  Layouts = Layouts[, c("Product.Name", "Target")]
  Layouts$Product.Name = as.character(Layouts$Product.Name)
  Layouts$Target = as.character(Layouts$Target)
  Layouts = unique(Layouts)
  
  # Some entries might be duplicated. For those that are, remove any without pathway
  # information
  dup_drugs = names(which(table(Layouts$Product.Name) == 2))
  remove_entries = which(
    is.na(Layouts$Target) &
      Layouts$Product.Name %in% dup_drugs)
  if(length(remove_entries) > 0) Layouts = Layouts[-remove_entries, ]
  
  # Merge the Targets of duplicated entries
  merged_layouts = list()
  for(prod in unique(Layouts$Product.Name)) {
    if(sum(Layouts$Product.Name == prod) == 1) {
      merged_layouts[[prod]] = as.character(Layouts[Layouts$Product.Name == prod, ])
    } else {
      full_targets = paste0(Layouts[Layouts$Product.Name == prod, "Target"], collapse = ", ")
      merged_layouts[[prod]] = c(prod, full_targets)
    }
  }
  merged_layouts = data.frame(do.call(rbind, merged_layouts), stringsAsFactors = FALSE)
  colnames(merged_layouts) = c("Product.Name", "Target")
  rownames(merged_layouts) = merged_layouts$Product.Name
  merged_layouts$Product.Name = NULL
  return(merged_layouts[drugs, "Target"])
}
