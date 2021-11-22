#' These functions perform an enrichment analysis for hierarchical clustering
#' based on pre-defined labels. The idea of this comes from
#'

#' @title Hierarchical Cluster Enrichment Analysis
#'
#' @description Find label enrichments in a hierarchical cluster object. This
#'     function is based on:
#'     Freudenberg et al., "CLEAN: CLustering Enrichment ANalysis",
#'     BMC Bioinformatics, 2009, 10:234
#'
#' @param clustering An hclust object.
#' @param labels A list of labels, each entry being a vector of integers
#'        indicating the IDs of the group members. The IDs should match
#'        'clustering$labels'
#' @param min_cluster_size The minimum cluster size to consider.
#' @param max_cluster_size THe maximum cluster size to consider. Defaults to
#'        1/4 of the total number of elements or the size of the largest label
#'         group, depending on what is larger.
#' @param min_odds_ratio A lower limit for the odds ratio for a significant
#'        cluster. This can be set to 0 if no limitation is desired.
#' @param max_annotated_labels Optionally filter the fisher test results and
#'        use only the N most significant labels when generating the annotation
#'        vector.
#' @param keep_only_best_cluster A boolean that determines if only the most
#'        significant cluster for each pathway is kept. Not implemented yet.
#' @return A list with two elements: a vector of significant annotations,
#'        corresponding to the clustered object IDs and the results of the
#'        fisher test
#'
#' @author Jan Sauer
#'
#' @usage get_cluster_enrichment()
#'
#' @examples print(get_cluster_enrichment)
#' @export
get_cluster_enrichment = function(
  clustering, labels, min_cluster_size, max_cluster_size=NULL,
  min_odds_ratio=0, max_annotated_labels=NULL,
  keep_only_best_cluster=FALSE) {

  # Remove NA labels
  labels = labels[!is.na(names(labels))]

  # Remove labels with fewer than min_cluster_size
  labels = labels[lengths(labels) >= min_cluster_size]

  # Set clusters
  clusters = list()
  for(ii in seq_len(nrow(clustering$merge))) {
    hcrow = clustering$merge[ii, ]
    if(hcrow[1] < 0) {
      set1 = -1 * hcrow[1]
    } else {
      set1 = clusters[[hcrow[1]]]
    }
    if(hcrow[2] < 0) {
      set2 = -1 * hcrow[2]
    } else {
      set2 = clusters[[hcrow[2]]]
    }
    clusters[[ii]] = c(set1, set2)
  }

  # Keep only clusters with at least N drugs
  clusters = clusters[sapply(clusters, length) >= min_cluster_size]

  # Keep only clusters with at most M drugs
  if(is.null(max_cluster_size)) {
    max_cluster_size = floor(max(lengths(clusters) / 4, max(lengths(labels))))}
  clusters = clusters[sapply(clusters, length) <= max_cluster_size]

  # Perform Fisher test for each cluster and each label represented in each
  # cluster
  ftest_results = expand.grid(
    "Cluster.ID" = seq_along(clusters),
    "Label" = names(labels),
    stringsAsFactors = FALSE)
  ftest_results$PValue = NA
  ftest_results$OddsRatio = NA
  for(ii in seq_along(clusters)) {
    cluster = clusters[[ii]]
    for(label in names(labels)) {
      # If fewer than 'min_cluster_size' entries in the cluster are
      # annotated with 'label' then don't process it.
      # This marks it for removal further below.
      if(length(intersect(labels[[label]], cluster)) < min_cluster_size)
        next
      fisher_matrix = matrix(c(
        # Label in Cluster
        length(intersect(labels[[label]], cluster)),
        # Label not in Cluster
        length(setdiff(labels[[label]], cluster)),
        # Other labels in Cluster
        length(setdiff(cluster, labels[[label]])),
        # Other labels not in Cluster
        length(setdiff(unlist(labels[
          -which(names(labels) == label)]), cluster))),
        nrow = 2)
      ftest = fisher.test(fisher_matrix, alternative = "greater")
      ftest_results[
        ftest_results$Cluster.ID == ii &
          ftest_results$Label == label,
        c("PValue", "OddsRatio")] = c(
          ftest$p.value, ftest$estimate)
    }
  }

  # Remove NA-entries, i.e. entries without any overlaps between label and
  # cluster
  ftest_results = ftest_results[!is.na(ftest_results$PValue), ]

  # For each label, remove superclusters in the hierarchy with no additional
  # entries for that label ("true superclusters")
  entries_to_keep = rep(TRUE, nrow(ftest_results))
  for(label in unique(ftest_results$Label)) {
    label_results = ftest_results[ftest_results$Label == label, ]
    if(nrow(label_results) == 1) next
    for(ii in seq(1, nrow(label_results), 1)) {
      cluster_1_ID = label_results[ii, "Cluster.ID"]
      cluster_1 = clusters[[cluster_1_ID]]
      cluster_1_labeled = intersect(labels[[label]], cluster_1)
      for(jj in seq(1, nrow(label_results), 1)) {
        if(ii == jj) next
        cluster_2_ID = label_results[jj, "Cluster.ID"]
        cluster_2 = clusters[[cluster_2_ID]]
        cluster_2_labeled = intersect(labels[[label]], cluster_2)
        # Cluster 2 must be a super cluster, i.e. contain all elements of
        # cluster 1, to be flagged for removal
        if(length(setdiff(cluster_1, cluster_2)) != 0) next
        # Cluster 2's elements annotated with 'label' must be entirely
        # inside cluster 1's elements annotated with 'label' for cluster 2
        # to be flagged for removal.
        if(length(setdiff(cluster_2_labeled, cluster_1_labeled)) != 0) next
        # If neither of the previous conditions were met and the loop reaches
        # this point, cluster 2 is a 'true supercluster' with regards to
        # 'label' and the corresponding entry in ftest_results should be removed
        flag_index = which(
          ftest_results$Cluster.ID == cluster_2_ID &
          ftest_results$Label == label)
        # This is a sanity check to ensure that the supercluster removal only
        # flags each true supercluster once
        if(!entries_to_keep[flag_index]) stop(
          "Critical Error while removing super clusters")
        entries_to_keep[flag_index] = FALSE
        # I only want to flag the immediate supercluster to avoid repeats
        break
      }
    }
  }
  ftest_results = ftest_results[entries_to_keep, ]

  # Perform multiple testing correction
  ftest_results$PAdj = p.adjust(ftest_results$PValue, "BH")
  ftest_results = ftest_results[ftest_results$PAdj <= 0.1, ]

  # For each label, remove superclusters in the hierarchy with a higher
  # adjusted p-value
  ftest_results = ftest_results[order(ftest_results$PAdj), ]
  entries_to_keep = rep(TRUE, nrow(ftest_results))
  for(label in unique(ftest_results$Label)) {
    label_results = ftest_results[ftest_results$Label == label, ]
    if(nrow(label_results) == 1) next
    for(ii in seq(1, nrow(label_results)-1, 1)) {
      cluster_1_ID = label_results[ii, "Cluster.ID"]
      cluster_1 = clusters[[cluster_1_ID]]
      for(jj in seq(ii+1, nrow(label_results), 1)) {
        cluster_2_ID = label_results[jj, "Cluster.ID"]
        cluster_2 = clusters[[cluster_2_ID]]
        # If cluster 2 is a supercluster of cluster 1 then, due to the sorting
        # above, it will have a higher p-value and should be removed
        if(length(setdiff(cluster_1, cluster_2)) == 0) {
          entries_to_keep[which(
            ftest_results$Cluster.ID == cluster_2_ID &
              ftest_results$Label == label)] = FALSE
        }
        # if cluster 2 is a subcluster of cluster 1 then, due to the sorting
        # above, it will have a higher p-value and can be removed as it contains
        # only redundant information contained in cluster 1.
        if(length(setdiff(cluster_2, cluster_1)) == 0) {
          entries_to_keep[which(
            ftest_results$Cluster.ID == cluster_2_ID &
              ftest_results$Label == label)] = FALSE
        }
      }
    }
  }
  ftest_results = ftest_results[entries_to_keep, ]

  # Remove insufficient odds ratios
  ftest_results = ftest_results[ftest_results$OddsRatio >= min_odds_ratio, ]

  # Generate a new labels list with only the significantly enriched labels
  # If an object is enriched for more than one label, keep the one with the
  # min adjusted pvalue
  num_objects = -min(clustering$merge)
  new_labels = rep(NA, num_objects)
  label_pvals = rep(1, num_objects)
  for(ii in seq_len(nrow(ftest_results))) {
    label = ftest_results[ii, "Label"]
    cluster_id = ftest_results[ii, "Cluster.ID"]
    pval = ftest_results[ii, "PAdj"]
    cluster_entries = intersect(clusters[[cluster_id]], labels[[label]])
    old_pvals = label_pvals[cluster_entries]
    new_labels[cluster_entries[old_pvals > pval]] = label
    label_pvals[cluster_entries[old_pvals > pval]] = pval
  }

  # If any of the annotated labels in new_labels are present less than
  # min_cluster_size times, due to duplicate annotation, then remove
  # these from the results altogether
  # to_remove = names(which(table(new_labels) < min_cluster_size))
  # ftest_results = ftest_results[!ftest_results$Label %in% to_remove, ]
  # new_labels[new_labels %in% to_remove] = NA

  # # This version generates compound labels if more than one target is
  # # significantly enriched for a given entry.
  # num_objects = -min(clustering$merge)
  # new_labels = rep(NA, num_objects)
  # for(ii in seq_len(nrow(ftest_results))) {
  #   label = ftest_results[ii, "Label"]
  #   cluster_id = ftest_results[ii, "Cluster.ID"]
  #   cluster_entries = intersect(clusters[[cluster_id]], labels[[label]])
  #   for(ce in cluster_entries) {
  #     new_labels[ce] = ifelse(
  #       test = is.na(new_labels[ce]),
  #       yes = label,
  #       no = paste0(new_labels[ce], ", ", label))
  #   }
  # }

  return(list("Labels" = new_labels, "Fisher.Test" = ftest_results))
}
