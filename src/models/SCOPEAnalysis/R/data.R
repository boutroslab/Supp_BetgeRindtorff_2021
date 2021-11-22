#' Viability Classifier Accuracies
#'
#' The full accuracy matrix of classifiers for individual organoid lines (rows)
#' applied to features of all other lines (columns)
#'
#' @format A matrix with 15 rows and columns
#' \describe{
#'   \item{rows}{Viability classifiers trained on each organoid line}
#'   \item{columns}{Validation data for each organoid line}
#' }
"acc_matrix"

#' Viability Classifier Accuracies for Two-Channel Classifiers
#'
#' A list of the full accuracy matrices of classifiers for individual organoid
#' lines (rows) applied to features of all other lines (columns). The list has
#' four entries, one for each reduced classifier: only actin, only cell nuclei,
#' only dead cells, or actin + cell nuclei.
#'
#' @format A list with 4 entries. Each entry is a 15 x 15 matrix.
"acc_matrix_reduced"

#' CellTiter-Glo Data
#'
#' CellTiter-Glo viability scores and corresponding dose-response AUCs
#'
#' @format Data frame with 12288 rows and 29 columns
#' \describe{
#'   \item{screen}{}
#'   \item{date}{}
#'   \item{operator}{}
#'   \item{mithras}{}
#'   \item{full_barcode}{}
#'   \item{Line}{}
#'   \item{plate}{}
#'   \item{library}{}
#'   \item{tag1}{}
#'   \item{tag2}{}
#'   \item{well}{}
#'   \item{row}{}
#'   \item{col}{}
#'   \item{pcount}{}
#'   \item{batch}{}
#'   \item{drug}{}
#'   \item{substance}{}
#'   \item{Concentration}{}
#'   \item{rack.position}{}
#'   \item{volume}{}
#'   \item{norm_fac}{}
#'   \item{pcount_norm}{}
#'   \item{med_ctrl}{}
#'   \item{viability}{}
#'   \item{rep}{}
#'   \item{s_AUC_fit}{}
#'   \item{s_AUC_actual}{}
#'   \item{l_AUC_fit}{}
#'   \item{l_AUC_actual}{}
#' }
"auc_ctg"

#' Dose-Response Data for Viability Classifier
#'
#' Dose-Response AUCs for Clinical Cancer Panel subset of image-based
#' viabilities.
#'
#' @format Data frame with 10464 rows and 16 columns
#' \describe{
#'   \item{Product.Name}{}
#'   \item{Concentration}{}
#'   \item{Plate.ID}{}
#'   \item{Well.ID}{}
#'   \item{Num.Objects}{}
#'   \item{Percent.Dead}{Percent of organoids that are dead in the well}
#'   \item{Percent.Live}{Percent of organoids that are alive in the well}
#'   \item{Mean.Certainty.Live}{The mean classification accuracy of all organoids classified as alive}
#'   \item{Mean.Certainty.Dead}{The mean classification accuracy of all organoids classified as dead}
#'   \item{Line}{}
#'   \item{Replicate}{}
#'   \item{Layout}{}
#'   \item{s_AUC_fit}{}
#'   \item{s_AUC_actual}{}
#'   \item{l_AUC_fit}{}
#'   \item{l_AUC_actual}{}
#' }
"auc_img_full"

#' Dose-Response Data for Viability Classifier
#'
#' Dose-Response AUCs for Clinical Cancer Panel subset of image-based
#' viabilities as determined by the two-channel (actin + cell nuclei)
#' classifier.
#'
#' @format Data frame with 10464 rows and 16 columns
#' \describe{
#'   \item{Product.Name}{}
#'   \item{Concentration}{}
#'   \item{Plate.ID}{}
#'   \item{Well.ID}{}
#'   \item{Num.Objects}{}
#'   \item{Percent.Dead}{Percent of organoids that are dead in the well}
#'   \item{Percent.Live}{Percent of organoids that are alive in the well}
#'   \item{Mean.Certainty.Live}{The mean classification accuracy of all organoids classified as alive}
#'   \item{Mean.Certainty.Dead}{The mean classification accuracy of all organoids classified as dead}
#'   \item{Line}{}
#'   \item{Replicate}{}
#'   \item{Layout}{}
#'   \item{s_AUC_fit}{}
#'   \item{s_AUC_actual}{}
#'   \item{l_AUC_fit}{}
#'   \item{l_AUC_actual}{}
#' }
"auc_img_twoChannel"

#' Drug Target Annotations
#'
#' Pathways and targets for all drugs in the screen.
#'
#' @format Data frame with 528 rows and 3 columns
#' \describe{
#'   \item{Drug}{}
#'   \item{Pathways}{}
#'   \item{Targets}{}
#' }
"drug_annotations"

#' Viabilities of wells
#'
#' Classification results for wells.
#'
#' @format Data frame with 32866 rows and 12 columns
#' \describe{
#'   \item{Product.Name}{}
#'   \item{Concentration}{}
#'   \item{Plate.ID}{}
#'   \item{Well.ID}{}
#'   \item{Num.Objects}{}
#'   \item{Percent.Dead}{Percent of organoids that are dead in the well}
#'   \item{Percent.Live}{Percent of organoids that are alive in the well}
#'   \item{Mean.Certainty.Live}{The mean classification accuracy of all organoids classified as alive}
#'   \item{Mean.Certainty.Dead}{The mean classification accuracy of all organoids classified as dead}
#'   \item{Line}{}
#'   \item{Replicate}{}
#'   \item{Layout}{}
#' }
"mortality"

#' Viabilities of wells
#'
#' A list of classification results for wells for all reduced-channel
#' classifiers. The list has four entries, one for each reduced classifier:
#' only actin, only cell nuclei, only dead cells, or actin + cell nuclei.
#'
#' @format A list with 4 entries. Each entry is a data frame with 32866 rows
#' and 12 columns.
#' \describe{
#'   \item{Product.Name}{}
#'   \item{Concentration}{}
#'   \item{Plate.ID}{}
#'   \item{Well.ID}{}
#'   \item{Num.Objects}{}
#'   \item{Percent.Dead}{Percent of organoids that are dead in the well}
#'   \item{Percent.Live}{Percent of organoids that are alive in the well}
#'   \item{Mean.Certainty.Live}{The mean classification accuracy of all organoids classified as alive}
#'   \item{Mean.Certainty.Dead}{The mean classification accuracy of all organoids classified as dead}
#'   \item{Line}{}
#'   \item{Replicate}{}
#'   \item{Layout}{}
#' }
"mortality_reduced"

#' ROC curves for viability classifiers
#'
#' ROC curves for viability classifiers for all combinations of classifiers
#' applied to data of other lines.
#'
#' @format Data frame with 5 columns
#' \describe{
#'   \item{Threshold}{The classification threshold being tested. Sometimes, a
#'   threshold of 2.0 is selected by python. This is only a technical measure
#'   to ensure that the ROC curve goes from (0, 0) to (1, 1).}
#'   \item{FalsePosRate}{False positive rate at given threshold}
#'   \item{TruePosRate}{True positive rate at given threshold}
#'   \item{clf_line}{The source line of the classifier being tested}
#'   \item{data_line}{The source line for the data being tested on}
#' }
"roc_data"

#' ROC curve AUCs
#'
#' AUC values for ROC curves found in 'roc_data'
#'
#' @format Data frame with 3 columns
#' \describe{
#'   \item{clf_line}{The source line of the classifier being tested}
#'   \item{data_line}{The source line for the data being tested on}
#'   \item{AUC}{Area under receiver operating curve}
#' }
"roc_aucs"

#' Drug-induced phenotype profiles
#'
#' Drug-induced phenotype vector as calculated by the SVMs. Each row in the
#' matrix corresponds to one treatment. A treatment consists of a line, drug,
#' and concentration (when applicable). The treatment ID can be read from the
#' row names with the format 'LINE.DRUG__CONCENTRATION' (for Clinical Cancer
#' Panel drugs) or 'LINE.DRUG' (for KiStem drugs). 'drug_effect_metadata' also
#' has this information in a more easily manageable table format.
#'
#' Profiles are averaged over 10 cross validation iterations.
#'
#' @format Matrix with 11620 rows and 25 columns
"drug_effect_profiles"

#' Drug-induced phenotype profile metadata
#'
#' Metadata belonging to 'drug_effect_profiles'. Row X in the metadata
#' corresponds to row X in the profiles (rownames are identical)
#'
#' @format Data frame with 11620 rows and 6 columns
#' \describe{
#'   \item{AUC_Mean}{The mean SVM AUROC, i.e. classification performance,
#'   averaged over 10 cross validation iterations.}
#'   \item{AUC_Std}{The standard deviation of the SVM AUROC over 10 cross
#'   validation iterations.}
#'   \item{Distance}{The euclidean distance between the median features of
#'   negative controls and treated drug (after PCA-transformation).}
#'   \item{Line}{}
#'   \item{Drug}{}
#'   \item{Concentration}{}
#' }
"drug_effect_metadata"

#' Well-averaged features
#'
#' Well-averages of processed single-organoid features. Note that this still
#' contains all features, even those that were too constant, i.e. resulted
#' in a value of NA/NaN after preprocessing. Row names indicate plate and
#' well ID. Note that some wells are missing, either because data was corrupted
#' (D013T01) or because imaging quality was too low for correct segmentation
#' (D010T01, D046T01, D054T01, and D055T01). Also note that D054T01 was only
#' treated with the clinical cancer panel library (L08).
#'
#' @format Matrix with 32939 rows and 3143 columns
"well_features"

#' Well-averaged features metadata
#'
#' Metadata belonging to 'well_features'.
#'
#' @format Matrix with 32939 rows and 3143 columns
"well_metadata"

