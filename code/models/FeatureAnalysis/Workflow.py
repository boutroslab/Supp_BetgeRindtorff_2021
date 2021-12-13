from OrganoidFeatures import OrganoidFeatures
import DrugEffects
import LineDifferences
import OrganoidViabilityClassifier
import Utils

# Features are saved one file per well. This function call collates them into
# a single file for quicker and easier loading. The 'plates' argument must
# always be an iterable, even if only one plate is passed.
OrganoidFeatures.combine_wells(
    plates=["D004T01P003L02", "D004T01P004L03", "D004T01P005L08",
            "D004T01P007L02", "D004T01P008L03", "D004T01P009L08"],
    verbose=True)

# Preprocess features. This consists of size outlier removal, removing
# organoids at the edges of images, blurry organoid removal, and plate
# effect correction. The combined feature file for each plate of a line
# must be present. The relevant plates are read in from
# Utils.py::get_rep_ids(). A normalized feature file is saved in the same
# folder as the raw feature file. Likewise, a "*_stageOne.h5" file is saved.
# These are the features after organoid quality control but before plate
# effect removal.
# IMPORTANT: It is HIGHLY recommended to preprocess all features for all lines
# in the screen before performing downstream analysis. Some of the downstream
# functions rely on the features for all lines to be present, in particular
# viability analysis comparisons and line clustering. Lines can be excluded or
# included in the analysis by editing Utils.py::get_rep_ids().
lines = Utils.get_all_lines("human")
for line in lines:
    OrganoidFeatures.preprocess_features(line=line)

# Calculate organoid viabilities and various diagnostics
# A separate classifier is trained for various channel combinations, i.e.
# only Actin, only DAPI, only FITC, and Actin + DAPI. The training is
# quick enough to make this cheap.
# Output is saved into Config.DEADORGANOIDCLASSIFIERDIR
OrganoidViabilityClassifier.run_on_all_lines()

# Calculate line differences.
# This function generates PCA-transformed features for negative controls
# of all lines that can be further analyzed and visualized as desired, e.g.
# hierarchical clustering.
# Output is saved into Config.LINEDIFFERENCESDIR
LineDifferences.process_dmso_organoids("human")

# Calculate drug-induced phenotypes
DrugEffects.train("human", 25)
