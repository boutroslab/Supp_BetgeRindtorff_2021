from OrganoidFeatures import OrganoidFeatures
import DrugEffects
import LineDifferences
import Utils
import Config
import CreateOrganoidCutouts

# Calculate drug-induced phenotypes

# list all lines that are currently included within the analysis
lines = Utils.get_all_lines("human")
print(lines)

# list all plates that are currently included
plates = Utils.get_all_plates("human")
print(plates)

# run incremental PCA across all included features
#LineDifferences.process_dmso_organoids_no_filter("human")
LineDifferences.process_all_organoids("human")

