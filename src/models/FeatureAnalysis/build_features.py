from OrganoidFeatures import OrganoidFeatures
import DrugEffects
import LineDifferences
import OrganoidViabilityClassifier
import CreateOrganoidCutouts
import Utils

lines = Utils.get_all_lines("human")
for line in lines:
  #CreateOrganoidCutouts.get_typical_organoids_for_dmso_detailed(line, 'data/processed/FeatureAnalysis', sizethresh=0, num_organoids=300)
  CreateOrganoidCutouts.get_small_organoids_for_dmso_detailed(line, 'data/processed/FeatureAnalysis', sizethresh=10, num_organoids=300)
