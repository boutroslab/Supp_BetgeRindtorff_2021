import os

# DATA PARAMETERS
# These are data-derived parameters that will vary between data types
SIZETHRESHOLD = 300

# CONFIGURATION PARAMETERS
IMAGESIZE = (2048, 2048)
NUMFIELDS = 4

# CUSTOM DIRECTORIES
# These directories point to the two base directories needed for this module,
# the directory for the analysis output and the directory for the data input
BASEDIR = "/omics/groups/OE0049/b110_data/B110-Isilon2/promise/data/interim/FeatureAnalysis"
PROMISEDIR = "/omics/groups/OE0049/b110_data/B110-Isilon2/promise/data/raw/PROMISE/data-10x-4t-c-16z"

# PROMISE
# These directories describe the PROMISE folder structure.
FEATUREDIR = os.path.join(PROMISEDIR, "features")
INCEPTIONDIR = os.path.join(PROMISEDIR, "features", "inception")
LAYOUTDIR = os.path.join(PROMISEDIR, "layouts", "python_friendly")
SEGMENTATIONDIR = os.path.join(PROMISEDIR, "segmentation")
IMAGEDIR = os.path.join(PROMISEDIR, "hdf5projection")

# RESULTS
# These directories (and files) describe the folder structure of the analysis
# results
PYTHONDIR = os.path.join(BASEDIR, "python")
SIZESTATISTICSDIR = os.path.join(BASEDIR, "size_statistics")
if not os.path.isdir(SIZESTATISTICSDIR):
    os.makedirs(SIZESTATISTICSDIR)
BLURRYORGANOIDCLSDIR = os.path.join(BASEDIR, "blurry_organoids")
if not os.path.isdir(BLURRYORGANOIDCLSDIR):
    os.makedirs(BLURRYORGANOIDCLSDIR)
DEADORGANOIDCLASSIFIERDIR = os.path.join(BASEDIR, "organoid_viability", "human")
if not os.path.isdir(DEADORGANOIDCLASSIFIERDIR):
    os.makedirs(DEADORGANOIDCLASSIFIERDIR)
LINEDIFFERENCESDIR = os.path.join(BASEDIR, "line_differences")
if not os.path.isdir(LINEDIFFERENCESDIR):
    os.makedirs(LINEDIFFERENCESDIR)
DRUGEFFECTSDIR = os.path.join(BASEDIR, "drug_effects")
if not os.path.isdir(DRUGEFFECTSDIR):
    os.makedirs(DRUGEFFECTSDIR)
CLUSTERDIR = os.path.join(BASEDIR, "cluster_scripts")
if not os.path.isdir(CLUSTERDIR):
    os.makedirs(CLUSTERDIR)
BATCHEFFECTDIR = os.path.join(BASEDIR, "batch_effects")
if not os.path.isdir(BATCHEFFECTDIR):
    os.makedirs(BATCHEFFECTDIR)

# LEGACY PARAMETERS
# This parameter comes from when I compared features for individual
# organoids and entire connected clumps. There was no notable difference in
# early analysis, so I abandoned the comparison. The code is still in place
# to make reimplementation easier. This means an HDF5 file containing features
# can contain 'features_organoids, metadata_organoids, etc.' as well as
# 'features_XYZ, metadata_XYZ, etc.'
# That means, do NOT change this unless you are absolutely certain of what
# you're doing.
FEATURETYPE = "organoids"
