"""
This contains global parameters for the package that should be defined in a
single place. Directories are defined as variables, but everything else
should be a function for more flexibility.

Anything that might change from project to project or run to run should be
placed here. Note that the function names should stay the same or be
thoroughly refactored in the codebase.
"""
import os
import numpy as np

# === DIRECTORIES === #
# Describes the directory structure and the base folder for input images
# =================== #
# Image Input
INPUT_DIR = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/hdf5projection"
# Base directory for all results
BASE_DIR = "/collab-ag-fischer/Segmentation"
# Subdirectories
INTENSITY_SEGMENTATION = os.path.join(BASE_DIR, "intensity_segmentation")
DNN_SEGMENTATION = os.path.join(BASE_DIR, "dnn_segmentation")
DNN_FEATURES = os.path.join(BASE_DIR, "dnn_features")
PREVIEW_DIR = os.path.join(BASE_DIR, "previews")

# === FILENAMES === #
# Various filenames kept here for consistency
# ================= #
DNN_TRAINING_SEGMENTS = os.path.join(
    DNN_SEGMENTATION, "DNN_TrainingSegments.h5")
DNN_WEIGHTS = os.path.join(DNN_SEGMENTATION, "DNN_weights.h5")
DNN_ACCURACY = os.path.join(DNN_SEGMENTATION, "DNN_accuracy.csv")
LOGFILE = os.path.join(BASE_DIR, "log.txt")


# === CONFIGURATION === #
# Various global configuration parameters
# ===================== #
def get_image_names():
    """
    Generates the names for input image names.

    If your project contains arbitrarily named images then this function can
    just return a hard-coded tuple of those.

    The images should be returned in the format "[image_dir]/[image_name]",
    i.e. it should be the relative path starting from Config.INPUT_DIR

    :return: A tuple of strings
    """

    human_directories = (
        "D004T01P003L02", "D004T01P004L03", "D004T01P005L08",
        "D004T01P007L02", "D004T01P008L03", "D004T01P009L08",
        "D007T01P003L02", "D007T01P004L03", "D007T01P005L08",
        "D007T01P007L02", "D007T01P008L03", "D007T01P009L08",
        "D010T01P017L02", "D010T01P018L03", "D010T01P019L08",
        "D010T01P021L02", "D010T01P022L03", "D010T01P023L08",
        "D013T01P001L02", "D013T01P002L03", "D013T01P003L08",
        "D013T01P905L02", "D013T01P906L03", "D013T01P907L08",
        "D018T01P901L02", "D018T01P902L03", "D018T01P003L08",
        "D018T01P905L02", "D018T01P906L03", "D018T01P907L08",
        "D019T01P001L02", "D019T01P002L03", "D019T01P003L08",
        "D019T01P005L02", "D019T01P006L03", "D019T01P007L08",
        "D020T01P001L02", "D020T01P002L03", "D020T01P003L08",
        "D020T01P905L02", "D020T01P906L03", "D020T01P907L08",
        "D020T02P009L02", "D020T02P010L03", "D020T02P011L08",
        "D020T02P013L02", "D020T02P014L03", "D020T02P015L08",
        "D021T01P901L02", "D021T01P902L03", "D021T01P003L08",
        "D021T01P905L02", "D021T01P906L03", "D021T01P907L08",
        "D022T01P001L02", "D022T01P002L03", "D022T01P003L08",
        "D022T01P005L02", "D022T01P906L03", "D022T01P907L08",
        "D027T01P001L02", "D027T01P002L03", "D027T01P003L08",
        "D027T01P905L02", "D027T01P906L03", "D027T01P907L08",
        "D030T01P001L02", "D030T01P002L03", "D030T01P003L08",
        "D030T01P905L02", "D030T01P906L03", "D030T01P907L08",
        "D046T01P001L02", "D046T01P002L03", "D046T01P003L08",
        "D046T01P005L02", "D046T01P006L03", "D046T01P007L08",
        "D054T01P004L08", "D054T01P006L08",
        "D055T01P007L02", "D055T01P008L03", "D055T01P009L08",
        "D055T01P011L02", "D055T01P012L03", "D055T01P013L08")
    mouse_directories = (
        "M001A03P003L02", "M001A03P004L03", "M001A03P005L04",
        "M001A03P006L05", "M001A03P007L06", "M001A03P008L07",
        "M001A03P009L02", "M001A03P010L03", "M001A03P011L04",
        "M001A03P012L05", "M001A03P013L06", "M001A03P014L07",
        "M001B04P003L02", "M001B04P004L03", "M001B04P005L04",
        "M001B04P006L05", "M001B04P007L06", "M001B04P008L07",
        "M001B04P009L02", "M001B04P010L03", "M001B04P011L04",
        "M001B04P012L05", "M001B04P013L06", "M001B04P014L07",
        "M001K02P003L02", "M001K02P004L03", "M001K02P005L04",
        "M001K02P006L05", "M001K02P007L06", "M001K02P008L07",
        "M001K02P009L02", "M001K02P010L03", "M001K02P011L04",
        "M001K02P012L05", "M001K02P013L06", "M001K02P014L07",
        "M001W01P003L02", "M001W01P004L03", "M001W01P005L04",
        "M001W01P006L05", "M001W01P007L06", "M001W01P908L07",
        "M001W01P009L02", "M001W01P010L03", "M001W01P011L04",
        "M001W01P012L05", "M001W01P013L06", "M001W01P014L07")

    rows = (
        "A", "B", "C", "D", "E", "F", "G", "H",
        "I", "J", "K", "L", "M", "N", "O", "P")
    cols = (
        "01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12",
        "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24")

    img_dirs = human_directories + mouse_directories
    img_names = []
    for img_dir in img_dirs:
        for row in rows:
            for col in cols:
                img_names.append(os.path.join(
                    img_dir, img_dir + "_" + row + "_" +
                    col + "_contrastProjections.h5"))

    # HACK: These are corrupted wells and they'll lead to inevitable errors
    missing_wells = np.array((
        "D013T01P001L02/D013T01P001L02_A_01_contrastProjections.h5",
        "D013T01P001L02/D013T01P001L02_A_02_contrastProjections.h5",
        "D013T01P001L02/D013T01P001L02_A_03_contrastProjections.h5",
        "D013T01P001L02/D013T01P001L02_A_06_contrastProjections.h5",
        "D013T01P001L02/D013T01P001L02_A_07_contrastProjections.h5",
        "D013T01P001L02/D013T01P001L02_A_10_contrastProjections.h5",
        "D013T01P001L02/D013T01P001L02_A_12_contrastProjections.h5",
        "D013T01P001L02/D013T01P001L02_A_15_contrastProjections.h5",
        "D013T01P001L02/D013T01P001L02_A_18_contrastProjections.h5",
        "D013T01P001L02/D013T01P001L02_A_21_contrastProjections.h5",
        "D013T01P001L02/D013T01P001L02_A_22_contrastProjections.h5",
        "D013T01P001L02/D013T01P001L02_A_23_contrastProjections.h5",
        "D013T01P001L02/D013T01P001L02_A_24_contrastProjections.h5",
        "D013T01P001L02/D013T01P001L02_B_04_contrastProjections.h5",
        "D013T01P001L02/D013T01P001L02_B_05_contrastProjections.h5",
        "D013T01P001L02/D013T01P001L02_B_07_contrastProjections.h5",
        "D013T01P001L02/D013T01P001L02_B_14_contrastProjections.h5",
        "D013T01P001L02/D013T01P001L02_I_11_contrastProjections.h5",
        "D013T01P001L02/D013T01P001L02_J_01_contrastProjections.h5",
        "D013T01P001L02/D013T01P001L02_K_09_contrastProjections.h5",
        "D013T01P001L02/D013T01P001L02_M_10_contrastProjections.h5",
        "D013T01P001L02/D013T01P001L02_M_17_contrastProjections.h5",
        "D013T01P001L02/D013T01P001L02_N_05_contrastProjections.h5"))
    img_names = list(np.array(img_names)[~np.in1d(img_names, missing_wells)])

    # Generate a (predictably) random subset of the images
    r = np.random.RandomState(1001)
    sel_indices = r.choice(range(len(img_names)), size=10000, replace=False)
    img_names = np.array(img_names)[sel_indices]

    return list(img_names)
    # return ["D026T11P004L11/D026T11P004L11_A_01_contrastProjections.h5"]
    # return ["D010T01P019L08/D010T01P019L08_D_04_contrastProjections.h5"]


def get_dnn_parameters():
    """
    Parameters for the DNN training
    :return: Dictionary
    """

    params = {
        # Total number of segments per category
        "num_segments": 1e6,
        # 'Radius' of square segments. Their total size is
        # (seg_radius * 2 + 1, seg_radius * 2 + 1)
        "seg_radius": (10, 10),
        # The number of images to simultaneously segment. Note that this does
        # not consider the number of fields per image. So if there are 4 fields
        # per image, the DNN will actually process n * 4 images simultaneously!
        "segmentation_batch_size": 1,
        # Training parameters
        "training_epochs": 10,
        "training_batch_size": 32
    }

    return params


def get_image_size():
    """
    Returns the image dimensions (x, y, channels). This COULD be read from the
    data, but is useful to have without reading data.
    :return:
    """
    return 2048, 2048, NUM_CHANNELS


# === EXPERIMENTAL PARAMETERS === #
# These are here primarily so that the code contains no 'magic numbers'.
# They might be customizable but that is probably unnecessary and could
# break the code
# =============================== #
NUM_CHANNELS = 2
FOREGROUND = 1
BACKGROUND = 0
NONCATEGORY = -1
