# This is the main module that calls the segmentation pipeline
import Config
import IntensitySegmenter
import DeepLearning
import os


def workflow():
    """
    This is the main segmentation workflow. It's a "push-and-forget" function
    that ultimately results in a set of DNN segmented images.
    :return:
    """

    # Set up results folder dir
    if not os.path.isdir(Config.BASE_DIR):
        os.makedirs(Config.BASE_DIR)
    if not os.path.isdir(Config.INTENSITY_SEGMENTATION):
        os.makedirs(Config.INTENSITY_SEGMENTATION)
    if not os.path.isdir(Config.DNN_SEGMENTATION):
        os.makedirs(Config.DNN_SEGMENTATION)
    if not os.path.isdir(Config.DNN_FEATURES):
        os.makedirs(Config.DNN_FEATURES)
    if not os.path.isdir(Config.PREVIEW_DIR):
        os.makedirs(Config.PREVIEW_DIR)

    IntensitySegmenter.run_intensity_segmentation(logfile=Config.LOGFILE)
    DeepLearning.run_training_segment_creation(logfile=Config.LOGFILE)
    DeepLearning.run_dnn_training(logfile=Config.LOGFILE)
    DeepLearning.run_dnn_segmentation(logfile=Config.LOGFILE)

