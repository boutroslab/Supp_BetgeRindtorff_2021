import Config
import h5py
import os
import numpy as np
import ImagePreprocessor
import pandas as pd


def load_projection(well, field_index=None):
    """
    Load a projection
    :param well: String. E.g. "D001T01P001L01_A_01"
    :param field_index: Integer. Starts at 0.
    :return:
    """
    plate, row, col = well.split("_")[0:3]
    img_fn = os.path.join(
        Config.IMAGEDIR, plate,
        "{}_contrastProjections.h5".format(well))
    if field_index is None:
        with h5py.File(img_fn, "r") as h5handle:
            proj = h5handle["images"][()]
    else:
        with h5py.File(img_fn, "r") as h5handle:
            proj = h5handle["images"][field_index, ...]
    return proj


def load_segmentation(well, field_index=None, probs=False):
    """
    Load a segmentation mask or probabilities
    :param well: String. E.g. "D001T01P001L01_A_01"
    :param field_index: Integer. Starts at 0.
    :param probs: Boolean. 'True' to load probability matrix,
        'False' to load binary mask. Probability matrix is
        not scaled to original image size but depends on the
        DNN segmentation network structure
    :return:
    """
    plate, row, col = well.split("_")[0:3]
    seg_fn = os.path.join(
        Config.SEGMENTATIONDIR, plate,
        "{}_DNNsegmentation.h5".format(well))
    hdf5key = "probs" if probs else "mask"
    if field_index is None:
        with h5py.File(seg_fn, "r") as h5handle:
            seg = h5handle[hdf5key][()]
    else:
        with h5py.File(seg_fn, "r") as h5handle:
            seg = h5handle[hdf5key][field_index, ...]
    seg[seg > 0] = 1

    return seg


def save_well_to_file(well, file, outsize=None):
    proj = [
        ImagePreprocessor.preprocess_image(img)
        for img in load_projection(well)]
    xdim = proj[0].shape[1]
    ydim = proj[0].shape[0]
    zdim = proj[0].shape[2]
    border = 10
    num_fields = len(proj)
    if num_fields == 4:
        mosaic = np.zeros(
            shape=(2 * xdim + border, 2 * ydim + border, zdim),
            dtype=proj[0].dtyle)
        mosaic[0:xdim, 0:ydim, :] = proj[0]
        mosaic[xdim + border, 0:ydim, :] = proj[1]
        mosaic[0:xdim, ydim + border, :] = proj[3]
        mosaic[xdim + border, ydim + border, :] = proj[2]


def get_organoid_image(metadata_entry, radius=None):
    """
    Get a cutout of an individual organoid.

    The radius of the image is estimated from the area as:
    (radius = sqrt(Size / pi) * 1.5)
    but can be overridden as segmentation errors and irregular
    shapes may make this incorrect.

    The pixel values can optionally be preprocessed by ImagePreprocessor
    :param metadata_entry: A pandas Series
    :param radius:
    :return:
    """

    if not isinstance(metadata_entry, pd.Series):
        raise ValueError("'metadata_entry' must be a pandas Series instance")

    well_id = "{}_{}_{}".format(
        metadata_entry.Plate, metadata_entry.Well[0],
        metadata_entry.Well[1:3])
    img = load_projection(well=well_id, field_index=metadata_entry.Field - 1)
    img = ImagePreprocessor.preprocess_image(img)

    y = int(np.round(metadata_entry.OriginalX))
    x = int(np.round(metadata_entry.OriginalY))
    if radius is None:
        radius = 1.5 * np.sqrt(metadata_entry.Size / np.pi)

    low_x = max(int(x - radius), 0)
    low_y = max(int(y - radius), 0)
    high_x = min(int(x + radius), img.shape[0])
    high_y = min(int(y + radius), img.shape[1])
    img_seg = img[low_x:high_x, low_y:high_y, :]

    return img_seg
