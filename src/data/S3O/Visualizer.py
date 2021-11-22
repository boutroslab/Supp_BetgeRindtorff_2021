"""
This module visualizes images and segmentations
"""
import os
import numpy as np
import scipy.misc
import scipy.ndimage
import ImageIO
import ImagePreprocessor
import Config
import Utils


def generate_preview_image(well, zoom=(0.125, 0.125), fmt="png"):
    """
    Generate a composite image that shows, from left to right, the original
    image, the intensity segmentation, and the DNN segmentation.

    This function is by no means "scientific", i.e. it treats pixels that
    were classified as neither foreground or background (non-category)
    as background. It's heavily biased to show off how awesome the DNN
    segmentation is.

    The scaling is also potentially buggy. If there are some issues with
    the preview images, check the actual data before panicking.

    :param well: An image file name in the format returned by
                 Config.get_image_names()
    :param zoom: The factor by which to scale the image
    :param fmt: The image format. This is simply appended to the filename
                and scipy.misc.imsave attempts to identify the correct
                format, i.e. this should be interpretable by scipy.
    :return:
    """

    image = ImageIO.load_image(os.path.join(Config.INPUT_DIR, well))
    image = np.stack([
        ImagePreprocessor.preprocess_image(field)
        for field in image])
    image_mosaic = np.concatenate((
        np.concatenate((image[0], image[1]), axis=1),
        np.concatenate((image[3], image[2]), axis=1)),
        axis=0)
    image_mosaic = scipy.ndimage.zoom(
        input=image_mosaic, zoom=zoom + (1,), order=1)
    image_mosaic = (image_mosaic * 65535).astype(np.uint16)

    # If there are only two channels then the FITC (green) channel is missing
    if image_mosaic.shape[-1] == 2:
        image_mosaic = np.stack((
            image_mosaic[..., 0],
            np.zeros(shape=image_mosaic.shape[0:2]),
            image_mosaic[..., 1]), axis=-1)

    int_seg = ImageIO.load_intensity_segmentation(
        os.path.join(Config.INTENSITY_SEGMENTATION, well))
    int_seg_mosaic = np.concatenate((
        np.concatenate((int_seg[0], int_seg[1]), axis=1),
        np.concatenate((int_seg[3], int_seg[2]), axis=1)),
        axis=0)
    int_seg_mosaic = scipy.ndimage.zoom(
        input=int_seg_mosaic, zoom=zoom, order=0)
    int_seg_mosaic = np.stack((
        int_seg_mosaic, int_seg_mosaic, int_seg_mosaic), axis=2)
    int_seg_mosaic[int_seg_mosaic == Config.NONCATEGORY] = Config.BACKGROUND
    int_seg_scaling = 65535 // int_seg_mosaic.max()
    int_seg_mosaic = int_seg_mosaic.astype(np.uint16) * int_seg_scaling

    dnn_seg = ImageIO.load_dnn_segmentation(
        os.path.join(Config.DNN_SEGMENTATION, well))
    dnn_seg_mosaic = np.concatenate((
        np.concatenate((dnn_seg[0], dnn_seg[1]), axis=1),
        np.concatenate((dnn_seg[3], dnn_seg[2]), axis=1)),
        axis=0)
    dnn_seg_mosaic = scipy.ndimage.zoom(
        input=dnn_seg_mosaic, zoom=zoom, order=0)
    dnn_seg_mosaic = np.stack((
        dnn_seg_mosaic, dnn_seg_mosaic, dnn_seg_mosaic), axis=2)
    dnn_seg_scaling = 65535 // dnn_seg_mosaic.max()
    dnn_seg_mosaic = dnn_seg_mosaic.astype(np.uint16) * dnn_seg_scaling

    border = np.ones((image_mosaic.shape[0], 10, 3))
    out = np.concatenate((
        image_mosaic, border,
        int_seg_mosaic, border,
        dnn_seg_mosaic), axis=1)

    out_fn = os.path.join(Config.PREVIEW_DIR, well)
    # This seems rather stable but odd filenames may cause issues here.
    out_fn = os.path.splitext(out_fn)[0] + "." + fmt
    if not os.path.isdir(os.path.dirname(out_fn)):
        os.makedirs(os.path.dirname(out_fn))

    scipy.misc.imsave(name=out_fn, arr=out)


def make_all_images():
    """
    Loops through all data and generates previews
    :return:
    """
    wells = Config.get_image_names()
    for well in wells:
        Utils.write_log(
            msg="Creating Preview for {} of {} ('{}')".format(
                wells.index(well)+1, len(wells), well),
            filename=Config.LOGFILE)
        generate_preview_image(well)


if __name__ == "__main__":
    make_all_images()
