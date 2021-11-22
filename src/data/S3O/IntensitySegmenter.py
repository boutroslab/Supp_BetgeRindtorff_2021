# This module performs intensity based segmentation.
import Config
import Utils
import ImageIO
import ImagePreprocessor
import os
import numpy as np
import scipy.optimize
import scipy.ndimage
import skimage.filters


def run_intensity_segmentation(low=None, high=None, logfile=None):
    """
    Performs intensity segmentation

    This is a cluster-friendly function in that it only requires a lower
    and an upper index limit and will then process the corresponding
    training images.

    'logfile' dictates where the output should be written to, i.e. if
    'logfile == None' then output is directed to STDOUT (convenient for
    writing output to a cluster log file), otherwise output is written
    to the file given by 'logfile'

    If low == None then it is assumed to be 0
    If high == None then it is assumed to be the length of the image list

    :param low: Low index of the images to process
    :param high: High index of the images to process (exclusive)
    :param logfile: Path of logfile. Directs to STDOUT if None
    :return:
    """

    image_names = Config.get_image_names()

    if low is None:
        low = 0
    if high is None:
        high = len(image_names)

    for image_name in image_names[low:high]:
        outf = os.path.join(Config.INTENSITY_SEGMENTATION, image_name)

        if os.path.isfile(outf):
            Utils.write_log(
                msg="File exists: '{}'".format(outf),
                filename=logfile)
            continue

        if not os.path.isdir(os.path.dirname(outf)):
            os.makedirs(os.path.dirname(outf))

        Utils.write_log(
            msg="Processing image {} of {} ('{}')".format(
                image_names.index(image_name) + 1,
                len(image_names), image_name),
            filename=logfile)

        # Load image
        in_image = os.path.join(Config.INPUT_DIR, image_name)
        fields = ImageIO.load_image(fn=in_image)
        fields = [ImagePreprocessor.preprocess_image(field) for field in fields]

        # Intensity segmentation
        intensity_seg = [segment_by_intensity(image=f) for f in fields]
        intensity_seg = np.stack(intensity_seg)

        ImageIO.save_intensity_segmentation(image=intensity_seg, fn=outf)


def segment_by_intensity(image):
    """
    Performs an intensity-based approximate segmentation

    I perform otsu segmentation on the grayscale image.

    To improve the segmentation, I demand that the background be contiguous.
    This should prevent invalid segmentation of parts of organoids.

    :param image: A 3D numpy array.
    :return: A segmentation map
    """

    # Make sure the shape is (x, y, channels)
    image = ImagePreprocessor.auto_transpose(image=image)

    # === FOREGROUND ===
    # The foreground is determined from the actin channel
    actin = image[..., 0]
    fg_otsu = skimage.filters.threshold_otsu(image=actin)

    foreground = actin >= fg_otsu

    # Close holes in the foreground mask
    foreground = scipy.ndimage.binary_dilation(
        input=foreground, structure=np.ones((5, 5)))
    foreground = scipy.ndimage.binary_fill_holes(input=foreground)
    foreground = scipy.ndimage.binary_erosion(
        input=foreground, structure=np.ones((5, 5)))

    # === BACKGROUND === #
    # Grayscale image
    gs = np.sum(a=image, axis=-1)
    gs = ImagePreprocessor.normalize_image(image=gs)

    bg_otsu = skimage.filters.threshold_otsu(image=gs)
    background = gs < bg_otsu

    # Perform a small opening on the background mask to separate off
    # smaller segments that aren't actual background
    background = scipy.ndimage.binary_opening(
        input=background, structure=np.ones((5, 5)))

    # The background must be contiguous, with the exception of elements
    # directly on the image boundary
    bg_labeled, n_labels = scipy.ndimage.label(background)
    region_sizes = [
        np.sum(bg_labeled == label)
        for label in range(1, n_labels + 1)]
    real_bg_value = [
        (i + 1) for i in range(len(region_sizes)) if
        region_sizes[i] == max(region_sizes)]

    # If there are multiple maximal background areas then
    # the image cannot be automatically parsed and is skipped
    if len(real_bg_value) != 1:
        return np.ones(image.shape[0:2], dtype=np.int8) * Config.NONCATEGORY

    real_bg_value = real_bg_value[0]

    # Remove all background segments that aren't the largest
    for label in range(1, n_labels + 1):
        if label == real_bg_value:
            continue
        # Uncomment this to also keep background segments on the edge.
        # These could be true background segments separated by an organoids
        # clump from the rest of the background but it is more likely that
        # they are partial organoids.
        # region_x, region_y = np.where(bg_labeled == label)
        # if np.min(region_x) == 0 or \
        #    np.max(region_x) == (gs.shape[0] - 1) or \
        #    np.min(region_y) == 0 or \
        #    np.max(region_y) == (gs.shape[1] - 1):
        #     continue
        else:
            background[bg_labeled == label] = 0

    # In addition, I only want to use background that is sufficiently distant
    # from the foreground to prevent incorrectly classified organoids from
    # being trained as background.
    background = np.bitwise_and(
        background,
        ~scipy.ndimage.maximum_filter(input=foreground, size=150))

    output = np.ones(image.shape[0:2], dtype=np.int8) * Config.NONCATEGORY
    output[background] = Config.BACKGROUND
    output[foreground] = Config.FOREGROUND

    return output


def segment_by_intensity_v1(image):
    """
    WARNING: THIS METHOD DIDN'T WORK. I MAY HAVE NOT USED ENOUGH TRAINING
    SAMPLES

    Performs an intensity-based approximate segmentation

    I perform otsu segmentation for the background on the grayscale image
    but also for the foreground on the actin channel.

    To improve the segmentation, I demand that the background be contiguous.
    This should prevent invalid segmentation of parts of organoids.

    :param image: A 3D numpy array.
    :return: A segmentation map
    """

    # Make sure the shape is (x, y, channels)
    image = ImagePreprocessor.auto_transpose(image=image)

    # === BACKGROUND === #
    # Grayscale image
    gs = np.sum(a=image, axis=-1)
    gs = ImagePreprocessor.normalize_image(image=gs)

    bg_otsu = skimage.filters.threshold_otsu(image=gs)
    background = gs < bg_otsu

    # Perform a small opening on the background mask to separate off
    # smaller segments that aren't actual background
    background = scipy.ndimage.binary_opening(
        input=background, structure=np.ones((5, 5)))

    # The background must be contiguous, with the exception of elements
    # directly on the image boundary
    bg_labeled, n_labels = scipy.ndimage.label(background)
    region_sizes = [
        np.sum(bg_labeled == label)
        for label in range(1, n_labels + 1)]
    real_bg_value = [
        (i + 1) for i in range(len(region_sizes)) if
        region_sizes[i] == max(region_sizes)]

    # If there are multiple maximal background areas then
    # the image cannot be automatically parsed and is skipped
    if len(real_bg_value) != 1:
        return np.ones(image.shape[0:2], dtype=np.int8) * -1

    real_bg_value = real_bg_value[0]

    # Remove all background segments that aren't the largest
    for label in range(1, n_labels + 1):
        if label == real_bg_value:
            continue
        # Uncomment this to also keep background segments on the edge.
        # These could be true background segments separated by an organoids
        # clump from the rest of the background but it is more likely that
        # they are partial organoids.
        # region_x, region_y = np.where(bg_labeled == label)
        # if np.min(region_x) == 0 or \
        #    np.max(region_x) == (gs.shape[0] - 1) or \
        #    np.min(region_y) == 0 or \
        #    np.max(region_y) == (gs.shape[1] - 1):
        #     continue
        else:
            background[bg_labeled == label] = 0

    # === FOREGROUND ===
    # The foreground is determined from the actin channel
    actin = image[..., 0]
    fg_otsu = skimage.filters.threshold_otsu(image=actin)

    foreground = actin >= fg_otsu

    # Close holes in the foreground mask
    foreground = scipy.ndimage.binary_dilation(
        input=foreground, structure=np.ones((5, 5)))
    foreground = scipy.ndimage.binary_fill_holes(input=foreground)
    foreground = scipy.ndimage.binary_erosion(
        input=foreground, structure=np.ones((5, 5)))

    output = np.ones(image.shape[0:2], dtype=np.int8) * Config.NONCATEGORY
    output[background] = Config.BACKGROUND
    output[foreground] = Config.FOREGROUND

    return output
