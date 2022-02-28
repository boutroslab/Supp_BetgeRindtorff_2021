"""
This module contains the deep learning functionality, including the training
data generation
"""

import os
import logging
import h5py
import numpy as np
import scipy.ndimage
from tensorflow import keras
import keras.layers
import keras.models
from image_preprocessor import preprocess_image

logging.basicConfig(
    format="[%(levelname)s] %(asctime)s :: %(message)s",
    datefmt="%m/%d/%Y %H:%M:%S",
    level=logging.DEBUG,
)
LOGGER = logging.getLogger("PROMISE/OrganoidSegmentation")


# Various magic numbers and file paths
NUM_CHANNELS = 3
FOREGROUND = 1
BACKGROUND = 0
NONCATEGORY = -1
DNN_WEIGHTS = "DNN_weights.h5"
IO_DIR = "/io"


def load_images_from_h5(fn):
    """
    Loads an image stored on disk.

    This function should return an iterable of numpy arrays with the shape
    (x, y, channels). This iterable could, of course, also be a numpy array
    with the shape (num_images, x, y, channels).

    :param fn: Full path of image
    :return: An iterable of numpy arrays
    """
    with h5py.File(fn, "r") as h5handle:
        proj = h5handle["images"][()]

    return proj


def save_dnn_segmentation(mask, probs, fn):
    """
    Saves the results of the DNN segmentation

    :param mask: A numpy array to save
    :param probs: A numpy array to save
    :param fn: Full path of image
    :return:
    """
    with h5py.File(fn, "w-") as h5handle:
        h5handle.create_dataset(name="images", data=probs, compression=3)
        h5handle.create_dataset(name="mask", data=mask, compression=3)


def run_inference(image_filename):
    """Run inference on a single image.
    Args:
        image_filename (str): Path to the image file
    """
    LOGGER.info(f"Running inference on image at '{image_filename}'.")

    # Load image
    LOGGER.debug(f"({image_filename}) Loading image from disk ...")
    img = load_images_from_h5(fn=image_filename)
    if img.shape[1:] != (3, 2048, 2048):
        raise IOError(
            (
                "Expected input images to have the dimensions 2048 x 2048 with 3 channels but got "
                f"images with dimensions {img.shape[2]} x {img.shape[1]} with {img.shape[3]} "
                "channels."
            )
        )

    # Preprocess image and bring it into the expected [batches, height, width, channels] format
    LOGGER.debug(f"({image_filename}) Running preprocessing ...")
    img = np.stack([preprocess_image(field) for field in img])

    # Set up model and get predictions
    LOGGER.debug(f"({image_filename}) Loading model and generating predictions ...")
    model = build_model(image_size=img.shape[1:], weights=DNN_WEIGHTS)
    res_raw = model.predict(
        x={
            "input_ch0": img[..., 0:1],
            "input_ch1": img[..., 1:2],
            "input_ch2": img[..., 2:3],
        }
    )

    LOGGER.debug(f"({image_filename}) Postprocessing results ...")
    # Softmax to get probabilities
    # Subtracting the maximum value from each label has no impact on the
    # softmax-ing, i.e. exp(x+C)/sum(exp(x+C)) == exp(x)/sum(exp(x)),
    # but prevents float-overflow errors.
    res = np.exp(res_raw - np.max(res_raw, axis=-1, keepdims=True))
    res /= np.sum(res, axis=-1, keepdims=True)
    res = res[..., 1]

    # Rescale back to the original image size
    scaling = [img.shape[ii] / res.shape[ii] for ii in range(len(res.shape))]
    res = scipy.ndimage.zoom(input=res, zoom=scaling, order=1)
    res_bin = (res > 0.5).astype(np.uint8)

    def close_holes(img):
        img = scipy.ndimage.binary_dilation(input=img, structure=np.ones((5, 5)))
        img = scipy.ndimage.binary_fill_holes(input=img)
        img = scipy.ndimage.binary_erosion(input=img, structure=np.ones((5, 5)))
        return img

    res_bin = np.stack([close_holes(field) for field in res_bin]).astype(np.uint8)

    # Cast to float16 to conserve space
    res = np.round(res, 4).astype(np.float16)

    # Save binary prediction and probabilities
    out_path = image_filename + "_segmentation.h5"
    LOGGER.debug(f"({image_filename}) Saving output to {out_path} ...")
    save_dnn_segmentation(mask=res_bin, probs=res, fn=out_path)

    LOGGER.debug(f"({image_filename}) Inference successful! ")


def build_model(image_size, weights=None):
    """
    Builds a fully convolutional keras model that can be trained on segments
    and applied to entire images.

    If 'weights' is None then the network is assumed to be in training mode.
    This means:
    - The kernel size of the layers 'fc_ch*' will be generated dynamically so
      that their output shape is (batch, 1, 1, filters).
    - A flatten layer is added to the output for the loss function to work.

    Else 'weights' must be the full path of the keras weights file. The sizes
    of the weights are extracted and a FCN built accordingly.

    This model has one branch per channel of 'input_shape'

    :param image_size: Numpy array (x, y, channels)
    :param weights: Full path to weights file
    :return:
    """

    w = {}
    with h5py.File(weights, "r") as h5handle:
        groups = list(h5handle.keys())
        for group in groups:
            keys = h5handle[group].keys()
            if len(keys) == 0:
                continue
            if len(keys) != 1:
                raise AssertionError("Expected one key per group")
            w[group] = h5handle[group][list(keys)[0]]["kernel:0"].shape

    image_size = tuple(image_size)
    xy = image_size[0:2]

    inputs = []
    branch_outputs = []
    for ii in range(image_size[-1]):
        input_layer = keras.layers.Input(shape=xy + (1,), name="input_ch{}".format(ii))
        inputs.append(input_layer)
        conv1 = keras.layers.Conv2D(
            filters=32,
            kernel_size=(3, 3),
            padding="valid",
            strides=(1, 1),
            name="conv1_ch{}".format(ii),
        )(input_layer)
        conv1 = keras.layers.LeakyReLU(alpha=0.001, name="conv1_ch{}_lrelu".format(ii))(conv1)
        maxpool1 = keras.layers.MaxPool2D(
            pool_size=(3, 3), strides=(2, 2), padding="valid", name="maxpool1_ch{}".format(ii)
        )(conv1)

        conv2 = keras.layers.Conv2D(
            filters=64,
            kernel_size=(3, 3),
            padding="valid",
            strides=(1, 1),
            name="conv2_ch{}".format(ii),
        )(maxpool1)
        conv2 = keras.layers.LeakyReLU(alpha=0.001, name="conv2_ch{}_lrelu".format(ii))(conv2)
        maxpool2 = keras.layers.MaxPool2D(
            pool_size=(3, 3), strides=(2, 2), padding="valid", name="maxpool2_ch{}".format(ii)
        )(conv2)

        if weights is None:
            conv3 = keras.layers.Conv2D(
                filters=128,
                kernel_size=maxpool2.shape.as_list()[1:3],
                padding="valid",
                strides=(1, 1),
                name="conv3_ch{}".format(ii),
            )(maxpool2)
        else:
            conv3 = keras.layers.Conv2D(
                filters=128,
                kernel_size=w["conv3_ch{}".format(ii)][0:2],
                padding="valid",
                strides=(1, 1),
                name="conv3_ch{}".format(ii),
            )(maxpool2)
        conv3 = keras.layers.LeakyReLU(alpha=0.001, name="conv3_ch{}_lrelu".format(ii))(conv3)
        branch_outputs.append(conv3)

    concat = keras.layers.Concatenate(name="concat_branches")(branch_outputs)

    output = keras.layers.Conv2D(
        filters=2,
        kernel_size=(1, 1),
        padding="valid",
        strides=(1, 1),
        activation="linear",
        name="output_conv",
    )(concat)

    model = keras.models.Model(inputs=inputs, outputs=output)
    model.load_weights(weights)
    return model


if __name__ == "__main__":
    # Look for images in the I/O directory
    all_image_files = os.listdir(IO_DIR)
    LOGGER.info(f"Running inference on {len(all_image_files)} images.")
    for image_filename in all_image_files:
        run_inference(image_filename=os.path.join(IO_DIR, image_filename))
