from __future__ import division
import os
import ImagePreprocessor
import h5py
import numpy as np
import scipy.ndimage
import sys
import keras.layers
import keras.models


def segment_plate(in_file, weights_file, out_file=None):
    """
    Load the hdf5projection and apply the DNN segmentation to it. Saves to
    outfile, otherwise returns the segmentation results.

    The outfile is an HDF5 file with two datasets, the segmentation
    probabilities and the resized mask:
    FILE
        DATASET probs
        DATASET mask

    This function assumes the PROMISE folder and file structure.

    :param in_file: The full path of the hdf5projection to load. Should have
        the shape (num_fields, num_channels, spatial, spatial)
    :param weights_file: The keras model weights file
    :param out_file: The full path of the file to write to. If 'None' then
        return the segmentation.
    :return: If outfile is None: A tuple (segmentation probabilities,
        resized segmentation mask). If outfile is not None: no return value.
    """

    # Load file
    with h5py.File(in_file, "r") as h5:
        images = h5["images"][()]
        metadata = h5["metadata"][()]

    images = np.stack([
        ImagePreprocessor.preprocess_image(image)
        for image in images])

    # Set up model
    imshape = images.shape[1:]
    model = build_model(image_size=imshape, weights=weights_file)

    # Generate input dictionary
    input_dict = {
        "input_ch{}".format(ii): images[..., ii:ii+1] 
        for ii in range(images.shape[-1])}
    res_raw = model.predict(x=input_dict)
    # Subtracting the maximum value from each label has no impact on the
    # softmax-ing, i.e. exp(x+C)/sum(exp(x+C)) == exp(x)/sum(exp(x)),
    # but prevents float-overflow errors.
    res = np.exp(res_raw - np.max(res_raw, axis=-1, keepdims=True))
    res /= np.sum(res, axis=-1)[..., None]
    res_bin = (res[..., 1] > 0.5).astype(np.uint8)

    # Rescale back to the original image size
    scaling = [
        images.shape[ii] / res_bin.shape[ii]
        for ii in range(len(res_bin.shape))]
    res_bin = scipy.ndimage.zoom(input=res_bin, zoom=scaling, order=0)

    def close_holes(img):
        img = scipy.ndimage.binary_dilation(
            input=img, structure=np.ones((5, 5)))
        img = scipy.ndimage.binary_fill_holes(input=img)
        img = scipy.ndimage.binary_erosion(
            input=img, structure=np.ones((5, 5)))
        return img

    res_bin = np.stack([close_holes(field) for field in res_bin]).astype(
        np.uint8)

    # Cast to float16 to conserve space
    res = np.round(res, 4).astype(np.float16)

    if out_file is None:
        return res, res_bin
    else:
        # Make the outfile directory if necessary
        outdir = os.path.dirname(out_file)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)

        with h5py.File(out_file, "w-") as h5:
            h5.create_dataset(name="probs", data=res, compression=3)
            h5.create_dataset(name="mask", data=res_bin, compression=3)
            h5.create_dataset(name="metadata", data=metadata)
        return None


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
    if weights is not None:
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
        input_layer = keras.layers.Input(
            shape=xy + (1,),
            name="input_ch{}".format(ii))
        inputs.append(input_layer)
        conv1 = keras.layers.Conv2D(
            filters=32, kernel_size=(3, 3), padding="valid",
            strides=(1, 1), name="conv1_ch{}".format(ii))(input_layer)
        conv1 = keras.layers.LeakyReLU(
            alpha=0.001, name="conv1_ch{}_lrelu".format(ii))(conv1)
        maxpool1 = keras.layers.MaxPool2D(
            pool_size=(3, 3), strides=(2, 2),
            padding="valid", name="maxpool1_ch{}".format(ii))(conv1)

        conv2 = keras.layers.Conv2D(
            filters=64, kernel_size=(3, 3), padding="valid",
            strides=(1, 1), name="conv2_ch{}".format(ii))(maxpool1)
        conv2 = keras.layers.LeakyReLU(
            alpha=0.001, name="conv2_ch{}_lrelu".format(ii))(conv2)
        maxpool2 = keras.layers.MaxPool2D(
            pool_size=(3, 3), strides=(2, 2),
            padding="valid", name="maxpool2_ch{}".format(ii))(conv2)

        if weights is None:
            conv3 = keras.layers.Conv2D(
                filters=128, kernel_size=maxpool2.shape.as_list()[1:3],
                padding="valid", strides=(1, 1),
                name="conv3_ch{}".format(ii))(maxpool2)
        else:
            conv3 = keras.layers.Conv2D(
                filters=128, kernel_size=w["conv3_ch{}".format(ii)][0:2],
                padding="valid", strides=(1, 1),
                name="conv3_ch{}".format(ii))(maxpool2)
        conv3 = keras.layers.LeakyReLU(
            alpha=0.001, name="conv3_ch{}_lrelu".format(ii))(conv3)
        branch_outputs.append(conv3)

    concat = keras.layers.Concatenate(name="concat_branches")(branch_outputs)

    output = keras.layers.Conv2D(
        filters=2, kernel_size=(1, 1),
        padding="valid", strides=(1, 1),
        activation="linear", name="output_conv")(concat)

    if weights is None:
        output = keras.layers.Flatten()(output)

    model = keras.models.Model(
        inputs=inputs,
        outputs=output)

    if weights is not None:
        model.load_weights(weights)

    return model


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python %s <PLATE> <INFILE> <OUTFILE> <PYDIR>" % sys.argv[0])

    prog, plate, infile, outfile, pydir = sys.argv

    weightsfile = os.path.join(pydir, "DNN_weights.pkl")

    # If the weights file doesn't exist yet, then skip this entire procedure
    if os.path.isfile(weightsfile):
        segment_plate(
            in_file=infile, out_file=outfile,
            weights_file=weightsfile)
