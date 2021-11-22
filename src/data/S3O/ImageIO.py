# Contains functions to read and write images
import h5py


def load_image(fn):
    """
    Loads an image stored on disk.

    This function should return an iterable of numpy arrays with the shape
    (channels, x, y). This iterable could, of course, also be a numpy array
    with the shape (num_images, channels, x, y).

    :param fn: Full path of image
    :return: An iterable of numpy arrays
    """
    with h5py.File(fn, "r") as h5handle:
        proj = h5handle["images"][()]

    # This is a hack to ensure that only the Actin+Cell Nucelus channels are
    # used
    return proj[:, (0, 2), :, :]


def load_intensity_segmentation(fn):
    """
    Loads an intensity segmentation stored on disk.

    This function should return an iterable of numpy arrays with the shape
    (x, y). This iterable could, of course, also be a numpy array
    with the shape (num_images, x, y).

    :param fn: Full path of image
    :return: An iterable of numpy arrays
    """
    with h5py.File(fn, "r") as h5handle:
        proj = h5handle["images"][()]

    return proj


def save_intensity_segmentation(image, fn):
    """
    Saves the results of the intensity segmentation

    :param image: A numpy array to save
    :param fn: Full path of image
    :return:
    """

    with h5py.File(fn, "w-") as h5handle:
        h5handle.create_dataset(name="images", data=image, compression=3)


def load_dnn_segmentation(fn, mask=True):
    """
    Loads an intensity segmentation stored on disk.

    This function should return an iterable of numpy arrays with the shape
    (x, y). This iterable could, of course, also be a numpy array
    with the shape (num_images, x, y).

    :param fn: Full path of image
    :param mask: A boolean indicating if either the binary mask (True) or the
                 foreground probability map (False) should be returnd.
    :return: An iterable of numpy arrays
    """

    key = "mask" if mask else "images"

    with h5py.File(fn, "r") as h5handle:
        image = h5handle[key][()]

    return image


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


def save_dnn_features(features, labels, fn):
    """
    Saves the pre-classifier features of the DNN segmentation

    :param features: A numpy array to save
    :param labels: A numpy array to save
    :param fn: Full path of image
    :return:
    """

    with h5py.File(fn, "w-") as h5handle:
        h5handle.create_dataset(name="features", data=features, compression=3)
        h5handle.create_dataset(name="labels", data=labels, compression=3)
