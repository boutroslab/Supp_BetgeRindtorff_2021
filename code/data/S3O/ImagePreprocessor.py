import numpy as np
import scipy.ndimage
import scipy.ndimage.interpolation

# This is a config dictionary to store all magic numbers
CONFIG = {
    "truncation_low": 0,
    "truncation_high": 99.95,
    "illumination_correction_p": 0.5,
    "background_removal_numbins": 5000,
    "image_features_filter_sizes": (2, 4, 8, 16, 32, 64)
}


def preprocess_image(image):
    """
    Prepares an image for further handling. This involves transposing it so it
    has the shape (spatial, spatial, channels), truncating outliers,
    correcting illumination, and normalizing the channels to lie between 0 and
    1.
    :param image: A 3D numpy array or an object
    :return:
    """

    image = auto_transpose(image)
    image = remove_background(image)
    image = truncate_outliers(
        image, low=CONFIG["truncation_low"], high=CONFIG["truncation_high"])
    image = correct_illumination(
        image=image, p=CONFIG["illumination_correction_p"])
    image = normalize_image(image)

    return image


def auto_transpose(image):
    """
    Automatically transposes an image into the format (X, Y, channels) by
    assuming the smallest dimension is the number of channels
    :param image:
    :return:
    """

    if len(image.shape) != 3:
        raise Exception("'auto_transpose()' only works on 3D numpy arrays")
    channel_axis = [
        s for s in range(len(image.shape)) if
        image.shape[s] == min(image.shape)]
    if len(channel_axis) > 1:
        raise Exception("Image '%s' may only have one axis "
                        "with the smalles value")
    if channel_axis[0] == 2:
        return image
    new_axes = tuple(
        [s for s in range(len(image.shape)) if
         s != channel_axis[0]] + channel_axis)
    return np.transpose(a=image, axes=new_axes)


def truncate_outliers(image, low, high):
    """
    Adjust pixel intensities in image. Any pixel value smaller than the
    'low'-th quantile is set to the 'low'-th quantile and any pixel value
    larger than the 'high'-th quantile is set to the 'high'-th quantile. This
    function only takes 2D/3D arrays or objects that can be coerced to this.
    :param image: A 2D or 3D numpy array or object that can be coerced to
    this. If 3D, then the correction is performed channel-wise. The smallest
    axis is assumed to be the channel axis.
    :param low: The bottom percentile below which to truncate the image. Value
    must lie between (0, 100)
    :param high: The top percentile below which to truncate the image. Value
    must lie between (0, 100)
    :return: A numpy array of the same dimension as the input array (if 3D,
    then the image is transposed into the shape (spatial, spatial, channel))
    """

    image = np.array(image)

    if len(image.shape) == 2:
        q_low = np.percentile(a=image, q=low)
        q_high = np.percentile(a=image, q=high)
        image[image < q_low] = q_low
        image[image > q_high] = q_high
    elif len(image.shape) == 3:
        image = auto_transpose(image)
        for c in range(image.shape[2]):
            q_low = np.percentile(a=image[:, :, c], q=low)
            q_high = np.percentile(a=image[:, :, c], q=high)
            image[:, :, c][image[:, :, c] < q_low] = q_low
            image[:, :, c][image[:, :, c] > q_high] = q_high
    else:
        raise Exception("'image' must be a 2D or 3D numpy array")

    return image


def correct_illumination(image, p):
    """
    Performs illumination correction by looking at the lowest 'p' percentile
    of pixels in each row and column. It then subtracts the row and column
    percentiles from the original image, taking care to retain the input
    dtype.
    :param image:
    :param p:
    :return:
    """
    image = np.array(image)
    dtype = image.dtype

    if len(image.shape) == 2:
        image = np.expand_dims(a=image, axis=2)

    image = auto_transpose(image=image)

    output = []
    for ch in range(image.shape[2]):
        fld = np.array(image[..., ch], dtype=np.float_)
        # Axis 0
        percentiles = np.percentile(a=fld, q=p, axis=0).astype(fld.dtype)
        bg = np.reshape(
            a=np.repeat(a=percentiles, repeats=len(percentiles), axis=0),
            newshape=fld.shape)
        bg = np.transpose(bg)
        fld -= bg

        # Axis 1
        percentiles = np.percentile(a=fld, q=p, axis=1).astype(fld.dtype)
        bg = np.reshape(
            a=np.repeat(a=percentiles, repeats=len(percentiles), axis=0),
            newshape=fld.shape)
        fld -= bg

        # Convert back to original dtype and adjust values to be valid
        try:
            info = np.iinfo(dtype)
        except ValueError:
            info = np.finfo(dtype)
        maxval = info.max
        minval = info.min
        if fld.min() < minval:
            fld -= (fld.min() - minval)
        if fld.max() > maxval:
            oldmin = fld.min()
            fld -= fld.min()
            fld /= fld.max()
            fld *= (maxval - oldmin)
            fld += oldmin

        output.append(fld.astype(dtype=dtype))

    if len(output) == 1:
        return output[0]
    else:
        return np.stack(output, axis=2)


def remove_background(image):
    """
    This function attempts to darken the background and improve contrast. It
    assumes that there will be more background than foreground and chooses the
    intensity histogram maximum as the "background value" to normalize to 0
    :param image: A 2D or 3D numpy array. If 3D, the smallest axis is assumed
    to represent the channels.
    :return:
    """
    dtype = image.dtype
    dims = len(image.shape)
    image = np.array(image, dtype=np.float32)

    # Adjust dimensions
    if dims == 2:
        image = np.expand_dims(a=image, axis=2)
    image = auto_transpose(image=image)

    # Go through each channel and perform background removal
    nbins = CONFIG["background_removal_numbins"]
    for channel in range(image.shape[-1]):
        field = image[..., channel]
        vals = np.reshape(field, (-1,))
        hist = np.histogram(a=vals, bins=nbins)
        modeindex = np.argmax(hist[0])
        modeval = np.mean(hist[1][modeindex:modeindex + 2])

        field -= modeval
        field[field < 0] = 0

    output = np.round(a=image).astype(dtype=dtype)

    if len(output) == 1:
        return output[0, ...]
    else:
        return output


def normalize_image(image, dtype=np.float32):
    """
    Normalizes an image so that its values are between 0 and 1 (casts it to
    'dtype' as well)
    :param image: A 2D or 3D numpy array. If 3D, the smallest axis is assumed
    to represent the channels.
    :param dtype: Optionally a dtype to cast the values to. Defaults to
    float32.
    :return:
    """

    image = np.array(image, dtype=dtype)

    if len(image.shape) == 3:
        image = auto_transpose(image)

    image -= np.min(image, axis=(0, 1))
    image /= np.max(image, axis=(0, 1))
    return image


def elastic_transform(image, alpha, sigma, random_state=None):
    """Elastic deformation of images as described in [Simard2003]_.
    .. [Simard2003] Simard, Steinkraus and Platt, "Best Practices for
       Convolutional Neural Networks applied to Visual Document Analysis", in
       Proc. of the International Conference on Document Analysis and
       Recognition, 2003.
    """
    assert len(image.shape) == 2

    if random_state is None:
        random_state = np.random.RandomState(None)

    shape = image.shape

    dx = scipy.ndimage.gaussian_filter(
        (random_state.rand(*shape) * 2 - 1), sigma,
        mode="constant", cval=0) * alpha
    dy = scipy.ndimage.gaussian_filter(
        (random_state.rand(*shape) * 2 - 1), sigma,
        mode="constant", cval=0) * alpha

    x, y = np.meshgrid(np.arange(shape[0]), np.arange(shape[1]), indexing='ij')
    indices = np.reshape(x + dx, (-1, 1)), np.reshape(y + dy, (-1, 1))

    return scipy.ndimage.interpolation.map_coordinates(image, indices, order=1).reshape(shape)
