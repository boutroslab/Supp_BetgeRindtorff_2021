# TODO: Refactor this to make every function output a float array
import numpy as np
import scipy.ndimage
import skimage.morphology
import skimage.filters
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

    image = np.array(image, dtype=np.float32)
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

    image = np.array(image)

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
    return image.transpose(new_axes)


def truncate_outliers(image, low, high):
    """
    Adjust pixel intensities in image. Any pixel value smaller than the 
    'low'-th quantile is set to the 'low'-th quantile and any pixel value 
    larger than the 'high'-th quantile is set to the 'high'-th quantile. This 
    function only takes 2D arrays or objects that can be coerced to this.
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
    to represent the channels. Must be a floating-point dtype
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


def calc_image_features(image):
    """
    Calculates predefined features of an image
    :param image: A 3D numpy array. The smallest axis is assumed to represent 
    the channels. The values of 'image' must lie between 0 and 1 (i.e. it 
    must be normalized) 
    :return:
    """

    if image.min() != 0 or image.max() != 1:
        raise Exception("min(image) must be 0 and max(image) must be 1")

    image = auto_transpose(image)

    features = []

    # WARNING: Entropy takes a VERY long time to calculate
    image_binned = np.cast[np.uint8](np.round(np.array(image) * 255))
    for size in tuple(CONFIG["image_features_filter_sizes"]):
        features.append(calc_feature(
            image=image_binned, func=skimage.filters.rank.entropy,
            selem=skimage.morphology.disk(size)))

    # Contrast
    def calc_contrast(im, s):
        mean_intensity_squared = scipy.ndimage.uniform_filter(
            im, size=s) ** 2
        mean_squared_intensity = scipy.ndimage.uniform_filter(
            im ** 2, size=s)
        return mean_squared_intensity - mean_intensity_squared
    for size in tuple(CONFIG["image_features_filter_sizes"]):
        features.append(calc_feature(
            image=image, func=calc_contrast,
            s=size))

    # Blurring
    for sigma in tuple(CONFIG["image_features_filter_sizes"]):
        features.append(calc_feature(
            image=image, func=scipy.ndimage.filters.gaussian_filter,
            sigma=sigma))

    # Edges
    grad_filter = [[[1, 2, 1], [0, 0, 0], [-1, -2, -1]],
                   [[-1, -2, -1], [0, 0, 0], [1, 2, 1]],
                   [[2, 1, 0], [1, 0, -1], [0, -1, -2]],
                   [[-2, -1, 0], [-1, 0, 1], [0, 1, 2]],
                   [[1, 0, -1], [2, 0, -2], [1, 0, -1]],
                   [[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]],
                   [[0, 1, 2], [-1, 0, 1], [-2, -1, 0]],
                   [[0, -1, -2], [1, 0, -1], [2, 1, 0]]]

    for filt in grad_filter:
        features.append(calc_feature(
            image=image, func=scipy.ndimage.filters.convolve,
            weights=np.array(filt)))

    # Morphological gradient
    for size in tuple(CONFIG["image_features_filter_sizes"]):
        features.append(calc_feature(
            image=image, func=scipy.ndimage.morphology.morphological_gradient,
            size=size))

    # Morphological laplace
    for size in tuple(CONFIG["image_features_filter_sizes"]):
        features.append(calc_feature(
            image=image, func=scipy.ndimage.morphology.morphological_laplace,
            size=size))

    features = np.concatenate([image] + features, axis=2)

    return features


def calc_feature(image, func, **kwargs):
    """
    A helper function to generate features for a multichannel image for each 
    channel individually.
    :param image: A 3D numpy array. The smallest axis is assumed to represent 
    the channels
    :param func: A function that can be applied to a 2D numpy array. The first 
    argument of 'func' must be the array
    :param kwargs: Further arguments to be passed to func
    :return: 
    """

    image = auto_transpose(image)

    output = []
    for c in range(image.shape[2]):
        output.append(func(image[:, :, c], **kwargs))
    output = np.stack(output, axis=2)
    return output


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
