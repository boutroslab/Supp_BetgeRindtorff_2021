"""
This module performs data augmentation on DNN training data
"""
import numpy as np
import scipy.ndimage
import scipy.ndimage.interpolation
import scipy.ndimage.filters
import ImagePreprocessor
import h5py


# Main function
def augment_data(input_image):
    """
    Performs data augmentation and returns the original data and augmented
    sets. It reshapes the input image, if necessary, to have the shape
    (spatial, spatial, channels).

    It expects the data to be normalized to the range [0, 1].
    :param input_image: A 2D or 3D numpy array. If 3D, it expects the smallest
    dimension to be the channel
    :return:
    """

    if input_image.max() > 1 or input_image.min() < 0:
        raise Exception("'input_image' must be in the range [0, 1]")

    input_dim = len(input_image.shape)

    if input_dim != 2 and input_dim != 3:
        raise Exception("'image' must be a 2D or 3D array")

    if input_dim == 3:
        input_image = ImagePreprocessor.auto_transpose(input_image)

    data_aug = np.expand_dims(input_image, axis=0)
    instructions = augment_config()
    for entry in instructions:
        func = entry["func"]
        args = {k: v for k, v in entry.items() if k != "func"}
        step_output = []
        for i in range(data_aug.shape[0]):
            out_img = func(image=data_aug[i, ...], **args)
            # Transformations may occasionally lead to values that are
            # outside [0, 1], so a manual fix is needed here
            out_img[out_img > 1] = 1
            out_img[out_img < 0] = 0
            step_output.append(out_img)
        step_output = np.concatenate(step_output)
        data_aug = np.concatenate((data_aug, step_output), axis=0)

    return data_aug


def augment_config():
    """
    This is a placeholder function that contains the instructions for the
    augmentations to perform and indicates how many total images are created
    for each input image.

    The instruction list should be a list of dictionaries with the structure
    [{"func": function_name, "arg1": arg1, ..., "argN": argN}, ... ]. The keys of
    the entries should be 'func' and the names of the arguments that the
    functions take. The value for 'func' should be an actual reference to the
    function while the values for the arguments should be the data types that
    the function expects.

    Note that the order is important! Instructions are applied to the FULL
    output of the previous instruction, i.e. if 'mirror_and_rotate' is
    followed by 'elastic_transform' then an independent elastic transformation
    is applied to each rotated image whereas if the order is reversed then the
    elastically transformed image(s) is/are rotated

    NOTE: Data augmentation has been shown to improve DNN training, however
    we have SO MANY organoid images, it is completely unnecessary here.
    :return:
    """

    # instructions = [
    #     {"func": add_noise, "nt": 3},
    #     {"func": elastic_transform, "nt": 2}]
    # instructions = [
    #     {"func": mirror_and_rotate}]
    instructions = []
    return instructions


def elastic_transform(image, nt=1, alpha=34, sigma=4):
    """
    Elastic deformation of images as described in:
    Simard, Steinkraus and Platt, "Best Practices for Convolutional Neural
    Networks applied to Visual Document Analysis", in Proc. of the
    International Conference on Document Analysis and Recognition, 2003.

    Implemented by https://gist.github.com/chsasank (elastic_transform.py)
    and adapted by me for 3D, channel-independent input

    The default values for alpha and sigma are the optimal values found by
    Simard et al. (2003) for the MNIST data set.

    :param image: A 2D or 3D numpy array. If 3D then each channel is
    transformed with the same parameters. If 3D, then it expects the
    smallest dimension to be the number of channels.
    :param nt: The number of (independent) transformations to perform
    :param alpha:
    :param sigma:
    :return A 3D numpy array with the shape (spatial, spatial, channels)
    """
    input_dim = len(image.shape)

    if input_dim != 2 and input_dim != 3:
        raise Exception("'image' must be a 2D or 3D array")

    if input_dim == 2:
        image = np.expand_dims(a=image, axis=2)

    image = ImagePreprocessor.auto_transpose(image)

    random_state = np.random.RandomState(None)

    output_img = []
    for i in range(nt):
        dx = scipy.ndimage.filters.gaussian_filter(
            (random_state.rand(image.shape[0], image.shape[1]) * 2 - 1), sigma,
            mode="constant", cval=0) * alpha
        dy = scipy.ndimage.filters.gaussian_filter(
            (random_state.rand(image.shape[0], image.shape[1]) * 2 - 1), sigma,
            mode="constant", cval=0) * alpha

        x, y = np.meshgrid(
            np.arange(image.shape[0]),
            np.arange(image.shape[1]),
            indexing='ij')
        indices = np.reshape(x + dx, (-1, 1)), np.reshape(y + dy, (-1, 1))

        output_img.append(np.stack([
            scipy.ndimage.interpolation.map_coordinates(
                image[..., i], indices, order=1).reshape(image.shape[0:2])
            for i in range(image.shape[2])], axis=2))

    output_img = np.stack(output_img)

    # Ensure that output_img has the same shape as image
    if input_dim == len(output_img.shape[1:]):
        return output_img
    elif input_dim == 2 and len(output_img.shape[1:]) == 3:
        return np.squeeze(a=output_img, axis=-1)
    elif input_dim == 3 and len(output_img.shape[1:]) == 2:
        return np.expand_dims(a=output_img, axis=-1)

    return None


def mirror_and_rotate(image):
    """
    Mirrors and rotates an image in multiples of 90 degrees (non-right
    rotations would lead to black pixels that might influence the importance
    of pixels at the edge of the image). It creates a total of 7 new images
    (the flipped image and the rotations of the original and the flipped
    images)
    :param image: A 2D or 3D numpy array. If 3D then each channel is
    transformed with the same parameters. If 3D, then it expects the
    smallest dimension to be the number of channels.
    :return: A numpy array with the dimensions (7, spatial, spatial,
    channels)
    """

    if len(image.shape) != 2 and len(image.shape) != 3:
        raise Exception("'image' must be a 2D or 3D array")

    if len(image.shape) == 3:
        image = ImagePreprocessor.auto_transpose(image)

    flipped = np.flip(m=image, axis=1)

    output_img = [
        flipped,
        np.rot90(m=image, k=1, axes=(0, 1)),
        np.rot90(m=image, k=2, axes=(0, 1)),
        np.rot90(m=image, k=3, axes=(0, 1)),
        np.rot90(m=flipped, k=1, axes=(0, 1)),
        np.rot90(m=flipped, k=2, axes=(0, 1)),
        np.rot90(m=flipped, k=3, axes=(0, 1))]

    return np.stack(output_img)


def add_noise(image, nt=1, scale=0.1):
    """
    Adds salt & pepper noise to an image

    The input image is assumed to be normalized to [0,1] and the output
    intensities are truncated to stay in this range.
    :param image:
    :param nt: The number of (independent) transformations to perform
    :param scale: The intensity change is chosen from a truncated normal
           distribution with a standard deviation of 'scale'.
    :return:
    """
    if image.max() > 1 or image.min() < 0:
        raise Exception("'image' must be in the range [0, 1]")

    outimg = []
    for i in range(nt):
        noise = np.random.normal(loc=0, scale=scale, size=image.shape)
        outimg.append(image + noise)
    outimg = np.stack(outimg)
    outimg[outimg < 0] = 0
    outimg[outimg > 1] = 1

    return outimg


def make_test_image(radius, demo_image="circle", noise=0):
    """
    This function creates a demo image
    :param radius: The radius of the image
    :param demo_image: The type of demo image (currently implemented are:
    circle, circle_3d)
    :param noise: The standard deviation of the normally distributed noise
           that is added to the image. A value of 0 means no noise.
    :return:
    """
    if demo_image == "circle":
        circle = np.zeros((2 * radius + 1, 2 * radius + 1))
        circle_r = np.sum(np.square(np.indices(circle.shape) - radius), axis=0)
        width = 2
        circle[
            (circle_r >= (radius - width) ** 2) *
            (circle_r <= (radius + width) ** 2)] = 1
        outimg = circle
    elif demo_image == "circle_3d":
        circle = np.zeros((2 * radius + 1, 2 * radius + 1))
        circle_r = np.sum(np.square(np.indices(circle.shape) - radius), axis=0)
        width = 2
        circle[
            (circle_r >= (radius - width) ** 2) *
            (circle_r <= (radius + width) ** 2)] = 1
        outimg = np.stack([circle] * 3, axis=2)
    elif demo_image == "organoid":
        with h5py.File(
                "/Users/jansauer/Thesis/Projects/MultiStageSegmentation/"
                "PROMISE/hdf5projection/plate1_2525_15_35/"
                "plate1_2525_15_35_A_03_contrastProjections.h5",
                "r") as f:
            full_img = f["images"][0, ...]
        full_img = ImagePreprocessor.preprocess_image(full_img)

        with h5py.File(
                "/Volumes/B210-Projekte/users/sauerja/MultiStageSegmentation/"
                "Results_4x_Mouse/dotmaps/"
                "plate1_2525_15_35_A_03_contrastProjections.h5",
                "r") as f:
            dotmap = f["images"][()]

        _, px, py = np.transpose(np.where(dotmap == 1))[0]
        outimg = full_img[px - 10:px + 11, py - 10: py + 11, :]
    else:
        return None

    if noise > 0:
        outimg += np.random.normal(loc=0, scale=noise, size=outimg.shape)
        outimg[outimg > 1] = 1
        outimg[outimg < 0] = 0

    return outimg
