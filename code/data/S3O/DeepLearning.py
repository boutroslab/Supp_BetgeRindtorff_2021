"""
This module contains the deep learning functionality, including the training
data generation
"""
import Config
import ImageIO
import ImagePreprocessor
import Utils
import DataAugmenter
import numpy as np
import h5py
import os
import sklearn.model_selection
import sklearn.metrics
import scipy.ndimage
import keras.layers
import keras.models


def run_training_segment_creation(low=None, high=None, logfile=None):
    """
    Loads the intensity segmentation, creates training segments with a radius
    of Config.DNN_SEGMENT_RADIUS, and stores them in an HDF5 file with the
    format:

    FILE:
        Dataset "images" (num_images, x, y, channels)
        Dataset "labels" (num_labels,)

    It can be run on a subset of all training images, i.e. will create multiple
    hdf5 files that must be merged (see Utils.py)

    The function extracts an equal number of segments from each of the full
    images.

    'logfile' dictates where the output should be written to, i.e. if
    'logfile == None' then output is directed to STDOUT (convenient for
    writing output to a cluster log file), otherwise output is written
    to the file given by 'logfile'

    If low == None then it is assumed to be 0
    If high == None then it is assumed to be the length of the image list

    :param low: Low index of the images to process
    :param high: High index of the images to process (exclusive)
    :param logfile: (Optional) Path of logfile. Directs to STDOUT if None
    :return:
    """
    Utils.write_log(msg="Creating Training Segments", filename=logfile)

    # Load image names and remember the fraction of images used
    image_names = Config.get_image_names()
    total_num = len(image_names)
    if low is None:
        low = 0
    if high is None:
        high = len(image_names)
    image_names = image_names[low:high]

    # Determine the degree of data augmentation
    seg_radius = Config.get_dnn_parameters()["seg_radius"]
    black_square = np.ones(
        shape=(seg_radius[0] * 2 + 1, seg_radius[1] * 2 + 1,
               Config.NUM_CHANNELS))
    black_square_aug = DataAugmenter.augment_data(black_square)
    aug_degree = black_square_aug.shape[0]

    # Determine how many points to sample from each image for each category
    num_segments = Config.get_dnn_parameters()["num_segments"]
    samples_per_image = int(np.ceil(
        num_segments / (total_num * aug_degree)))

    # Create the hdf5 file
    # hdf5_fn = Config.DNN_TRAINING_SEGMENTS + "_{}_{}".format(low, high)
    hdf5_fn = Config.DNN_TRAINING_SEGMENTS
    if not os.path.isdir(os.path.dirname(hdf5_fn)):
        os.makedirs(os.path.dirname(hdf5_fn))
    seg_size = [s * 2 + 1 for s in seg_radius]
    with h5py.File(hdf5_fn, "w-") as h5handle:
        h5handle.create_dataset(
            name="image_names",
            data=np.array(image_names).astype(np.bytes_),
            compression=3)
        h5handle.create_dataset(
            name="images",
            shape=(samples_per_image * len(image_names) * 2,
                   seg_size[0], seg_size[1], Config.NUM_CHANNELS),
            chunks=(1, seg_size[0], seg_size[1], Config.NUM_CHANNELS),
            compression=3, dtype=np.float32)
        h5handle.create_dataset(
            name="labels",
            shape=(samples_per_image * len(image_names) * 2,),
            compression=3, dtype=np.uint8)

    # Go through each image and extract + save segments
    # It is inefficient to write every single segment individually to disk,
    # so they are packaged into batches and then written in a single operation
    current_index = 0
    write_batch_size = 10000
    segment_batch = []
    label_batch = []
    for image_name in image_names:
        Utils.write_log(
            msg="Processing image {} of {} ('{}')".format(
                image_names.index(image_name) + 1,
                len(image_names), image_name),
            filename=logfile)
        # Load projection
        proj_fn = os.path.join(Config.INPUT_DIR, image_name)
        proj = ImageIO.load_image(proj_fn)
        proj = np.stack([
            ImagePreprocessor.preprocess_image(field) for field in proj])

        # Load intensity segmentation
        int_seg_fn = os.path.join(Config.INTENSITY_SEGMENTATION, image_name)
        int_seg = ImageIO.load_intensity_segmentation(int_seg_fn)

        # Extract segments
        points = np.transpose(np.where(int_seg != Config.NONCATEGORY))
        # Keep only points that are not on the border
        points = points[
            (points[:, 1] >= seg_radius[0]) *
            (points[:, 1] < proj.shape[1] - seg_radius[0]) *
            (points[:, 2] >= seg_radius[1]) *
            (points[:, 2] < proj.shape[2] - seg_radius[1]), :]
        labels = int_seg[tuple(np.transpose(points))]
        label_set = tuple(set(labels))
        for label in label_set:
            sel_indices = np.random.choice(
                a=np.where(labels == label)[0],
                size=samples_per_image, replace=False)
            sel_points = points[sel_indices]
            for point in sel_points:
                fld, x, y = point
                x_low = x - seg_radius[0]
                x_high = x + seg_radius[0] + 1
                y_low = y - seg_radius[1]
                y_high = y + seg_radius[1] + 1

                img_seg = proj[fld, x_low:x_high, y_low:y_high, :]
                # Even if there is no data augmentation, this call ensures
                # that the segment has the shape (num_images, x, y, channels)
                img_aug = DataAugmenter.augment_data(img_seg)
                aug_length = img_aug.shape[0]

                segment_batch.append(img_aug)
                label_batch.append(np.repeat(a=label, repeats=aug_length))

                write_batch = False
                # If this is the last segment then definitely write the batch
                last_seg = (np.all(point == sel_points[-1])) * \
                           (label == label_set[-1]) * \
                           (image_name == image_names[-1])
                if last_seg:
                    write_batch = True
                if len(segment_batch) >= write_batch_size:
                    write_batch = True

                if write_batch:
                    segment_batch = np.concatenate(segment_batch, axis=0)
                    label_batch = np.concatenate(label_batch, axis=0)
                    batch_length = segment_batch.shape[0]

                    with h5py.File(hdf5_fn, "r+") as h5handle:
                        h5handle["images"][
                            current_index:(current_index + batch_length),
                            ...] = segment_batch
                        h5handle["labels"][
                            current_index:(current_index + batch_length),
                            ...] = label_batch

                    # Reset the segment and label batches
                    segment_batch = []
                    label_batch = []

                    current_index += batch_length


def run_dnn_training(logfile=None):
    """
    Train a DNN
    :param logfile:
    :return:
    """
    Utils.write_log(msg="Train DNN", filename=logfile)

    # Load data
    with h5py.File(Config.DNN_TRAINING_SEGMENTS, "r") as h5handle:
        images = h5handle["images"][()]
        labels = h5handle["labels"][()]

    # Turn labels into one-hot vector
    labels = keras.utils.to_categorical(labels)

    # Scramble data
    scramble_index = np.random.choice(
        a=np.arange(0, len(images)),
        size=len(images), replace=False)
    images = images[scramble_index, ]
    labels = labels[scramble_index]

    # Split into training, validation, and testing dataset
    x_train, xx, y_train, yy = sklearn.model_selection.train_test_split(
        images, labels, test_size=0.5)

    x_val, x_test, y_val, y_test = sklearn.model_selection.train_test_split(
        xx, yy, test_size=0.5)

    # Train model
    model = build_model(image_size=images.shape[1:])
    model.compile(
        optimizer="rmsprop", loss="binary_crossentropy",
        metrics=["binary_accuracy"])
    checkpointer = keras.callbacks.ModelCheckpoint(
        filepath=os.path.join(
            Config.DNN_SEGMENTATION, "model_checkpoint_{epoch:02d}.h5"),
        verbose=1)
    input_dict_train = {
        "input_ch{}".format(ii): x_train[..., ii:ii + 1]
        for ii in range(x_train.shape[-1])}
    input_dict_val = {
        "input_ch{}".format(ii): x_val[..., ii:ii + 1]
        for ii in range(x_val.shape[-1])}
    model.fit(
        x=input_dict_train,
        y=y_train,
        epochs=Config.get_dnn_parameters()["training_epochs"],
        batch_size=Config.get_dnn_parameters()["training_batch_size"],
        validation_data=(input_dict_val, y_val),
        callbacks=[checkpointer])

    # Determine the best epoch to use based on validation accuracy
    # TODO(jan) Merge this into the *.fit call with a callback
    val_accs = []
    for ii in range(Config.get_dnn_parameters()["training_epochs"]):
        chk_fn = os.path.join(
            Config.DNN_SEGMENTATION,
            "model_checkpoint_{:02d}.h5".format(ii+1))
        val_model = keras.models.load_model(filepath=chk_fn)
        val_pred = val_model.predict(x=input_dict_val)
        val_pred_bin = np.argmax(val_pred, axis=-1)
        val_true_bin = np.argmax(y_val, axis=-1)
        val_accs.append(np.mean(val_pred_bin == val_true_bin))

    best_epoch = np.argmax(val_accs)
    best_fn = os.path.join(
        Config.DNN_SEGMENTATION,
        "model_checkpoint_{:02d}.h5".format(best_epoch+1))
    model = keras.models.load_model(filepath=best_fn)

    # Save weights
    model.save_weights(Config.DNN_WEIGHTS)

    # Get test statistics of model
    input_dict_test = {
        "input_ch{}".format(ii): x_test[..., ii:ii + 1]
        for ii in range(x_test.shape[-1])}
    y_pred = model.predict(x=input_dict_test)
    y_pred_bin = np.argmax(y_pred, axis=-1)
    y_true_bin = np.argmax(y_test, axis=-1)

    acc = np.mean(y_pred_bin == y_true_bin)
    cm = sklearn.metrics.confusion_matrix(
        y_true=y_true_bin, y_pred=y_pred_bin)
    with open(Config.DNN_ACCURACY, "w") as csv:
        csv.write("# Accuracy = {}\n".format(acc))
        csv.write("# Confusion Matrix\n")
        csv.write("#            | predicted is bg | predicted is fg\n")
        csv.write("# true is bg | {:>15} | {:>15}\n".format(cm[0, 0], cm[0, 1]))
        csv.write("# true is fg | {:>15} | {:>15}\n".format(cm[1, 0], cm[1, 1]))


def run_dnn_segmentation(low=None, high=None, logfile=None):
    """
    Applies the DNN to the project images and saves the outputs (both the
    probability maps and the binary masks)

    'logfile' dictates where the output should be written to, i.e. if
    'logfile == None' then output is directed to STDOUT (convenient for
    writing output to a cluster log file), otherwise output is written
    to the file given by 'logfile'

    If low == None then it is assumed to be 0
    If high == None then it is assumed to be the length of the image list

    DNNs are extremely good at parallelizing tasks, so a single call to the
    DNN is done with multiple images

    :param low: Low index of the images to process
    :param high: High index of the images to process (exclusive)
    :param logfile: (Optional) Path of logfile. Directs to STDOUT if None
    :return:
    """
    Utils.write_log(msg="Starting DNN segmentation", filename=logfile)

    # Set up image batches
    image_names = Config.get_image_names()
    if low is None:
        low = 0
    if high is None:
        high = len(image_names)
    image_names = image_names[low:high]
    batch_size = Config.get_dnn_parameters()["segmentation_batch_size"]
    image_batches = [
        image_names[i:i+batch_size] for i in
        range(0, len(image_names), batch_size)]

    # Set up model
    image_size = Config.get_image_size()
    model = build_model(image_size=image_size, weights=Config.DNN_WEIGHTS)

    for image_batch in image_batches:
        Utils.write_log(
            msg="Processing image batch {}".format(
                image_batches.index(image_batch)),
            filename=logfile)
        # If there are multiple fields then remember the image IDs for later grouping
        image_ids = []
        images = []
        for image_fn in image_batch:
            full_fn = os.path.join(Config.INPUT_DIR, image_fn)
            image = ImageIO.load_image(full_fn)
            image = np.stack([
                ImagePreprocessor.preprocess_image(field)
                for field in image])
            images.append(image)
            image_ids += [image_fn] * len(image)
        images = np.concatenate(images, axis=0)

        input_dict = {
            "input_ch{}".format(ii): images[..., ii:ii + 1]
            for ii in range(images.shape[-1])}
        res_raw = model.predict(x=input_dict)
        # Subtracting the maximum value from each label has no impact on the
        # softmax-ing, i.e. exp(x+C)/sum(exp(x+C)) == exp(x)/sum(exp(x)),
        # but prevents float-overflow errors.
        res = np.exp(res_raw - np.max(res_raw, axis=-1, keepdims=True))
        res /= np.sum(res, axis=-1)[..., None]
        res = res[..., 1]

        # Rescale back to the original image size
        scaling = [
            images.shape[ii] / res.shape[ii]
            for ii in range(len(res.shape))]
        res = scipy.ndimage.zoom(input=res, zoom=scaling, order=1)
        res_bin = (res > 0.5).astype(np.uint8)

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

        # Separate results back into individual fields
        image_ids = np.array(image_ids)
        res_by_images = []
        res_bin_by_images = []
        unique_image_ids = []
        for image_id in np.unique(image_ids):
            res_by_images.append(res[image_ids == image_id, ...])
            res_bin_by_images.append(res_bin[image_ids == image_id, ...])
            unique_image_ids += [image_id]

        # Save images in batch
        for ii in range(len(unique_image_ids)):
            fn = os.path.join(Config.DNN_SEGMENTATION, unique_image_ids[ii])
            if not os.path.isdir(os.path.dirname(fn)):
                os.makedirs(os.path.dirname(fn))
            ImageIO.save_dnn_segmentation(
                mask=res_bin_by_images[ii],
                probs=res_by_images[ii], fn=fn)


def run_dnn_featureextraction(low=None, high=None, logfile=None):
    """
    Applies the trained DNN to the project images but only computes the
    pre-classifier features and saves these to disk along with the
    input and predicted labels.

    'logfile' dictates where the output should be written to, i.e. if
    'logfile == None' then output is directed to STDOUT (convenient for
    writing output to a cluster log file), otherwise output is written
    to the file given by 'logfile'

    If low == None then it is assumed to be 0
    If high == None then it is assumed to be the length of the image list

    DNNs are extremely good at parallelizing tasks, so a single call to the
    DNN is done with multiple images

    :param low: Low index of the images to process
    :param high: High index of the images to process (exclusive)
    :param logfile: (Optional) Path of logfile. Directs to STDOUT if None
    :return:
    """
    Utils.write_log(msg="Starting DNN segmentation", filename=logfile)

    # Set up image batches
    image_names = Config.get_image_names()
    if low is None:
        low = 0
    if high is None:
        high = len(image_names)
    image_names = image_names[low:high]
    batch_size = Config.get_dnn_parameters()["segmentation_batch_size"]
    image_batches = [
        image_names[i:i+batch_size] for i in
        range(0, len(image_names), batch_size)]

    # Set up model
    image_size = Config.get_image_size()
    model = build_model(image_size=image_size, weights=Config.DNN_WEIGHTS)

    for image_batch in image_batches:
        Utils.write_log(
            msg="Processing image batch {}".format(
                image_batches.index(image_batch)),
            filename=logfile)
        # If there are multiple fields then remember the image IDs for later grouping
        image_ids = []
        images = []
        for image_fn in image_batch:
            full_fn = os.path.join(Config.INPUT_DIR, image_fn)
            image = ImageIO.load_image(full_fn)
            image = np.stack([
                ImagePreprocessor.preprocess_image(field)
                for field in image])
            images.append(image)
            image_ids += [image_fn] * len(image)
        images = np.concatenate(images, axis=0)

        # Calculate and process the label
        input_dict = {
            "input_ch{}".format(ii): images[..., ii:ii + 1]
            for ii in range(images.shape[-1])}
        res_raw = model.predict(x=input_dict)

        # Create intermediate model to capture features
        intermediate_model = keras.models.Model(
            inputs=model.inputs,
            outputs=model.layers[-2].output)

        res_features = intermediate_model.predict(x=input_dict)

        # Subtracting the maximum value from each label has no impact on the
        # softmax-ing, i.e. exp(x+C)/sum(exp(x+C)) == exp(x)/sum(exp(x)),
        # but prevents float-overflow errors.
        res = np.exp(res_raw - np.max(res_raw, axis=-1, keepdims=True))
        res /= np.sum(res, axis=-1)[..., None]
        res = res[..., 1]

        res_unscaled = (res > 0.5).astype(np.uint8)

        # Rescale back to the original image size
        scaling = [
            images.shape[ii] / res.shape[ii]
            for ii in range(len(res.shape))]
        res = scipy.ndimage.zoom(input=res, zoom=scaling, order=1)
        res_bin = (res > 0.5).astype(np.uint8)

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
        res_features = np.round(res_features, 4).astype(np.float16)

        # Separate results back into individual fields
        image_ids = np.array(image_ids)
        res_by_images = []
        res_bin_by_images = []
        res_unscaled_by_images = []
        res_features_by_images = []
        unique_image_ids = []
        for image_id in np.unique(image_ids):
            res_by_images.append(res[image_ids == image_id, ...])
            res_bin_by_images.append(res_bin[image_ids == image_id, ...])
            res_unscaled_by_images.append(
                res_unscaled[image_ids == image_id, ...])
            res_features_by_images.append(
                res_features[image_ids == image_id, ...])
            unique_image_ids += [image_id]

        # Save images and features in batch
        for ii in range(len(unique_image_ids)):
            # fn = os.path.join(Config.DNN_SEGMENTATION, unique_image_ids[ii])
            # if not os.path.isdir(os.path.dirname(fn)):
            #     os.makedirs(os.path.dirname(fn))
            # ImageIO.save_dnn_segmentation(
            #     mask=res_bin_by_images[ii],
            #     probs=res_by_images[ii], fn=fn)

            img_features = np.reshape(
                a=res_features_by_images[ii],
                newshape=(-1, res_features_by_images[ii].shape[-1]))
            img_labels = np.reshape(
                a=res_unscaled_by_images[ii],
                newshape=(-1,))

            fn = os.path.join(Config.DNN_FEATURES, unique_image_ids[ii])
            if not os.path.isdir(os.path.dirname(fn)):
                os.makedirs(os.path.dirname(fn))
            ImageIO.save_dnn_features(
                features=img_features,
                labels=img_labels, fn=fn)


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


def build_model_static(image_size, weights=None):
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

    :param image_size: Numpy array (x, y)
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

    input_ch1 = keras.layers.Input(
        shape=image_size + (1,), name="input_ch1")
    conv1_ch1 = keras.layers.Conv2D(
        filters=32, kernel_size=(3, 3), padding="valid",
        strides=(1, 1), name="conv1_ch1")(input_ch1)
    conv1_ch1 = keras.layers.LeakyReLU(
        alpha=0.001, name="conv1_ch1_lrelu")(conv1_ch1)
    maxpool1_ch1 = keras.layers.MaxPool2D(
        pool_size=(3, 3), strides=(2, 2),
        padding="valid", name="maxpool1_ch1")(conv1_ch1)
    conv2_ch1 = keras.layers.Conv2D(
        filters=64, kernel_size=(3, 3), padding="valid",
        strides=(1, 1), name="conv2_ch1")(maxpool1_ch1)
    conv2_ch1 = keras.layers.LeakyReLU(
        alpha=0.001, name="conv2_ch1_lrelu")(conv2_ch1)
    maxpool2_ch1 = keras.layers.MaxPool2D(
        pool_size=(3, 3), strides=(2, 2),
        padding="valid", name="maxpool2_ch1")(conv2_ch1)
    if weights is None:
        conv3_ch1 = keras.layers.Conv2D(
            filters=128, kernel_size=maxpool2_ch1.shape.as_list()[1:3],
            padding="valid", strides=(1, 1),
            name="conv3_ch1")(maxpool2_ch1)
    else:
        conv3_ch1 = keras.layers.Conv2D(
            filters=128, kernel_size=w["conv3_ch1"][0:2],
            padding="valid", strides=(1, 1),
            name="conv3_ch1")(maxpool2_ch1)
    conv3_ch1 = keras.layers.LeakyReLU(
        alpha=0.001, name="conv3_ch1_lrelu")(conv3_ch1)

    input_ch2 = keras.layers.Input(
        shape=image_size + (1,), name="input_ch2")
    conv1_ch2 = keras.layers.Conv2D(
        filters=32, kernel_size=(3, 3), padding="valid",
        strides=(1, 1), name="conv1_ch2")(input_ch2)
    conv1_ch2 = keras.layers.LeakyReLU(
        alpha=0.001, name="conv1_ch2_lrelu")(conv1_ch2)
    maxpool1_ch2 = keras.layers.MaxPool2D(
        pool_size=(3, 3), strides=(2, 2),
        padding="valid", name="maxpool1_ch2")(conv1_ch2)
    conv2_ch2 = keras.layers.Conv2D(
        filters=64, kernel_size=(3, 3), padding="valid",
        strides=(1, 1), name="conv2_ch2")(maxpool1_ch2)
    conv2_ch2 = keras.layers.LeakyReLU(
        alpha=0.001, name="conv2_ch2_lrelu")(conv2_ch2)
    maxpool2_ch2 = keras.layers.MaxPool2D(
        pool_size=(3, 3), strides=(2, 2),
        padding="valid", name="maxpool2_ch2")(conv2_ch2)
    if weights is None:
        conv3_ch2 = keras.layers.Conv2D(
            filters=128, kernel_size=maxpool2_ch2.shape.as_list()[1:3],
            padding="valid", strides=(1, 1),
            name="conv3_ch2")(maxpool2_ch2)
    else:
        conv3_ch2 = keras.layers.Conv2D(
            filters=128, kernel_size=w["conv3_ch2"][0:2],
            padding="valid", strides=(1, 1),
            name="conv3_ch2")(maxpool2_ch2)
    conv3_ch2 = keras.layers.LeakyReLU(
        alpha=0.001, name="conv3_ch2_lrelu")(conv3_ch2)

    input_ch3 = keras.layers.Input(
        shape=image_size + (1,), name="input_ch3")
    conv1_ch3 = keras.layers.Conv2D(
        filters=32, kernel_size=(3, 3), padding="valid",
        strides=(1, 1), name="conv1_ch3")(input_ch3)
    conv1_ch3 = keras.layers.LeakyReLU(
        alpha=0.001, name="conv1_ch3_lrelu")(conv1_ch3)
    maxpool1_ch3 = keras.layers.MaxPool2D(
        pool_size=(3, 3), strides=(2, 2),
        padding="valid", name="maxpool1_ch3")(conv1_ch3)
    conv2_ch3 = keras.layers.Conv2D(
        filters=64, kernel_size=(3, 3), padding="valid",
        strides=(1, 1), name="conv2_ch3")(maxpool1_ch3)
    conv2_ch3 = keras.layers.LeakyReLU(
        alpha=0.001, name="conv2_ch3_lrelu")(conv2_ch3)
    maxpool2_ch3 = keras.layers.MaxPool2D(
        pool_size=(3, 3), strides=(2, 2),
        padding="valid", name="maxpool2_ch3")(conv2_ch3)
    if weights is None:
        conv3_ch3 = keras.layers.Conv2D(
            filters=128, kernel_size=maxpool2_ch3.shape.as_list()[1:3],
            padding="valid", strides=(1, 1),
            name="conv3_ch3")(maxpool2_ch3)
    else:
        conv3_ch3 = keras.layers.Conv2D(
            filters=128, kernel_size=w["conv3_ch3"][0:2],
            padding="valid", strides=(1, 1),
            name="conv3_ch3")(maxpool2_ch3)
    conv3_ch3 = keras.layers.LeakyReLU(
        alpha=0.001, name="conv3_ch3_lrelu")(conv3_ch3)

    concat = keras.layers.Concatenate(
        name="concat_branches")([conv3_ch1, conv3_ch2, conv3_ch3])

    output = keras.layers.Conv2D(
        filters=2, kernel_size=(1, 1),
        padding="valid", strides=(1, 1),
        activation="linear", name="output_conv")(concat)

    if weights is None:
        output = keras.layers.Flatten()(output)

    model = keras.models.Model(
        inputs=[input_ch1, input_ch2, input_ch3],
        outputs=output)

    if weights is not None:
        model.load_weights(weights)

    return model
