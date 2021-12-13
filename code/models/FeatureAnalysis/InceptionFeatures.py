# This module calculates inception features using the Keras package.

import keras
import ImageLoader
import ImagePreprocessor
import Config
import scipy.misc
import numpy as np
import os
import pandas as pd
import sys


QUEUE_TEMPLATE = """#!/bin/sh
#PBS -N {job_name}
#PBS -l walltime={wall_time}
#PBS -l select=1:ppn=1
#PBS -V

#PBS -M jan.sauer@dkfz-heidelberg.de
#PBS -m bea
#PBS -e {err_file}.err
#PBS -o {log_file}.log
{script}
"""


def get_model():
    """
    A wrapper function to return a keras model. This could later be altered
    to allow different types of pretrained networks or different parameter
    sets.
    :return:
    """
    model = keras.applications.inception_resnet_v2.InceptionResNetV2(
        include_top=False, weights='imagenet', pooling="max")
    return model


def load_image(well, field):
    """
    Load, preprocess, and resize an image.
    :param well:
    :param field:
    :return:
    """
    img = ImageLoader.load_projection(well=well, field_index=field)
    img = ImagePreprocessor.preprocess_image(image=img)
    img = scipy.misc.imresize(arr=img, size=(299, 299))

    # Renormalize img to [0,1] channel-wise
    img = img.astype(np.float_)
    img -= np.min(img, axis=(0, 1))
    img /= np.max(img, axis=(0, 1))
    return img


def split_image(img):
    """
    Microscope images aren't natural RGB images and the relative intensities
    are irrelevant. This function splits a single RGB image into three pseudo-
    RGB images (one for each channel) and adds a "batch" dimension.
    :param img:
    :return:
    """
    min_axis = np.argmin(img.shape)
    split_img = np.split(
        ary=img, indices_or_sections=img.shape[min_axis], axis=int(min_axis))
    split_img = [np.expand_dims(
        a=np.repeat(a=split_img[ii], repeats=3, axis=int(min_axis)),
        axis=0) for ii in range(len(split_img))]
    return split_img


def load_plate_images(plate, batch=None, num_batches=None):
    """
    Load all wells for a plate and prepare them for analysis with Keras. This
    way, Keras and Inception only need to be instantiated once. An index list
    is also returned to indicate the well, field, and channel of each array
    :param plate:
    :param batch:
    :param num_batches:
    :return:
    """

    img_dir = os.path.join(Config.IMAGEDIR, plate)
    wells = [
        "_".join(well.split("_")[0:3]) for
        well in sorted(os.listdir(img_dir))]
    batch_size = int(np.ceil(len(wells) / num_batches))
    batch_low = batch_size * batch
    batch_high = min(batch_size * (batch + 1), len(wells))
    wells = wells[batch_low:batch_high]
    # Load images
    imgs = []
    metadata = []
    for well in wells:
        for field in range(Config.NUMFIELDS):
            img = load_image(well=well, field=field)
            img = split_image(img=img)
            imgs += img
            metadata += tuple([
                (well, field, channel) for channel in range(len(img))])

    return np.concatenate(imgs), metadata


def calc_inception_features(plate, batch=0, num_batches=1):
    """
    Calculate the inception features.

    :param plate:
    :param batch: Default 0 for single-batch execution
    :param num_batches: Default 1 for single-batch execution
    :return:
    """

    if batch >= num_batches:
        raise ValueError("'batch' must be strictly less than 'num_batches'")

    out_fn = os.path.join(
        Config.FEATUREDIR, "inception", plate,
        "{}_features_inception_{}.csv".format(plate, batch))
    if os.path.isfile(out_fn):
        return None

    if not os.path.isdir(os.path.dirname(out_fn)):
        os.makedirs(os.path.dirname(out_fn))

    model = get_model()

    # Memory constraints make loading in batches easier
    imgs, labels = load_plate_images(
        plate=plate, batch=batch, num_batches=num_batches)

    features = pd.DataFrame(data=model.predict(imgs))
    metadata = pd.DataFrame(
        data=labels,
        columns=("Well", "Field", "Channel"))

    # Average the features across fields
    features = features.groupby((metadata.Well, metadata.Channel)).mean()

    # Concatenate the channels
    combined_features = []
    well_names = []
    for well in features.index.get_level_values(0).unique():
        combined_features.append(features.loc[well, :].values.reshape((-1, )))
        well_names.append(well)
    combined_features = np.stack(combined_features)
    combined_features = pd.DataFrame(
        data=combined_features,
        index=well_names,
        columns=[
            "F{}".format(ii+1) for ii
            in range(combined_features.shape[1])])

    combined_features.to_csv(out_fn)


def create_cluster_script(plate, num_batches=48):
    """
    A helper script to generate PBS/Torque-compatible scripts that can be
    submitted to a queue
    :param plate:
    :param num_batches: '48' is the magic number of cores currently active on
           my openstack cluster
    :return:
    """

    for job_id in range(num_batches):
        pbs_fn = os.path.join(
            Config.CLUSTERDIR, "PBS_INCEPTION_{}.sh".format(job_id))
        pyscript_fn = os.path.join(
            Config.CLUSTERDIR, "pyscript_INCEPTION_{}.sh".format(job_id))

        # Create python script
        with open(pyscript_fn, "w") as f:
            f.writelines([
                "cd %s\n" % Config.PYTHONDIR,
                "python InceptionFeatures.py {} {}".format(plate, job_id)])

        with open(pbs_fn, "w") as f:
            f.write(QUEUE_TEMPLATE.format(
                script=pyscript_fn,
                job_name="INCEPTION_{}".format(job_id),
                wall_time="48:00:00",
                err_file=os.path.join(
                    Config.CLUSTERDIR, "INCEPTION_{}_err".format(job_id)),
                log_file=os.path.join(
                    Config.CLUSTERDIR, "INCEPTION_{}_out".format(job_id)),
                working_dir=Config.PYTHONDIR))


if __name__ == "__main__":
    if len(sys.argv) == 2:
        calc_inception_features(plate=sys.argv[1])
    elif len(sys.argv) == 4:
        calc_inception_features(
            plate=sys.argv[1], batch=int(sys.argv[2]),
            num_batches=int(sys.argv[3]))
    else:
        print(
            "Usage:\n"
            "python InceptionFeatures.py PLATE\n"
            "-- or --\n"
            "python InceptionFeatures.py PLATE BATCH NUM_BATCHES\n")
