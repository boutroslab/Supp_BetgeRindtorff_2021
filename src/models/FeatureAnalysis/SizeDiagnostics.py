"""
This module contains methods for size diagnostics
"""
from OrganoidFeatures import OrganoidFeatures
import Utils
import Config
import os
import numpy as np
import h5py
import pandas as pd
import ImageLoader
import ImagePreprocessor
import scipy.misc


def get_object_sizes(species):
    """
    Return, or calculate, the object size of an entire data set
    :param species:
    :return:
    """

    out_fn = os.path.join(
        Config.SIZESTATISTICSDIR, "all_sizes_{}.h5".format(species))
    if os.path.isfile(out_fn):
        with h5py.File(out_fn, "r") as h5handle:
            size_features = pd.DataFrame(
                data=h5handle["data"][()],
                index=h5handle["wells"][()].astype(np.str),
                columns=h5handle["feature"][()].astype(np.str))
    else:
        plates = Utils.get_all_plates(species=species)
        size_features = OrganoidFeatures.load_single_feature(
            "x.0.s.area", wells=None, plates=plates, normalized=False)
        with h5py.File(out_fn, "w-") as h5handle:
            h5handle.create_dataset(
                name="feature", data=size_features.columns.values.astype(np.string_))
            h5handle.create_dataset(name="data", data=size_features.values)
            h5handle.create_dataset(
                name="wells", data=size_features.index.values.astype(np.string_))
    return size_features


def get_size_frequencies(sizes=None, species=None, delta=5, normalize=False):
    """
    Calculates the number of objects in a well within a certain size
    bandwidth, i.e. a value of (40, 200) for a well means that there are 200
    objects across the entire dataset with an area between 40**2 and
    (40 + delta)**2.

    :param sizes: A vector of sizes to analyze
    :param species: 'human' or 'mouse'
    :param delta: The bin size
    :param normalize: Boolean
    :return:
    """
    # 'sizes' takes precedence over 'species'
    if sizes is None and species is None:
        return None
    elif sizes is None and species is not None:
        sizes = np.sqrt(get_object_sizes(species))
    elif sizes is not None and species is None:
        pass
    elif sizes is not None and species is not None:
        pass

    # Getting a count table of multiples of delta is faster than looping and
    # summing over thresholds
    sizes_table = np.floor(sizes / delta)["x.0.s.area"].value_counts(
        sort=False).sort_index()
    x = sizes_table.index.values.astype(np.int) * delta
    y = sizes_table.values

    if normalize:
        y = y.astype(np.float)
        y /= np.sum(y)
    return x, y


def demo_function():
    """
    This function contains a demo and an explanation for choices. This should
    be added to the thesis and possibly a vignette later
    :return:
    """
    import BlurryOrganoidClassifier
    import ImageLoader
    import ImagePreprocessor
    import pickle
    import matplotlib.pyplot as plt

    """
    Objects with an area of < 300 pixels are predominately uninteresting and 
    can be filtered away from further analysis. These objects make up 
    approximately 10% of the data
    """
    sizes = get_object_sizes("human")

    # Do this with one cell line for demo purposes
    sel_indices = [
        ii for ii in range(len(sizes)) if
        sizes.index.values[ii].startswith("D004T01")]
    sizes = sizes.iloc[sel_indices, :]

    # Load focus stats
    plates = [s[0:14] for s in sizes.index.values]
    focus_stats = []
    processed_plates = []
    for plate in plates:
        if plate in processed_plates:
            continue
        processed_plates.append(plate)
        fn = os.path.join(
            Config.BLURRYORGANOIDCLSDIR,
            "BlurryOrganoidResults_{}_{}.pkl".format(
                plate, Config.FEATURETYPE))
        with open(fn, "rb") as f:
            focus_stats.append(pickle.load(f))
    focus_stats = pd.concat(focus_stats)

    # Look only at the sizes of in-focus objects
    sizes = sizes.loc[focus_stats["Focus"] == "good", :]
    x, y = get_size_frequencies(sizes=sizes, delta=1, normalize=True)
    plt.figure()
    plt.semilogx(x, np.cumsum(y))
    plt.axvline(x=300)
    plt.show()
    # plt.violinplot(np.log10(np.sqrt(sizes["x.0.s.area"].values)))

    """
    Here are some example images from two wells that show that these objects
    are largely irrelevant
    """

    # well = "D021T01P905L02_B_09"
    well = "D021T01P905L02_H_04"
    field = 4

    # Load image
    proj = ImageLoader.load_projection(well=well, field_index=field - 1)
    proj = ImagePreprocessor.preprocess_image(proj)
    mask = ImageLoader.load_segmentation(
        well=well, field_index=field - 1, probs=False)

    features_object = OrganoidFeatures.load_wells(
        wells=[well], normalize=False)
    features = features_object.get_features()
    metadata = features_object.get_metadata()

    focus_class = BlurryOrganoidClassifier.apply_classifier_to_features(
        features_object)

    index_vec = (metadata["Field"] == field).values
    features = features.loc[index_vec, :]
    metadata = metadata.loc[index_vec, :]
    focus_class = focus_class[index_vec]

    size_thresh = 300
    # Get the maximal radius
    max_radius = 0
    for organoid_id in range(len(features)):
        if focus_class.iloc[organoid_id, 0] != "good":
            continue
        y, x, rad, size = features.iloc[organoid_id][[
            "x.0.m.cx", "x.0.m.cy", "x.0.s.radius.max", "x.0.s.area"]].astype(int)
        if size < size_thresh:
            max_radius = max(max_radius, rad)

    # Get segments
    segs = []
    for organoid_id in range(len(features)):
        if focus_class.iloc[organoid_id, 0] != "good":
            continue
        y, x, size = features.iloc[organoid_id][[
            "x.0.m.cx", "x.0.m.cy", "x.0.s.area"]].astype(int)
        # if size <= 2000 or size >= 3000:
        if size > size_thresh:
            continue
        if x < max_radius or x >= (2048 - max_radius) or \
                y < max_radius or y >= (2048 - max_radius):
            continue

        xmin = max(0, x - max_radius)
        xmax = min(2048, x + max_radius)
        ymin = max(0, y - max_radius)
        ymax = min(2048, y + max_radius)
        s = proj[xmin:xmax, ymin:ymax, :]
        m = np.stack((
            mask[xmin:xmax, ymin:ymax],
            mask[xmin:xmax, ymin:ymax],
            mask[xmin:xmax, ymin:ymax]),
            axis=2).astype(np.float)

        c = s + 0.1*m
        c[c > 1] = 1
        segs.append(c)

    # Combine segments
    n_cols = np.ceil(np.sqrt(len(segs)))
    n_rows = np.ceil(len(segs) / n_cols)
    border = 5
    mosaic_matrix = np.ones(shape=(
        int(max_radius * 2 * n_rows + border * (n_rows - 1)),
        int(max_radius * 2* n_cols + border * (n_cols - 1)), 3))
    for ii in range(len(segs)):
        nr = int(np.floor(ii / n_cols))
        nc = int(ii - nr * n_cols)
        sx = (max_radius*2 + border) * nr
        sy = (max_radius*2 + border) * nc
        mosaic_matrix[sx:sx + max_radius*2, sy:sy + max_radius*2, :] = segs[ii]

    # plt.imshow(mosaic_matrix)
    return mosaic_matrix


def generate_example_images():
    """
    This function generates example flashcards for small objects in various
    size ranges
    :return:
    """

    cell_line = "D018T01"

    # Load data for cell line
    data = OrganoidFeatures.load_line(line=cell_line, normalize="raw")

    # Find all objects in certain size ranges and select random sample
    for min_size, max_size in (
            (0, 250), (250, 500), (500, 750), (750, 1000),
            (1000, 1500), (1500, 2000)):
        out_dir = os.path.join(
            Config.SIZESTATISTICSDIR,
            "segments_{}_{}".format(min_size, max_size))
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        sel_indices = np.logical_and(
            (data.get_metadata().loc[:, "Size"].values > min_size),
            (data.get_metadata().loc[:, "Size"].values <= max_size))
        relevant_objects = data.get_metadata().loc[sel_indices, :]
        max_radii = data.get_features().loc[sel_indices, "x.0.s.radius.max"]
        relevant_objects = pd.concat((relevant_objects, max_radii), axis=1)
        relevant_objects = relevant_objects.sample(50)

        # Extract image segments
        segments = []
        for ii in range(len(relevant_objects)):
            well = relevant_objects.index.values[ii]
            tmp = relevant_objects.iloc[ii, :]
            x, y, field, rad = tmp[
                ["OriginalX", "OriginalY", "Field", "x.0.s.radius.max"]]
            x = np.ceil(x).astype(np.int)
            y = np.ceil(y).astype(np.int)
            rad = np.ceil(rad + 5).astype(np.int)

            proj = ImageLoader.load_projection(
                well=well, field_index=field - 1)
            proj = ImagePreprocessor.preprocess_image(proj)
            # mask = ImageLoader.load_segmentation(
            #     well=well, field_index=field - 1, probs=False)

            xmin = max(0, x - rad)
            xmax = min(proj.shape[0], x + rad)
            ymin = max(0, y - rad)
            ymax = min(proj.shape[1], y + rad)
            # The x and y coordinate have to be flipped
            s = proj[ymin:ymax, xmin:xmax, :]
            # m = mask[ymin:ymax, xmin:xmax]
            # m = np.stack((m, m, m), axis=2)
            # segments.append(s + 0.2*m)
            scipy.misc.imsave(
                name=os.path.join(out_dir, "seg_{}.png".format(ii)), arr=s)
