"""
This is a utilities module that encapsulates all project-specific tasks as
well as commonly used functions, e.g.
 - Naming schemes of plates and wells and assignment to replicates
 - Extraction of images for specific wells/treatments
 - Robust z-score calculation
"""
import scipy.stats
import os
import re
import pandas as pd
import Config
import numpy as np
import h5py
import sklearn.decomposition
import scipy.spatial.distance


def get_rep_ids():
    rep_ids = {
        "D004T01P003L02": 1, "D004T01P004L03": 1, "D004T01P005L08": 1,
        "D004T01P007L02": 2, "D004T01P008L03": 2, "D004T01P009L08": 2,
        "D007T01P003L02": 1, "D007T01P004L03": 1, "D007T01P005L08": 1,
        "D007T01P007L02": 2, "D007T01P008L03": 2, "D007T01P009L08": 2,
        "D010T01P017L02": 1, "D010T01P018L03": 1, "D010T01P019L08": 1,
        "D010T01P021L02": 2, "D010T01P022L03": 2, "D010T01P023L08": 2,
        #"D013T01P001L02": 1, 
        "D013T01P002L03": 1, 
        "D013T01P003L08": 1,
        "D013T01P905L02": 2, "D013T01P906L03": 2, "D013T01P907L08": 2,
        "D018T01P901L02": 1, "D018T01P902L03": 1, "D018T01P003L08": 1,
        "D018T01P905L02": 2, "D018T01P906L03": 2, "D018T01P907L08": 2,
        "D019T01P001L02": 1, "D019T01P002L03": 1, "D019T01P003L08": 1,
        "D019T01P005L02": 2, "D019T01P006L03": 2, "D019T01P007L08": 2,
        "D020T01P001L02": 1, "D020T01P002L03": 1, "D020T01P003L08": 1,
        "D020T01P905L02": 2, 
        #"D020T01P906L03": 2,
        "D020T01P907L08": 2,
        "D020T02P009L02": 1, "D020T02P010L03": 1, "D020T02P011L08": 1,
        "D020T02P013L02": 2, "D020T02P014L03": 2, "D020T02P015L08": 2,
        #"D021T01P901L02": 1, "D021T01P902L03": 1, "D021T01P003L08": 1,
        #"D021T01P905L02": 2, "D021T01P906L03": 2, "D021T01P907L08": 2,
        "D022T01P001L02": 1, "D022T01P002L03": 1, "D022T01P003L08": 1,
        "D022T01P005L02": 2, "D022T01P906L03": 2, "D022T01P907L08": 2,
        "D027T01P001L02": 1, "D027T01P002L03": 1, "D027T01P003L08": 1,
        "D027T01P905L02": 2, 
        #"D027T01P906L03": 2, 
        "D027T01P907L08": 2,
        "D030T01P001L02": 1, "D030T01P002L03": 1, "D030T01P003L08": 1,
        "D030T01P905L02": 2, "D030T01P906L03": 2, "D030T01P907L08": 2 ,
        "D046T01P001L02": 1, "D046T01P002L03": 1, "D046T01P003L08": 1,
        "D046T01P005L02": 2, "D046T01P006L03": 2, "D046T01P007L08": 2 #,
        #"D052T01P001L08": 1, "D052T01P003L08": 2,
        #"D054T01P004L08": 1, "D054T01P006L08": 2,
        #"D055T01P007L02": 1, "D055T01P008L03": 1, "D055T01P009L08": 1,
        #"D055T01P011L02": 2, "D055T01P012L03": 2, "D055T01P013L08": 2
        }
    return rep_ids


def get_rep_ids_real():
    rep_ids = {
        "D004T01P003L02": 1, "D004T01P004L03": 1, "D004T01P005L08": 1,
        "D004T01P007L02": 2, "D004T01P008L03": 2, "D004T01P009L08": 2,
        "D007T01P003L02": 1, "D007T01P004L03": 1, "D007T01P005L08": 1,
        "D007T01P007L02": 2, "D007T01P008L03": 2, "D007T01P009L08": 2,
        "D010T01P017L02": 1, "D010T01P018L03": 1, "D010T01P019L08": 1,
        "D010T01P021L02": 2, "D010T01P022L03": 2, "D010T01P023L08": 2,
        #"D013T01P001L02": 1,
        "D013T01P002L03": 1, "D013T01P003L08": 1,
        "D013T01P905L02": 2, "D013T01P906L03": 2, "D013T01P907L08": 2,
        "D018T01P901L02": 1, "D018T01P902L03": 1, "D018T01P003L08": 1,
        "D018T01P905L02": 2, "D018T01P906L03": 2, "D018T01P907L08": 2,
        "D019T01P001L02": 1, "D019T01P002L03": 1, "D019T01P003L08": 1,
        "D019T01P005L02": 2, "D019T01P006L03": 2, "D019T01P007L08": 2,
        "D020T01P001L02": 1, "D020T01P002L03": 1, "D020T01P003L08": 1,
        "D020T01P905L02": 2, 
        #"D020T01P906L03": 2, 
        "D020T01P907L08": 2,
        "D020T02P009L02": 1, "D020T02P010L03": 1, "D020T02P011L08": 1,
        "D020T02P013L02": 2, "D020T02P014L03": 2, "D020T02P015L08": 2,
        #"D021T01P901L02": 1, "D021T01P902L03": 1, "D021T01P003L08": 1,
        #"D021T01P905L02": 2, "D021T01P906L03": 2, "D021T01P907L08": 2,
        "D022T01P001L02": 1, "D022T01P002L03": 1, "D022T01P003L08": 1,
        "D022T01P005L02": 2, "D022T01P906L03": 2, "D022T01P907L08": 2,
        "D027T01P001L02": 1, "D027T01P002L03": 1, "D027T01P003L08": 1,
        "D027T01P905L02": 2, 
        #"D027T01P906L03": 2, 
        "D027T01P907L08": 2,
        "D030T01P001L02": 1, "D030T01P002L03": 1, "D030T01P003L08": 1,
        "D030T01P905L02": 2, "D030T01P906L03": 2, "D030T01P907L08": 2,
        "D046T01P001L02": 1, "D046T01P002L03": 1, "D046T01P003L08": 1,
        "D046T01P005L02": 2, "D046T01P006L03": 2, "D046T01P007L08": 2,
        # "D052T01P001L08": 1, "D052T01P003L08": 2,
        # "D054T01P004L08": 1, "D054T01P006L08": 2,
        # "D055T01P007L02": 1, "D055T01P008L03": 1, "D055T01P009L08": 1,
        # "D055T01P011L02": 2, "D055T01P012L03": 2, "D055T01P013L08": 2,
        "M001A03P003L02": 1, "M001A03P004L03": 1, "M001A03P005L04": 1,
        "M001A03P006L05": 1, "M001A03P007L06": 1, "M001A03P008L07": 1,
        "M001A03P009L02": 2, "M001A03P010L03": 2, "M001A03P011L04": 2,
        "M001A03P012L05": 2, "M001A03P013L06": 2, "M001A03P014L07": 2,
        "M001B04P003L02": 1, "M001B04P004L03": 1, "M001B04P005L04": 1,
        "M001B04P006L05": 1, "M001B04P007L06": 1, "M001B04P008L07": 1,
        "M001B04P009L02": 2, "M001B04P010L03": 2, "M001B04P011L04": 2,
        "M001B04P012L05": 2, "M001B04P013L06": 2, "M001B04P014L07": 2,
        "M001K02P003L02": 1, "M001K02P004L03": 1, "M001K02P005L04": 1,
        "M001K02P006L05": 1, "M001K02P007L06": 1, "M001K02P008L07": 1,
        "M001K02P009L02": 2, "M001K02P010L03": 2, "M001K02P011L04": 2,
        "M001K02P012L05": 2, "M001K02P013L06": 2, "M001K02P014L07": 2,
        "M001W01P003L02": 1, "M001W01P004L03": 1, "M001W01P005L04": 1,
        "M001W01P006L05": 1, "M001W01P007L06": 1, "M001W01P908L07": 1,
        "M001W01P009L02": 2, "M001W01P010L03": 2, "M001W01P011L04": 2,
        "M001W01P012L05": 2, "M001W01P013L06": 2, "M001W01P014L07": 2}
    return rep_ids


def get_batch_ids():
    batch_ids = {
        "D004T01": 1, "D007T01": 1, "D010T01": 2, "D013T01": 1, "D018T01": 1,
        "D019T01": 1, "D020T01": 1, "D020T02": 2, "D021T01": 1, "D022T01": 1,
        "D027T01": 1, "D030T01": 1, "D046T01": 2, "D052T01": 2, "D054T01": 2,
        "D055T01": 2,
        "M001A03": 1, "M001B04": 1, "M001K02": 1, "M001W01": 1}
    return batch_ids


def get_species(line_or_plate_or_well):
    """
    A helper function to retrieve the species of well, plate, or line
    :param line_or_plate_or_well:
    :return:
    """
    if line_or_plate_or_well.startswith("D"):
        return "human"
    elif line_or_plate_or_well.startswith("M"):
        return "mouse"
    else:
        raise ValueError("Species not recognized!")


def get_replicate_id(plate):
    """
    Hard coded replicate IDs
    :param plate:
    :return:
    """

    rep_ids = get_rep_ids()

    if plate not in rep_ids.keys():
        raise ValueError(
            "No hard-coded replicate ID for plate {}".format(plate))

    return rep_ids[plate]


def get_batch_id(line):
    batch_ids = get_batch_ids()

    if line not in batch_ids.keys():
        raise ValueError(
            "No hard-coded replicate ID for line {}".format(line))

    return batch_ids[line]


def get_all_plates(species='all'):
    """
    Retrieve all plates, optionally only 'human' or 'mouse' plates
    :param species: String. One of ('all', 'human', 'mouse')
    :return:
    """

    plates = sorted(get_rep_ids().keys())

    if species == 'all':
        plates = [
            plate for plate in plates if
            plate.startswith("D0") or plate.startswith("M001")]
    elif species == 'human':
        plates = [plate for plate in plates if plate.startswith("D0")]
    elif species == 'mouse':
        plates = [plate for plate in plates if plate.startswith("M001")]
    else:
        raise ValueError("'species' must be one of ('all', 'human', 'mouse')")

    return plates


def get_all_lines(species='all'):
    """
    Retrieve all plates, optionally only 'human' or 'mouse' plates
    :param species: String. One of ('all', 'human', 'mouse')
    :return:
    """
    return sorted(set([
        s[0:7] for s in
        get_all_plates(species=species)]))


def get_plates_for_line(line):
    plates = sorted([
        plate for plate in sorted(get_rep_ids().keys())
        if plate.startswith(line)])
    return plates


def robust_z(xx, trim=0.1):
    """
    Calculate the robust z-score:
    z = (x - trim_mean(x)) / trim_sd(x)

    Ignores non-finite values.

    Use a trimmed mean here as the median absolute deviation would result
    in sparse features being thrown out even though they may be relevant
    to certain cell lines.

    :param xx:
    :param trim:
    :return:
    """
    xxt = scipy.stats.trimboth(xx, trim)
    xxm = np.nanmean(xxt)
    xxs = np.nanstd(xxt)
    return (xx - xxm) / xxs


def load_layout(layout_or_plate=None):
    """
    Loads the layout from file. It allows either a plate or a layout ID and
    only requires that the last three letters are in the format LXX
    :param layout_or_plate: A string.
    :return:
    """
    layout_id = layout_or_plate[-3:]
    if re.match("^L[0-9]{2}$", layout_id) is None:
        raise ValueError("'layout_or_plate' must end in LXX")
    layout_fn = os.path.join(Config.LAYOUTDIR, "{}.xlsx".format(layout_id))
    return pd.read_excel(layout_fn)


def get_target_and_pathway_for_drugs(metadata):
    """
    Adds the target and pathway to the metadata dataframe
    :param metadata:
    :return:
    """

    plates = metadata.Plate.unique().astype(np.str)
    lids = sorted(set([s[-3:] for s in plates]))

    layouts = []
    for lid in lids:
        layout = load_layout(lid[0:3])
        layout = layout.loc[:, ("Product.Name", "Target", "Pathway")]
        layouts.append(layout)
    layouts = pd.DataFrame(pd.concat(layouts))
    layouts = layouts.drop_duplicates()
    layouts.index = layouts["Product.Name"]
    del layouts["Product.Name"]

    # Remove entries with no Target or Pathway
    layouts = layouts.loc[layouts.isnull().sum(axis=1) == 0, :]
    metadata.loc[:, "Target"] = layouts.loc[
        metadata.Drug.values.astype(np.str), "Target"].values
    metadata.loc[:, "Pathway"] = layouts.loc[
        metadata.Drug.values.astype(np.str), "Pathway"].values

    return metadata


def make_drug_filename_safe(drugname):
    """
    This function ensures that the irregular drug names can be saved as file
    or folder names (i.e. removes illegal characters and spaces). Note that
    this function is Unix-specific and WILL NOT work on Windows.
    :param drugname:
    :return:
    """
    drugname = drugname.replace("/", "|")
    drugname = drugname.replace(" ", "_")
    return drugname


def get_characteristic_organoids(features, metadata, num_organoids=1):
    """
    This function calculates the mean across the given features and returns
    the most characteristic organoid(s), i.e. the organoids most similar to
    the sample average. The similarity is measured with the mahalanobis
    distance to the median/mad as the euclidean distance in such high-
    dimensional spaces becomes unintuitive.

    Features with ANY NaN entries are ignored.

    Location features are ignored so as to preserve location invariance.
    :param features:
    :param metadata:
    :param num_organoids:
    :return:
    """

    if num_organoids < 1:
        print("WARNING: Setting 'num_organoids' to 1")
        num_organoids = 1

    orig_features = features.copy()

    features = features.loc[:, np.sum(~np.isfinite(features), axis=0) == 0]
    location_features = [
        fname for fname in features.columns if
        re.search(pattern="cx|cy", string=fname)
        is not None]
    features = features.drop(labels=location_features, axis=1)

    # Use z-scaling so that each feature has the same weight
    features = features.apply(robust_z, trim=0)

    # Repeat removal of NaN in case any were generated by the robust z-calculation
    features = features.loc[:, np.sum(~np.isfinite(features), axis=0) == 0]

    # Perform PCA dimensionality reduction and select the top features
    pca = sklearn.decomposition.PCA()
    pca_trafo = pca.fit_transform(X=features.values)
    features = pd.DataFrame(
        data=pca_trafo,
        columns=["PC{}".format(ii+1) for ii in range(pca_trafo.shape[1])],
        index=features.index)
    features = features.loc[:, np.cumsum(pca.explained_variance_ratio_) < 0.9]

    median = features.median()

    # Inverse Covariance Matrix
    invcov = np.linalg.inv(features.cov())

    # Mahalanobis Distance
    dist = features.apply(
        axis=1,
        func=lambda xx: scipy.spatial.distance.mahalanobis(
            u=xx, v=median, VI=invcov),
        raw=True)

    # Euclidean Distance
    # dist = features.apply(
    #     axis=1,
    #     func=lambda xx: np.linalg.norm(xx - median),
    #     raw=True)

    min_organoid_features = orig_features.iloc[
        np.argsort(dist.values)[0:num_organoids], :]
    min_organoid_metadata = metadata.iloc[
        np.argsort(dist.values)[0:num_organoids], :]

    return min_organoid_features, min_organoid_metadata


def rename_h5_datasets(fname, name_mapping):
    """
    Renames an hdf5 dataset name. 'name_mapping' should be a dictionary that
    maps old names (keys) to new names (values)
    :param fname:
    :param name_mapping:
    :return:
    """
    with h5py.File(fname, "r+") as h5handle:
        for old, new in name_mapping.items():
            if old in h5handle.keys():
                h5handle[new] = h5handle[old]
                del h5handle[old]


def corr_coeff(a, b):
    """
    Calculate the correlation between a and b.
    :param a: pandas dataframe or numpy array with the shape
              (n_samples, n_features)
    :param b: pandas dataframe or numpy array with the shape
              (n_samples, n_features)
    :return:
    """
    if hasattr(a, "values"):
        a = a.values
    if hasattr(b, "values"):
        b = b.values

    a = a.T
    b = b.T
    a_ma = a - a.mean(1)[:, None]
    b_mb = b - b.mean(1)[:, None]

    # Sum of squares across rows
    ssa = (a_ma ** 2).sum(1)
    ssb = (b_mb ** 2).sum(1)

    # Finally get corr coeff
    return np.dot(a_ma, b_mb.T) / np.sqrt(np.dot(ssa[:, None], ssb[None]))
