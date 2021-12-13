# This script looks at line differences of DMSO wells
from OrganoidFeatures import OrganoidFeatures
import os
import Config
import Utils
import numpy as np
import pandas as pd
import sklearn.discriminant_analysis
import sklearn.ensemble
import sklearn.model_selection
import sklearn.svm
import sklearn.manifold
import pickle
#import keras.layers
#import keras.models


def get_drug_features(species, drug):
    """
    Either generates the drug features for all cell lines
    or loads them from file
    :param species:
    :param drug:
    :return:
    """

    drug_dirname = Utils.make_drug_filename_safe(drug)

    in_fn = os.path.join(
        Config.LINEDIFFERENCESDIR, species, drug_dirname,
        "Features_{}_{}.h5".format(species, drug_dirname))
    if os.path.isfile(in_fn):
        drug_feats, drug_md = OrganoidFeatures.load_features(in_fn)
    else:
        print("Generate {} drug features ...".format(species))
        lines = sorted(set([
            s[0:7] for s in
            Utils.get_all_plates(species=species)]))

        drug_feats = []
        drug_md = []
        for line in lines:
            print(" - {} ...".format(line))
            dat = OrganoidFeatures.load_line(
                line=line, normalize="normalized", drug=drug)
            drug_feats.append(dat.get_features().loc[
                dat.get_metadata().loc[:, "Drug"] == drug])
            drug_md.append(dat.get_metadata().loc[
                dat.get_metadata().loc[:, "Drug"] == drug])

        print("Combine arrays ...")
        drug_feats = pd.DataFrame(pd.concat(drug_feats))
        drug_md = pd.DataFrame(pd.concat(drug_md))

        print("Remove invalid features ...")
        drug_feats = drug_feats.loc[:, np.sum(~np.isfinite(drug_feats)) == 0]

        print("Save to disk ...")
        OrganoidFeatures.save_features(drug_feats, drug_md, in_fn)

    # # Calculate the z scores for all features
    # print("Calculating z-score ...")
    # dmso_feats = dmso_feats.loc[:, np.sum(~np.isfinite(dmso_feats)) == 0]
    # dmso_feats = dmso_feats.apply(robust_z)

    return drug_feats, drug_md
    


# def detect_outliers_with_autoencoder(features, encoding_dim):
#     """
#     Trains and applies an autoencoder to detect outliers. Returns the
#     (euclidean) distances of the reconstructed organoids to the
#     original features.
#     :param features:
#     :param encoding_dim: The size of the hidden layer
#     :return:
#     """
# 
#     # z-scale features
#     features = (features - features.mean()) / features.mad()
# 
#     # Set up network
#     input_layer = keras.layers.Input(shape=(features.shape[1],))
#     encoded = keras.layers.Dense(
#         units=encoding_dim, activation="relu")(input_layer)
#     decoded = keras.layers.Dense(
#         units=features.shape[1], activation="relu")(encoded)
#     autoencoder = keras.models.Model(inputs=input_layer, outputs=decoded)
#     encoder = keras.models.Model(inputs=input_layer, outputs=encoded)
#     encoded_input = keras.layers.Input(shape=(encoding_dim,))
#     decoder_layer = autoencoder.layers[-1]
#     decoder = keras.models.Model(
#         inputs=encoded_input, outputs=decoder_layer(encoded_input))
# 
#     autoencoder.compile(
#         optimizer="adam", loss="mean_squared_error")
# 
#     autoencoder.fit(
#         x=features.values, y=features.values, epochs=50,
#         batch_size=128, shuffle=True, verbose=0)
# 
#     xx_pred = decoder.predict(encoder.predict(features.values))
# 
#     dist = np.zeros(len(features.values))
#     for i, x in enumerate(features.values):
#         dist[i] = np.linalg.norm(x - xx_pred[i])
# 
#     return dist


def get_ranked_features(
        feats, label, key_features=None, cor_thresh=0.85):
    """
    Either calculate or load correlations from file and return them.

    Iteratively select features with the lowest correlation to the previously
    selected features, starting with the key features. Then, perform a PCA to
    rank features based on their importance

    :param feats:
    :param label:
    :param key_features:
    :param cor_thresh:
    :return:
    """

    if key_features is None:
        key_features = [
            'x.0.s.area', 'x.a.b.q05',
            'x.b.b.q05', 'x.c.b.q05']

    # Reduce the correlation matrix by iteratively selecting the least
    # correlating features with previously selected ones
    selected_features = np.array(key_features)
    remaining_features = feats.columns.drop(key_features).values
    # correlations = [0] * len(selected_features)
    allcors = Utils.corr_coeff(a=feats, b=feats)
    while len(remaining_features) > 0:
        cors = allcors[
            np.in1d(feats.columns.values, selected_features), :][
            :, np.in1d(feats.columns.values, remaining_features)]
        cors_max = np.max(np.abs(cors), axis=0)
        feature_to_keep = remaining_features[cors_max == cors_max.min()][0]
        # correlations.append(cors_max.min())
        selected_features = np.append(selected_features, feature_to_keep)
        # Remove selected feature as well as features that correlate too
        # strongly
        sel_indices = np.logical_and(
            remaining_features != feature_to_keep,
            cors_max < cor_thresh)
        remaining_features = remaining_features[sel_indices]

    selected_features = np.array(selected_features)
    # correlations = np.array(correlations)

    feats = feats.loc[:, selected_features]
    correlations = Utils.corr_coeff(a=feats, b=feats)
    np.fill_diagonal(correlations, np.nan)
    max_cors = np.nanmax(correlations, axis=0)

    # I haven't decided whether an unsupervised or a supervised method is best
    # for ranking feature importances, so I'm just gonna use three methods

    pca = sklearn.decomposition.PCA()
    pca.fit(X=feats)
    # The "importance" of each feature is calculated as its relative
    # contribution to each PC times the explained variance of that PC,
    # summed over all PCs.
    pca_components = np.abs(pca.components_)
    pca_importances = np.sum(
        pca_components /
        np.sum(pca_components, axis=1)[:, None] *
        pca.explained_variance_ratio_[:, None], axis=0)

    rf = sklearn.ensemble.RandomForestClassifier(n_estimators=25)
    rf.fit(X=feats.values, y=label)
    rf_importances = rf.feature_importances_

    lda = sklearn.discriminant_analysis.LinearDiscriminantAnalysis()
    lda.fit(X=feats.values, y=label)
    lda_coef = np.abs(lda.coef_)
    lda_importances = np.sum(
        lda_coef / np.sum(lda_coef, axis=1)[:, None],
        axis=0)

    importances = pd.DataFrame(
        data=np.stack(
            arrays=(pca_importances, lda_importances,
                    rf_importances, max_cors),
            axis=1),
        columns=["PCA_Importance", "LDA_Importance",
                 "RF_Importance", "Correlations"],
        index=feats.columns.values)

    return importances


def process_dmso_organoids(species):
    """
    This function first removes outlier organoids with an autoencoder and then
    calculates the correlation of features between
    :param species:
    :return:
    """
    results_dir = os.path.join(
        Config.LINEDIFFERENCESDIR, species,
        "DMSO", "results")
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)

    print("Load Features ...")
    dmso_feats, dmso_md = get_drug_features(species=species, drug="DMSO")

    print("Trim outliers with autoencoder ...")
    dists_fn = os.path.join(
        results_dir, "AutoEncoder_Organoid_Distances_{}.pkl".format(species))
    if os.path.isfile(dists_fn):
        with open(dists_fn, "rb") as f:
            dists = pickle.load(f)
    else:
        dists = {}
        for line in dmso_md.Line.unique():
            print(" - {}".format(line))
            dists[line] = detect_outliers_with_autoencoder(
                features=dmso_feats.loc[dmso_md.Line == line, :],
                encoding_dim=30)
        with open(dists_fn, "wb") as f:
            pickle.dump(dists, f, pickle.HIGHEST_PROTOCOL)

    trim_thresh = 80
    sel_indices = []
    # The somewhat bizarre loop is to ensure the selection indices are
    # generated in the right order.
    for line in dmso_md.groupby(dmso_md.Line).first().index.values:
        sel_indices.append(
            dists[line] <= np.percentile(dists[line], trim_thresh))
    sel_indices = np.concatenate(sel_indices)
    dmso_feats = dmso_feats.loc[sel_indices, :]
    dmso_md = dmso_md.loc[sel_indices, :]

    print("Calculate Z-Scores of Features ...")
    dmso_feats = (dmso_feats - dmso_feats.mean()) / dmso_feats.std()

    print("Perform PCA")
    pca = sklearn.decomposition.PCA()
    dmso_feats_pca = pd.DataFrame(
        data=pca.fit_transform(X=dmso_feats),
        columns=["PC{}".format(ii+1) for ii in
                 range(dmso_feats.shape[1])],
        index=dmso_feats.index)

    # Save original and PCA features
    feats_fn = os.path.join(
        results_dir, "ReducedFeaturesPCA_{}_{}.h5".format(species, "DMSO"))
    OrganoidFeatures.save_features(
        features=dmso_feats_pca.iloc[:, 0:50],
        metadata=dmso_md, filename=feats_fn)

    feats_fn = os.path.join(
        results_dir, "ReducedFeatures_{}_{}.h5".format(species, "DMSO"))
    OrganoidFeatures.save_features(
        features=dmso_feats,
        metadata=dmso_md, filename=feats_fn)

## Below are two modifications without filters or with a LDC based filter.

def process_dmso_organoids_no_filter(species):
    """
    This function does not remove outlier organoids. 

    :param species:
    :return:
    """
    results_dir = os.path.join(
        Config.LINEDIFFERENCESDIR, species,
        "DMSO", "results")
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)

    print("Load Features ...")
    dmso_feats, dmso_md = get_drug_features(species=species, drug="DMSO")

    # print("Trim outliers with autoencoder ...")
    # dists_fn = os.path.join(
    #     results_dir, "AutoEncoder_Organoid_Distances_{}.pkl".format(species))
    # if os.path.isfile(dists_fn):
    #     with open(dists_fn, "rb") as f:
    #         dists = pickle.load(f)
    # else:
    #     dists = {}
    #     for line in dmso_md.Line.unique():
    #         print(" - {}".format(line))
    #         dists[line] = detect_outliers_with_autoencoder(
    #             features=dmso_feats.loc[dmso_md.Line == line, :],
    #             encoding_dim=30)
    #     with open(dists_fn, "wb") as f:
    #         pickle.dump(dists, f, pickle.HIGHEST_PROTOCOL)
    # 
    # trim_thresh = 80
    # sel_indices = []
    # # The somewhat bizarre loop is to ensure the selection indices are
    # # generated in the right order.
    # for line in dmso_md.groupby(dmso_md.Line).first().index.values:
    #     sel_indices.append(
    #         dists[line] <= np.percentile(dists[line], trim_thresh))
    # sel_indices = np.concatenate(sel_indices)
    # dmso_feats = dmso_feats.loc[sel_indices, :]
    # dmso_md = dmso_md.loc[sel_indices, :]

    print("Calculate Z-Scores of Features ...")
    dmso_feats = (dmso_feats - dmso_feats.mean()) / dmso_feats.std()

    print("Perform PCA")
    pca = sklearn.decomposition.PCA()
    dmso_feats_pca = pd.DataFrame(
        data=pca.fit_transform(X=dmso_feats),
        columns=["PC{}".format(ii+1) for ii in
                 range(dmso_feats.shape[1])],
        index=dmso_feats.index)

    # Save original and PCA features
    feats_fn = os.path.join(
        results_dir, "ReducedFeaturesPCA_noaa_{}_{}.h5".format(species, "DMSO"))
    OrganoidFeatures.save_features(
        features=dmso_feats_pca.iloc[:, 0:50],
        metadata=dmso_md, filename=feats_fn)

    feats_fn = os.path.join(
        results_dir, "ReducedFeatures_noaa_{}_{}.h5".format(species, "DMSO"))
    OrganoidFeatures.save_features(
        features=dmso_feats,
        metadata=dmso_md, filename=feats_fn)

# helper function to for process_dmso_organoids_ldc
def aggregate_ldc_logs(species):
    """
    This function aggregates all LDC logs for a species for viability-based filtering
    
    :param species:
    :return:
    """
    # aggregating logs from directory
    results_dir = os.path.join(
        Config.DEADORGANOIDCLASSIFIERDIR, species, "results")
    
    #for f in glob('{}/somefile*.csv'.format(results_dir)):
    #df = pd.read_csv(f)


def process_dmso_organoids_ldc(species):
    """
    This function does not remove outlier organoids. 

    :param species:
    :return:
    """
    results_dir = os.path.join(
        Config.LINEDIFFERENCESDIR, species,
        "DMSO", "results")
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)

    print("Load Features ...")
    dmso_feats, dmso_md = get_drug_features(species=species, drug="DMSO")

    print("Filter outliers with LDC")
    # dists_fn = os.path.join(
    #     results_dir, "AutoEncoder_Organoid_Distances_{}.pkl".format(species))
    # if os.path.isfile(dists_fn):
    #     with open(dists_fn, "rb") as f:
    #         dists = pickle.load(f)
    # else:
    #     dists = {}
    #     for line in dmso_md.Line.unique():
    #         print(" - {}".format(line))
    #         dists[line] = detect_outliers_with_autoencoder(
    #             features=dmso_feats.loc[dmso_md.Line == line, :],
    #             encoding_dim=30)
    #     with open(dists_fn, "wb") as f:
    #         pickle.dump(dists, f, pickle.HIGHEST_PROTOCOL)
    # 
    # trim_thresh = 80
    # sel_indices = []
    # # The somewhat bizarre loop is to ensure the selection indices are
    # # generated in the right order.
    # for line in dmso_md.groupby(dmso_md.Line).first().index.values:
    #     sel_indices.append(
    #         dists[line] <= np.percentile(dists[line], trim_thresh))
    # sel_indices = np.concatenate(sel_indices)
    # dmso_feats = dmso_feats.loc[sel_indices, :]
    # dmso_md = dmso_md.loc[sel_indices, :]

    print("Calculate Z-Scores of Features ...")
    dmso_feats = (dmso_feats - dmso_feats.mean()) / dmso_feats.std()

    print("Perform PCA")
    pca = sklearn.decomposition.PCA()
    dmso_feats_pca = pd.DataFrame(
        data=pca.fit_transform(X=dmso_feats),
        columns=["PC{}".format(ii+1) for ii in
                 range(dmso_feats.shape[1])],
        index=dmso_feats.index)

    # Save original and PCA features
    feats_fn = os.path.join(
        results_dir, "ReducedFeaturesPCA_noaa_{}_{}.h5".format(species, "DMSO"))
    OrganoidFeatures.save_features(
        features=dmso_feats_pca.iloc[:, 0:50],
        metadata=dmso_md, filename=feats_fn)

    feats_fn = os.path.join(
        results_dir, "ReducedFeatures_noaa_{}_{}.h5".format(species, "DMSO"))
    OrganoidFeatures.save_features(
        features=dmso_feats,
        metadata=dmso_md, filename=feats_fn)


def get_all_features(species):
    """
    Either generates the features for all cell lines
    or loads them from file
    :param species:
    :return:
    """

    drug_dirname = Utils.make_drug_filename_safe("all_drugs")

    in_fn = os.path.join(
        Config.LINEDIFFERENCESDIR, species, drug_dirname,
        "Features_{}_{}.h5".format(species, drug_dirname))
    if os.path.isfile(in_fn):
        drug_feats, drug_md = OrganoidFeatures.load_features(in_fn)
    else:
        print("Generate {} drug features ...".format(species))
        lines = sorted(set([
            s[0:7] for s in
            Utils.get_all_plates(species=species)]))

        drug_feats = []
        drug_md = []
        for line in lines:
            print(" - {} ...".format(line))
            dat = OrganoidFeatures.load_line(
                line=line, normalize="normalized")
            drug_feats.append(dat.get_features())
            drug_md.append(dat.get_metadata())

        print("Combine arrays ...")
        drug_feats = pd.DataFrame(pd.concat(drug_feats))
        drug_md = pd.DataFrame(pd.concat(drug_md))

        print("Remove invalid features ...")
        drug_feats = drug_feats.loc[:, np.sum(~np.isfinite(drug_feats)) == 0]

        print("Save to disk ...")
        OrganoidFeatures.save_features(drug_feats, drug_md, in_fn)

    # # Calculate the z scores for all features
    # print("Calculating z-score ...")
    # dmso_feats = dmso_feats.loc[:, np.sum(~np.isfinite(dmso_feats)) == 0]
    # dmso_feats = dmso_feats.apply(robust_z)

    return drug_feats, drug_md


def process_all_organoids(species):
    """
    This function does **not** remove outlier organoids with an autoencoder. 
    It runs this step of feature processing using a high-mem instance.
    :param species:
    :return:
    """
    results_dir = os.path.join(
        Config.LINEDIFFERENCESDIR, species,
        "all_drugs", "results")
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)

    print("Load Features ...")
    dmso_feats, dmso_md = get_all_features(species=species)

    print("Calculate Z-Scores of Features ...")
    dmso_feats = (dmso_feats - dmso_feats.mean()) / dmso_feats.std()

    print("Perform PCA")
    pca = sklearn.decomposition.PCA()
    dmso_feats_pca = pd.DataFrame(
        data=pca.fit_transform(X=dmso_feats),
        columns=["PC{}".format(ii+1) for ii in
                 range(dmso_feats.shape[1])],
        index=dmso_feats.index)
        
    print("preparing eigenvalues")
    eigenvector = pca.components_.T 
    eigenvalue = pca.explained_variance_
    
    # Hacking my way to export the PCA components
    feats_fn = os.path.join(
        results_dir, "ReducedFeaturesPCA_all_drugs_components{}.csv".format(species))
    pd.DataFrame(eigenvector).to_csv(feats_fn)
        
    feats_fn = os.path.join(
        results_dir, "ReducedFeaturesPCA_all_drugs_variance{}.csv".format(species))
    pd.DataFrame(eigenvalue).to_csv(feats_fn)


    # Save original and PCA features
    feats_fn = os.path.join(
        results_dir, "ReducedFeaturesPCA_all_drugs_{}.h5".format(species))
    OrganoidFeatures.save_features(
        features=dmso_feats_pca.iloc[:, 0:50],
        metadata=dmso_md, filename=feats_fn)

    feats_fn = os.path.join(
        results_dir, "ReducedFeatures_all_drugs_{}.h5".format(species))
    OrganoidFeatures.save_features(
        features=dmso_feats,
        metadata=dmso_md, filename=feats_fn)
