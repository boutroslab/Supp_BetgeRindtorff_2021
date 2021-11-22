# These functions train a classifier to differentiate between drugs and DMSO
# organoids
# ---
# 1. Perform an incremental PCA on all organoids across all cell lines
# 2. Train Classifier to differentiate DMSO from Drug organoids
#    (drug/line-wise)
# 3. Remove any inactive drugs
# 4. Perform an incremental PCA on only organoids of active drugs across all
#    cell lines
# 5. Train classifier to differentiate DMSO from Drug organoids for only the
#    active drugs (drug/line-wise) and extract "cleaner" effect vectors
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
import sklearn.decomposition
import sklearn.calibration
import pickle


def normalize_features(features, metadata, zscore=True):
    """
    Normalize features to DMSO values and optionally z-scale them.

    NOTE: The division by the DMSO standard deviation has no effect
    if z-scaling occurs as well. In principle, 'zscore' determines
    if the scaling occurs by the entire dataset or just the DMSO
    controls.
    :param features:
    :param metadata:
    :param zscore:
    :return:
    """

    dmso_avg = np.stack(
        [features.loc[metadata.Drug == "DMSO", :].mean().values] *
        len(features))
    features_diff = features - dmso_avg

    if zscore:
        features_diff = (features_diff - features_diff.mean()) / \
            features_diff.std()
    else:
        dmso_var = np.stack(
            [features.loc[metadata.Drug == "DMSO", :].std().values] *
            len(features))
        features_diff /= dmso_var

    return features_diff


def scale_features(features, metadata, zscore=True):
    """
    Normalize features to DMSO values and optionally z-scale them.

    NOTE: The division by the DMSO standard deviation has no effect
    if z-scaling occurs as well. In principle, 'zscore' determines
    if the scaling occurs by the entire dataset or just the DMSO
    controls.
    :param features:
    :param metadata:
    :param zscore:
    :return:
    """

    avg = np.stack(
        [features.mean().values] *
        len(features))
    features_diff = features - avg

    if zscore:
        features_diff = (features_diff - features_diff.mean()) / \
            features_diff.std()
    else:
        var = np.stack(
            [features.std().values] *
            len(features))
        features_diff /= var

    return features_diff


def get_features(line):
    """
    Loads features
    :param line:

    :return:
    """
    features_object = OrganoidFeatures.load_line(
        line=line, normalize="normalized")
    features = features_object.get_features()
    metadata = features_object.get_metadata()

    # Rename drugs column to take concentrations into account. L08 drugs,
    # except for DMSO and Staurosporine, should have the format DRUG__CONC,
    # e.g. Bortezomib__0.02. This does not apply to DMSO and Staurosporine
    sel_indices = np.logical_and(
        metadata.Concentration != "nan",
        ~np.in1d(metadata.Drug.values, ("DMSO", "Staurosporine_500nM")))
    metadata.loc[sel_indices, "Drug"] = metadata.loc[sel_indices, "Drug"].str.cat(
        others=metadata.loc[sel_indices, "Concentration"],
        sep="__")

    return features, metadata


def get_common_features(species):
    """
    Goes through the whole dataset and records which features are not NaN
    for all lines
    :return:
    """
    out_fn = os.path.join(
        Config.DRUGEFFECTSDIR,
        "Common_Features_{}.txt".format(species))
    if os.path.isfile(out_fn):
        with open(out_fn, "r") as f:
            used_features = np.array([s.strip() for s in f.readlines()])
        return used_features

    print("Determining common features ...")

    lines = Utils.get_all_lines(species)
    used_features = None

    for line in lines:
        features, metadata = get_features(line)

        if used_features is None:
            used_features = features.columns.values[
                np.sum(~np.isfinite(features)) == 0]
        else:
            used_features = np.intersect1d(
                used_features, features.columns.values[
                    np.sum(~np.isfinite(features)) == 0])

        print(line, used_features.shape[0])

    with open(out_fn, "w") as f:
        f.writelines("\n".join(used_features))

    return used_features


def get_incremental_pca(species, n_components):
    """
    Incrementally performs a PCA on the entire dataset to find a common,
    reduced set of features.

    It explicitly uses only features shared by all lines

    :param species:
    :param n_components:
    :return:
    """

    out_fn = os.path.join(
        Config.DRUGEFFECTSDIR,
        "IncrementalPCA_{}_{}components.pkl".format(
            species, n_components))
    if os.path.isfile(out_fn):
        with open(out_fn, "rb") as f:
            return pickle.load(f)

    print("Calculating incremental PCA ...")

    lines = Utils.get_all_lines(species=species)
    used_features = get_common_features(species)

    pca = sklearn.decomposition.IncrementalPCA(
        n_components=n_components)

    for line in lines:
        print(line)
        features, metadata = get_features(line)

        # Use only the common features
        features = features.loc[:, used_features]

        # Normalize to DMSO
        features = normalize_features(
            features=features, metadata=metadata)

        # Run PCA
        pca.partial_fit(X=features)

    with open(out_fn, "wb") as f:
        pickle.dump(pca, f, pickle.HIGHEST_PROTOCOL)

    return pca


def get_incremental_pca_scaled(species, n_components):
    """
    Incrementally performs a PCA on the entire dataset to find a common,
    reduced set of features.

    It explicitly uses only features shared by all lines

    :param species:
    :param n_components:
    :return:
    """

    out_fn = os.path.join(
        Config.DRUGEFFECTSDIR,
        "IncrementalPCA_scaled_{}_{}components.pkl".format(
            species, n_components))
    if os.path.isfile(out_fn):
        with open(out_fn, "rb") as f:
            return pickle.load(f)

    print("Calculating scaled incremental PCA ...")

    lines = Utils.get_all_lines(species=species)
    used_features = get_common_features(species)

    pca = sklearn.decomposition.IncrementalPCA(
        n_components=n_components)

    for line in lines:
        print(line)
        features, metadata = get_features(line)

        # Use only the common features
        features = features.loc[:, used_features]

        # Scale Features instead of DMSO Normalization
        features = scale_features(
            features=features, metadata=metadata)

        # Run PCA
        pca.partial_fit(X=features)

    with open(out_fn, "wb") as f:
        pickle.dump(pca, f, pickle.HIGHEST_PROTOCOL)

    return pca


def get_transformed_features(species, n_components):
    """
    Apply a PCA to transform organoid features
    :param species:
    :param n_components:
    :return:
    """

    out_fn = os.path.join(
        Config.DRUGEFFECTSDIR, species,
        "TransformedFeatures_{}_{}components.h5".format(
            species, n_components))
    if os.path.isfile(out_fn):
        features, metadata = OrganoidFeatures.load_features(out_fn)
        return features, metadata

    print("Calculating transformed features ...")

    pca = get_incremental_pca(
        species=species, n_components=n_components)

    lines = Utils.get_all_lines(species=species)
    used_features = get_common_features(species)

    all_features = []
    all_metadata = []
    for line in lines:
        print(line)
        features, metadata = get_features(line)

        # Use only the common features
        features = features.loc[:, used_features]

        # Normalize to DMSO
        features = normalize_features(
            features=features, metadata=metadata)

        all_features.append(pca.transform(X=features))
        all_metadata.append(metadata)

    all_features = np.concatenate(all_features)
    all_metadata = pd.DataFrame(pd.concat(all_metadata))

    all_features = pd.DataFrame(
        data=all_features,
        columns=["PC%d" % (ii+1) for ii in range(all_features.shape[1])],
        index=all_metadata.index)

    OrganoidFeatures.save_features(
        features=all_features,
        metadata=all_metadata,
        filename=out_fn)

    return all_features, all_metadata

def get_transformed_features_scaled(species, n_components):
    """
    Apply a PCA to transform organoid features
    :param species:
    :param n_components:
    :return:
    """

    out_fn = os.path.join(
        Config.DRUGEFFECTSDIR, species,
        "TransformedFeatures_scaled_{}_{}components.h5".format(
            species, n_components))
    if os.path.isfile(out_fn):
        features, metadata = OrganoidFeatures.load_features(out_fn)
        return features, metadata

    print("Calculating transformed features ...")

    pca = get_incremental_pca_scaled(
        species=species, n_components=n_components)

    lines = Utils.get_all_lines(species=species)
    used_features = get_common_features(species)
    
    print(lines)

    all_features = []
    all_metadata = []
    for line in lines:
        print(line)
        features, metadata = get_features(line)

        # Use only the common features
        features = features.loc[:, used_features]

        # Normalize to DMSO
        features = scale_features(
            features=features, metadata=metadata)

        all_features.append(pca.transform(X=features))
        all_metadata.append(metadata)

    all_features = np.concatenate(all_features)
    all_metadata = pd.DataFrame(pd.concat(all_metadata))

    all_features = pd.DataFrame(
        data=all_features,
        columns=["PC%d" % (ii+1) for ii in range(all_features.shape[1])],
        index=all_metadata.index)

    OrganoidFeatures.save_features(
        features=all_features,
        metadata=all_metadata,
        filename=out_fn)

    return all_features, all_metadata



def train_classifier(species, n_components):
    """
    Trains a classifier to differentiate drugs from DMSO.

    This version uses the reduced features of the entire dataset. The dimension
    reduction is performed by the IncrementalPCA module

    :param species:
    :param n_components:
    :return:
    """

    features, metadata = get_transformed_features(
        species=species, n_components=n_components)

    lines = metadata.Line.unique()
    for line in lines:
        acc_out_fn = os.path.join(
            Config.DRUGEFFECTSDIR, species, line,
            "SVM_Accuracies_PCA_{}_{}components.csv".format(line, n_components))
        if os.path.isfile(acc_out_fn):
            print("Accuracies for {} already exists".format(line))
            continue

        profiles_mean_out_fn = os.path.join(
            Config.DRUGEFFECTSDIR, species, line,
            "SVM_Profiles_mean_PCA_{}_{}components.csv".format(line, n_components))
        if os.path.isfile(profiles_mean_out_fn):
            print("Profiles for {} already exists".format(line))
            continue

        profiles_std_out_fn = os.path.join(
            Config.DRUGEFFECTSDIR, species, line,
            "SVM_Profiles_std_PCA_{}_{}components.csv".format(line, n_components))
        if os.path.isfile(profiles_std_out_fn):
            print("Profiles for {} already exists".format(line))
            continue

        if not os.path.isdir(os.path.join(Config.DRUGEFFECTSDIR, species, line)):
            os.makedirs(os.path.join(Config.DRUGEFFECTSDIR, species, line))

        features_line = features.loc[metadata.Line == line, :]
        metadata_line = metadata.loc[metadata.Line == line, :]

        features_dmso = features_line.loc[metadata_line.Drug == "DMSO", :]

        drugs = metadata_line.Drug.unique()
        drug_accuracies = pd.DataFrame(
            data=None,
            columns=["AUC_Mean", "AUC_Std", "Distance"],
            index=drugs)
        drug_profiles_mean = pd.DataFrame(
            data=None,
            columns=features.columns.values,
            index=drugs)
        drug_profiles_std = pd.DataFrame(
            data=None,
            columns=features.columns.values,
            index=drugs)

        for drug in drugs:
            print(line, "|", np.where(drugs == drug)[0][0]+1, "of", len(drugs))
            features_drug = features_line.loc[metadata_line.Drug == drug, :]

            drug_dist = np.linalg.norm(
                features_drug.mean().values - features_dmso.mean().values)

            x = np.concatenate((features_dmso, features_drug), axis=0)
            y = np.repeat(("DMSO", "DRUG"), repeats=(
                len(features_dmso), len(features_drug)))

            # This version of the training uses a time-intensive parameter search
            # Train SVM
            x_param, x_train, y_param, y_train = \
                sklearn.model_selection.train_test_split(
                    x, y, test_size=0.5)
            param_grid = [{"C": [5e-4, 1e-4, 5e-3, 1e-3, 5e-2, 1e-2]}]
            svm_l1 = sklearn.model_selection.GridSearchCV(
                estimator=sklearn.linear_model.LogisticRegression(
                    penalty="l2", dual=False),
                param_grid=param_grid, cv=5, scoring="roc_auc")
            svm_l1.fit(x_param, y_param == "DRUG")
            cv = list(sklearn.model_selection.KFold(
                n_splits=10, shuffle=True).split(X=x_train, y=y_train))
            clf = sklearn.calibration.CalibratedClassifierCV(
                base_estimator=sklearn.linear_model.LogisticRegression(
                    penalty="l2", dual=False,
                    C=svm_l1.best_estimator_.get_params()["C"]), cv=cv)
            clf.fit(X=x_train, y=y_train)

            aucs = []
            profiles = []
            for ii in range(len(cv)):
                try:
                    train_i, test_i = cv[ii]
                    estimator = clf.calibrated_classifiers_[ii]
                    y_pred = estimator.predict_proba(X=x_train[test_i])
                    aucs.append(sklearn.metrics.roc_auc_score(
                        y_true=(y_train[test_i] == "DRUG").astype(np.int),
                        y_score=y_pred[:, 1]))
                    profiles.append(estimator.base_estimator.coef_)
                except ValueError:
                    # continue
                    pass

            drug_accuracies.loc[drug, :] = {
                "AUC_Mean": np.mean(aucs),
                "AUC_Std": np.std(aucs),
                "Distance": drug_dist}

            drug_profiles_mean.loc[drug, :] = np.mean(
                np.concatenate(profiles), axis=0)[None]

            drug_profiles_std.loc[drug, :] = np.std(
                np.concatenate(profiles), axis=0)[None]

        drug_accuracies.to_csv(acc_out_fn)
        drug_profiles_mean.to_csv(profiles_mean_out_fn)
        drug_profiles_std.to_csv(profiles_std_out_fn)


def train_classifier_by_replicate(species, n_components):
    """
    Trains a classifier to differentiate drugs from DMSO replicate-wise.

    This version uses the reduced features of the entire dataset. The dimension
    reduction is performed by the IncrementalPCA module

    :param species:
    :param n_components:
    :return:
    """

    features, metadata = get_transformed_features(
        species=species, n_components=n_components)

    lines = metadata.Line.unique()
    for line in lines:
        acc_out_fn = os.path.join(
            Config.DRUGEFFECTSDIR, species, line,
            "SVM_Accuracies_PCA_replicatewise_{}_{}components.csv".format(line, n_components))
        if os.path.isfile(acc_out_fn):
            print("Accuracies for {} already exists".format(line))
            continue

        profiles_mean_out_fn = os.path.join(
            Config.DRUGEFFECTSDIR, species, line,
            "SVM_Profiles_mean_PCA_replicatewise_{}_{}components.csv".format(line, n_components))
        if os.path.isfile(profiles_mean_out_fn):
            print("Profiles for {} already exists".format(line))
            continue

        profiles_std_out_fn = os.path.join(
            Config.DRUGEFFECTSDIR, species, line,
            "SVM_Profiles_std_PCA_replicatewise_{}_{}components.csv".format(line, n_components))
        if os.path.isfile(profiles_std_out_fn):
            print("Profiles for {} already exists".format(line))
            continue

        if not os.path.isdir(os.path.join(Config.DRUGEFFECTSDIR, species, line)):
            os.makedirs(os.path.join(Config.DRUGEFFECTSDIR, species, line))

        features_line = features.loc[metadata.Line == line, :]
        metadata_line = metadata.loc[metadata.Line == line, :]

        features_dmso = features_line.loc[metadata_line.Drug == "DMSO", :]
        metadata_dmso = metadata_line.loc[metadata_line.Drug == "DMSO", :]

        drugs = metadata_line.Drug.unique()
        replicates = metadata_line.Replicate.unique()
        drug_reps = np.array(
            [d + "___" + str(r) for d in drugs for r in replicates])
        drug_accuracies = pd.DataFrame(
            data=None,
            columns=["AUC_Mean", "AUC_Std", "Distance"],
            index=drug_reps)
        drug_profiles_mean = pd.DataFrame(
            data=None,
            columns=features.columns.values,
            index=drug_reps)
        drug_profiles_std = pd.DataFrame(
            data=None,
            columns=features.columns.values,
            index=drug_reps)

        for drug_rep in drug_reps:
            drug, rep = drug_rep.split("___")
            rep = int(rep)
            print(
                line, "|", np.where(drug_rep == drug_reps)[0][0]+1,
                "of", len(drug_reps))
            features_drug = features_line.loc[np.logical_and(
                metadata_line.Drug == drug, metadata_line.Replicate == rep), :]
            features_dmso_rep = features_dmso.loc[metadata_dmso.Replicate == rep, :]

            drug_dist = np.linalg.norm(
                features_drug.mean().values - features_dmso_rep.mean().values)

            x = np.concatenate((features_dmso_rep, features_drug), axis=0)
            y = np.repeat(("DMSO", "DRUG"), repeats=(
                len(features_dmso_rep), len(features_drug)))

            # This version of the training uses a time-intensive parameter search
            # Train SVM
            x_param, x_train, y_param, y_train = \
                sklearn.model_selection.train_test_split(
                    x, y, test_size=0.5)
            param_grid = [{"C": [5e-4, 1e-4, 5e-3, 1e-3, 5e-2, 1e-2]}]
            svm_l1 = sklearn.model_selection.GridSearchCV(
                estimator=sklearn.linear_model.LogisticRegression(
                    penalty="l2", dual=False),
                param_grid=param_grid, cv=5, scoring="roc_auc")
            svm_l1.fit(x_param, y_param == "DRUG")
            cv = list(sklearn.model_selection.KFold(
                n_splits=10, shuffle=True).split(X=x_train, y=y_train))
            clf = sklearn.calibration.CalibratedClassifierCV(
                base_estimator=sklearn.linear_model.LogisticRegression(
                    penalty="l2", dual=False,
                    C=svm_l1.best_estimator_.get_params()["C"]), cv=cv)
            clf.fit(X=x_train, y=y_train)

            aucs = []
            profiles = []
            for ii in range(len(cv)):
                try:
                    train_i, test_i = cv[ii]
                    estimator = clf.calibrated_classifiers_[ii]
                    y_pred = estimator.predict_proba(X=x_train[test_i])
                    aucs.append(sklearn.metrics.roc_auc_score(
                        y_true=(y_train[test_i] == "DRUG").astype(np.int),
                        y_score=y_pred[:, 1]))
                    profiles.append(estimator.base_estimator.coef_)
                except ValueError:
                    # continue
                    pass

            drug_accuracies.loc[drug_rep, :] = {
                "AUC_Mean": np.mean(aucs),
                "AUC_Std": np.std(aucs),
                "Distance": drug_dist}

            drug_profiles_mean.loc[drug_rep, :] = np.mean(
                np.concatenate(profiles), axis=0)[None]

            drug_profiles_std.loc[drug_rep, :] = np.std(
                np.concatenate(profiles), axis=0)[None]

        drug_accuracies.to_csv(acc_out_fn)
        drug_profiles_mean.to_csv(profiles_mean_out_fn)
        drug_profiles_std.to_csv(profiles_std_out_fn)


def train_classifier_with_pvalues(species, n_components):
    """
    Trains a classifier to differentiate drugs from DMSO.

    This function attempts to determine a p-value for the separation but this
    doesn't work so well due to computation time.

    This version uses the reduced features of the entire dataset. The dimension
    reduction is performed by the IncrementalPCA module

    :param species:
    :param n_components:
    :return:
    """

    features, metadata = get_transformed_features(
        species=species, n_components=n_components)

    lines = metadata.Line.unique()
    for line in lines:
        acc_out_fn = os.path.join(
            Config.DRUGEFFECTSDIR, species, line,
            "SVM_Accuracies_PCA_{}_{}components_pvalues.csv".format(
                line, n_components))
        if os.path.isfile(acc_out_fn):
            print("Accuracies for {} already exists".format(line))
            continue

        if not os.path.isdir(os.path.join(Config.DRUGEFFECTSDIR, species, line)):
            os.makedirs(os.path.join(Config.DRUGEFFECTSDIR, species, line))

        features_line = features.loc[metadata.Line == line, :]
        metadata_line = metadata.loc[metadata.Line == line, :]

        full_features_dmso = features_line.loc[metadata_line.Drug == "DMSO", :]

        drugs = metadata_line.Drug.unique()
        drug_accuracies = pd.DataFrame(
            data=None,
            columns=["AUC_Mean", "AUC_Std", "Distance", "PValue"],
            index=drugs)

        for drug in drugs:
            print(line, "|", np.where(drugs == drug)[0][0]+1, "of", len(drugs))
            features_drug = features_line.loc[metadata_line.Drug == drug, :]
            # features_dmso = full_features_dmso.sample(n=len(features_drug))
            features_dmso = full_features_dmso

            drug_dist = np.linalg.norm(
                features_drug.mean().values - full_features_dmso.mean().values)

            x = np.concatenate((features_dmso, features_drug), axis=0)
            y = np.repeat(("DMSO", "DRUG"), repeats=(
                len(features_dmso), len(features_drug)))

            # Perform parameter search before permutation
            x_param, x_train, y_param, y_train = \
                sklearn.model_selection.train_test_split(
                    x, y, test_size=0.5)
            param_grid = [{"C": [5e-4, 1e-4, 5e-3, 1e-3, 5e-2, 1e-2]}]
            svm_l1 = sklearn.model_selection.GridSearchCV(
                estimator=sklearn.linear_model.LogisticRegression(
                    penalty="l2", dual=False),
                param_grid=param_grid, cv=5, scoring="roc_auc")
            svm_l1.fit(x_param, y_param == "DRUG")
            cv = list(sklearn.model_selection.KFold(
                n_splits=10, shuffle=True).split(X=x_train, y=y_train))
            clf = sklearn.calibration.CalibratedClassifierCV(
                base_estimator=sklearn.linear_model.LogisticRegression(
                    penalty="l2", dual=False,
                    C=svm_l1.best_estimator_.get_params()["C"]), cv=cv)
            clf.fit(X=x_train, y=y_train)

            aucs = []
            for ii in range(len(cv)):
                try:
                    train_i, test_i = cv[ii]
                    estimator = clf.calibrated_classifiers_[ii]
                    y_pred = estimator.predict_proba(X=x_train[test_i])
                    aucs.append(sklearn.metrics.roc_auc_score(
                        y_true=(y_train[test_i] == "DRUG").astype(np.int),
                        y_score=y_pred[:, 1]))
                except ValueError:
                    continue

            real_auc = np.mean(aucs)
            real_auc_std = np.std(aucs)

            # Permute label and test the separation accuracy
            permutation_aucs = []
            while len(permutation_aucs) < 1000:
                y_permute = np.random.permutation(y_train)
                cv = list(sklearn.model_selection.KFold(
                    n_splits=10, shuffle=True).split(X=x_train, y=y_permute))
                clf = sklearn.calibration.CalibratedClassifierCV(
                    base_estimator=sklearn.linear_model.LogisticRegression(
                        penalty="l2", dual=False,
                        C=svm_l1.best_estimator_.get_params()["C"]), cv=cv)
                clf.fit(X=x_train, y=y_permute)

                aucs = []
                for ii in range(len(cv)):
                    try:
                        train_i, test_i = cv[ii]
                        estimator = clf.calibrated_classifiers_[ii]
                        y_pred = estimator.predict_proba(X=x_train[test_i])
                        aucs.append(sklearn.metrics.roc_auc_score(
                            y_true=(y_permute[test_i] == "DRUG").astype(np.int),
                            y_score=y_pred[:, 1]))
                    except ValueError:
                        continue

                permutation_aucs.append(np.mean(aucs))

            p_value = (np.sum(np.abs(permutation_aucs) >= real_auc) + 1) / \
                      (len(permutation_aucs) + 1)

            drug_accuracies.loc[drug, :] = {
                "AUC_Mean": real_auc,
                "AUC_Std": real_auc_std,
                "Distance": drug_dist,
                "PValue": p_value}

        drug_accuracies.to_csv(acc_out_fn)


def get_dmso_separation_accuracies(species, n_components):
    """
    Randomly splits the DMSO data of each line and attempts to separate them
    in order to determine the distriubtion of probabilities for false
    positives.
    :param species:
    :param n_components:
    :return:
    """
    dmso_acc_fn = os.path.join(
        Config.DRUGEFFECTSDIR, species,
        "SVM_Accuracies_DMSO_PCA_{}_{}components.csv".format(
            species, n_components))
    if os.path.isfile(dmso_acc_fn):
        print("DMSO Accuracies for {} already exists".format(species))
        return None

    if not os.path.isdir(os.path.join(Config.DRUGEFFECTSDIR, species)):
        os.makedirs(os.path.join(Config.DRUGEFFECTSDIR, species))

    n_splits = 100

    features, metadata = get_transformed_features(
        species=species, n_components=n_components)

    # Keep only DMSO features
    features = features.loc[metadata.Drug == "DMSO", :]
    metadata = metadata.loc[metadata.Drug == "DMSO", :]

    lines = metadata.Line.unique()
    accuracies = []
    for line in lines:
        features_line = features.loc[metadata.Line == line, :]

        y = np.repeat(a=("DMSO1", "DMSO2"), repeats=(
            len(features_line) // 2, len(features_line) - len(features_line) // 2))

        for counter in range(n_splits):
            print(line, "|", counter+1, "of", n_splits+1)
            x = features_line.sample(frac=1).values

            x_param, x_train, y_param, y_train = \
                sklearn.model_selection.train_test_split(
                    x, y, test_size=0.5)
            param_grid = [{"C": [5e-4, 1e-4, 5e-3, 1e-3, 5e-2, 1e-2]}]
            svm_l1 = sklearn.model_selection.GridSearchCV(
                estimator=sklearn.linear_model.LogisticRegression(
                    penalty="l2", dual=False),
                param_grid=param_grid, cv=5, scoring="roc_auc")
            svm_l1.fit(x_param, y_param == "DMSO2")
            cv = list(sklearn.model_selection.KFold(
                n_splits=10, shuffle=True).split(X=x_train, y=y_train))
            clf = sklearn.calibration.CalibratedClassifierCV(
                base_estimator=sklearn.linear_model.LogisticRegression(
                    penalty="l2", dual=False,
                    C=svm_l1.best_estimator_.get_params()["C"]), cv=cv)
            clf.fit(X=x_train, y=y_train)

            aucs = []
            for ii in range(len(cv)):
                train_i, test_i = cv[ii]
                estimator = clf.calibrated_classifiers_[ii]
                y_pred = estimator.predict_proba(X=x_train[test_i])
                aucs.append(sklearn.metrics.roc_auc_score(
                    y_true=(y_train[test_i] == "DMSO2").astype(np.int),
                    y_score=y_pred[:, 1]))

            accuracies.append((np.mean(aucs), np.std(aucs), line))

    accuracies = pd.DataFrame(
        data=np.stack(accuracies),
        columns=["AUC_Mean", "AUC_Std", "Line"])
    accuracies.to_csv(dmso_acc_fn)


def get_noise_separation_accuracies(n_components):
    """
    Generates random data to test the "false positive" separation AUC
    :param n_components:
    :return:
    """

    random_acc_fn = os.path.join(
        Config.DRUGEFFECTSDIR, "random_noise",
        "SVM_Accuracies_RandomPoints_{}components.csv".format(n_components))
    if os.path.isfile(random_acc_fn):
        print("Accuracies for random points already exist")
        return None

    if not os.path.isdir(os.path.join(Config.DRUGEFFECTSDIR, "random_noise")):
        os.makedirs(os.path.join(Config.DRUGEFFECTSDIR, "random_noise"))

    n_splits = 100
    n_samples = 10000

    y = np.repeat(a=("A", "B"), repeats=n_samples // 2)

    accuracies = []
    for counter in range(n_splits):
        print(counter)
        x = np.random.normal(size=(n_samples, n_components))

        x_param, x_train, y_param, y_train = \
            sklearn.model_selection.train_test_split(
                x, y, test_size=0.5)
        param_grid = [{"C": [5e-4, 1e-4, 5e-3, 1e-3, 5e-2, 1e-2]}]
        svm_l1 = sklearn.model_selection.GridSearchCV(
            estimator=sklearn.linear_model.LogisticRegression(
                penalty="l2", dual=False),
            param_grid=param_grid, cv=5, scoring="roc_auc")
        svm_l1.fit(x_param, y_param == "B")
        cv = list(sklearn.model_selection.KFold(
            n_splits=10, shuffle=True).split(X=x_train, y=y_train))
        clf = sklearn.calibration.CalibratedClassifierCV(
            base_estimator=sklearn.linear_model.LogisticRegression(
                penalty="l2", dual=False,
                C=svm_l1.best_estimator_.get_params()["C"]), cv=cv)
        clf.fit(X=x_train, y=y_train)

        aucs = []
        for ii in range(len(cv)):
            train_i, test_i = cv[ii]
            estimator = clf.calibrated_classifiers_[ii]
            y_pred = estimator.predict_proba(X=x_train[test_i])
            aucs.append(sklearn.metrics.roc_auc_score(
                y_true=(y_train[test_i] == "B").astype(np.int),
                y_score=y_pred[:, 1]))

        accuracies.append((np.mean(aucs), np.std(aucs)))

    accuracies = pd.DataFrame(
        data=np.stack(accuracies),
        columns=["AUC_Mean", "AUC_Std"])
    accuracies.to_csv(random_acc_fn)


def train(species, n_components):
    """
    The functions of this package all call each other in a pipeline fashion.
    This helper function exists as an API function so that external code
    doesn't need to understand the internal pipeline.

    :param species:
    :param n_components:
    :return:
    """

    train_classifier(species=species, n_components=n_components)
    # train_classifier_by_replicate(species=species, n_components=n_components)
