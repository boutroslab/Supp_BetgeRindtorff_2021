# This script learns to differentiate between dead (positive controls) and
# live (negative controls) organoids.
import os
import Config
import pickle
import Utils
import numpy as np
import sklearn.model_selection
import sklearn.ensemble
import sklearn.svm
import scipy.stats
import pandas as pd
import sklearn.linear_model
import sklearn.metrics
from OrganoidFeatures import OrganoidFeatures


CHANNELCOMBOS = (None, ("Actin",), ("DAPI",), ("FITC",), ("Actin", "DAPI"))


def trim_func(func1d, x, p):
    xx = scipy.stats.trimboth(x, p)
    return func1d(xx)


def get_trimmed_summary(x, p=0.1):
    """
    :param x:
    :param p:
    :return:
    """
    scale_mean = np.apply_along_axis(
        func1d=lambda xx: trim_func(np.nanmean, xx, p),
        axis=0, arr=x)
    scale_std = np.apply_along_axis(
        func1d=lambda xx: trim_func(np.nanstd, xx, p),
        axis=0, arr=x)
    return scale_mean, scale_std


def load_features(cell_line, channels):
    """
    Load features for a cell line, filter channels, and z-normalize

    I load the entire line here to ensure the normalization is performed
    over the entire cell line's data

    :param cell_line:
    :param channels:
    :return:
    """
    features_object = OrganoidFeatures.load_line(
        line=cell_line, normalize="normalized")
    features = features_object.get_features()
    feature_names = features.columns.values
    metadata = features_object.get_metadata()

    # Keep only features related to the selected channels
    if channels is not None:
        channel_dict = {"Actin": "a", "FITC": "b", "DAPI": "c"}
        valid_fch = np.array(
            [channel_dict[ch] for ch in channels] + ["B"])
        feature_names = feature_names[np.array([
            np.all(np.in1d(list(s.split(".")[1]), valid_fch))
            for s in features.columns.values])]
        # Additionally, keep only texture/intensity features
        # feature_names = feature_names[np.in1d(np.array(
        #     [s.split(".")[2] for s in feature_names]), ("h", "b"))]
        features = features.loc[:, feature_names]

    # Normalize by DMSO values
    dmso_features = features.loc[metadata.Drug == "DMSO", :]
    dmso_mean, dmso_std = get_trimmed_summary(dmso_features.values)
    features = (features - dmso_mean) / dmso_std

    # To allow for transfer learning and keep the number of features consistent,
    # I set the features that have NaN entries to 0 so that the RFC isn't trained
    # on them
    features.loc[:, np.sum(~np.isfinite(features), axis=0) != 0] = 0

    return features, metadata, feature_names


def train_classifier(cell_line, channels=None):
    """
    Train a classifier to learn the difference between a live and a dead
    organoid. The live and dead organoids are taken from "control" wells.

    Positive controls are, as determined visually and via initial feature
    analysis, Bortezomib and Irinotecan / SN-38 at concentrations of 0.2
    and 1.0. Negative controls are DMSO wells.

    The classifier can optionally be trained to limit itself to features that
    only relate to specific channels.

    :param cell_line:
    :param channels: An iterable containing any combination of
           ("Actin", "DAPI", "FITC")
    :return: A tuple (classifier, accuracy, features used) or
        (classifier, accuracy, features used, data_x, data_y) if
        return_data == True
    """
    if cell_line not in Utils.get_all_lines("human"):
        raise ValueError(
            "'cell_line' must be one of {}".format(Utils.get_all_lines()))

    if channels is not None:
        if not np.all(np.in1d(channels, ("Actin", "DAPI", "FITC"))):
            raise ValueError(
                "'channels' must be an iterable containing any combination "
                "of the values ('Actin', 'DAPI', 'FITC')")

    if np.sum(np.in1d(channels, ("Actin", "DAPI", "FITC"))) == 3:
        channels = None

    if channels is None:
        clf_fn = os.path.join(
            Config.DEADORGANOIDCLASSIFIERDIR, "classifiers",
            "OrganoidClassifier_%s.pkl" % cell_line)
    else:
        clf_fn = os.path.join(
            Config.DEADORGANOIDCLASSIFIERDIR, "classifiers",
            "OrganoidClassifier_{}__{}.pkl".format(
                cell_line, "_".join(sorted(channels))))
    if not os.path.isdir(os.path.dirname(clf_fn)):
        os.makedirs(os.path.dirname(clf_fn))

    if os.path.isfile(clf_fn):
        with open(clf_fn, "rb") as f:
            return pickle.load(f)

    # Select plates
    plates = sorted(
        [plate for plate in Utils.get_plates_for_line(cell_line)
         if plate.startswith(cell_line) if plate.endswith("L08")],
        key=lambda xx: xx[9:14])

    # Load features
    features, metadata, feature_names = load_features(
        cell_line=cell_line, channels=channels)

    # From here, keep only the L08 organoids for training
    sel_indices = np.in1d(metadata.Plate.values, plates)
    features = features.loc[sel_indices, :]
    metadata = metadata.loc[sel_indices, :]

    # Select negative (DMSO) controls
    neg_ctrls = features.loc[metadata.Drug == "DMSO", :]

    # Select positive controls:
    # Bortezomib and Irinotecan / SN-38 with a concentration of 0.2 or 1.0
    pos_ctrls = features.loc[np.logical_and(
        np.in1d(metadata.Drug, ("Bortezomib", "Irinotecan / SN-38")),
        np.in1d(metadata.Concentration, ("1.0", "0.2"))), :]

    # Create training data (sample to equalize groups)
    num_samples = min(pos_ctrls.shape[0], neg_ctrls.shape[0])
    x = np.concatenate((
        pos_ctrls.sample(num_samples).values,
        neg_ctrls.sample(num_samples).values), axis=0)
    y = np.repeat(("POS", "NEG"), num_samples)

    x_train, x_test, y_train, y_test = sklearn.model_selection.train_test_split(
        x, y, test_size=0.40)

    # Run classifier
    clf = sklearn.ensemble.RandomForestClassifier(n_estimators=10)
    clf.fit(x_train, y_train)
    clf_acc = clf.score(x_test, y_test)

    with open(clf_fn, "wb") as f:
        pickle.dump((
            clf, clf_acc, feature_names, x_test, y_test),
            f, pickle.HIGHEST_PROTOCOL)

    return clf, clf_acc, feature_names, x_test, y_test


def run_classifier(cell_line, channels=None):
    """
    Run the trained classifier on a cell line

    The classifier can optionally be run on features that
    only relate to specific channels.
    :param cell_line:
    :param channels:
    :return:
    """
    if cell_line not in Utils.get_all_lines("human"):
        raise ValueError(
            "'cell_line' must be one of {}".format(Utils.get_all_lines()))

    if channels is not None:
        if not np.all(np.in1d(channels, ("Actin", "DAPI", "FITC"))):
            raise ValueError(
                "'channels' must be an iterable containing any combination "
                "of the values ('Actin', 'DAPI', 'FITC')")

    if np.sum(np.in1d(channels, ("Actin", "DAPI", "FITC"))) == 3:
        channels = None

    clf, clf_acc, feature_names, x_test, y_test = train_classifier(
        cell_line=cell_line, channels=channels)

    features, metadata, feature_names = load_features(
        cell_line=cell_line, channels=channels)

    prediction = clf.predict(features.values) == "POS"
    prediction_prob = clf.predict_proba(features.values)

    results = pd.DataFrame(
        data=np.concatenate((prediction[:, None], prediction_prob), axis=1),
        columns=["Prediction", "ProbLive", "ProbDead"])

    entry = metadata.loc[:, ["Drug", "Concentration"]].groupby(
        features.index.values).first()
    entry.loc[entry.Drug == "DMSO", "Concentration"] = "nan"
    entry.columns = ("Product.Name", "Concentration")
    entry["Product.Name"] = entry["Product.Name"].astype(str)
    entry["Concentration"] = entry["Concentration"].astype(str)
    entry["Plate.ID"] = metadata.Plate.groupby(features.index.values).first()
    entry["Well.ID"] = metadata.Well.groupby(features.index.values).first()
    entry["Num.Objects"] = results.Prediction.groupby(
        features.index.values).count()
    entry["Percent.Dead"] = results.Prediction.groupby(
        features.index.values).mean()
    entry["Percent.Live"] = 1 - entry["Percent.Dead"]
    entry["Mean.Certainty.Live"] = np.nan
    mcl = results.loc[results.Prediction == 0, "ProbLive"].groupby(
        features.index.values[results.Prediction == 0]).mean()
    entry.loc[mcl.index, "Mean.Certainty.Live"] = mcl
    entry["Mean.Certainty.Dead"] = np.nan
    mcd = results.loc[results.Prediction == 1, "ProbDead"].groupby(
        features.index.values[results.Prediction == 1]).mean()
    entry.loc[mcd.index, "Mean.Certainty.Dead"] = mcd

    if channels is None:
        out_fn = os.path.join(
            Config.DEADORGANOIDCLASSIFIERDIR, "results",
            "{}_organoid_classification.csv".format(cell_line))
    else:
        out_fn = os.path.join(
            Config.DEADORGANOIDCLASSIFIERDIR, "results",
            "{}_organoid_classification__{}.csv".format(
                cell_line, "_".join(sorted(channels))))
    if not os.path.isdir(os.path.dirname(out_fn)):
        os.makedirs(os.path.dirname(out_fn))

    entry.to_csv(out_fn, index=False)
    
def run_classifier_log(cell_line, channels=None):
    """
    Run the trained classifier on a cell line with a detailed output log

    The classifier can optionally be run on features that
    only relate to specific channels.
    :param cell_line:
    :param channels:
    :return:
    """
    if cell_line not in Utils.get_all_lines("human"):
        raise ValueError(
            "'cell_line' must be one of {}".format(Utils.get_all_lines()))

    if channels is not None:
        if not np.all(np.in1d(channels, ("Actin", "DAPI", "FITC"))):
            raise ValueError(
                "'channels' must be an iterable containing any combination "
                "of the values ('Actin', 'DAPI', 'FITC')")

    if np.sum(np.in1d(channels, ("Actin", "DAPI", "FITC"))) == 3:
        channels = None

    clf, clf_acc, feature_names, x_test, y_test = train_classifier(
        cell_line=cell_line, channels=channels)

    features, metadata, feature_names = load_features(
        cell_line=cell_line, channels=channels)

    prediction = clf.predict(features.values) == "POS"
    prediction_prob = clf.predict_proba(features.values)

    results = pd.DataFrame(
        data=np.concatenate((prediction[:, None], prediction_prob), axis=1),
        columns=["Prediction", "ProbLive", "ProbDead"])

    entry = results

    if channels is None:
        out_fn = os.path.join(
            Config.DEADORGANOIDCLASSIFIERDIR, "results",
            "{}_organoid_classification_log.csv".format(cell_line))
    else:
        out_fn = os.path.join(
            Config.DEADORGANOIDCLASSIFIERDIR, "results",
            "{}_organoid_classification__{}_log.csv".format(
                cell_line, "_".join(sorted(channels))))
    if not os.path.isdir(os.path.dirname(out_fn)):
        os.makedirs(os.path.dirname(out_fn))

    entry.to_csv(out_fn, index=False)


def run_for_all_channel_combinations(line):
    """
    Run classifier for all channel combinations
    :param line:
    :return:
    """
    for ch in CHANNELCOMBOS:
        print(line, ch)
        run_classifier(cell_line=line, channels=ch)

def run_on_all_lines_without_dx():
    """
    This is a convenience function to run the classifier over all plates for
    all lines while skipping diagnostics and returning a very detailed log of organoid classification results.
    :return:
    """
    lines = Utils.get_all_lines("human")
    for line in lines:
        print(line)
        run_classifier_log(cell_line=line, channels=None)
      

def run_on_all_lines():
    """
    This is a convenience function to run the classifier over all plates for
    all lines and for multiple combinations of channels
    :return:
    """
    lines = Utils.get_all_lines("human")
    for line in lines:
        for ch in CHANNELCOMBOS:
            print(line, ch)
            run_classifier(cell_line=line, channels=ch)

    get_best_acc_matrices()
    create_diagnostic_data()


def create_roc_data(clf_line, data_line=None, channels=None, plot=False):
    """
    Create ROC data and either plot it or save it to disk

    :param clf_line:
    :param data_line:
    :param channels:
    :param plot:
    :return:
    """

    if data_line is None:
        data_line = clf_line

    clf = train_classifier(cell_line=clf_line, channels=channels)
    classifier = clf[0]
    val_clf = train_classifier(data_line, channels=channels)
    x_val = val_clf[3]
    y_val = val_clf[4]

    y_true = sklearn.preprocessing.label_binarize(
        y_val, classes=classifier.classes_)[..., 0]
    y_pred = classifier.predict_proba(x_val)[..., 1]
    fpr, tpr, thresh = sklearn.metrics.roc_curve(y_true, y_pred)
    roc_auc = sklearn.metrics.auc(fpr, tpr)

    if plot:
        import matplotlib.pyplot as plt
        plt.figure()
        lw = 2
        plt.plot(fpr, tpr, color='darkorange',
                 lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
        plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        if channels is None:
            plt.title('ROC for {} on {}'.format(clf_line, data_line))
        else:
            plt.title('ROC for {} on {} for {}'.format(
                clf_line, data_line, "_".join(sorted(channels))))
        plt.legend(loc="lower right")
        plt.show()
    else:
        if channels is None:
            out_fn = os.path.join(
                Config.DEADORGANOIDCLASSIFIERDIR, "diagnostics",
                "roc_data_{}_on_{}.csv".format(
                    clf_line, data_line))
        else:
            out_fn = os.path.join(
                Config.DEADORGANOIDCLASSIFIERDIR, "diagnostics",
                "roc_data_{}_on_{}__{}.csv".format(
                    clf_line, data_line, "_".join(sorted(channels))))
        if not os.path.isdir(os.path.dirname(out_fn)):
            os.makedirs(os.path.dirname(out_fn))

        if os.path.isfile(out_fn):
            return None

        out_dat = np.stack((thresh, fpr, tpr), axis=1)
        with open(out_fn, "w") as f:
            if channels is None:
                f.write(
                    "# ROC Data for {} classifier applied to {} "
                    "data\n".format(clf_line, data_line))
            else:
                f.write(
                    "# ROC Data for {} classifier applied to {} "
                    "data for channels {}\n".format(
                        clf_line, data_line, "_".join(sorted(channels))))
            f.write("# AUC: %s\n" % roc_auc)
            f.write("Threshold,FalsePosRate,TruePosRate\n")
            for ii in range(out_dat.shape[0]):
                f.write(",".join(out_dat[ii, ...].astype(str)) + "\n")


def create_diagnostic_data():
    """
    This function generates diagnostic data for the classifiers. This includes:

      - ROC curves for each individual cell line on the test data
      - Transfer-learning accuracy grid
      - Feature importances

    This generates diagnostic data for all channel combinations

    :return:
    """

    lines = Utils.get_all_lines("human")

    # ROC curves
    for clf_line in lines:
        for data_line in lines:
            for ch_combo in CHANNELCOMBOS:
                create_roc_data(
                    clf_line=clf_line, data_line=data_line, channels=ch_combo)

    # Transfer Learning Accuracy Grid
    # Load validation data and classifiers
    for ch_combo in CHANNELCOMBOS:
        val_data = []
        clfs = []
        for line in lines:
            clf = train_classifier(cell_line=line, channels=ch_combo)
            val_data.append((clf[3], clf[4]))
            clfs.append((clf[0], clf[2]))

        # Calculate accuracy matrix
        acc_matrix = []
        for clf in clfs:
            entry = []
            for ii in range(len(val_data)):
                vd = val_data[ii]
                entry.append(clf[0].score(*vd))
            acc_matrix.append(entry)

        if ch_combo is None:
            out_fn = os.path.join(
                Config.DEADORGANOIDCLASSIFIERDIR, "diagnostics",
                "TransferLearningMatrix.csv")
        else:
            out_fn = os.path.join(
                Config.DEADORGANOIDCLASSIFIERDIR, "diagnostics",
                "TransferLearningMatrix__{}.csv".format(
                    "_".join(sorted(ch_combo))))

        if not os.path.isdir(os.path.dirname(out_fn)):
            os.makedirs(os.path.dirname(out_fn))

        with open(out_fn, "w") as f:
            if ch_combo is None:
                f.write(
                    "# Accuracy of classifier for cell line (ROW) applied to "
                    "validation data for cell line (COL)\n")
            else:
                f.write(
                    "# Accuracy of classifier for cell line (ROW) applied to "
                    "validation data for cell line (COL) if trained only on " 
                    "{}\n".format("_".join(sorted(ch_combo))))
            f.write("CLASSIFIER_DATA," + ",".join(lines) + "\n")
            for ii in range(len(acc_matrix)):
                f.write(lines[ii] + "," +
                        ",".join([str(s) for s in acc_matrix[ii]]) + "\n")

    # Feature accuracies
    for ch_combo in CHANNELCOMBOS:
        for line in lines:
            clf = train_classifier(cell_line=line, channels=ch_combo)
            feature_names = clf[2]
            importances = pd.DataFrame(
                data=clf[0].feature_importances_,
                index=feature_names,
                columns=["Importance"])
            if ch_combo is None:
                out_fn = os.path.join(
                    Config.DEADORGANOIDCLASSIFIERDIR, "diagnostics",
                    "FeatureImportances_{}.csv".format(line))
            else:
                out_fn = os.path.join(
                    Config.DEADORGANOIDCLASSIFIERDIR, "diagnostics",
                    "FeatureImportances_{}__{}.csv".format(
                        line, "_".join(sorted(ch_combo))))

            if not os.path.isdir(os.path.dirname(out_fn)):
                os.makedirs(os.path.dirname(out_fn))

            importances.to_csv(out_fn)

    # Adjusted accuracies
    get_best_acc_matrices()


def get_acc_for_threshs(clf_line, data_line=None, channels=None, threshs=None):
    """
    Applies a classifier trained on 'clf_line' to the validation data of
    'data_line' and returns the accuracy for various thresholds.

    Note that this is NOT a statistically robust method! In principle, the
    dataset needs to be separated into three categories, one for training
    via cross validation, one for choosing the optimal thresholds and other
    diagnostic parameter tuning, and one to assess the final metrics.
    :param clf_line:
    :param data_line:
    :param channels:
    :param threshs:
    :return:
    """

    if data_line is None:
        data_line = clf_line

    if threshs is None:
        threshs = np.arange(0, 1, 0.1)
    else:
        threshs = np.array(threshs)

    clf = train_classifier(cell_line=clf_line, channels=channels)
    classifier = clf[0]
    val_clf = train_classifier(data_line, channels=channels)
    x_val = val_clf[3]
    y_val = val_clf[4]

    y_true = sklearn.preprocessing.label_binarize(
        y_val, classes=classifier.classes_).ravel()
    y_pred = classifier.predict_proba(x_val)[
        ..., classifier.classes_ == "POS"].ravel()

    accs = []
    for thresh in threshs:
        accs.append(np.mean(
            y_true == (y_pred >= thresh)))

    return pd.DataFrame(
        data=accs,
        columns=["Accuracy"],
        index=threshs)


def get_best_acc_matrices(threshs=None):
    """
    Get the optimal accuracy matrix using variable thresholds
    :param threshs:
    :return:
    """

    lines = Utils.get_all_lines("human")
    for ch_combo in CHANNELCOMBOS:
        acc_matrix = []
        thresh_matrix = []
        for clf_line in lines:
            acc_row = []
            thresh_row = []
            for data_line in lines:
                accs = get_acc_for_threshs(
                    clf_line=clf_line, data_line=data_line,
                    channels=ch_combo, threshs=threshs)
                acc_row.append(accs.Accuracy.max())
                thresh_row.append(accs.index.values[
                    accs.Accuracy == accs.Accuracy.max()][0])
            acc_matrix.append(acc_row)
            thresh_matrix.append(thresh_row)

        if ch_combo is None:
            acc_out_fn = os.path.join(
                Config.DEADORGANOIDCLASSIFIERDIR, "diagnostics",
                "TransferLearningMatrix_variableThresholds.csv")
            thresh_out_fn = os.path.join(
                Config.DEADORGANOIDCLASSIFIERDIR, "diagnostics",
                "TransferLearningMatrix_variableThresholds_thresholds.csv")
        else:
            acc_out_fn = os.path.join(
                Config.DEADORGANOIDCLASSIFIERDIR, "diagnostics",
                "TransferLearningMatrix_variableThresholds__{}.csv".format(
                    "_".join(sorted(ch_combo))))
            thresh_out_fn = os.path.join(
                Config.DEADORGANOIDCLASSIFIERDIR, "diagnostics",
                "TransferLearningMatrix_variableThresholds_thresholds__{}.csv".format(
                    "_".join(sorted(ch_combo))))

        if not os.path.isdir(os.path.dirname(acc_out_fn)):
            os.makedirs(os.path.dirname(acc_out_fn))

        with open(acc_out_fn, "w") as f:
            if ch_combo is None:
                f.write(
                    "# Accuracy of classifier for cell line (ROW) applied to "
                    "validation data for cell line (COL) with variable "
                    "thresholds based on ROC optimization\n")
            else:
                f.write(
                    "# Accuracy of classifier for cell line (ROW) applied to "
                    "validation data for cell line (COL) if trained only on " 
                    "{} with variable thresholds based on ROC "
                    "optimization\n".format("_".join(sorted(ch_combo))))
            f.write("CLASSIFIER_DATA," + ",".join(lines) + "\n")
            for ii in range(len(acc_matrix)):
                f.write(lines[ii] + "," + ",".join(
                    [str(s) for s in acc_matrix[ii]]) + "\n")

        with open(thresh_out_fn, "w") as f:
            if ch_combo is None:
                f.write(
                    "# Thresholds for the optimal accuracy of classifier for "
                    "cell line (ROW) applied to validation data for cell line "
                    "(COL) with variable thresholds based on ROC "
                    "optimization\n")
            else:
                f.write(
                    "# Thresholds for the optimal accuracy of classifier for "
                    "cell line (ROW) applied to validation data for cell line "
                    "(COL) if trained only on {} with variable thresholds "
                    "based on ROC optimization\n".format(
                        "_".join(sorted(ch_combo))))
            f.write("CLASSIFIER_DATA," + ",".join(lines) + "\n")
            for ii in range(len(acc_matrix)):
                f.write(lines[ii] + "," + ",".join(
                    [str(s) for s in thresh_matrix[ii]]) + "\n")
