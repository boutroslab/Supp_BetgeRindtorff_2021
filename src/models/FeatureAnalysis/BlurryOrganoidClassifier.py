"""
Trains and applies a classifier to detect blurry organoids. The
classification is performed plate-wise to avoid batch effects as these
cannot be normalized until after blurry organoids have been filtered
away
"""
import Config
import os
import pickle
import numpy as np
import sklearn.model_selection
import sklearn.ensemble


def learn_classifier(features):
    """
    This is a generalized, static function to train a classifier to
    differentiate between blurry and in-focus organoids.

    The function consists of two steps:
    1. Initial annotation with a high-precision/low-recall segmentation
       based on the standard deviation of intensities. The idea is that
       blurry organoids will have a lower standard deviation due to the
       superimposed guassian blur. I use organoids in the top and bottom
       20%-ile as training samples for in-focus and blurry organoids,
       respectively.
    2. Train classifier on all unrelated features to prevent
       overfitting. "Unrelated" in this context means that all correlating
       features are removed

    :param features: A pandas object containing the features
    :return:
    """

    # selection_features = np.mean(
    #     features[["x.a.b.sd", "x.c.b.sd"]].values,
    #     axis=1)
    selection_features = np.max(
        features[["x.a.b.mean", "x.c.b.mean"]].values,
        axis=1)

    pos_samples = np.where(
        selection_features >= np.percentile(selection_features, 50))[0]
    neg_samples = np.where(
        selection_features <= np.percentile(selection_features, 30))[0]

    # Downsample to match sizes
    sample_size = min(len(pos_samples), len(neg_samples))
    pos_samples = np.random.choice(a=pos_samples, size=sample_size)
    neg_samples = np.random.choice(a=neg_samples, size=sample_size)

    x = features.iloc[np.concatenate((pos_samples, neg_samples))]
    y = np.repeat(("good", "blurry"), (len(pos_samples), len(neg_samples)))

    # Remove constant features
    x = x.iloc[:, np.std(x.values, axis=0) != 0]

    # Remove all features that correlated with the standard deviations on
    # the training set
    def maxcorr(mat_a, mat_b):
        """
        Calculates the maximum correlation of each column of matrix mat_a
        with all columns of mat_b
        Taken from
        https://stackoverflow.com/questions/33650188/
            efficient-pairwise-correlation-for-two-matrices-of-features
        :param mat_a:
        :param mat_b:
        :return:
        """
        n = mat_a.shape[0]
        sa = mat_a.sum(0)
        sb = mat_b.sum(0)
        p1 = n * np.dot(mat_b.T, mat_a)
        p2 = sa * sb[:, None]
        p3 = n * ((mat_b ** 2).sum(0)) - (sb ** 2)
        p4 = n * ((mat_a ** 2).sum(0)) - (sa ** 2)
        pcorr = ((p1 - p2) / np.sqrt(p4 * p3[:, None]))
        return np.max(np.abs(pcorr), axis=0)

    x_sel = x[["x.a.b.mean", "x.c.b.mean"]]
    x = x.iloc[:, maxcorr(x.values, x_sel.values) < 0.9]

    # Train classifier
    x_train, x_val, y_train, y_val = sklearn.model_selection.train_test_split(
        x, y, test_size=0.25)

    rfclf = sklearn.ensemble.RandomForestClassifier(n_estimators=100)
    # rfclf = sklearn.linear_model.LogisticRegression()
    rfclf.fit(x_train, y_train)
    val_acc = rfclf.score(x_val, y_val)

    return {
        "classifier": rfclf, "feature_names": x.columns.values,
        "accuracy": val_acc}


def get_classifier(plate, features=None, metadata=None):
    """
    Retrieves the classifier for the plate if it exists, otherwise trains it
    and then returns it. If 'train' is set to False then no training is ever
    performed.

    'features_object' can optionally be passed to avoid loading features from
    disk twice.
    :param plate:
    :param features:
    :param metadata:
    :return:
    """
    cls_fn = os.path.join(
        Config.BLURRYORGANOIDCLSDIR,
        "Classifier_%s_%s.pkl" % (plate, Config.FEATURETYPE))
    if os.path.isfile(cls_fn):
        with open(cls_fn, "rb") as f:
            clf = pickle.load(f)
    else:
        print(
            "Training blurry organoid classifier for '{}_{}'".format(
                plate, Config.FEATURETYPE))

        if features is None or metadata is None:
            raise ValueError(
                "'features' and 'metadata' cannot be 'None' if a classifier "
                "should be trained.")
        fo_plates = metadata["Plate"].unique()
        if len(fo_plates) != 1 or fo_plates[0] != plate:
            raise ValueError(
                "'features_object' may only contain data from 'plate'")
        clf = learn_classifier(features=features)
        clf["feature_type"] = Config.FEATURETYPE
        with open(cls_fn, "wb") as f:
            pickle.dump(clf, f)
    return clf
