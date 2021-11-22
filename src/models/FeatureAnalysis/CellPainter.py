# This module contains functions to draw overlays onto well images, e.g. to
# visualize which organoids are identified as blurry or dead versus alive
import BlurryOrganoidClassifier
from OrganoidFeatures import OrganoidFeatures
import ImageLoader
import ImagePreprocessor
import numpy as np
# import matplotlib.pyplot as plt
import scipy.ndimage.measurements
import skimage.transform
import OrganoidViabilityClassifier


def paint_blurry_organoids(well):
    """
    Loads features for a well, runs the first three steps of the feature
    normalization, shrapnel removal, edge object removal, and blurry organoid
    identification, and visualizes which organoids are identified as blurry and
    which are identified as focused.
    :param well: String, e.g. 'D004T01P003L02_A_01'
    :return:
    """

    border = 5
    scale_factor = 1/4

    # Load well feature
    feature_object = OrganoidFeatures.load_wells(wells=[well], normalize="raw")
    feature_object.remove_size_outliers()
    # feature_object.remove_organoids_at_boundary()
    features = feature_object.get_features()
    metadata = feature_object.get_metadata()

    # Run blurry organoid classifier over data
    plate = well.split("_")[0]
    cls = BlurryOrganoidClassifier.get_classifier(
        plate=plate, features=features, metadata=metadata)
    cls_labels = [s[0].upper() for s in cls["classifier"].classes_]
    values = np.eye(3)
    val_dict = {cls_labels[ii]: values[ii] for ii in range(len(cls_labels))}
    focus_class = cls["classifier"].predict(
        features.loc[:, cls["feature_names"]])
    # # For Keras classifier
    # features_cls = features.apply(func=lambda xx: (xx - xx.min()) / (xx.max() - xx.min()))
    # features_cls = features_cls.fillna(0)
    # focus_class = cls["classifier"].predict(features_cls.loc[:, cls["feature_names"]])
    # focus_class = (focus_class < 0.5).astype(np.str)
    # val_dict = {"TRUE": (0, 1, 0), "FALSE": (1, 0, 0)}

    all_fields = []
    for field in metadata.Field.unique():
        field_features = features.loc[metadata.Field == field, :]
        field_predicted_class = focus_class[metadata.Field == field]

        # Load projection and segmentation
        proj = ImageLoader.load_projection(well=well, field_index=field - 1)
        proj = ImagePreprocessor.preprocess_image(proj)
        mask = ImageLoader.load_segmentation(
            well=well, field_index=field - 1, probs=False)

        labeled_mask = scipy.ndimage.measurements.label(
            -(mask.astype(np.int16)) < 0)[0]
        painted_mask = np.zeros(shape=mask.shape + (3,))
        for organoid_id in range(len(field_features)):
            y, x = field_features.iloc[organoid_id, :][
                ["x.0.m.cx", "x.0.m.cy"]].astype(np.int32)
            label = field_predicted_class[organoid_id][0].upper()
            value = val_dict[label]
            region_id = labeled_mask[x, y]
            if region_id == 0:
                continue
            painted_mask[np.where(labeled_mask == region_id)] = value

        # Fill up the other organoids with a third color
        tmp = np.max(painted_mask, axis=2)
        unpainted_points = np.where((tmp == 0) * (mask != 0))
        unpainted_points = (
            unpainted_points[0], unpainted_points[1],
            np.repeat(2, len(unpainted_points[0])))
        painted_mask[unpainted_points] = 1

        proj_resize = skimage.transform.rescale(
            image=proj, scale=scale_factor, mode="reflect")
        painted_mask_resize = skimage.transform.rescale(
            image=painted_mask, scale=scale_factor, mode="reflect")

        combined = np.concatenate((
            proj_resize, np.ones(shape=(proj_resize.shape[0], border, 3)),
            painted_mask_resize), axis=1)

        if field < metadata.Field.max():
            combined = np.concatenate(
                (combined, np.ones(shape=(border, combined.shape[1], 3))),
                axis=0)

        all_fields.append(combined)

    all_fields = np.concatenate(all_fields, axis=0)
    return all_fields


def paint_organoid_viability(well):
    """
    Loads features for a well, runs the organoid viability classifier, and
    visualizes the output
    :param well: String, e.g. 'D004T01P003L02_A_01'
    :return:
    """

    border = 5
    scale_factor = 1/4

    # Load well feature
    line = well[0:7]
    all_features, all_metadata, feature_names = OrganoidViabilityClassifier.load_features(
        cell_line=line, channels=None)
    features = all_features.loc[all_features.index == well, :]
    metadata = all_metadata.loc[all_metadata.index == well, :]

    # Load classifier
    cls, cls_acc, fnames_cls, x_test, y_test = OrganoidViabilityClassifier.train_classifier(
        cell_line=line, channels=None)

    # Run classifier
    prediction = cls.predict(features.values)
    cls_labels = [s.upper() for s in cls.classes_][::-1]
    values = np.eye(3)
    val_dict = {cls_labels[ii]: values[ii] for ii in range(len(cls_labels))}

    all_fields = []
    for field in metadata.Field.unique():
        field_features = features.loc[metadata.Field == field, :]
        field_metadata = metadata.loc[metadata.Field == field, :]
        field_predicted_class = prediction[metadata.Field == field]

        # Load projection and segmentation
        proj = ImageLoader.load_projection(well=well, field_index=field - 1)
        proj = ImagePreprocessor.preprocess_image(proj)
        mask = ImageLoader.load_segmentation(
            well=well, field_index=field - 1, probs=False)

        labeled_mask = scipy.ndimage.measurements.label(
            -(mask.astype(np.int16)) < 0)[0]
        painted_mask = np.zeros(shape=mask.shape + (3,))
        for organoid_id in range(len(field_features)):
            y, x = field_metadata.iloc[organoid_id, :][
                ["OriginalX", "OriginalY"]].astype(np.int32)
            label = field_predicted_class[organoid_id].upper()
            value = val_dict[label]
            region_id = labeled_mask[x, y]
            if region_id == 0:
                continue
            painted_mask[np.where(labeled_mask == region_id)] = value

        proj_resize = skimage.transform.rescale(
            image=proj, scale=scale_factor, mode="reflect")
        painted_mask_resize = skimage.transform.rescale(
            image=painted_mask, scale=scale_factor, mode="reflect")

        combined = np.concatenate((
            proj_resize, np.ones(shape=(proj_resize.shape[0], border, 3)),
            painted_mask_resize), axis=1)

        if field < metadata.Field.max():
            combined = np.concatenate(
                (combined, np.ones(shape=(border, combined.shape[1], 3))),
                axis=0)

        all_fields.append(combined)

    all_fields = np.concatenate(all_fields, axis=0)
    return all_fields
