import os
import h5py
import scipy.stats
import numpy as np
import pandas as pd
import statsmodels.robust
import Config
import Utils
import BlurryOrganoidClassifier


class OrganoidFeatures(object):
    """
    A class to lead features and metadata associated with the features. This
    class also contains functions to modify the values, e.g. filtering,
    transforming, or normalizing them.
    """

    def __init__(self, **kwargs):
        """
        Initializer function
        :param args:
        :param kwargs:
        """

        self._features = None
        self._metadata = None
        self._verbose = False
        # This flag makes sure that feature preprocessing steps can only be
        # called on full cell line data
        self._loaded_line = False
        self.__dict__.update(kwargs)

    @classmethod
    def load_wells(cls, wells, normalize, **kwargs):
        """
        Loads features for specific wells
        :param wells:
        :param normalize:
        :param kwargs: Any additional members to pass to the OrganoidFeatures
               object.
        :return:
        """
        obj = cls(wells=wells, normalize=normalize, **kwargs)
        if normalize == "raw":
            obj._features, obj._metadata = obj.load_raw_features_from_wells(
                wells=wells)
        else:
            plates = sorted(set([well[0:14] for well in wells]))
            drug = kwargs["drug"] if "drug" in kwargs.keys() else None
            obj._features, obj._metadata = obj.load_features_from_file(
                plates=plates, normalize=normalize, drug=drug)
            sel_indices = np.in1d(obj._features.index, wells)
            obj._features = obj._features.loc[sel_indices, :]
            obj._metadata = obj._metadata.loc[sel_indices, :]

        return obj

    @classmethod
    def load_line(cls, line, normalize, **kwargs):
        """
        Loads features for an entire line
        :param line:
        :param normalize:
        :param kwargs: Any additional members to pass to the OrganoidFeatures
               object.
        :return:
        """
        plates = Utils.get_plates_for_line(line=line)
        drug = kwargs["drug"] if "drug" in kwargs.keys() else None
        obj = cls(plates=plates, normalize=normalize, **kwargs)
        obj._loaded_line = True
        obj._features, obj._metadata = obj.load_features_from_file(
            plates=plates, normalize=normalize, drug=drug)
        return obj

    @classmethod
    def load_plates(cls, plates, normalize, **kwargs):
        """
        Loads features for arbitrary plates
        :param plates:
        :param normalize:
        :param kwargs: Any additional members to pass to the OrganoidFeatures
               object.
        :return:
        """
        drug = kwargs["drug"] if "drug" in kwargs.keys() else None
        obj = cls(plates=plates, normalize=normalize, **kwargs)
        obj._features, obj._metadata = obj.load_features_from_file(
            plates=plates, normalize=normalize, drug=drug)
        return obj

    @classmethod
    def load_inception_features(cls, line, normalize, **kwargs):
        """
        Load the well-based inception features.

        WARNING: These are well-based features and the metadata will have
        a different structure. Some functions of this class are inherently
        incompatible with the inception features. Use with care.
        :param line:
        :param normalize:
        :param kwargs:
        :return:
        """
        plates = Utils.get_plates_for_line(line=line)
        fn_appendix = "processedFeatures_inception.csv" if \
            normalize else "features_inception.csv"
        features = []
        for plate in plates:
            feature_fn = os.path.join(
                Config.INCEPTIONDIR, plate,
                "{}_{}".format(plate, fn_appendix))
            if not os.path.isfile(feature_fn):
                continue
            features.append(pd.DataFrame.from_csv(path=feature_fn))
        features = pd.DataFrame(pd.concat(features))

        # Load libraries and create metadata
        well_names = features.index.values
        layout_ids = np.array([wn[11:14] for wn in well_names])
        plates = np.array([wn[0:14] for wn in well_names])
        rep_ids = np.array([
            Utils.get_replicate_id(plate) for plate in plates])
        well_ids = np.array([wn[15:16] + wn[17:19] for wn in well_names])
        well_layout_ids = np.array([
            "{}_{}".format(s[0], s[1]) for s in zip(layout_ids, well_ids)])

        layouts = []
        for layout_id in sorted(set(layout_ids)):
            layout = pd.read_excel(io=os.path.join(
                Config.LAYOUTDIR, "%s.xlsx" % layout_id))
            if "concentration" not in layout.columns:
                layout.loc[:, "concentration"] = None
            layout.loc[:, "Library_ID"] = layout_id
            layout.index = layout["Library_ID"].str.cat(
                layout["Well_ID_384"], "_").values
            layouts.append(layout.loc[
                           :, ["Product.Name", "Well_ID_384",
                               "concentration", "Library_ID"]])
        layouts = pd.DataFrame(pd.concat(layouts))

        drugs = layouts.loc[
            well_layout_ids, "Product.Name"].values.astype(np.str)
        concentrations = layouts.loc[
            well_layout_ids, "concentration"].values.astype(np.float)

        metadata = pd.DataFrame(
            data={
                "Drug": drugs,
                "Concentration": concentrations,
                "Line": line,
                "Plate": plates,
                "Replicate": rep_ids,
                "Well": well_ids,
                "Layout": layout_ids},
            index=features.index.values)

        obj = cls(plates=plates, normalize=normalize, **kwargs)
        obj._loaded_line = True
        obj._features = features
        obj._metadata = metadata
        return obj

    @classmethod
    def load_metadata(cls, line, normalize, **kwargs):
        """
        Load only the metadata for a line
        :param line:
        :param normalize
        :return:
        """
        plates = Utils.get_plates_for_line(line=line)
        obj = cls(**kwargs)
        obj._metadata = obj.load_metadata_from_file(
            plates=plates, normalize=normalize)
        return obj

    @classmethod
    def load_well_aggregated_features(cls, line):
        """
        Loads well-aggregated features
        :param line:
        :return:
        """
        features = []
        metadata = []
        plates = Utils.get_plates_for_line(line)
        for plate in plates:
            in_fn = os.path.join(
                Config.FEATUREDIR, plate,
                "{}_averaged_features.h5".format(plate))
            f, m = OrganoidFeatures.load_features(in_fn)
            features.append(f)
            metadata.append(m)
        features = pd.DataFrame(pd.concat(features))
        metadata = pd.DataFrame(pd.concat(metadata))
        obj = cls(plates=plates)
        obj._features = features
        obj._metadata = metadata

    @staticmethod
    def combine_wells(plates, verbose=False):
        """
        Combine the individual feature files for wells into single files for
        each plate
        :param plates: An iterable of plate ids. Even single plates must be
               passed as iterables
        :param verbose:
        :return:
        """
        for plate in plates:
            if verbose:
                print("Plate {} | #{} of {}".format(
                    plate, plates.index(plate)+1, len(plates)))
            fdir = os.path.join(Config.FEATUREDIR, plate, "wells")
            outfn = os.path.join(
                Config.FEATUREDIR, plate, "{}_features.h5".format(plate))
            if os.path.isfile(outfn):
                print("Feature File for {} exists".format(plate))
                continue
            well_files = [
                "_".join(wf.split("_")[0:3]) for
                wf in sorted(os.listdir(fdir))]
            feats, md = OrganoidFeatures.load_raw_features_from_wells(
                wells=well_files)
            OrganoidFeatures.save_features(
                features=feats, metadata=md, filename=outfn)

    def get_features(self):
        return self._features

    def get_metadata(self):
        return self._metadata

    @staticmethod
    def load_raw_features_from_wells(wells):
        """
        Loads the raw organoid features from the well files. Ignores 'plates'
        and 'line' parameters passed to constructor.

        :param wells
        :return:
        """

        if wells is not None:
            if not isinstance(wells, (list, tuple)):
                raise ValueError(
                    "'wells' must be either None or a list/tuple")

        # Check feature type
        keymap = {
            "organoids": {
                "features": "features_organoids",
                "feature_names": "feature_names_organoids"},
            "clumps": {
                "features": "features_clumps",
                "feature_names": "feature_names_clumps"}}
        hdf5_keys = keymap[Config.FEATURETYPE]

        # Set well folder
        well_folder = "wells"
        fname_appendix = "_features.h5"

        # Set wells to load
        if wells is None:
            raise ValueError("'wells' must be set")

        # Trim any trailing text from the wells and set the file names
        wells = ["_".join(well.split("_")[0:3]) for well in wells]
        wells = sorted([s + fname_appendix for s in wells])

        # Load features for all wells
        features = []
        feature_names = []
        well_names = []
        object_ids = []
        for well in wells:
            feature_fn = os.path.join(
                Config.FEATUREDIR, well.split("_")[0], well_folder, well)
            if not os.path.isfile(feature_fn):
                continue

            Utils.rename_h5_datasets(
                fname=feature_fn, name_mapping={
                    'feature_names': 'feature_names_organoids',
                    'features': 'features_organoids'})

            try:
                with h5py.File(feature_fn, "r") as h5handle:
                    well_features = h5handle[hdf5_keys["features"]][()]
                    well_feature_names = h5handle[
                        hdf5_keys["feature_names"]][()]
                well_names.append(np.repeat(
                    "_".join(well.split("_")[0:3]),
                    well_features.shape[1]))
                features.append(well_features)
                feature_names.append(well_feature_names)
                object_ids.append(np.arange(well_features.shape[1]))
            except KeyError:
                continue

        if len(features) == 0:
            return None

        # Check feature names consistency
        if len(feature_names) > 1:
            fname_iter = iter(feature_names)
            if not all(np.array_equal(next(fname_iter), rest)
                       for rest in fname_iter):
                raise ValueError(
                    "Not all feature names were identical between files")
        feature_names = feature_names[0]

        # A quirk in Python3 means that the strings are read in as binary
        # blobs. This converts the feature names to proper strings
        feature_names = feature_names.astype(np.str)

        # Combine features and well names
        features = np.concatenate(features, axis=1)
        well_names = np.concatenate(well_names)
        object_ids = np.concatenate(object_ids)

        # Remove FIELD field and put it into a separate object
        # NOTE: I don't think I've ever actually used this
        if "FIELD" in feature_names:
            fields = features[np.where(
                feature_names == "FIELD")[0][0], ...]
            features = np.delete(
                features, np.where(feature_names == "FIELD")[0], axis=0)
            feature_names = np.delete(
                feature_names, np.where(feature_names == "FIELD")[0])
            fields = fields.astype(np.int)
        else:
            fields = np.repeat(a=(np.nan,), repeats=features.shape[1])

        features = pd.DataFrame(
            data=features.transpose(),
            index=well_names,
            columns=feature_names)

        # Load libraries and create metadata
        layout_ids = np.array([wn[11:14] for wn in well_names])
        cell_lines = np.array([wn[0:7] for wn in well_names])
        plates = np.array([wn[0:14] for wn in well_names])
        rep_ids = np.array([
            Utils.get_replicate_id(plate) for plate in plates])
        well_ids = np.array([wn[15:16] + wn[17:19] for wn in well_names])
        well_layout_ids = np.array([
            "{}_{}".format(s[0], s[1]) for s in zip(layout_ids, well_ids)])

        layouts = []
        for layout_id in sorted(set(layout_ids)):
            layout = pd.read_excel(io=os.path.join(
                Config.LAYOUTDIR, "%s.xlsx" % layout_id))
            if "concentration" not in layout.columns:
                layout.loc[:, "concentration"] = None
            layout.loc[:, "Library_ID"] = layout_id
            layout.index = layout["Library_ID"].str.cat(
                layout["Well_ID_384"], "_").values
            layouts.append(layout.loc[
                :, ["Product.Name", "Well_ID_384",
                    "concentration", "Library_ID"]])
        layouts = pd.DataFrame(pd.concat(layouts))

        drugs = layouts.loc[
            well_layout_ids, "Product.Name"].values.astype(np.str)
        concentrations = layouts.loc[
            well_layout_ids, "concentration"].values.astype(np.float)

        metadata = pd.DataFrame(
            data=drugs, columns=["Drug"], index=features.index.values)
        metadata.loc[:, "Concentration"] = concentrations
        metadata.loc[:, "Line"] = cell_lines
        metadata.loc[:, "Plate"] = plates
        metadata.loc[:, "Well"] = well_ids
        metadata.loc[:, "Replicate"] = rep_ids
        metadata.loc[:, "Field"] = fields
        metadata.loc[:, "ObjectID"] = object_ids
        metadata.loc[:, "OriginalX"] = features.loc[:, "x.0.m.cx"]
        metadata.loc[:, "OriginalY"] = features.loc[:, "x.0.m.cy"]
        metadata.loc[:, "Size"] = features.loc[:, "x.0.s.area"]

        return features, metadata

    @staticmethod
    def load_features_from_file(plates, normalize, drug=None):
        """
        Load features and metadata from the collated well file
        :param plates:
        :param normalize: ("raw", "stage_one", "normalized")
        :param drug:
        :return:
        """
        if plates is not None:
            if not isinstance(plates, (list, tuple)):
                raise ValueError(
                    "'plates' must be either None or a list/tuple")

        if isinstance(normalize, bool):
            raise ValueError(
                "Deprecation Error: 'normalize' must now be a string value")

        if normalize.lower() not in ("raw", "stage_one", "normalized"):
            raise ValueError(
                "'normalize' must be one of "
                "('raw', 'stage_one','normalized')")

        # Check feature type
        keymap = {
            "organoids": {
                "features": "features_organoids",
                "feature_names": "feature_names_organoids",
                "metadata": "metadata_organoids",
                "metadata_names": "metadata_names_organoids",
                "well_names": "well_names_organoids"},
            "clumps": {
                "features": "features_clumps",
                "feature_names": "feature_names_clumps",
                "metadata": "metadata_clumps",
                "metadata_names": "metadata_names_clumps",
                "well_names": "well_names_clumps"}}
        hdf5_keys = keymap[Config.FEATURETYPE]

        fn_appendix = {
            "raw": "features.h5",
            "stage_one": "processedFeatures_stageOne.h5",
            "normalized": "processedFeatures.h5"}[normalize]

        features = []
        metadata = []
        for plate in plates:
            feature_fn = os.path.join(
                Config.FEATUREDIR, plate,
                "{}_{}".format(plate, fn_appendix))
            if not os.path.isfile(feature_fn):
                continue
            with h5py.File(feature_fn, "r") as h5handle:
                # Column names
                plate_feature_names = h5handle[
                    hdf5_keys["feature_names"]][()].astype(np.str)
                plate_md_names = h5handle[
                    hdf5_keys["metadata_names"]][()].astype(np.str)

                # Metadata
                plate_md = h5handle[hdf5_keys["metadata"]][()]

                # If drug is set then extract only those rows
                if drug is not None:
                    sel_indices = plate_md[:, np.where(
                        plate_md_names == "Drug")[0][0]].astype(
                        np.str) == drug
                    if np.sum(sel_indices) == 0:
                        continue
                    plate_md = plate_md[sel_indices, :]
                    well_names = h5handle[hdf5_keys["well_names"]][sel_indices]
                    plate_features = h5handle[
                        hdf5_keys["features"]][:, sel_indices].transpose()
                # If not then select everything
                else:
                    well_names = h5handle[hdf5_keys["well_names"]][()]
                    plate_features = h5handle[
                        hdf5_keys["features"]][()].transpose()

                # Add to global list
                features.append(pd.DataFrame(
                    data=plate_features,
                    index=well_names.astype(np.str),
                    columns=plate_feature_names))
                metadata.append(pd.DataFrame(
                    data=plate_md,
                    index=well_names.astype(np.str),
                    columns=plate_md_names))

        features = pd.DataFrame(pd.concat(features))
        metadata = pd.DataFrame(pd.concat(metadata))

        # Adjust metadata column formats
        metadata = pd.DataFrame(
            data={
                "Drug": metadata["Drug"].values.astype(np.str),
                "Concentration": metadata["Concentration"].values.astype(np.str),
                "Line": metadata["Line"].values.astype(np.str),
                "Plate": metadata["Plate"].values.astype(np.str),
                "Well": metadata["Well"].values.astype(np.str),
                "Replicate": metadata["Replicate"].values.astype(np.float).astype(np.int),
                "Field": metadata["Field"].values.astype(np.float).astype(np.int),
                "ObjectID": metadata["ObjectID"].values.astype(np.float).astype(np.int),
                "OriginalX": metadata["OriginalX"].values.astype(np.float),
                "OriginalY": metadata["OriginalY"].values.astype(np.float),
                "Size": metadata["Size"].values.astype(np.float).astype(np.int)},
            index=features.index)

        return features, metadata

    @staticmethod
    def load_metadata_from_file(plates, normalize, drug=None):
        """
        Load only metadata from the collated well file
        :param plates:
        :param normalize: ("raw", "stage_one", "normalized")
        :param drug:
        :return:
        """
        if plates is not None:
            if not isinstance(plates, (list, tuple)):
                raise ValueError(
                    "'plates' must be either None or a list/tuple")

        if isinstance(normalize, bool):
            raise ValueError(
                "Deprecation Error: 'normalize' must now be a string value")

        if normalize.lower() not in ("raw", "stage_one", "normalized"):
            raise ValueError(
                "'normalize' must be one of "
                "('raw', 'stage_one','normalized')")

        # Check feature type
        keymap = {
            "organoids": {
                "features": "features_organoids",
                "feature_names": "feature_names_organoids",
                "metadata": "metadata_organoids",
                "metadata_names": "metadata_names_organoids",
                "well_names": "well_names_organoids"},
            "clumps": {
                "features": "features_clumps",
                "feature_names": "feature_names_clumps",
                "metadata": "metadata_clumps",
                "metadata_names": "metadata_names_clumps",
                "well_names": "well_names_clumps"}}
        hdf5_keys = keymap[Config.FEATURETYPE]

        fn_appendix = {
            "raw": "features.h5",
            "stage_one": "processedFeatures_stageOne.h5",
            "normalized": "processedFeatures.h5"}[normalize]

        metadata = []
        for plate in plates:
            feature_fn = os.path.join(
                Config.FEATUREDIR, plate,
                "{}_{}".format(plate, fn_appendix))
            if not os.path.isfile(feature_fn):
                continue
            with h5py.File(feature_fn, "r") as h5handle:
                # Column names
                plate_md_names = h5handle[
                    hdf5_keys["metadata_names"]][()].astype(np.str)

                # Metadata
                plate_md = h5handle[hdf5_keys["metadata"]][()]

                # If drug is set then extract only those rows
                if drug is not None:
                    sel_indices = plate_md[:, np.where(
                        plate_md_names == "Drug")[0][0]].astype(
                        np.str) == drug
                    if np.sum(sel_indices) == 0:
                        continue
                    plate_md = plate_md[sel_indices, :]
                    well_names = h5handle[hdf5_keys["well_names"]][sel_indices]
                # If not then select everything
                else:
                    well_names = h5handle[hdf5_keys["well_names"]][()]

                # Add to global list
                metadata.append(pd.DataFrame(
                    data=plate_md,
                    index=well_names.astype(np.str),
                    columns=plate_md_names))

        metadata = pd.DataFrame(pd.concat(metadata))

        # Adjust metadata column formats
        metadata = pd.DataFrame(
            data={
                "Drug": metadata["Drug"].values.astype(np.str),
                "Concentration": metadata["Concentration"].values.astype(np.str),
                "Line": metadata["Line"].values.astype(np.str),
                "Plate": metadata["Plate"].values.astype(np.str),
                "Well": metadata["Well"].values.astype(np.str),
                "Replicate": metadata["Replicate"].values.astype(np.float).astype(np.int),
                "Field": metadata["Field"].values.astype(np.float).astype(np.int),
                "ObjectID": metadata["ObjectID"].values.astype(np.float).astype(np.int),
                "OriginalX": metadata["OriginalX"].values.astype(np.float),
                "OriginalY": metadata["OriginalY"].values.astype(np.float),
                "Size": metadata["Size"].values.astype(np.float).astype(np.int)},
            index=metadata.index)

        return metadata

    @staticmethod
    def load_single_feature_from_file(feature, plates, normalize):
        """
        Loads a single feature and no metadata

        Because this is fast enough, there is no well-specific implementation.
        Filtering out specific wells or drugs should be done outside the scope
        of this method.
        :param feature:
        :param plates:
        :param normalize: ("raw", "stage_one", "normalized")
        :return:
        """
        if plates is not None:
            if not isinstance(plates, (list, tuple)):
                raise ValueError(
                    "'plates' must be either None or a list/tuple")

        if isinstance(normalize, bool):
            raise ValueError(
                "Deprecation Error: 'normalize' must now be a string value")

        if normalize.lower() not in ("raw", "stage_one", "normalized"):
            raise ValueError(
                "'normalize' must be one of "
                "('raw', 'stage_one','normalized')")

        # Check feature type
        keymap = {
            "organoids": {
                "features": "features_organoids",
                "feature_names": "feature_names_organoids",
                "well_names": "well_names_organoids"},
            "clumps": {
                "features": "features_clumps",
                "feature_names": "feature_names_clumps",
                "well_names": "well_names_clumps"}}

        hdf5_keys = keymap[Config.FEATURETYPE]

        # Set well folder
        fn_appendix = {
            "raw": "features.h5",
            "stage_one": "processedFeatures_stageOne.h5",
            "normalized": "processedFeatures.h5"}[normalize]

        features = []
        for plate in plates:
            feature_fn = os.path.join(
                Config.FEATUREDIR, plate,
                "{}_{}".format(plate, fn_appendix))
            if not os.path.isfile(feature_fn):
                continue
            with h5py.File(feature_fn, "r") as h5handle:
                # Column names
                plate_feature_names = h5handle[
                    hdf5_keys["feature_names"]][()].astype(np.str)
                feature_id = np.where(plate_feature_names == feature)[0]
                if len(feature_id) > 1:
                    raise ValueError("Duplicate feature names found!")
                if len(feature_id) == 0:
                    raise ValueError("Invalid Feature Name!")

                well_names = h5handle[hdf5_keys["well_names"]][()]
                plate_features = h5handle[
                                     hdf5_keys["features"]][feature_id, :].transpose()

                # Add to global list
                features.append(pd.DataFrame(
                    data=plate_features,
                    index=well_names.astype(np.str),
                    columns=[feature]))

        return pd.DataFrame(pd.concat(features))

    # ======================= ORGANOID PREPROCESSING ======================== #
    # These methods remove bad organoids, i.e. shrapnel, incorrectly
    # segmented objects, and blurry organoids
    # ======================================================================= #

    def remove_size_outliers(self):
        """
        Remove small objects (shrapnel) as well as excessively large objects
        (top 5% of size outliers).
        :return:
        """
        if self._verbose:
            print("Removing shrapnel ...")
        sizes = self._metadata.loc[:, "Size"].values
        sel_indices = np.where(np.logical_and(
            sizes >= Config.SIZETHRESHOLD,
            sizes < np.percentile(sizes, 95)))[0]
        self._features = self._features.iloc[sel_indices, :]
        self._metadata = self._metadata.iloc[sel_indices, :]

    def remove_eccentricity_outliers(self):
        """
        Organoids should have reasonably circular eccentricies, hence outliers
        in this regard are removed.

        Note that EBImage eccentricities are defined as:
        ecc = sqrt(1 - (minoraxis**2 / majoraxis**2)),
        i.e. 0 for circles and 1 for a line.
        :return:
        """
        if self._verbose:
            print("Removing debris based on eccentricity ...")
        eccentricities = self._features.loc[:, "x.0.m.eccentricity"].values
        self._features = self._features.loc[eccentricities < 0.95, :]
        self._metadata = self._metadata.loc[eccentricities < 0.95, :]

    def remove_organoids_at_boundary(self):
        """
        Partial organoids may skew results, so keep only organoids that are
        fully in the image
        :return:
        """
        # Keep only organoids whose maximum radius is smaller than their
        # distance to the image boundaries
        if self._verbose:
            print("Removing organoids at image edges ...")
        min_dist = np.min(
            np.stack(
                (self._features.loc[:, "x.0.m.cx"],
                 Config.IMAGESIZE[0] - self._features.loc[:, "x.0.m.cx"],
                 self._features.loc[:, "x.0.m.cy"],
                 Config.IMAGESIZE[0] - self._features.loc[:, "x.0.m.cy"]),
                axis=1), axis=1)
        max_radii = self._features.loc[:, "x.0.s.radius.max"]
        self._features = self._features.loc[max_radii < min_dist, :]
        self._metadata = self._metadata.loc[max_radii < min_dist, :]

    def remove_blurry_organoids(self):
        """
        Removes blurry objects. Autonomously trains and loads classifiers from
        disk if necessary.
        :return:
        """
        if self._verbose:
            print("Removing blurry organoids ...")
        plates = sorted(self._metadata["Plate"].unique())
        focus_class = np.repeat([None], len(self._features))
        for plate in plates:
            if self._verbose:
                print("  - Plate {}: {} of {}".format(
                    plate, plates.index(plate) + 1, len(plates)))
            sel_indices = self._metadata["Plate"] == plate
            cls = BlurryOrganoidClassifier.get_classifier(
                plate=plate,
                features=self.get_features().loc[sel_indices, :],
                metadata=self.get_metadata().loc[sel_indices, :])
            x = self._features.loc[sel_indices][cls["feature_names"]]
            focus_class[sel_indices] = cls["classifier"].predict(x)
        sel_indices = np.where(focus_class == "good")[0]
        self._features = self._features.iloc[sel_indices, :]
        self._metadata = self._metadata.iloc[sel_indices, :]

    def remove_bad_organoids(self, save=True):
        """
        The first preprocessing step: removing organoids that are deemed "bad"
        :param save:
        :return:
        """
        self.remove_size_outliers()
        self.remove_organoids_at_boundary()
        self.remove_blurry_organoids()

        # Save features
        if save:
            if self._verbose:
                print("Saving normalized features ...")
            for plate in self._metadata.loc[:, "Plate"].unique():
                if self._verbose:
                    print(" - Plate {}".format(plate))
                plate_indices = self._metadata.loc[:, "Plate"] == plate
                feature_fn = os.path.join(
                    Config.FEATUREDIR, plate,
                    "{}_processedFeatures_stageOne.h5".format(plate))
                if os.path.isfile(feature_fn):
                    print(
                        "Processed feature file (Stage 1) for {} exists and "
                        "cannot be overwritten by this method".format(plate))
                    continue
                OrganoidFeatures.save_features(
                    features=self._features.loc[plate_indices, :],
                    metadata=self._metadata.loc[plate_indices, :],
                    filename=feature_fn)

    # ======================== FEATURE PREPROCESSING ======================== #
    # These methods fix batch and plate effects and transform features
    # ======================================================================= #

    def remove_plate_effects(self):
        """
        Remove plate effects by using an additive ratio based approach.

        The DMSO controls of each plate within a cell line should come from
        the same distribution. This function rescales the features of a plate
        accordingly.

        This correction occurs WITHIN lines and batches

        :return:
        """

        if not self._loaded_line:
            print("---\n"
                  "Plate effect correction\n"
                  "must be performed on an entire line. Please use:\n"
                  "OrganoidFeatures.load_line() to load the data or set the\n"
                  "variable '_loaded_line' to True manually to bypass this\n"
                  "warning message\n"
                  "---\n")
            return None

        if self._verbose:
            print("Removing batch effects ...")
        for line in self._metadata.loc[:, "Line"].unique():
            line_indices = self._metadata.loc[:, "Line"] == line
            dmso_indices = np.logical_and(
                line_indices,
                self._metadata.loc[:, "Drug"] == "DMSO")

            # Normalize middle
            # plate_middle = self._features.loc[dmso_indices.values, :].groupby(
            #     self._metadata.loc[dmso_indices.values, "Plate"]).apply(
            #     lambda df: df.apply(
            #         lambda col: np.nanmean(scipy.stats.trimboth(col, 0.1))))
            plate_middle = self._features.loc[dmso_indices.values, :].groupby(
                self._metadata.loc[dmso_indices.values, "Plate"]).median()
            for plate in self._metadata.loc[line_indices, "Plate"].unique():
                plate_indices = np.logical_and(
                    line_indices,
                    self._metadata.loc[:, "Plate"] == plate)
                self._features.loc[plate_indices, :] += \
                    (plate_middle.median() - plate_middle.loc[plate, :])

            # Normalize deviation
            # plate_dev = self._features.loc[dmso_indices.values, :].groupby(
            #     self._metadata.loc[dmso_indices.values, "Plate"]).apply(
            #     lambda df: df.apply(
            #         lambda col: np.nanstd(scipy.stats.trimboth(col, 0.1))))
            plate_middle = self._features.loc[dmso_indices.values, :].groupby(
                self._metadata.loc[dmso_indices.values, "Plate"]).median()
            plate_dev = self._features.loc[dmso_indices.values, :].groupby(
                self._metadata.loc[dmso_indices.values, "Plate"]).apply(
                lambda df: np.abs(df - df.median()).median())
            for plate in self._metadata.loc[line_indices, "Plate"].unique():
                plate_indices = np.logical_and(
                    line_indices,
                    self._metadata.loc[:, "Plate"] == plate)
                f = self._features.loc[plate_indices, :]
                self._features.loc[plate_indices, :] = \
                    (f - plate_middle.median()) * \
                    (plate_dev.median() / plate_dev.loc[plate, :]) + \
                    plate_middle.median()

    def remove_batch_effects(self):
        """
        There was a notable batch effect between the first set of lines and
        the second set. The current version of batch correction involves
        calculating the per feature shift and stretch so that the full
        distributions across all wells of D020T01 and D020T02 are identical
        and then applying those transformations to each line of the second
        batch.

        This function MUST come after plate effects are removed!

        Note that this batch correction method effectively removes features
        not properly represented in both D020T0X lines as there is no way to
        estimate the batch effect on these features.

        The batch correction consists of of the median-shift and the MAD
        (median absolute deviation) stretch required to match the
        distributions. These metrics were chosen to be as robust as possible
        towards outliers.

        Note that this correction is direction-dependent, i.e. the resulting
        parameters must be applied to batch 2.

        This correction is only required for the human lines.
        :return:
        """
        if not self._loaded_line:
            print("---\n"
                  "Batch effect correction\n"
                  "must be performed on an entire line. Please use:\n"
                  "OrganoidFeatures.load_line() to load the data or set the\n"
                  "variable '_loaded_line' to True manually to bypass this\n"
                  "warning message\n"
                  "---\n")
            return None

        shift_fn = os.path.join(
            Config.BATCHEFFECTDIR, "BatchEffectCorrection_human.csv")
        if not os.path.isdir(os.path.dirname(shift_fn)):
            os.makedirs(os.path.dirname(shift_fn))
        if not os.path.isfile(shift_fn):
            batch1 = OrganoidFeatures.load_line(line="D020T01", normalize="stage_one")
            batch1.remove_plate_effects()
            batch2 = OrganoidFeatures.load_line(line="D020T02", normalize="stage_one")
            batch2.remove_plate_effects()

            batch1median = batch1.get_features().median()
            batch2median = batch2.get_features().median()
            median_shift = batch1median - batch2median

            batch1mad = batch1.get_features().apply(
                lambda col: np.abs(col - col.median()).median())
            batch2mad = batch2.get_features().apply(
                lambda col: np.abs(col - col.median()).median())
            mad_shift = batch1mad / batch2mad

            shift = pd.DataFrame(pd.concat((median_shift, mad_shift), axis=1))
            shift.columns = ("median_shift", "mad_shift")
            shift.to_csv(shift_fn)
        shift_correction = pd.read_csv(shift_fn, index_col=0)

        lines = self._metadata.Line.unique()
        for line in lines:
            if Utils.get_batch_id(line) != 2:
                continue
            if Utils.get_species(line) != "human":
                continue
            line_indices = self._metadata.Line == line

            # Stretch the feature distributions by the mad_shift factor
            # Shift the features by the median_shift factor
            median_shift = shift_correction.loc[
                self._features.columns, "median_shift"]
            self._features.loc[line_indices, :] += median_shift

            mad_shift = shift_correction.loc[
                self._features.columns, "mad_shift"]
            f = self._features.loc[line_indices, :]
            m = f.median()
            self._features.loc[line_indices, :] = ((f - m) * mad_shift) + m

    def remove_inception_batch_effects(self):
        """
        There was a notable batch effect between the first set of lines and
        the second set. The current version of batch correction involves
        calculating the per feature shift and stretch so that the full
        distributions across all wells of D020T01 and D020T02 are identical
        and then applying those transformations to each line of the second
        batch.

        This function MUST come after plate effects are removed!

        Note that this batch correction method effectively removes features
        not properly represented in both D020T0X lines as there is no way to
        estimate the batch effect on these features.

        The batch correction consists of of the median-shift and the MAD
        (median absolute deviation) stretch required to match the
        distributions. These metrics were chosen to be as robust as possible
        towards outliers.

        Note that this correction is direction-dependent, i.e. the resulting
        parameters must be applied to batch 2.

        This correction is only required for the human lines.
        :return:
        """
        if not self._loaded_line:
            print("---\n"
                  "Batch effect correction\n"
                  "must be performed on an entire line. Please use:\n"
                  "OrganoidFeatures.load_line() to load the data or set the\n"
                  "variable '_loaded_line' to True manually to bypass this\n"
                  "warning message\n"
                  "---\n")
            return None

        shift_fn = os.path.join(
            Config.BATCHEFFECTDIR, "BatchEffectCorrection_human_inception.csv")
        if not os.path.isdir(os.path.dirname(shift_fn)):
            os.makedirs(os.path.dirname(shift_fn))
        if not os.path.isfile(shift_fn):
            batch1 = OrganoidFeatures.load_inception_features(
                line="D020T01", normalize=False)
            batch1.remove_plate_effects()
            batch2 = OrganoidFeatures.load_inception_features(
                line="D020T02", normalize=False)
            batch2.remove_plate_effects()

            batch1median = batch1.get_features().median()
            batch2median = batch2.get_features().median()
            median_shift = batch1median - batch2median

            batch1mad = batch1.get_features().apply(
                lambda col: np.abs(col - col.median()).median())
            batch2mad = batch2.get_features().apply(
                lambda col: np.abs(col - col.median()).median())
            mad_shift = batch1mad / batch2mad

            shift = pd.DataFrame(pd.concat((median_shift, mad_shift), axis=1))
            shift.columns = ("median_shift", "mad_shift")
            shift.to_csv(shift_fn)
        shift_correction = pd.read_csv(shift_fn, index_col=0)

        lines = self._metadata.Line.unique()
        for line in lines:
            if Utils.get_batch_id(line) != 2:
                continue
            if Utils.get_species(line) != "human":
                continue
            line_indices = self._metadata.Line == line

            # Stretch the feature distributions by the mad_shift factor
            # Shift the features by the median_shift factor
            median_shift = shift_correction.loc[
                self._features.columns, "median_shift"]
            self._features.loc[line_indices, :] += median_shift

            mad_shift = shift_correction.loc[
                self._features.columns, "mad_shift"]
            f = self._features.loc[line_indices, :]
            m = f.median()
            self._features.loc[line_indices, :] = ((f - m) * mad_shift) + m

    def transform_features(self, max_lambda=5):
        """
        Performs a grid search for a box-cox transformation with lambda
        between 0 and N and selects the best transformation based on the
        skew, i.e. the skew closest to 0 with a pvalue > 0.05 is chosen.
        If no transformation leads to sufficient normalization then the
        feature is not transformed.

        The optimal lambda is chosen only on the DMSO controls. This is
        because there is no guarantee that the entire screen should be
        symmetric across all features, e.g. if a majority of drugs has
        no effect then the distribution for feature values should be
        heavy around DMSO/no effect.
        :param max_lambda: The maximum box-cox exponent.
        :return:
        """
        if self._verbose:
            print("Performing Box-Cox transformation ...")
        plates = sorted(set(self._metadata.Plate))
        for plate in plates:
            if self._verbose:
                print("  - Plate {}: {} of {}".format(
                    plate, plates.index(plate) + 1,
                    len(plates)))
            plate_indices = np.where(self._metadata.Plate == plate)[0]
            plate_indices_dmso = np.where(np.logical_and(
                self._metadata.Plate == plate,
                self._metadata.Drug == "DMSO"))[0]
            # Perform box-cox transform
            for ii in range(len(self._features.columns)):
                # Only the DMSO organoids, without outliers, are used to
                # determine the exponent
                plate_feat = np.array(
                    self._features.iloc[plate_indices_dmso, ii])
                plate_feat = scipy.stats.trimboth(plate_feat, 0.05)
                plate_feat -= plate_feat.min()
                plate_feat_max = plate_feat.max()
                # No need to transform if the features are constant
                if plate_feat_max == 0:
                    continue
                plate_feat /= plate_feat_max
                plate_feat += 1e-6
                # Search through the exponent space for the lambda value that
                # leads to the smallest skew (skew ~= 0)
                start_sign = np.sign(scipy.stats.skewtest(plate_feat)[0])
                if start_sign == 1:
                    min_lambda_ii = 0
                    max_lambda_ii = 1
                    delta_lambda_sign = -1
                else:
                    min_lambda_ii = 1
                    max_lambda_ii = max_lambda
                    delta_lambda_sign = 1

                cur_lambda = (max_lambda_ii - min_lambda_ii) / 2
                delta_lambda = delta_lambda_sign * (cur_lambda / 2)
                skews = []
                skews_pval = []
                lambdas = []
                while True:
                    t = scipy.stats.boxcox(
                        x=plate_feat,
                        lmbda=cur_lambda)
                    skew = scipy.stats.skewtest(t)
                    skews.append(skew[0])
                    skews_pval.append(skew[1])
                    lambdas.append(cur_lambda)
                    # If the sign doesn't change, move the exponent further
                    # from 1. If it does change, refine the search grid.
                    sign = np.sign(skew[0])
                    if start_sign == sign:
                        cur_lambda += delta_lambda
                    else:
                        cur_lambda -= delta_lambda
                    delta_lambda /= 2
                    if abs(delta_lambda) < 1e-8:
                        break

                skews_abs = [abs(s) for s in skews]
                best_i = skews_abs.index(min(skews_abs))
                # if skews_pval[best_i] > 0.05:
                plate_feat = np.array(
                    self._features.iloc[plate_indices, ii])
                plate_median = np.median(plate_feat)
                plate_mad = statsmodels.robust.mad(plate_feat)
                plate_feat -= plate_feat.min()
                plate_feat_max = plate_feat.max()
                if plate_feat_max == 0:
                    raise Exception(
                        "Critical Error @ BoxCox transform for plate"
                        "== {} and ii == {}".format(plate, ii))
                plate_feat /= plate_feat_max
                plate_feat += 1e-6
                transformed = scipy.stats.boxcox(
                    x=plate_feat + 1e-6, lmbda=lambdas[best_i])
                transformed -= np.median(transformed)
                transformed *= (plate_mad / statsmodels.robust.mad(transformed))
                transformed += plate_median
                self._features.iloc[plate_indices, ii] = transformed

    def correct_batch_effects(self, save=True):
        """
        The second preprocessing step: Fix plate and batch effects by adjusting
        features.

        This function should only be called if 'remove_bad_organoids' was
        called first OR if the 'stage_one' features are loaded
        :param save:
        :return:
        """
        if not self._loaded_line:
            print("---\n"
                  "Feature preprocessing and batch/plate effect correction\n"
                  "must be performed on an entire line. Please use:\n"
                  "OrganoidFeatures.load_line() to load the data or set the\n"
                  "variable '_loaded_line' to True manually to bypass this\n"
                  "warning message\n"
                  "---\n")
            return None

        self.remove_plate_effects()
        # self.remove_batch_effects()
        # self.transform_features()

        # Save features
        if save:
            if self._verbose:
                print("Saving normalized features ...")
            for plate in self._metadata.loc[:, "Plate"].unique():
                if self._verbose:
                    print(" - Plate {}".format(plate))
                plate_indices = self._metadata.loc[:, "Plate"] == plate
                feature_fn = os.path.join(
                    Config.FEATUREDIR, plate,
                    "{}_processedFeatures.h5".format(plate))
                if os.path.isfile(feature_fn):
                    print(
                        "Processed feature file for {} exists and "
                        "cannot be overwritten by this method".format(plate))
                    continue
                OrganoidFeatures.save_features(
                    features=self._features.loc[plate_indices, :],
                    metadata=self._metadata.loc[plate_indices, :],
                    filename=feature_fn)

    @staticmethod
    def preprocess_features(line, save=True):
        """
        The full preprocessing workflow for a given line
        :param line:
        :param save:
        :return:
        """
        obj = OrganoidFeatures.load_line(line=line, normalize="raw")
        obj.remove_bad_organoids(save=save)
        obj.correct_batch_effects(save=save)

    def preprocess_inception_features(self):
        # self.transform_features()
        self.remove_plate_effects()
        self.remove_inception_batch_effects()
        self.save_inception_features()

    def calc_well_averages(self, fn_template=None):
        """
        Aggregate the wells
        :param fn_template: A custom template for filenames. Should contain
               a single '{}' element for ''.format(...) to insert the plate
               name. Defaults to '{}_averaged_features.h5'
        :return:
        """

        # Average values of each well
        well_expected = self._features.groupby(
            by=self._features.index.values).median()
        well_expected.columns = [
            "{}_expected".format(ii) for ii in well_expected.columns]
        # well_variation = self._features.groupby(
        #     by=self._features.index.values).std()
        well_variation = self._features.groupby(
            by=self._features.index.values).apply(
            lambda df: (np.abs(df - df.median())).median())
        # well_variation = self._features.groupby(
        #     by=self._features.index.values).apply(
        #     lambda df: pd.DataFrame(
        #         data=statsmodels.robust.mad(df), index=df.columns).T)
        # well_variation.index = well_variation.index.droplevel(1)
        well_variation.columns = [
            "{}_variation".format(ii) for ii in well_variation.columns]
        well_features = pd.DataFrame(pd.concat(
            objs=(well_expected, well_variation),
            axis=1))
        well_features.loc[:, "num_objects"] = self._features.groupby(
            by=self._features.index.values).count().iloc[:, 0]

        # Create metadata
        well_id = np.array([
            "{}_{}_{}".format(
                self._metadata.loc[:, "Plate"].values[ii],
                self._metadata.loc[:, "Well"].values[ii][0],
                self._metadata.loc[:, "Well"].values[ii][1:3])
            for ii in range(len(self._metadata))])
        metadata = self._metadata.groupby(
            by=well_id).first()
        del metadata["ObjectID"]
        del metadata["OriginalX"]
        del metadata["OriginalY"]
        del metadata["Field"]
        del metadata["Size"]

        # Save platewise
        for plate in np.unique(metadata.loc[:, "Plate"].values):
            if fn_template is None:
                out_fn = os.path.join(
                    Config.FEATUREDIR, plate,
                    "{}_averaged_features.h5".format(plate))
            else:
                out_fn = os.path.join(
                    Config.FEATUREDIR, plate, fn_template.format(plate))
            if os.path.isfile(out_fn):
                print("Averaged features for {}_{} already exist".format(
                    plate, Config.FEATURETYPE))
                continue

            sel_indices = np.where(metadata.loc[:, "Plate"] == plate)[0]
            OrganoidFeatures.save_features(
                features=well_features.iloc[sel_indices, :],
                metadata=metadata.iloc[sel_indices, :],
                filename=out_fn)

    @DeprecationWarning
    def calc_feature_correlation(self):
        """
        Calculate the correlation of features
        :return:
        """

        val = self._features.T.values
        val_m = val - val.mean(1)[:, None]
        ss_val = (val_m ** 2).sum(1)
        cors = np.dot(val_m, val_m.T) / np.sqrt(
            np.dot(ss_val[:, None], ss_val[None]))
        return cors

    @DeprecationWarning
    def extract_drug_features(self, drug, outfile):
        """
        Extracts features for a given drug and saves them to a file
        :return:
        """
        feats = self._features.loc[
            self._metadata.loc[:, "Drug"] == drug, :]
        md = self._metadata.loc[
            self._metadata.loc[:, "Drug"] == drug, :]

        with h5py.File(outfile, "w-") as h5handle:
            h5handle.create_dataset(
                name="features_{}".format(Config.FEATURETYPE),
                data=feats.values.T)
            h5handle.create_dataset(
                name="feature_names_{}".format(Config.FEATURETYPE),
                data=feats.columns.values.astype(np.string_))
            h5handle.create_dataset(
                name="well_names_{}".format(Config.FEATURETYPE),
                data=feats.index.values.astype(np.string_))
            h5handle.create_dataset(
                name="metadata_{}".format(Config.FEATURETYPE),
                data=md.values.astype(np.string_))
            h5handle.create_dataset(
                name="metadata_names_{}".format(Config.FEATURETYPE),
                data=md.columns.values.astype(np.string_))

    def save_features_by_plates(self, fn_template):
        """
        Save the features to disk. Each plate is saved to a separate file.

        :param fn_template: A custom template for filenames. Should contain
               a single '{}' element for ''.format(...) to insert the plate
               name.
        :return:
        """
        if self._verbose:
            print("Saving normalized features ...")
        for plate in self._metadata.loc[:, "Plate"].unique():
            if self._verbose:
                print("- Saving plate {}".format(plate))
            plate_indices = self._metadata.loc[:, "Plate"] == plate
            feature_fn = os.path.join(
                Config.FEATUREDIR, plate,
                fn_template.format(plate))
            if os.path.isfile(feature_fn):
                print(
                    "Processed feature file for {} exists and "
                    "cannot be overwritten by this method".format(plate))
                continue
            OrganoidFeatures.save_features(
                features=self._features.loc[plate_indices, :],
                metadata=self._metadata.loc[plate_indices, :],
                filename=feature_fn)

    @DeprecationWarning
    def save_raw_features_by_plates(self):
        """
        Save the features to disk. Each plate is saved to a separate file.
        :return:
        """
        if self._verbose:
            print("Saving raw features ...")
        for plate in self._metadata.loc[:, "Plate"].unique():
            if self._verbose:
                print("- Saving plate {}".format(plate))
            plate_indices = self._metadata.loc[:, "Plate"] == plate
            feature_fn = os.path.join(
                Config.FEATUREDIR, plate,
                "{}_features.h5".format(plate))
            if os.path.isfile(feature_fn):
                print(
                    "Processed feature file for {} exists and "
                    "cannot be overwritten by this method".format(plate))
                continue
            with h5py.File(feature_fn, "w-") as h5handle:
                h5handle.create_dataset(
                    name="features_{}".format(Config.FEATURETYPE),
                    data=self._features.loc[plate_indices, :].values.T)
                h5handle.create_dataset(
                    name="feature_names_{}".format(Config.FEATURETYPE),
                    data=self._features.columns.values.astype(np.string_))
                h5handle.create_dataset(
                    name="well_names_{}".format(Config.FEATURETYPE),
                    data=self._features.loc[
                         plate_indices, :].index.values.astype(np.string_))
                h5handle.create_dataset(
                    name="metadata_{}".format(Config.FEATURETYPE),
                    data=self._metadata.loc[
                         plate_indices, :].values.astype(np.string_))
                h5handle.create_dataset(
                    name="metadata_names_{}".format(Config.FEATURETYPE),
                    data=self._metadata.columns.values.astype(np.string_))

    def save_inception_features(self):
        """
        Save the processed inception features
        :return:
        """
        plates = self._metadata.Plate.unique()
        for plate in plates:
            out_fn = os.path.join(
                Config.INCEPTIONDIR, plate,
                "{}_featuresProcessed_inception.csv".format(plate))
            md_fn = os.path.join(
                Config.INCEPTIONDIR, plate,
                "{}_metadata_inception.csv".format(plate))
            self._features.loc[self._metadata.Plate == plate, :].to_csv(out_fn)
            self._metadata.loc[self._metadata.Plate == plate, :].to_csv(md_fn)

    @staticmethod
    def save_features(features, metadata, filename):
        """
        Save processed features in a uniform format
        :param features:
        :param metadata:
        :param filename:
        :return:
        """
        if not os.path.isdir(os.path.dirname(filename)):
            os.makedirs(os.path.dirname(filename))

        with h5py.File(filename, "w-") as h5handle:
            h5handle.create_dataset(
                name="features_{}".format(Config.FEATURETYPE),
                data=features.values.T)
            h5handle.create_dataset(
                name="feature_names_{}".format(Config.FEATURETYPE),
                data=features.columns.values.astype(np.string_))
            h5handle.create_dataset(
                name="well_names_{}".format(Config.FEATURETYPE),
                data=features.index.values.astype(np.string_))
            h5handle.create_dataset(
                name="metadata_{}".format(Config.FEATURETYPE),
                data=metadata.values.astype(np.string_))
            h5handle.create_dataset(
                name="metadata_names_{}".format(Config.FEATURETYPE),
                data=metadata.columns.values.astype(np.string_))

    @staticmethod
    def load_features(filename):
        """
        Load features saved with 'save_features'

        Yes, the *.astype(np.float).astype(np.int) is on purpose.
        :param filename:
        :return:
        """
        with h5py.File(filename, "r") as h5handle:
            features = pd.DataFrame(
                data=h5handle["features_{}".format(Config.FEATURETYPE)][()].T,
                columns=h5handle["feature_names_{}".format(
                    Config.FEATURETYPE)][()].astype(np.str),
                index=h5handle["well_names_{}".format(
                    Config.FEATURETYPE)][()].astype(np.str))
            metadata = pd.DataFrame(
                data=h5handle["metadata_{}".format(Config.FEATURETYPE)][()],
                columns=h5handle["metadata_names_{}".format(
                    Config.FEATURETYPE)][()].astype(np.str),
                index=h5handle["well_names_{}".format(
                    Config.FEATURETYPE)][()].astype(np.str))

        metadata = pd.DataFrame(
            data={
                "Drug": metadata["Drug"].values.astype(np.str),
                "Concentration": metadata["Concentration"].values.astype(
                    np.str),
                "Line": metadata["Line"].values.astype(np.str),
                "Plate": metadata["Plate"].values.astype(np.str),
                "Well": metadata["Well"].values.astype(np.str),
                "Replicate": metadata["Replicate"].values.astype(
                    np.float).astype(np.int),
                "Field": metadata["Field"].values.astype(
                    np.float).astype(np.int),
                "ObjectID": metadata["ObjectID"].values.astype(
                    np.float).astype(np.int),
                "OriginalX": metadata["OriginalX"].values.astype(np.float),
                "OriginalY": metadata["OriginalY"].values.astype(np.float),
                "Size": metadata["Size"].values.astype(np.float).astype(np.int)},
            index=features.index)
        return features, metadata
