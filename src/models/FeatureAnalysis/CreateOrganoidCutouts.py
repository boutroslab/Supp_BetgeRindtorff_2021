# These are functions I've used for my "workflow" that I save here as a log.

from OrganoidFeatures import OrganoidFeatures
import Utils
import ImageLoader
import numpy as np
import pandas as pd
import os
import scipy.misc
import scipy.stats
import Config
import imageio


def get_typical_organoids(line, drug, concentration=None, sizethresh=75, num_organoids=10):
    """
    Get typical organoids for a given line+drug+concentration.

    If 'concentration' is None, then there is no filtering for concentration

    This is primarily for organoid visualization. Smaller organoids aren't as
    "pretty" as large ones, so only the top (100-sizethresh)% organoids with
    regards to size are used for the detection of typical organoids.
    'sizethresh' is passed directly to np.percentile(..., sizethresh).

    :param line:
    :param drug:
    :param concentration
    :param sizethresh:
    :param num_organoids: The number of organoids to extract.
    :return:
    """

    # Load features
    obj = OrganoidFeatures.load_line(
        line=line,
        normalize="normalized",
        drug=drug)

    # Keep only the entries with the corresponding concentration
    if concentration is None:
        features = obj.get_features()
        metadata = obj.get_metadata()
    else:
        sel_indices = obj.get_metadata().Concentration == concentration
        features = obj.get_features().loc[sel_indices, :]
        metadata = obj.get_metadata().loc[sel_indices, :]

    # Keep only the top n% of organoids for the typical organoid extraction
    q = np.percentile(metadata.Size.values, sizethresh)
    sel_indices = metadata.Size >= q
    filtered_features = features.loc[sel_indices, :]
    filtered_metadata = metadata.loc[sel_indices, :]

    # Get typical features
    feature_entry, metadata_entry = Utils.get_characteristic_organoids(
        features=filtered_features,
        metadata=filtered_metadata,
        num_organoids=num_organoids)

    return feature_entry, metadata_entry
    

def get_small_organoids(line, drug, concentration=None, sizethresh=75, num_organoids=10):
    """
    Get typical organoids for a given line+drug+concentration.

    If 'concentration' is None, then there is no filtering for concentration

    This is primarily for organoid visualization. Smaller organoids aren't as
    "pretty" as large ones, so only the top (100-sizethresh)% organoids with
    regards to size are used for the detection of typical organoids.
    'sizethresh' is passed directly to np.percentile(..., sizethresh).

    :param line:
    :param drug:
    :param concentration
    :param sizethresh:
    :param num_organoids: The number of organoids to extract.
    :return:
    """

    # Load features
    obj = OrganoidFeatures.load_line(
        line=line,
        normalize="normalized",
        drug=drug)

    # Keep only the entries with the corresponding concentration
    if concentration is None:
        features = obj.get_features()
        metadata = obj.get_metadata()
    else:
        sel_indices = obj.get_metadata().Concentration == concentration
        features = obj.get_features().loc[sel_indices, :]
        metadata = obj.get_metadata().loc[sel_indices, :]

    # Keep only the top n% of organoids for the typical organoid extraction
    q = np.percentile(metadata.Size.values, sizethresh)
    sel_indices = metadata.Size <= q
    filtered_features = features.loc[sel_indices, :]
    filtered_metadata = metadata.loc[sel_indices, :]

    # Get typical features
    feature_entry, metadata_entry = Utils.get_characteristic_organoids(
        features=filtered_features,
        metadata=filtered_metadata,
        num_organoids=num_organoids)

    return feature_entry, metadata_entry


def get_typical_organoids_for_all_drugs(
        line, outdir, sizethresh=75, num_organoids=10):
    """
    Generate example organoids for all drug treatments.

    This is primarily for organoid visualization. Smaller organoids aren't as
    "pretty" as large ones, so only the top 1-'sizethresh'/100 organoids are
    used for the detection of typical organoids. 'sizethresh' is passed
    directly to np.percentile(..., sizethresh)

    :param line:
    :param outdir:
    :param sizethresh:
    :param num_organoids: The number of organoids to extract.
    :return:
    """

    line_metadata = OrganoidFeatures.load_metadata(
        line=line, normalize="normalized").get_metadata()

    drugs = line_metadata.loc[:, ("Drug", "Concentration")]
    drugs = drugs.loc[~drugs.duplicated(), :]

    for ii in range(len(drugs)):
        print(ii+1, "of", len(drugs))
        drug, concentration = drugs.iloc[ii, :]

        if drug in ("DMSO", "Staurosporine_500nM"):
            concentration = None

        feature_entry, metadata_entry = get_typical_organoids(
            line=line, drug=drug, concentration=concentration,
            sizethresh=sizethresh, num_organoids=num_organoids)

        # Create composite directory name if there is a concentration
        if concentration in ("nan", None):
            dirname = drug
        else:
            dirname = drug + "__" + concentration

        for jj in range(len(metadata_entry)):
            out_fn = os.path.join(
                outdir, line, dirname, "{}.tif".format(jj+1))
            if not os.path.isdir(os.path.dirname(out_fn)):
                os.makedirs(os.path.dirname(out_fn))
            entry = metadata_entry.iloc[jj, :]
            organoid_img = ImageLoader.get_organoid_image(metadata_entry=entry)
            organoid_img = (organoid_img * 255).astype(np.uint8)
            imageio.imwrite(out_fn, organoid_img)


def get_typical_organoids_for_dmso(
        line, outdir, sizethresh=75, num_organoids=10):
    """
    Generate example organoids for all drug treatments.

    This is primarily for organoid visualization. Smaller organoids aren't as
    "pretty" as large ones, so only the top 1-'sizethresh'/100 organoids are
    used for the detection of typical organoids. 'sizethresh' is passed
    directly to np.percentile(..., sizethresh)

    :param line:
    :param outdir:
    :param sizethresh:
    :param num_organoids: The number of organoids to extract.
    :return:
    """

    feature_entry, metadata_entry = get_typical_organoids(
        line=line, drug="DMSO", concentration=None,
        sizethresh=sizethresh, num_organoids=num_organoids)

    for jj in range(len(metadata_entry)):
        out_fn = os.path.join(
            outdir, line, "DMSO", "{}.tif".format(jj+1))
        if not os.path.isdir(os.path.dirname(out_fn)):
            os.makedirs(os.path.dirname(out_fn))
        entry = metadata_entry.iloc[jj, :]
        organoid_img = ImageLoader.get_organoid_image(metadata_entry=entry)
        organoid_img = (organoid_img * 255).astype(np.uint8)
        imageio.imwrite(out_fn, organoid_img)

def get_typical_organoids_for_dmso_detailed(
        line, outdir, sizethresh=0, num_organoids=100):
    """
    Generate example organoids for all drug treatments.

    This is primarily for organoid visualization. Smaller organoids aren't as
    "pretty" as large ones, so only the top 1-'sizethresh'/100 organoids are
    used for the detection of typical organoids. 'sizethresh' is passed
    directly to np.percentile(..., sizethresh)

    :param line:
    :param outdir:
    :param sizethresh:
    :param num_organoids: The number of organoids to extract.
    :return:
    """

    feature_entry, metadata_entry = get_typical_organoids(
        line=line, drug="DMSO", concentration=None,
        sizethresh=sizethresh, num_organoids=num_organoids)

    for jj in range(len(metadata_entry)):
        entry = metadata_entry.iloc[jj, :]
        plate = metadata_entry.iloc[jj, 3]
        well = metadata_entry.iloc[jj, 4]
        field = metadata_entry.iloc[jj, 6]
        id = metadata_entry.iloc[jj, 7]
        
        print(plate)
        print(well)
        print(field)
        print(id)
        
        out_fn = os.path.join(
            outdir, line, "DMSO","{}_{}_{}_{}.tif".format(str(plate), str(well), str(field), str(id)))
        
        if not os.path.isdir(os.path.dirname(out_fn)):
            os.makedirs(os.path.dirname(out_fn))
            
        organoid_img = ImageLoader.get_organoid_image(metadata_entry=entry)
        organoid_img = (organoid_img * 255).astype(np.uint8)
        imageio.imwrite(out_fn, organoid_img)
        

def get_small_organoids_for_dmso_detailed(
        line, outdir, sizethresh=10, num_organoids=100):
    """
    Generate example organoids for all drug treatments.

    This is primarily for organoid visualization. Smaller organoids aren't as
    "pretty" as large ones, so only the top 1-'sizethresh'/100 organoids are
    used for the detection of typical organoids. 'sizethresh' is passed
    directly to np.percentile(..., sizethresh)

    :param line:
    :param outdir:
    :param sizethresh:
    :param num_organoids: The number of organoids to extract.
    :return:
    """

    feature_entry, metadata_entry = get_small_organoids(
        line=line, drug="DMSO", concentration=None,
        sizethresh=sizethresh, num_organoids=num_organoids)

    for jj in range(len(metadata_entry)):
        entry = metadata_entry.iloc[jj, :]
        plate = metadata_entry.iloc[jj, 3]
        well = metadata_entry.iloc[jj, 4]
        field = metadata_entry.iloc[jj, 6]
        id = metadata_entry.iloc[jj, 7]
        
        print(plate)
        print(well)
        print(field)
        print(id)
        
        out_fn = os.path.join(
            outdir, line, "DMSO","{}_{}_{}_{}.tif".format(str(plate), str(well), str(field), str(id)))
        
        if not os.path.isdir(os.path.dirname(out_fn)):
            os.makedirs(os.path.dirname(out_fn))
            
        organoid_img = ImageLoader.get_organoid_image(metadata_entry=entry)
        organoid_img = (organoid_img * 255).astype(np.uint8)
        imageio.imwrite(out_fn, organoid_img)
        

def get_typical_organoids_for_gsk_drugs_in_kistem(
        outdir, sizethresh=75, num_organoids=10):
    """
    Generate example organoids for all drugs that target GSK for all lines

    This is primarily for organoid visualization. Smaller organoids aren't as
    "pretty" as large ones, so only the top 1-'sizethresh'/100 organoids are
    used for the detection of typical organoids. 'sizethresh' is passed
    directly to np.percentile(..., sizethresh)

    :param outdir:
    :param sizethresh:
    :param num_organoids: The number of organoids to extract.
    :return:
    """

    # These are the hard-coded drugs that target GSK according to the
    # library annotation
    drugs = (
        "1-Azakenpaullone", "AR-A014418", "AZD1080", "AZD2858", "Bikinin",
        "BIO", "CHIR-98014", "CHIR-99021 (CT99021)", "IM-12", "Indirubin",
        "LY2090314", "SB216763", "SB415286", "TDZD-8", "Tideglusib",
        "TWS119")

    for line in Utils.get_all_lines("human"):
        for drug in drugs:
            feature_entry, metadata_entry = get_typical_organoids(
                line=line, drug=drug, concentration="nan",
                sizethresh=sizethresh, num_organoids=num_organoids)

            for jj in range(len(metadata_entry)):
                out_fn = os.path.join(
                    outdir, line, drug, "{}.tif".format(jj+1))
                if not os.path.isdir(os.path.dirname(out_fn)):
                    os.makedirs(os.path.dirname(out_fn))
                entry = metadata_entry.iloc[jj, :]
                organoid_img = ImageLoader.get_organoid_image(metadata_entry=entry)
                organoid_img = (organoid_img * 255).astype(np.uint8)
                imageio.imwrite(out_fn, organoid_img)


def get_typical_organoids_for_pi3k(outdir, sizethresh=75, num_organoids=10):
    # These are the hard-coded drugs that target GSK according to the
    # library annotation
    drugs = (
        "Tyrphostin AG 879", "K02288", "PF-573228", "WIKI4", "AT7867",
        "WYE-354")

    for line in Utils.get_all_lines("human"):
        for drug in drugs:
            feature_entry, metadata_entry = get_typical_organoids(
                line=line, drug=drug, concentration="nan",
                sizethresh=sizethresh, num_organoids=num_organoids)

            for jj in range(len(metadata_entry)):
                out_fn = os.path.join(
                    outdir, line, drug, "{}.tif".format(jj + 1))
                if not os.path.isdir(os.path.dirname(out_fn)):
                    os.makedirs(os.path.dirname(out_fn))
                entry = metadata_entry.iloc[jj, :]
                organoid_img = ImageLoader.get_organoid_image(metadata_entry=entry)
                organoid_img = (organoid_img * 255).astype(np.uint8)
                imageio.imwrite(out_fn, organoid_img)


@PendingDeprecationWarning
# TODO: Alter this to call get_typical_organoids_for_drug(...)
def get_typical_organoids_for_enriched_clusters(outdir):
    """
    Using the results of the drug effects vignette, this function extracts
    typical organoids for every active drug.

    Note that it extracts sample organoids for ALL drugs, whether they're
    enriched or not
    :param outdir:
    :return:
    """

    enriched_drugs = pd.DataFrame.from_csv(
        path="/Users/jansauer/Thesis/Projects/PROMISE/FeatureAnalysis/"
             "drug_effects/human/EnrichedDrugs.csv", sep=";")
    # Add DMSO for each line
    dmso_df = pd.DataFrame(
        data={"Line": enriched_drugs.Line.unique(),
              "Pathway": np.nan,
              "Target": np.nan,
              "Enriched.Target": np.nan,
              "is.pathway.enriched": False,
              "is.target.enriched": False},
        index=("DMSO", )*len(enriched_drugs.Line.unique()))
    # Add pos. controls for each line
    borte_df = pd.DataFrame(
        data={"Line": enriched_drugs.Line.unique(),
              "Pathway": np.nan,
              "Target": np.nan,
              "Enriched.Target": np.nan,
              "is.pathway.enriched": False,
              "is.target.enriched": False},
        index=("Bortezomib",) * len(enriched_drugs.Line.unique()))
    irino_df = pd.DataFrame(
        data={"Line": enriched_drugs.Line.unique(),
              "Pathway": np.nan,
              "Target": np.nan,
              "Enriched.Target": np.nan,
              "is.pathway.enriched": False,
              "is.target.enriched": False},
        index=("Irinotecan / SN-38",) * len(enriched_drugs.Line.unique()))

    enriched_drugs = enriched_drugs.append(other=dmso_df)
    enriched_drugs = enriched_drugs.append(other=borte_df)
    enriched_drugs = enriched_drugs.append(other=irino_df)

    enriched_drugs = pd.DataFrame(
        data={"Line": "D021T01",
              "Pathway": np.nan,
              "Target": np.nan,
              "Enriched.Target": np.nan,
              "is.pathway.enriched": False,
              "is.target.enriched": False},
        index=("DMSO",))

    # enriched_drugs = dmso_df
    # enriched_drugs = pd.DataFrame(pd.concat((borte_df, irino_df)))
    for ii in range(len(enriched_drugs)):
        print(ii+1, "of", len(enriched_drugs))
        ed_row = enriched_drugs.iloc[ii, ]

        obj = OrganoidFeatures.load_line(
            line=ed_row.Line,
            normalize="normalized",
            drug=ed_row.name)
        q = np.percentile(obj.get_metadata().Size.values, 75)
        sel_indices = obj.get_metadata().Size >= q
        features = obj.get_features().loc[sel_indices, :]
        metadata = obj.get_metadata().loc[sel_indices, :]
        # If pos. controls use only top concentrations
        if ed_row.name in ("Bortezomib", "Irinotecan / SN-38"):
            sel_indices = np.in1d(metadata.Concentration, ("1.0", "0.2"))
            features = features.loc[sel_indices, :]
            metadata = metadata.loc[sel_indices, :]
        feature_entry, metadata_entry = Utils.get_characteristic_organoids(
            features=features, metadata=metadata, num_organoids=10)

        for jj in range(len(metadata_entry)):
            out_fn = os.path.join(
                outdir, ed_row.Line, ed_row.name, "{}.tif".format(jj+1))
            if not os.path.isdir(os.path.dirname(out_fn)):
                os.makedirs(os.path.dirname(out_fn))
            entry = metadata_entry.iloc[jj, :]
            organoid_img = ImageLoader.get_organoid_image(metadata_entry=entry)
            organoid_img = (organoid_img * 255).astype(np.uint8)
            imageio.imwrite(out_fn, organoid_img)


@PendingDeprecationWarning
# TODO: Alter this to call get_typical_organoids_for_drug(...)
def get_typical_organoids_for_drugs_at_all_concentrations(outdir):
    """
    Generate example organoid images for custom drugs at all concentrations
    from the CCP/L08 library
    :param outdir:
    :return:
    """
    entries = [
        # ("D046T01", "Bortezomib"),
        ("D020T02", "Bortezomib"),
        ("D021T01", "Panobinostat"),
        ("D021T01", "Doxorubicin"),
        ("D018T01", "Methotrexate")]

    for line, drug in entries:
        print(line, drug)

        obj = OrganoidFeatures.load_line(
            line=line,
            normalize="normalized",
            drug=drug)

        # Keep only L08
        layout = np.array([s[11:14] for s in obj.get_metadata().Plate])
        features = obj.get_features().loc[layout == "L08", :]
        metadata = obj.get_metadata().loc[layout == "L08", :]

        # Loop through concentrations
        concentrations = metadata.Concentration.unique()
        for c in concentrations:
            conc_features = features.loc[metadata.Concentration == c, :]
            conc_metadata = metadata.loc[metadata.Concentration == c, :]

            q = np.percentile(conc_metadata.Size.values, 75)
            sel_indices = conc_metadata.Size >= q
            conc_features = conc_features.loc[sel_indices, :]
            conc_metadata = conc_metadata.loc[sel_indices, :]

            feature_entry, metadata_entry = Utils.get_characteristic_organoids(
                features=conc_features, metadata=conc_metadata, num_organoids=5)

            for jj in range(len(metadata_entry)):
                out_fn = os.path.join(
                    outdir, line, "{}__{}".format(drug, c), "{}.tif".format(jj + 1))
                if not os.path.isdir(os.path.dirname(out_fn)):
                    os.makedirs(os.path.dirname(out_fn))
                entry = metadata_entry.iloc[jj, :]
                organoid_img = ImageLoader.get_organoid_image(
                    metadata_entry=entry)
                organoid_img = (organoid_img * 255).astype(np.uint8)
                imageio.imwrite(out_fn, organoid_img)
