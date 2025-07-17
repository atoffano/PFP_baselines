import pandas as pd
import os


def load_uniprot_mapping():
    """
    Map Uniprot IDs in the DataFrame using the provided id_mapping dictionary.
    """
    id_mapping = pd.read_csv(
        f"./2024_01/swissprot_2024_01_annotations.tsv",  # Most up to date mapping
        sep="\t",
        usecols=["EntryID", "Entry Name"],
    )
    id_mapping = id_mapping.set_index("Entry Name").to_dict()["EntryID"]
    return id_mapping


def load_pairwise_alignment(alignment_dir, id_mapping=None):
    """
    Load pairwise alignment data and map Query_id and Subject_id to EntryID using id_mapping
    """
    pairwise_alignment = pd.read_csv(
        alignment_dir,
        sep="\t",
        header=None,
        names=[
            "query_id",
            "subject_id",
            "perc_identity",
            "align_length",
            "mismatches",
            "gap_opens",
            "q_start",
            "q_end",
            "s_start",
            "s_end",
            "e_value",
            "bit_score",
        ],
    )

    ##### Load Uniprot ID mapping #####
    if id_mapping is not None:

        # Convert Query_id and Subject_id to EntryID using id_mapping
        pairwise_alignment["query_id"] = pairwise_alignment["query_id"].map(id_mapping)
        pairwise_alignment["subject_id"] = pairwise_alignment["subject_id"].map(
            id_mapping
        )

    # Drop rows with NaN values in Query_id or Subject_id
    pairwise_alignment = pairwise_alignment[
        pairwise_alignment["query_id"].notna()
        & pairwise_alignment["subject_id"].notna()
    ]

    # Remove self-alignments
    pairwise_alignment = pairwise_alignment[
        pairwise_alignment["query_id"] != pairwise_alignment["subject_id"]
    ]

    return pairwise_alignment


def load_swissprot(
    dataset,
    aspect,
    db_version,
    id_mapping=None,
):

    train = pd.read_csv(
        f"./{db_version}/swissprot_{db_version}_{aspect}_annotations.tsv",
        sep="\t",
    )

    # Uncomment to constrain the proteins from the full SwissProt to the original BeProf train proteins
    # train = pd.read_csv(
    #     f"./BeProf_D1/D1_{aspect}.tsv",
    #     sep="\t",
    # )

    # Proteins to annotate
    test = pd.read_csv(
        f"./{dataset}_test_annotations/{dataset}_{aspect}_test.tsv",
        sep="\t",
    )
    test = test[["EntryID"]]

    # Replace train EntryID with EntryName; Apply id_mapping dict
    if id_mapping is not None:
        # Apply id_mapping. if value not in id_mapping, leave as is
        train["EntryID"] = train["EntryID"].map(id_mapping).fillna(train["EntryID"])
        print("Mapping train SwissProt EntryID to EntryName...")
    train = train[~train["EntryID"].isin(test["EntryID"].unique())]  # Avoid leakage

    train["term"] = train["term"].str.split("; ")
    train = train.dropna().explode("term").drop_duplicates()

    return train, test


def load_atgo(
    aspect,
    db_version,
):
    """
    Load ATGO dataset and return train and test DataFrames.
    """
    base_dir = "./atgo/expanse/lustre/projects/mia174/yi930505/DeepTransformer/CAFA3/CAFA3/Benchmark_Dataset/Ours"

    # train = pd.read_csv(
    #     f"./{db_version}/swissprot_{db_version}_{aspect}_annotations.tsv",
    #     sep="\t",
    # )
    # train["term"] = train["term"].str.split("; ")

    # Uncomment to constrain the proteins from the full SwissProt to the original ATGO train proteins
    train = os.path.join(base_dir, aspect[:2], "train_protein_label")
    train = pd.read_csv(train, sep="  ", names=["EntryID", "term"])
    train["term"] = train["term"].str.split(",")

    test = os.path.join(base_dir, aspect[:2], "test_protein_label")
    test = pd.read_csv(test, sep="  ", header=None, names=["EntryID", "term"])
    test = test[["EntryID"]]  # Avoid leakage

    train = train.dropna().explode("term").drop_duplicates()
    train = train[~train["EntryID"].isin(test["EntryID"].unique())]  # Avoid leakage

    return train, test
