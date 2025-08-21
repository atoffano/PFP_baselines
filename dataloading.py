import pandas as pd


def load_uniprot_mapping():
    """
    Map Uniprot IDs in the DataFrame using the provided id_mapping dictionary.
    """
    id_mapping = pd.read_csv(
        f"./data/swissprot/2024_01/swissprot_2024_01_annotations.tsv",  # Most up to date mapping
        sep="\t",
        usecols=["EntryID", "Entry Name"],
    )
    id_mapping = id_mapping.set_index("Entry Name").to_dict()["EntryID"]
    return id_mapping


def load_pairwise_alignment(dataset, id_mapping=None):
    """
    Load pairwise alignment data and map Query_id and Subject_id to EntryID using id_mapping
    """
    # match method:
    #     case "stringdb":
    #         pairwise_alignment = pd.read_csv(
    #             "./data/stringdb/swissprot_stringdb.tsv",
    #             sep="\t",
    #             header=True,
    #             names=[
    #                 "protein1",
    #                 "protein2",
    #                 "neighborhood",
    #                 "fusion",
    #                 "cooccurence",
    #                 "coexpression",
    #                 "experimental",
    #                 "database",
    #                 "textmining",
    #                 "combined_score",
    #             ],
    #         )
    #         # Keep only relevant columns
    #         pairwise_alignment = pairwise_alignment[
    #             ["protein1", "protein2", "combined_score"]
    #         ].rename(
    #             columns={
    #                 "protein1": "query_id",
    #                 "protein2": "subject_id",
    #                 "combined_score": "score",
    #             }
    #         )

    #         if dataset in ["ATGO"]:
    #             # Apply reverse id_mapping, as ATGO uses Entry Name as protein IDs
    #             # StringDB uses EntryID (004R_FRG3G) as protein IDs
    #             reversed_id_mapping = {v: k for k, v in id_mapping.items()}
    #             pairwise_alignment["query_id"] = pairwise_alignment["query_id"].apply(
    #                 lambda x: reversed_id_mapping.get(x, x)
    #             )
    #             pairwise_alignment["subject_id"] = pairwise_alignment[
    #                 "subject_id"
    #             ].apply(lambda x: reversed_id_mapping.get(x, x))

    #     case "diamond":
    # Load pairwise alignment from Diamond output ran on SwissProt 2024_01
    pairwise_alignment = pd.read_csv(
        "./data/swissprot/2024_01/diamond_swissprot_2024_01_alignment.tsv",
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
    if dataset in ["D1", "CAFA3"]:
        # Diamond output uses EntryName (e.g. Q6GZX1) as protein IDs; Convert to EntryID
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
    # This is required to avoid self-annotation transfer
    pairwise_alignment = pairwise_alignment[
        pairwise_alignment["query_id"] != pairwise_alignment["subject_id"]
    ]

    return pairwise_alignment


def load_data(
    logger,
    dataset,
    aspect,
    db_version,
    annotations_2024_01=None,
    id_mapping=None,
    experimental_only=False,
    one_vs_all=False,
):
    logger.info(f"Loading SwissProt {db_version} annotations for training set")
    if db_version == "":
        # Load constrained dataset
        train = pd.read_csv(
            f"./data/{dataset}/{dataset}_{aspect}_train_annotations.tsv",
            sep="\t",
        )
    elif experimental_only:
        # Load experimental annotations only
        train = pd.read_csv(
            f"./data/swissprot/{db_version}/swissprot_{db_version}_{aspect}_exp_annotations.tsv",
            sep="\t",
        )
    else:
        train = pd.read_csv(
            f"./data/swissprot/{db_version}/swissprot_{db_version}_{aspect}_annotations.tsv",
            sep="\t",
        )

    # Proteins to annotate
    test = pd.read_csv(
        f"./data/{dataset}/{dataset}_{aspect}_test_annotations.tsv",
        sep="\t",
    )

    if not one_vs_all:
        test = test[["EntryID"]]  # Drop terms; Makes sure no leakage occurs
        train = train[
            ~train["EntryID"].isin(test["EntryID"].unique())
        ]  # To *really* make sure

    if annotations_2024_01:
        logger.info("Fixing train proteins' annotations to 2024 SwissProt version...")
        if experimental_only:
            annotations_2024_01 = pd.read_csv(
                f"./data/swissprot/2024_01/swissprot_2024_01_{aspect}_exp_annotations.tsv",
                sep="\t",
            )
        else:
            # Fix train proteins' annotations to 2024 SwissProt version
            annotations_2024_01 = pd.read_csv(
                f"./data/swissprot/2024_01/swissprot_2024_01_{aspect}_annotations.tsv",
                sep="\t",
            )
        # Get rows in annotations_2024_01 where EntryID is in train
        train = annotations_2024_01[
            annotations_2024_01["EntryID"].isin(train["EntryID"])
        ]

    if dataset in ["D1", "CAFA3"]:
        logger.info("Mapping train SwissProt EntryID to EntryName...")
        train["EntryID"] = train["EntryID"].map(id_mapping).fillna(train["EntryID"])

    train["term"] = train["term"].str.split("; ")
    train = train.dropna().explode("term").drop_duplicates()

    return train, test
