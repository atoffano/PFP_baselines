import pandas as pd
from collections import defaultdict
import os
import tqdm

DB_VERSIONS = [
    "2024_01",
    "2023_01",
    "2022_01",
    "2021_01",
    "2020_01",
    "2019_01",
    "2018_01",
    "2017_01",
    "2016_01",
    "2015_01",
    "2014_01",
    "2013_01",
    "2012_01",
    "2011_01",
    "2010_01",
    "15.0",
    "13.0",
    "10.0",
    "7.0",
    "4.0",
    "1.0",
]

def score(E):
    """
    Compute the Diamond Score for a query sequence against GO annotation terms.

    Parameters:
    E (dataframe): Dataframe of similar sequences according to Diamond blast results filtered by e-value of 0.001,
                    with columns 'subject_id', 'bit_score' and 'subject_annotations'.
    """
    go_terms = defaultdict(float)
    for _, row in E.iterrows():
        subject_annotations = row["subject_annotations"]
        bit_score = row["bit_score"]
        for annotation in subject_annotations:
            go_terms[annotation] += bit_score
    return go_terms


def alignment_score(E):
    """
    Compute the Diamond Score for a query sequence against GO annotation terms.

    Parameters:
    E (dataframe): Dataframe of similar sequences according to Diamond blast results filtered by e-value of 0.001,
                    with columns 'subject_id', 'bit_score' and 'subject_annotations'.
    """
    go_terms = score(E)
    # Normalize scores by the sum of bit scores
    for key in go_terms:
        go_terms[key] = go_terms[key] / E["bit_score"].sum()

    return go_terms


def alignment_knn(E, k=5):
    """
    Transfer annotations from the k most similar proteins based on bit score.

    Parameters:
    E (dataframe): Dataframe of similar sequences according to Diamond blast results,
                  with columns 'subject_id', 'bit_score' and 'subject_annotations'.
    k (int): Number of nearest neighbors to consider
    """

    # Sort by bit_score (descending) and get top k
    top_k = E.sort_values(by="bit_score", ascending=False).head(k)
    go_terms = score(top_k)

    # Normalize scores
    total_score = top_k["bit_score"].sum()
    if total_score > 0:  # Avoid division by zero
        for key in go_terms:
            go_terms[key] = go_terms[key] / total_score

    return go_terms


def naive_baseline(input_dir, train, val):
    # Compute the frequency of each GO term across all proteins in the training set
    go_term_counts = train["term"].value_counts()

    # Calculate the prediction score for each GO term
    go_term_scores = go_term_counts / train["EntryID"].nunique()

    # ----

    query_proteins = val["EntryID"].unique()

    # Assign the same scores to all query proteins
    predictions = pd.DataFrame(
        [
            (protein, term, score)
            for protein in query_proteins
            for term, score in go_term_scores.items()
        ],
        columns=["target_ID", "term_ID", "score"],
    )
    predictions.to_csv(
        f"{input_dir}/predictions/NaiveBaseline/predictions.tsv",
        sep="\t",
        index=False,
    )


def transfer_annotations(filtered_alignment, train, test, k_values):
    # Annotations for each sequence in the known protein set (used to transfer annotations)
    # Ts: A dictionary of annotations like {'protein_name': ['go_term1', 'go_term2', ...], ...}
    Ts = train.groupby("EntryID")["term"].apply(list).to_dict()

    filtered_alignment["subject_annotations"] = filtered_alignment["subject_id"].map(Ts)
    filtered_alignment = filtered_alignment[
        filtered_alignment["subject_annotations"].map(
            lambda x: len(x) if isinstance(x, list) else 0
        )
        > 0
    ]  # Drop rows without annotations (ie. val set proteins)

    grouped = filtered_alignment.groupby(
        "query_id"
    )  # Group by validation protein for faster processing

    ascore_pred = []
    blastknn_preds_dict = {k: [] for k in k_values}

    unaligned_proteins = 0
    unaligned_protein_ids = []  # Store unannotated protein IDs
    for protein in tqdm.tqdm(
        test["EntryID"].unique(), desc="Computing Diamond-based predictions"
    ):
        try:
            group = grouped.get_group(protein)
        except KeyError:
            unaligned_proteins += 1
            unaligned_protein_ids.append(protein)
            continue
        assert (
            not group["subject_id"].isin(test["EntryID"].unique()).any()
        ), "Annotation leakage has been found beetween validation proteins !"

        ascore_preds = alignment_score(group)  # Compute from all alignments
        for k in k_values:  # Compute from k closest alignments
            knn_preds = alignment_knn(group, k=k)

            blastknn_preds_dict[k].extend(
                [
                    {"target_ID": protein, "term_ID": term_id, "score": score}
                    for term_id, score in knn_preds.items()
                ]
            )

        ascore_pred.extend(
            [
                {"target_ID": protein, "term_ID": term_id, "score": score}
                for term_id, score in ascore_preds.items()
            ]
        )
    print(
        f"Number of unaligned proteins: {unaligned_proteins} out of {len(test['EntryID'].unique())} ({unaligned_proteins / test['EntryID'].nunique() * 100} %); No annotations have been transfered for alignment-based methods."
    )
    return unaligned_protein_ids, ascore_pred, blastknn_preds_dict


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    for db_version in tqdm.tqdm(DB_VERSIONS,
        desc="Processing databases",
    ):
        for aspect in ["BPO", "CCO", "MFO"]:
            print(f"Processing {aspect}...")
            input_dir = (
                f"./{db_version}/baselines_D1_{db_version}_{aspect}_2024_annotations"
            )
            alignment_dir = "./2024_01/diamond_swissprot_2024_01_alignment.tsv"

            os.makedirs(input_dir, exist_ok=True)
            os.makedirs(f"{input_dir}/predictions", exist_ok=True)
            os.makedirs(f"{input_dir}/predictions/AlignmentScore", exist_ok=True)
            os.makedirs(f"{input_dir}/predictions/NaiveBaseline", exist_ok=True)
            os.makedirs(f"{input_dir}/predictions/BlastKNN", exist_ok=True)

            k_values = [1, 3, 5, 10, 15, 20]  # K proteins to transfer annotations from

            # Uniprot ID mapping
            id_mapping = pd.read_csv(
                f"./2024_01/swissprot_2024_01_annotations.tsv",  # Most up to date mapping
                sep="\t",
                usecols=["EntryID", "Entry Name"],
            )
            id_mapping = id_mapping.set_index("Entry Name").to_dict()["EntryID"]

            # Proteins to annotate
            # test = pd.read_csv(
            #     f"/home/atoffano/these-antoine/data/uniprotkb/H30_{aspect}_test.tsv",
            #     sep="\t",
            # )
            test = pd.read_csv(
                f"/home/atoffano/PFP_baselines/D1_test_annotations/D1_{aspect}_test.tsv",
                sep="\t",
            )

            test = test[["EntryID"]]

            # Up to date annotations (2024_01 SwissProt version)
            train = pd.read_csv(
                f"./2024_01/swissprot_2024_01_{aspect}_annotations.tsv",
                sep="\t",
            )

            # Replace train EntryID with EntryName using id_mapping
            train["EntryID"] = train["EntryID"].map(id_mapping)
            train = train[train["EntryID"].notna()]
            train = train[
                ~train["EntryID"].isin(test["EntryID"].unique())
            ]  # Avoid leakage

            train["term"] = train["term"].str.split("; ")
            train = train.dropna().explode("term").drop_duplicates()

            print(f'Proteins in DB version {db_version}: {train["EntryID"].nunique()}')
            print(f'Proteins to annotate: {test["EntryID"].nunique()}')

            ##### Naive Baseline #####
            naive_baseline(input_dir, train, test)
            print("Naive Baseline predictions saved to", f"{input_dir}/predictions")

            ##### Alignment Methods #####
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
            # Convert Query_id and Subject_id to EntryID using id_mapping
            pairwise_alignment["query_id"] = pairwise_alignment["query_id"].map(
                id_mapping
            )
            pairwise_alignment["subject_id"] = pairwise_alignment["subject_id"].map(
                id_mapping
            )
            # Drop rows with NaN values in Query_id or Subject_id
            pairwise_alignment = pairwise_alignment[
                pairwise_alignment["query_id"].notna()
                & pairwise_alignment["subject_id"].notna()
            ]

            # Remove self-alignments
            filtered_alignment = pairwise_alignment[
                pairwise_alignment["query_id"] != pairwise_alignment["subject_id"]
            ]

            # Keep only alignments between proteins to annotate and protein that transfer annotations
            filtered_alignment = filtered_alignment[
                filtered_alignment["query_id"].isin(test["EntryID"].unique())
                & filtered_alignment["subject_id"].isin(train["EntryID"].unique())
            ]

            unaligned_protein_ids, ascore_pred, blastknn_preds_dict = (
                transfer_annotations(
                    filtered_alignment,
                    train,
                    test,
                    k_values,
                )
            )

            # Write unannotated proteins to a txt file
            unannotated_path = os.path.join(
                input_dir, f"unaligned_proteins_D1_{db_version}_{aspect}.txt"
            )
            with open(unannotated_path, "w") as f:
                for pid in unaligned_protein_ids:
                    f.write(f"{pid}\n")

            pd.DataFrame(ascore_pred).to_csv(
                f"{input_dir}/predictions/AlignmentScore/predictions.tsv",
                sep="\t",
                index=False,
            )
            for k in k_values:
                pd.DataFrame(blastknn_preds_dict[k]).to_csv(
                    f"{input_dir}/predictions/BlastKNN/k{k}_predictions.tsv",
                    sep="\t",
                    index=False,
                )

            print("Predictions saved to", f"{input_dir}/predictions")

        print("Done!")


if __name__ == "__main__":
    main()
