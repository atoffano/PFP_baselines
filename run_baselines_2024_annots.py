import pandas as pd
import os
import tqdm
from run_baselines import naive_baseline, transfer_annotations


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    for db_version in tqdm.tqdm(
        [
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
            "2009_03",
            "2008_01",
            "2007_03",
            "2006_02",
            "2005_01",
            "2003_12",
        ],
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

            ## --- 2024_annots
            # Fixed annotations to 2024 - Evaluate impact of alignment
            # Version-specific SwissProt annotations
            db_version_annots = pd.read_csv(
                f"./{db_version}/swissprot_{db_version}_{aspect}_annotations.tsv",
                sep="\t",
            )
            db_version_annots["EntryID"] = db_version_annots["EntryID"].map(id_mapping)
            db_version_proteins = db_version_annots["EntryID"].unique()
            # Filter train to only include proteins in db_version_annots
            train = train[
                train["EntryID"].isin(db_version_proteins)
            ]  # Filter to include only proteins in the current db_version_annots
            # ---

            train = train[train["term"].notna()]
            train["term"] = train["term"].apply(lambda x: ast.literal_eval(x))
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

            # Keep only alignments between protein to annotate and protein that transfer annotations
            filtered_alignment = filtered_alignment[
                filtered_alignment["query_id"].isin(test["EntryID"].unique())
                & filtered_alignment["subject_id"].isin(train["EntryID"].unique())
            ]

            ## --- 2024_annots
            # Filter alignments to keep subjects in db_version_proteins
            print(
                f"Number of proteins before 2024_annotations filtering: {filtered_alignment['subject_id'].nunique()}"
            )
            filtered_alignment = filtered_alignment[
                filtered_alignment["subject_id"].isin(db_version_proteins)
            ]  # Filter to include only proteins in the current db_version_annots
            print(
                f"Number of proteins after 2024_annotations filtering: {filtered_alignment['subject_id'].nunique()}"
            )
            ## ---

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
