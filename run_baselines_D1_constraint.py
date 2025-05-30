import pandas as pd
import os
from run_baselines import naive_baseline, transfer_annotations


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    # Print cwd
    print("Current working directory:", os.getcwd())
    for aspect in ["BPO", "CCO", "MFO"]:
        print(f"Processing {aspect}...")
        input_dir = f"./D1_constrained_experiment/baselines_D1_{aspect}"
        alignment_dir = "./2024_01/diamond_swissprot_2024_01_alignment.tsv"  # Most up to date alignment file

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
        test = pd.read_csv(
            f"./D1_test_annotations/D1_{aspect}_test.tsv",
            sep="\t",
        )
        test = test[["EntryID"]]

        # Annotations
        # - BeProf D1 training set for the constrained experiment
        train = pd.read_csv(
            f"./BeProf_D1/D1_{aspect}.tsv",
            sep="\t",
        )
        # Avoid leakage
        train = train[~train["EntryID"].isin(test["EntryID"].unique())]
        # Split col 'term' by '; ' and explode it
        train["term"] = train["term"].str.split("; ")
        train = train.dropna().explode("term").drop_duplicates()

        print(f'Proteins in BeProf D1 Train set: {train["EntryID"].nunique()}')
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
        filtered_alignment = pairwise_alignment[
            pairwise_alignment["query_id"] != pairwise_alignment["subject_id"]
        ]

        # Keep only alignments between protein to annotate and protein that transfer annotations
        filtered_alignment = filtered_alignment[
            filtered_alignment["query_id"].isin(test["EntryID"].unique())
            & filtered_alignment["subject_id"].isin(train["EntryID"].unique())
        ]

        unaligned_protein_ids, ascore_pred, blastknn_preds_dict = transfer_annotations(
            filtered_alignment,
            train,
            test,
            k_values,
        )

        # Write unannotated proteins to a txt file
        unannotated_path = os.path.join(
            input_dir, f"unaligned_proteins_D1_{aspect}.txt"
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
