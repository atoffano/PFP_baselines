import pandas as pd
import os
import tqdm
import argparse
import logging
from dataloading import *
import methods
import evaluation

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


def setup_logging(output_dir, aspect):
    """Setup logging for each aspect with separate log files."""
    log_dir = os.path.join(output_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    logger = logging.getLogger(f"{aspect}")
    logger.setLevel(logging.INFO)
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    log_file = os.path.join(log_dir, f"{aspect}.log")
    file_handler = logging.FileHandler(log_file, mode="w")
    file_handler.setLevel(logging.INFO)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger


def main():
    parser = argparse.ArgumentParser(
        description="Run baseline annotation transfer methods."
    )
    parser.add_argument("--dataset", type=str, required=True, help="Dataset name.")
    parser.add_argument(
        "--alignment_dir",
        type=str,
        required=True,
        help="Path to pairwise alignment file.",
    )
    parser.add_argument(
        "--k_values",
        type=int,
        nargs="+",
        default=[1, 3, 5, 10, 15, 20],
        help="List of k closest sequences to transfer annotations from.",
    )
    parser.add_argument(
        "--db_versions",
        type=str,
        nargs="+",
        default=DB_VERSIONS,
        help="SwissProt DB versions to process.",
    )
    parser.add_argument(
        "--aspects",
        type=str,
        nargs="+",
        default=["BPO", "CCO", "MFO"],
        help="Ontology aspects to process.",
    )
    parser.add_argument(
        "--id_mapping", action="store_true", help="Use Uniprot ID mapping."
    )

    args = parser.parse_args()

    if args.id_mapping:
        print("Loading Uniprot ID mapping...")
        id_mapping = load_uniprot_mapping()
    else:
        id_mapping = None

    for db_version in tqdm.tqdm(args.db_versions, desc="Processing databases"):
        for aspect in args.aspects:
            output_dir = f"./results/{args.dataset}/baselines_{args.dataset}_{db_version}_{aspect}"
            os.makedirs(output_dir, exist_ok=True)
            os.makedirs(f"{output_dir}/predictions", exist_ok=True)
            os.makedirs(f"{output_dir}/predictions/AlignmentScore", exist_ok=True)
            os.makedirs(f"{output_dir}/predictions/NaiveBaseline", exist_ok=True)
            os.makedirs(f"{output_dir}/predictions/BlastKNN", exist_ok=True)

            # Setup logging for this aspect
            logger = setup_logging(output_dir, aspect)
            logger.info(
                f"Starting processing for {aspect} with database version {db_version}"
            )

            # Load data
            logger.info(f"Loading data for {args.dataset} with aspect {aspect}")
            train, test = load_data(
                args.dataset,
                aspect,
                db_version,
                id_mapping=id_mapping,
            )
            logger.info(
                f"Loaded {train['EntryID'].nunique()} training proteins and {test['EntryID'].nunique()} test proteins"
            )
            logger.info(f"Train set:\n{train}")
            logger.info(f"Test set:\n{test}")

            # Comment to skip Naive Baseline
            logger.info("Running Naive Baseline...")
            methods.naive_baseline(output_dir, train, test)
            logger.info(f"Naive Baseline predictions saved to {output_dir}/predictions")
            # ---

            logger.info("Loading pairwise alignments...")
            pairwise_alignment = load_pairwise_alignment(
                args.alignment_dir, id_mapping=id_mapping
            )
            pairwise_alignment = pairwise_alignment[
                pairwise_alignment["query_id"].isin(test["EntryID"].unique())
                & pairwise_alignment["subject_id"].isin(train["EntryID"].unique())
            ]
            logger.info(f"Loaded {len(pairwise_alignment)} pairwise alignments")

            logger.info("Running alignment-based methods...")
            unaligned_protein_ids, ascore_pred, blastknn_preds_dict = (
                methods.transfer_annotations(
                    logger,
                    pairwise_alignment,
                    train,
                    test,
                    args.k_values,
                )
            )

            logger.info(f"Found {len(unaligned_protein_ids)} unaligned proteins")

            unannotated_path = os.path.join(
                output_dir,
                f"unaligned_proteins_{args.dataset}_{db_version}_{aspect}.txt",
            )
            with open(unannotated_path, "w") as f:
                for pid in unaligned_protein_ids:
                    f.write(f"{pid}\n")

            pd.DataFrame(ascore_pred).to_csv(
                f"{output_dir}/predictions/AlignmentScore/predictions.tsv",
                sep="\t",
                index=False,
            )
            logger.info(f"Saved {len(ascore_pred)} AlignmentScore predictions")

            for k in args.k_values:
                pred_count = len(blastknn_preds_dict[k])
                pd.DataFrame(blastknn_preds_dict[k]).to_csv(
                    f"{output_dir}/predictions/BlastKNN/k{k}_predictions.tsv",
                    sep="\t",
                    index=False,
                )
                logger.info(f"Saved {pred_count} BlastKNN predictions for k={k}")

            logger.info(f"All predictions saved to {output_dir}/predictions")
            logger.info(f"Completed processing for {aspect}")

            logger.info("Evaluating predictions...")

            evaluation.evaluate(
                logger, output_dir, args.dataset, aspect, k_values=args.k_values
            )
            logger.info(f"Evaluation completed for aspect {aspect}")

        print("Done!")


if __name__ == "__main__":
    main()

    # Example usage:
    # python main.py --dataset ATGO --output_dir ./atgo/baselines_unconstrained \
    # --alignment_dir ./2024_01/diamond_swissprot_2024_01_alignment.tsv --k_values 1 3 5 10 15 20 \
    # --db_versions --aspects BPO CCO MFO
