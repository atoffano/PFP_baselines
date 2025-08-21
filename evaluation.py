import pandas as pd
import os
import sys
import subprocess
import pickle
from collections import defaultdict
import tqdm
import argparse
import logging


def setup_logging(output_dir, aspect):
    """Setup logging for evaluation."""
    log_dir = os.path.join(output_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    logger = logging.getLogger(f"{aspect}")
    logger.setLevel(logging.INFO)
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    log_file = os.path.join(log_dir, f"{aspect}_eval.log")
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


def evaluate(logger, output_dir, dataset, aspect, k_values):
    """
    Evaluate the predictions using the ground truth (GT) annotations and the BeProf evaluation method.
    """
    background_pkl = f"./data/{dataset}/background_{dataset}.pkl"
    go_obo_file = "./data/go.obo"

    # Check if GT exists in pkl format. If not, convert GT TSV to pkl using gt_convert
    gt_pkl = f"./data/{dataset}/{dataset}_{aspect}_test_annotations.pkl"
    if not os.path.exists(gt_pkl):
        gt_tsv = f"./data/{dataset}/{dataset}_{aspect}_test_annotations.tsv"
        if os.path.exists(gt_tsv):
            logger.info(f"Converting Ground Truth TSV {gt_tsv} to pkl format")
            gt_convert(gt_tsv)
        else:
            logger.error(f"Ground Truth TSV file {gt_tsv} does not exist.")
            raise FileNotFoundError(f"Ground Truth TSV file {gt_tsv} does not exist.")

    # Evaluate NaiveBaseline predictions
    logger.info(f"Evaluating NaiveBaseline predictions")
    # Convert predictions to pkl
    pred_file = f"{output_dir}/predictions/NaiveBaseline/predictions.tsv"
    pred_pkl = f"{output_dir}/predictions/NaiveBaseline/predictions.pkl"
    if os.path.exists(pred_file):
        pred_dict = convert_predictions(pred_file, aspect)
        with open(pred_pkl, "wb") as f:
            pickle.dump(pred_dict, f)

        run_beprof_evaluation(
            logger,
            pred_pkl,
            gt_pkl,
            background_pkl,
            go_obo_file,
            f"{output_dir}/evaluation/NaiveBaseline",
        )
    else:
        logger.warning(f"NaiveBaseline predictions file {pred_file} does not exist.")

    # Evaluate IDScore predictions
    logger.info(f"Evaluating IDScore predictions")
    # Convert predictions to pkl
    pred_file = f"{output_dir}/predictions/IDScore/predictions.tsv"
    pred_pkl = f"{output_dir}/predictions/IDScore/predictions.pkl"
    if os.path.exists(pred_file):
        pred_dict = convert_predictions(pred_file, aspect)
        with open(pred_pkl, "wb") as f:
            pickle.dump(pred_dict, f)

        run_beprof_evaluation(
            logger,
            pred_pkl,
            gt_pkl,
            background_pkl,
            go_obo_file,
            f"{output_dir}/evaluation/IDScore",
        )
    else:
        logger.warning(f"IDScore predictions file {pred_file} does not exist.")

    # Evaluate AlignmentScore predictions
    logger.info(f"Evaluating AlignmentScore predictions")
    # Convert predictions to pkl
    pred_file = f"{output_dir}/predictions/AlignmentScore/predictions.tsv"
    pred_pkl = f"{output_dir}/predictions/AlignmentScore/predictions.pkl"
    if os.path.exists(pred_file):
        pred_dict = convert_predictions(pred_file, aspect)
        with open(pred_pkl, "wb") as f:
            pickle.dump(pred_dict, f)
        run_beprof_evaluation(
            logger,
            pred_pkl,
            gt_pkl,
            background_pkl,
            go_obo_file,
            f"{output_dir}/evaluation/AlignmentScore",
        )
    else:
        logger.warning(f"AlignmentScore predictions file {pred_file} does not exist.")

    # # Evaluate BlastKNN predictions for each k value
    for k in k_values:
        logger.info(f"Evaluating BlastKNN predictions for k={k}")
        pred_file = f"{output_dir}/predictions/BlastKNN/k{k}_predictions.tsv"
        pred_pkl = f"{output_dir}/predictions/BlastKNN/k{k}_predictions.pkl"
        if os.path.exists(pred_file):
            pred_dict = convert_predictions(pred_file, aspect)
            with open(pred_pkl, "wb") as f:
                pickle.dump(pred_dict, f)

            run_beprof_evaluation(
                logger,
                pred_pkl,
                gt_pkl,
                background_pkl,
                go_obo_file,
                f"{output_dir}/evaluation/BlastKNN_k{k}",
            )
        else:
            logger.warning(
                f"BlastKNN predictions file {pred_file} for k={k} does not exist."
            )


def run_beprof_evaluation(
    logger, pred_pkl, gt_pkl, background_pkl, go_obo_file, eval_output_dir
):
    """
    Run beprof_eval.py as a subprocess.
    """
    # Create evaluation output directory
    os.makedirs(eval_output_dir, exist_ok=True)

    # Prepare command arguments
    cmd = [
        sys.executable,  # Use the same Python interpreter
        "beprof_eval.py",
        "--predict",
        pred_pkl,
        "--output_path",
        eval_output_dir,
        "--true",
        gt_pkl,
        "--background",
        background_pkl,
        "--go",
        go_obo_file,
        "--metrics",
        "0,1,2,3,4,5",  # All metrics: F_max, Smin, Aupr, ICAupr, DPAupr, threshold
    ]

    logger.info(f"Running BeProf evaluation: {' '.join(cmd)}")

    try:
        # Run the subprocess and wait for completion
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
            cwd=os.path.dirname(os.path.abspath(__file__)),  # Run from script directory
        )

        logger.info(f"BeProf evaluation completed successfully")
        logger.info(f"Results saved to: {eval_output_dir}")

        # Log stdout if there's any output
        if result.stdout:
            logger.info(f"BeProf stdout: {result.stdout}")

    except subprocess.CalledProcessError as e:
        logger.error(f"BeProf evaluation failed with return code {e.returncode}")
        if e.stdout:
            logger.error(f"BeProf stdout:\n{e.stdout}")
        if e.stderr:
            logger.error(f"BeProf stderr:\n{e.stderr}")
        raise
    except Exception as e:
        logger.error(f"Error running BeProf evaluation: {str(e)}")
        raise


def gt_convert(gt_tsv):
    """
    Convert GT to pkl.
    """
    # Extract the ontology type (BPO, CCO, MFO) from the filename
    base_name = os.path.basename(gt_tsv)
    ontology_type = base_name.split("_")[1]  # Extract BPO, CCO, or MFO
    ontology_mapping = {"BPO": "all_bp", "CCO": "all_cc", "MFO": "all_mf"}
    ontology_key = ontology_mapping[ontology_type]

    # Output filename: same as input but with .pkl extension
    gt_pkl = os.path.splitext(gt_tsv)[0] + ".pkl"

    print(f"Processing {base_name}...")

    # Read the input file
    df = pd.read_csv(gt_tsv, sep="\t")
    df["term"] = df["term"].str.split("; ")
    df = df.explode("term")

    # Build the nested dictionary for pickling
    protein_go_terms = defaultdict(
        lambda: {"all_bp": set(), "all_cc": set(), "all_mf": set()}
    )
    for _, row in df.iterrows():
        protein_id = row["EntryID"]
        term_id = row["term"]
        protein_go_terms[protein_id][ontology_key].add(term_id)

    # Convert defaultdict to regular dict for pickling
    protein_go_terms_dict = dict(protein_go_terms)

    # Save to pickle file
    with open(gt_pkl, "wb") as f:
        pickle.dump(protein_go_terms_dict, f)
    print(f"Saved pickle file: {gt_pkl}")

    # # Print sample data for verification
    # if protein_go_terms_dict:
    #     sample_protein = next(iter(protein_go_terms_dict))
    #     print(f"Sample data for {sample_protein} in {ontology_type}:")
    #     print(protein_go_terms_dict[sample_protein])


def convert_predictions(pred_file, aspect):
    """
    Converts a TSV prediction file with columns: target_ID, term_ID
    into a dictionary where each protein gets a 'bp' dictionary of GO term predictions.
    """
    subontology = aspect[:2].lower()
    df = pd.read_csv(pred_file, sep="\t")
    # Split term column by '; ' and explode
    df["term_ID"] = df["term_ID"].str.split("; ")
    df = df.explode("term_ID")
    pred_dict = {}
    for prot, group in tqdm.tqdm(df.groupby("target_ID")):
        # if subontology == "all":
        #     for sub in ["cc", "mf", "bp"]:
        #         pred_dict[prot] = {
        #             f"{sub}": dict(zip(group["term_ID"], [1] * len(group["term_ID"])))
        #         }
        pred_dict[prot] = {
            f"{subontology}": dict(zip(group["term_ID"], group["score"]))
        }
    return pred_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate predictions using BeProf.")
    parser.add_argument(
        "--input_dir",
        required=True,
        help="Directory containing predictions and evaluation results. Will also serve as output",
    )
    parser.add_argument("--dataset", required=True, help="Dataset name.")
    parser.add_argument(
        "--aspect", required=True, help="Ontology aspect (BPO, CCO, MFO)."
    )
    parser.add_argument(
        "--k_values",
        nargs="+",
        type=int,
        default=[1, 3, 5, 10, 15, 20],
        help="List of k values for BlastKNN.",
    )
    args = parser.parse_args()

    logger = setup_logging(args.input_dir, args.aspect)
    evaluate(logger, args.input_dir, args.dataset, args.aspect, args.k_values)
