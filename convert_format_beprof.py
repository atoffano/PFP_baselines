#!/usr/bin/env python3
import pandas as pd
import pickle
import argparse
import tqdm
import time


def convert_predictions(pred_file, subontology):
    """
    Converts a TSV prediction file with columns: target_ID, term_ID
    into a dictionary where each protein gets a 'bp' dictionary of GO term predictions.
    """
    df = pd.read_csv(pred_file, sep="\t")
    pred_dict = {}
    for prot, group in tqdm.tqdm(df.groupby("target_ID")):
        if subontology == "all":
            for sub in ["cc", "mf", "bp"]:
                pred_dict[prot] = {
                    f"{sub}": dict(zip(group["term_ID"], [1] * len(group["term_ID"])))
                }
        pred_dict[prot] = {
            f"{subontology}": dict(zip(group["term_ID"], group["score"]))
        }
    return pred_dict


def main():
    parser = argparse.ArgumentParser(
        description="Convert prediction and ground truth TSV files to .pkl format required for BeProf evaluation."
    )
    parser.add_argument(
        "--pred_file", required=True, help="Path to prediction TSV file"
    )
    parser.add_argument("--pred_out", required=True, help="Output pickle file path")
    args = parser.parse_args()

    # Get subontology from file name
    if "CCO" in args.pred_file:
        subontology = "cc"
    elif "BPO" in args.pred_file:
        subontology = "bp"
    elif "MFO" in args.pred_file:
        subontology = "mf"
    else:
        subontology = "all"
    # Convert files
    pred_dict = convert_predictions(args.pred_file, subontology)

    # Save dictionaries as pickle files
    with open(args.pred_out, "wb") as f:
        pickle.dump(pred_dict, f)

    print(f"Saved prediction file to: {args.pred_out}")


if __name__ == "__main__":
    main()

# python convert_format_beprof.py --pred_file /home/atoffano/these-antoine/results_beprof/run_20250217_121916_CCO_667649D2_F_NET/predictions/GNN/predictions_test.tsv --pred_out /home/atoffano/these-antoine/results_beprof/run_20250217_121916_CCO_667649D2_F_NET/cco_pred_test.pkl
