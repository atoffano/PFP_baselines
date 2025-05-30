import pandas as pd
import ast
import os
import pickle
import argparse
from collections import defaultdict


def process_file(input_file):
    # Extract the ontology type (BPO, CCO, MFO) from the filename
    base_name = os.path.basename(input_file)
    ontology_type = base_name.split("_")[1]  # Extract BPO, CCO, or MFO
    ontology_mapping = {"BPO": "all_bp", "CCO": "all_cc", "MFO": "all_mf"}
    ontology_key = ontology_mapping[ontology_type]

    # Output filename: same as input but with .pkl extension
    output_file = os.path.splitext(input_file)[0] + ".pkl"

    print(f"Processing {base_name}...")

    # Read the input file
    df = pd.read_csv(input_file, sep="\t")
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
    with open(output_file, "wb") as f:
        pickle.dump(protein_go_terms_dict, f)
    print(f"Saved pickle file: {output_file}")

    # Print sample data for verification
    if protein_go_terms_dict:
        sample_protein = next(iter(protein_go_terms_dict))
        print(f"Sample data for {sample_protein} in {ontology_type}:")
        print(protein_go_terms_dict[sample_protein])


def main():
    parser = argparse.ArgumentParser(description="Convert GT format TSV to pickle.")
    parser.add_argument("--input", nargs="+", help="Input TSV files")
    args = parser.parse_args()

    for input_file in args.input:
        process_file(input_file)

    print("Done!")


if __name__ == "__main__":
    main()
