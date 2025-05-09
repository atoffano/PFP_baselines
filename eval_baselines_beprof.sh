#!/bin/bash
# filepath: /home/atoffano/these-antoine/results_beprof/beprof_eval.sh

# Usage:
#   ./beprof_eval.sh --pred_folder /path/to/base_folder
echo "Running evaluation script using BeProf script..."
# Parse input arguments to get --pred_folder value
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --pred_folder)
            PRED_FOLDER="$2"
            shift 2;;
        --pred_folder=*)
            PRED_FOLDER="${1#*=}"
            shift 1;;
        *)
            echo "Unknown parameter passed: $1"
            exit 1;;
    esac
done

echo "Pred folder: $PRED_FOLDER"

if [ -z "$PRED_FOLDER" ]; then
    echo "Usage: $0 --pred_folder /path/to/base_folder"
    exit 1
fi

## --- DB version study folders ---
DB_VERSION=(
    "2024_01"
    "2023_01"
    "2022_01"
    "2021_01"
    "2020_01"
    "2019_01"
    "2018_01"
    "2017_01"
    "2016_01"
    "2015_01"
    "2014_01"
    "2013_01"
    "2012_01"
    "2011_01"
    "2010_01"
    "2009_03"
    "2008_01"
    "2007_03"
    "2006_02"
    "2005_01"
    "2003_12"
)


ALL_FOLDERS=""
for DB_VERSION in "${DB_VERSION[@]}"; do
    FOLDERS=$(find "$PRED_FOLDER/$DB_VERSION" -maxdepth 1 -type d -name "baselines_H30*" 2>/dev/null)
    ALL_FOLDERS="$ALL_FOLDERS"$'\n'"$FOLDERS"
done
# ---

FOLDERS=$(echo "$ALL_FOLDERS" | grep -v '^$' | tac)
echo "Found the following folders:"
echo "$FOLDERS"

# # Alternative: limit to current level of folders in folder
# FOLDERS=$(find "$PRED_FOLDER" -maxdepth 1 -type d -name "baselines_H30_unconstrained_*")
# echo "Found the following folders:"
# echo "$FOLDERS"

# Revert folder order
FOLDERS=$(echo "$FOLDERS" | tac)
for PRED_FOLDER in $FOLDERS; do
    echo "Processing folder: $PRED_FOLDER"
    
        if [[ "$PRED_FOLDER" == *"D1"* ]]; then
            echo "Using D1 background"
            BACKGROUND="./background/background_D1_2024_01.pkl"
        fi
        if [[ "$PRED_FOLDER" == *"H30"* ]]; then
            BACKGROUND="./background/background_H30_2024_01.pkl"
            # Determine ontology based on the PRED_FOLDER path
            if [[ "$PRED_FOLDER" == *"CCO"* ]]; then
                ONTOLOGY="CCO"
                TRUE_FILE="./gt_H30_CCO_test.pkl"
                #TRUE_FILE="/home/atoffano/these-antoine/data/BeProf/Dataset1/prop_nofilt/gt_CCO_test.pkl"
            elif [[ "$PRED_FOLDER" == *"BPO"* ]]; then
                ONTOLOGY="BPO"
                TRUE_FILE="/home/atoffano/these-antoine/data/uniprotkb/gt_H30_BPO_test.pkl"
                #TRUE_FILE="/home/atoffano/these-antoine/data/BeProf/Dataset1/prop_nofilt/gt_BPO_test.pkl"
            elif [[ "$PRED_FOLDER" == *"MFO"* ]]; then
                ONTOLOGY="MFO"
                TRUE_FILE="/home/atoffano/these-antoine/data/uniprotkb/gt_H30_MFO_test.pkl"
                #TRUE_FILE="/home/atoffano/these-antoine/data/BeProf/Dataset1/prop_nofilt/gt_MFO_test.pkl"
            else
                echo "Unable to determine ontology (CCO/BPO/MFO) from the pred_file path: $PRED_FILE"
                continue
            fi
        fi


        echo "Ontology: $ONTOLOGY"

        # Get all prediction files tsv in the run directory. They are in subfolders AlignmentScore, BlastKNN and NeiveBaseline.
        PRED_FILE=$(find "$PRED_FOLDER" -type f -name "*.tsv" | grep -E "AlignmentScore|BlastKNN")

        # For all the prediction files, run the evaluation script
        for PRED_FILE in $PRED_FILE; do
            echo "Processing prediction file: $PRED_FILE"

            # Prediction output same as pred file but with .pkl extension
            PRED_OUT="${PRED_FILE%.tsv}.pkl"
            # Output path for the evaluation results = same folder as where the predictions are
            OUTPUT_PATH=$(dirname "$PRED_FILE")
            # if 'k*_' is in the path, add the k* to the output path
            if [[ "$PRED_FILE" == *"k"* ]]; then
                OUTPUT_PATH="${OUTPUT_PATH}/$(basename "$PRED_FILE" | cut -d'_' -f1)"
            fi

            echo "Run directory: $PRED_FOLDER"
            echo "Ontology: $ONTOLOGY"
            echo "True file: $TRUE_FILE"
            echo "Predictions file: $PRED_FILE"
            echo "Predictions output: $PRED_OUT"
            echo "Output path: $OUTPUT_PATH"
            echo "Background: $BACKGROUND"

            # Run the conversion script to create the .pkl prediction file
            echo "Converting predictions to .pkl format..."
            python convert_format_beprof.py --pred_file "$PRED_FILE" --pred_out "$PRED_OUT"

            # Run the evaluation script using the converted prediction file and true file
            echo "Running evaluation..."
            python beprof_eval.py --predict "$PRED_OUT" \
                                --output_path "$OUTPUT_PATH" \
                                --true "$TRUE_FILE" \
                                --background "$BACKGROUND" \
                                --go "/home/atoffano/these-antoine/data/ontologies/go.obo" \
                                --metrics "0,1,2,3,4"

            echo "-------------------------------------------------"
        done
done
