#!/bin/bash

PRED_FOLDER="."

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

# Get all folders in the prediction folder
# that start with "baselines_" and are in the DB_VERSION list
# The folders are in the format "baselines_<DB_VERSION>_<ONTOLOGY>_<DATASET>"
ALL_FOLDERS=""
for DB_VERSION in "${DB_VERSION[@]}"; do
    FOLDERS=$(find "$PRED_FOLDER/$DB_VERSION" -maxdepth 1 -type d -name "baselines_*" 2>/dev/null)
    ALL_FOLDERS="$ALL_FOLDERS"$'\n'"$FOLDERS"
done
# ---

FOLDERS=$(echo "$ALL_FOLDERS" | grep -v '^$' | tac)
echo "Found the following folders:"
echo "$FOLDERS"

# FOLDERS=$(echo "$FOLDERS" | tac)
for PRED_FOLDER in $FOLDERS; do
    echo "Processing folder: $PRED_FOLDER"        

    # Get all prediction files tsv in the run directory.
    # Remove the method name from grep if you want to skip its evaluation
    PRED_FILE=$(find "$PRED_FOLDER" -type f -name "*.tsv" | grep -E "AlignmentScore|BlastKNN|NaiveBaseline")

    for PRED_FILE in $PRED_FILE; do
        echo "Processing prediction file: $PRED_FILE"

        # Check if eval_beprof_evaluation_results_detailed.pkl already exists, and skip if it does
        if [[ -f "${PRED_FILE%.tsv}_eval_beprof_evaluation_results_detailed.pkl" ]]; then
            echo "Evaluation results already exist for $PRED_FILE. Skipping..."
            continue
        fi

        # Get the ontology name from the prediction file name
        if [[ "$PRED_FOLDER" == *"CCO"* ]]; then
            ONTOLOGY="CCO"
        elif [[ "$PRED_FOLDER" == *"BPO"* ]]; then
            ONTOLOGY="BPO"
        elif [[ "$PRED_FOLDER" == *"MFO"* ]]; then
            ONTOLOGY="MFO"
        fi

        # Get the dataset name from the prediction file name
        if [[ "$PRED_FOLDER" == *"D1"* ]]; then
            echo "Using D1 background and ground truth files"
            DATASET="D1"
        elif [[ "$PRED_FOLDER" == *"H30"* ]]; then
            echo "Using H30 background and ground truth files"
            DATASET="H30"
        fi
        
        BACKGROUND="./background/background_${DATASET}_2024_01.pkl"
        TRUE_FILE="./${DATASET}_test_annotations/${DATASET}_${ONTOLOGY}_test.pkl"
        PRED_OUT="${PRED_FILE%.tsv}.pkl"
        OUTPUT_PATH=$(dirname "$PRED_FILE")
        # if 'k*_' is in the path, add the k* to the output path
        if [[ "$PRED_FILE" == *"k"* ]]; then
            OUTPUT_PATH="${OUTPUT_PATH}/$(basename "$PRED_FILE" | cut -d'_' -f1)"
        fi

        echo "---"
        echo "Run directory: $PRED_FOLDER"
        echo "Ontology: $ONTOLOGY"
        echo "Ground Truth file: $TRUE_FILE"
        echo "Background: $BACKGROUND"
        echo "Predictions file: $PRED_FILE"
        echo "Predictions output: $PRED_OUT"
        echo "Output path: $OUTPUT_PATH"
        echo "---"

        # Convert predictions to .pkl format on the fly
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

        # Remove created pkl files once done
        rm "$PRED_OUT"

        echo "-------------------------------------------------"
    done
    echo "Finished processing folder: $PRED_FOLDER"
    echo "-------------------------------------------------"
done
echo "All folders processed."