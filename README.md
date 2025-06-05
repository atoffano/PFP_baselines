# PFP_baselines

Protein Function Prediction (PFP) is a critical task in bioinformatics, where the goal is to predict the function of proteins based on their sequences and annotations.
Current methods often rely on machine learning and deep learning techniques, with alignement-based methods being a common baseline.
This repository shows that their implementation in the litterature is often sub-optimal and does not represent real-world scenarios. This repository provides a pipeline to run these methods in a more realistic setting, showing that they can outperform state-of-the-art by a large margin.


## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Step-by-Step Instructions](#step-by-step-instructions)
  - [1. Data Downloading](#1-data-downloading)
  - [2. Data Parsing](#2-data-parsing)
  - [3. Propagation](#3-propagation)
  - [4. Alignment](#4-alignment)
  - [5. Running Baselines](#5-running-baselines)
  - [6. Format Conversion](#6-format-conversion)
  - [7. Evaluation](#7-evaluation)

---

## Overview

This repository provides a pipeline for running baseline alignment-based methods (e.g. Diamond) for Protein Function Prediction (PFP) tasks. It covers all steps from data downloading to evaluation, including parsing, annotation propagation, alignment, running baselines, format conversion, and final evaluation.

## Requirements

- Python 3.x
- Required Python packages (install with `pip install -r requirements.txt`)

## Step-by-Step Instructions

### 0. Diamond Installation
Make sure you have the Diamond software installed.
You can download it from the [Diamond GitHub repository](http://github.com/bbuchfink/diamond) or simply execute the following command on Linux-based systems:

```sh
wget http://github.com/bbuchfink/diamond/releases/download/v2.1.11/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
```

### 1. Data Downloading & Preparation

Download all necessary information (e.g., protein sequences, GO annotations) from UniProt (SwissProt) at different releases and parse them into the required formats.

```sh
python download_swissprot.py
```

### 2. Propagation
Propagate GO annotations using the ontology structure:

```sh
python scripts/propagate_annotations.py --input data/parsed_annotations.tsv --ontology data/go.obo --output data/propagated_annotations.tsv
```
This step uses the GO ontology to propagate annotations from parent to child terms, ensuring that all relevant annotations are included.
Note that the propagation step can take a while.

### 3. Alignment
Run sequence alignment using Diamond on the most up-to-date SwissProt database:

```sh
echo "Creating Diamond database..."
diamond makedb --in 2024_01/swissprot_2024_01.fasta -d 2024_01/swissprot_2024_01_proteins_set

echo "Running Diamond blast on protein sequences against themselves..."
diamond blastp --very-sensitive --db 2024_01/swissprot_2024_01_proteins_set.dmnd --query 2024_01/swissprot_2024_01.fasta --out 2024_01/diamond_swissprot_2024_01_alignment.tsv -e 0.001
```
This step creates a Diamond database from the SwissProt protein sequences and performs a sequence alignment to find similar proteins. The output will be stored in `2024_01/diamond_swissprot_2024_01_alignment.tsv`.
Note that as the 2024 release of SwissProt contains >500,000 proteins, the alignment step can take a while (around 1 hour).

### 4. Running Baselines
```sh
python baselines.py --dataset D1
```
Run the baseline methods (e.g., Naive, Diamond-KNN, AlignmentScore). If running on all Swissprot versions, this can take a while. Dataset can be set to `D1` (BeProf D1 dataset) or `H30` (Low homology dataset).
To greatly speed up the process, you can skip the Naive baseline by uncommenting the corresponding line in the `baselines.py` script.


### 5. Preparing the evaluation
Evaluate the baseline predictions using the BeProf evaluation script.
This step first requires the conversion of the ground truth targets (actual test annotations) to the format required by the BeProf evaluation script and the generation of the background file (which allows to compute weighted scores).
Those can be done with the following commands:

#### Convert the test annotations to the .pkl format required by the BeProf evaluation script:
```sh
# BeProf D1 dataset
python convert_test_annot.py --input ./D1_test_annotations/D1_BPO_test.tsv ./D1_test_annotations/D1_CCO_test.tsv ./D1_test_annotations/D1_MFO_test.tsv

# Low homology dataset
python convert_test_annot.py --input ./H30_test_annotations/H30_BPO_test.tsv ./H30_test_annotations/H30_CCO_test.tsv ./H30_test_annotations/H30_MFO_test.tsv
```
Files will be stored in the same directory as the input files.

#### Generate the background file
The background file is used to compute the term's Information Content (IC) to obtain weighted scores. It contains the distribution of GO terms in the dataset used to transfer annotations.
```sh
# BeProf D1 dataset
python background.py --cco ./BeProf_D1/D1_CCO.tsv --bpo ./BeProf_D1/D1_BPO.tsv --mfo ./BeProf_D1/D1_MFO.tsv --output ./background/background_D1.pkl --test_cco ./D1_test_annotations/D1_BPO_test.tsv --test_bpo ./D1_test_annotations/D1_BPO_test.tsv --test_mfo ./D1_test_annotations/D1_MFO_test.tsv

# Low homology dataset
python background.py --cco ./2024_01/swissprot_2024_01_CCO_annotations.tsv --bpo ./2024_01/swissprot_2024_01_BPO_annotations.tsv --mfo ./2024_01/swissprot_2024_01_MFO_annotations.tsv --test_cco ./H30_test_annotations/H30_BPO_test.tsv --test_bpo ./H30_test_annotations/H30_BPO_test.tsv --test_mfo ./H30_test_annotations/H30_MFO_test.tsv --output ./background/background_2024_01.pkl
```

### 6. Evaluation
Evaluating predictions can be done simply by running the `eval_baselines.sh` script. The script evaulates all predictions generated by the baselines and outputs the results in the same directory as the prediction file, automatically selecting the correct backgrounds and ground truth files based on folder's name.
For fine-grained control over the evaluation, you can run the following commands:

```sh 
# Convert the predictions to the .pkl format required by the BeProf evaluation script
python convert_format_beprof.py --pred_file "$PRED_FILE" --pred_out 'path/to/predictions.pkl'
# Run the actual evaluation
python beprof_eval.py --predict 'path/to/predictions.pkl' \
                      --output_path 'output/path' \
                      --true 'path/to/ground_truth_annotations.tsv' \
                      --background "$BACKGROUND" \
                      --go "go.obo" \
                      --metrics "0,1,2,3,4"
```
See the  [BeProf evaluation script github](https://github.com/CSUBioGroup/BeProf/tree/main) for details on the parameters, options and file formats.

---

### Notes
Adjust file paths and parameters as needed for your specific setup.
Additional performance could likely be further squeezed by augmenting the Diamond sensitivity to `--ultra-sensitive` and by tuning the scoring functions and parameters of the baseline according to the following paper: [A large-scale assessment of sequence database search tools for homology-based protein function prediction](https://doi.org/10.1093/bib/bbae349).