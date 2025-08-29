# PFP_baselines

Protein Function Prediction (PFP) is a critical task in bioinformatics, where the goal is to predict the function of proteins based on their sequences and annotations.

Current methods often rely on machine learning and deep learning techniques, with alignement-based methods being a common baseline.
This repository shows that their implementation in the literature is often sub-optimal and does not represent real-world scenarios.
It provides a pipeline to run these methods in a more realistic setting, showing that they can outperform state-of-the-art by a large margin.


## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Step-by-Step Instructions](#step-by-step-instructions)
  - [0. Diamond Installation](#0-diamond-installation)
  - [1. Data Downloading & Preparation](#1-data-downloading--preparation)
  - [2. Annotation Propagation](#2-annotation-propagation)
  - [3. Alignment](#3-alignment)
  - [4. Preparing the evaluation](#4-preparing-the-evaluation)
  - [5. Running Baselines](#5-running-baselines)
  - [6. Evaluation](#6-evaluation)
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

### 2. Annotation Propagation
Propagate GO annotations using the ontology structure:

```sh
python propagate_annotations.py
```
This step uses the GO ontology to propagate annotations from parent to child terms, ensuring that all relevant annotations are included.
Note that propagating terms can take a while.

### 3. Alignment
Run sequence alignment using Diamond on the most up-to-date SwissProt database:

```sh
echo "Creating Diamond database..."
diamond makedb --in data/2024_01/swissprot_2024_01.fasta -d data/2024_01/swissprot_2024_01_proteins_set

echo "Running Diamond blast on protein sequences against themselves..."
diamond blastp --very-sensitive --db data/swissprot/2024_01/swissprot_2024_01_proteins_set.dmnd --query data/swissprot/2024_01/swissprot_2024_01.fasta --out data/swissprot/2024_01/diamond_swissprot_2024_01_alignment.tsv -e 0.001
```
This step creates a Diamond database from the SwissProt protein sequences and performs a sequence alignment to find similar proteins. The output will be stored in `data/swissprot/2024_01/diamond_swissprot_2024_01_alignment.tsv`.
Note that as the 2024 release of SwissProt contains >500,000 proteins, the all-vs-all alignment step can be long (about 1 hour).

### 4. Preparing the evaluation
To evaluate the performance of a method, IC-weighted scores are used. These scores are computed based on the Information Content (IC) of the GO terms, which is derived from the background distribution of GO terms in the dataset.
Background files needs to be generated prior to running the baselines.
Example for the ATGO dataset:
```sh
python background.py --cco ./data/ATGO/ATGO_CCO_train_annotations.tsv --bpo ./data/ATGO/ATGO_BPO_train_annotations.tsv --mfo ./data/ATGO/ATGO_MFO_train_annotations.tsv --output ./data/ATGO/background_ATGO.pkl --test_cco ./data/ATGO/ATGO_MFO_test_annotations.tsv --test_bpo ./data/ATGO/ATGO_BPO_test_annotations.tsv --test_mfo ./data/ATGO/ATGO_CCO_test_annotations.tsv
```

### 5. Running Baselines
Run the baseline methods (e.g., Naive, DiamondKNN, AlignmentScore) using the prepared data and alignment results. The following command runs the baselines on the ATGO dataset, under the constrained setup:
```sh
python main.py \
--db_version '' \
--dataset ATGO \
--alignment_dir ./data/swissprot/2024_01/diamond_swissprot_2024_01_alignment.tsv \
--k_values 1 3 5 10 15 20 \
--aspects BPO CCO MFO \

```
Dataset can be set to `D1` (BeProf D1 dataset), `H30` (Low homology dataset), `ATGO` or `CAFA3`.
To greatly speed up the process, you can skip the Naive baseline by uncommenting the corresponding line in the `main.py` script.
`--db_version` can be set to the SwissProt version you want to use, e.g. `2024_01`, or a collection of versions, e.g. `2024_01 2021_01`. If not set, the script will use all available versions.
To apply annotation transfer from the train set to the test set (i.e. the usual setup in the litterature), set `--db_version` to an empty string `''`.
`--alignment_dir` specifies the path to the Diamond alignment file generated in step 3.
`--k_values` specifies the k values to use for the KNN baseline. You can adjust these values based on your needs.
`--aspects` specifies the GO subontologies to consider (BPO, CCO, MFO). Defaults to all three aspects.
Additional arguments can be passed to the script, such as `--experimental_only` to run only using experimental annotations. Leaving this flag unset will include all manually curated GO annotations present in SwissProt.
`--one_vs_all` will run the baselines in a one-vs-all setup, where each test protein can receive annotations from the rest of the proteins in the dataset (excluding themselves), regardless of their train/test split.
`--annotations_2024_01` will transfer annotations from the 2024_01 SwissProt release, using only proteins present in the specified `--db_version`.

Example usage on the ATGO dataset, applied to all SwissProt version, using experimental annotations and a one-vs-all setup:

```python
python main.py \
--dataset ATGO \
--alignment_dir ./data/swissprot/2024_01/diamond_swissprot_2024_01_alignment.tsv \
--k_values 1 3 5 10 15 20 \
--aspects BPO CCO MFO \
--experimental_only \
--one_vs_all
```

### 6. Evaluation
Evaluating predictions is ran automatically when running the `main.py` script.
However, this can also be done manually by running the `evaluation.py` script.
The script evaluates all predictions generated by the baselines and outputs the results in the same directory as the prediction file, automatically selecting the correct backgrounds and ground truth files based on folder's name.
Example usage:
```sh
python evaluation.py \
--input_dir './results/ATGO/baselines_ATGO_2024_01_BPO_exp_one_vs_all' \
--aspect 'BPO CCO MFO' \
--k_values 1 3 5 10 15 20
```

See the  [BeProf evaluation script github](https://github.com/CSUBioGroup/BeProf/tree/main) for details on the parameters, options and file formats.


---

### Notes
Adjust file paths and parameters as needed for your specific setup.
Additional performance could likely be further squeezed by tuning the scoring functions and parameters of the baseline according to the following paper: [A large-scale assessment of sequence database search tools for homology-based protein function prediction](https://doi.org/10.1093/bib/bbae349).