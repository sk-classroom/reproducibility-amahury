import numpy as np
from os.path import join as j

DATA_DIR = "data"
FIG_DIR = "figs"

# =============================
# Load the data
# =============================

# We will download the data from https://raw.githubusercontent.com/datasciencedojo/datasets/master/titanic.csv and save it in the data/raw directory

RAW_DATA_DIR = j(DATA_DIR, "raw")
RAW_DATA_FILE = j(RAW_DATA_DIR, "data~titanic.csv")

rule download_data:
    output:
        RAW_DATA_FILE
    shell: "curl -o {output}  https://raw.githubusercontent.com/datasciencedojo/datasets/master/titanic.csv"


# =============================
# Preprocess data
# =============================

# 1. We will preprocess the data (imputing missing values and encoding categorical variables into one-hot vector representations, etc)
# 2. We will then split the data into K folds for cross-validation

PREP_DATA_DIR = j(DATA_DIR, "preprocessed")
PROCESSED_DATA_FILE = j(PREP_DATA_DIR, "preprocessed_data~titanic.csv")

rule preprocess_data:
    input:
        ...
    output:
        ...
    shell:
        ...

N_FOLDS = 5
TRAIN_FILE = j(PREP_DATA_DIR, "train_data~titanic_fold~{fold}.csv")
TEST_FILE = j(PREP_DATA_DIR, "test_data~titanic_fold~{fold}.csv")

rule train_test_data_split:
    input:
        input_file = ...
    output:
        train_files = ...
        test_files = ...
    params:
        n_folds = N_FOLDS
    script:
        ...


# =============================
# Train
# =============================

MODEL_LIST = ["Logistic", "RandomForest", "SVC", "KNeighborsClassifier"]

DERIVED_DIR = j(DATA_DIR, "derived")
MODEL_FILE = j(DERIVED_DIR, "model~{model}_data~titanic_fold~{fold}.pkl")

rule train:
    input:
        train_file = ...
    output:
        model_file = ...
    params:
        model_name = ...
    script:
        ...


# =============================
# Test & Evaluation
# =============================

TEST_PRED_FILE = j(DERIVED_DIR, "test_prediction_model~{model}_data~titanic_fold~{fold}.csv")
EVALUATION_RESULTS_FILE = j(DERIVED_DIR, "results_model~{model}_data~titanic_fold~{fold}.csv")

rule test:
    input:
        test_file = ...,
        model_file = ...
    output:
        output_file = ...
    script:
        ...

rule evaluation:
    input:
        ...
    output:
        ...
    params:
        ...
    script:
        ...

# =============================
# Plot
# =============================

FIG_FILE = j(FIG_DIR, "aucroc_curve~data~titanic.png")
rule plot:
    input:
        ...
    output:
        output_file = FIG_FILE
    script:
        ...

# =============================
# All
# =============================
rule all:
    input:
        FIG_FILE
