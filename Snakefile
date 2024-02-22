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
        RAW_DATA_FILE
    output:
        PREPROCESSED_DATA_FILE
    shell:
        "python3 workflow/preprocessing.py {input} {output}"

N_FOLDS = 5
TRAIN_FILE = j(PREP_DATA_DIR, "train_data~titanic_fold~{fold}.csv")
TEST_FILE = j(PREP_DATA_DIR, "test_data~titanic_fold~{fold}.csv")

rule train_test_data_split:
    input:
        input_files = PREPROCESSED_DATA_FILE
    output:
        train_files = expand(TRAIN_FILE, fold=range(N_FOLDS))
        test_files = expand(TEST_FILE, fold=range(N_FOLDS))
    params:
        n_folds = N_FOLDS
    script:
        "workflow/train_test_data_split.py"


# =============================
# Train
# =============================

MODEL_LIST = ["Logistic", "RandomForest", "SVC", "KNeighborsClassifier"]

DERIVED_DIR = j(DATA_DIR, "derived")
MODEL_FILE = j(DERIVED_DIR, "model~{model}_data~titanic_fold~{fold}.pkl")

rule train:
    input:
        train_file = TRAIN_FILE
    output:
        model_file = MODEL_FILE
    params:
        model_name = lambda wildcards: wildcards.model
    script:
        "workflow/train.py"


# =============================
# Test & Evaluation
# =============================

TEST_PRED_FILE = j(DERIVED_DIR, "test_prediction_model~{model}_data~titanic_fold~{fold}.csv")
EVALUATION_RESULTS_FILE = j(DERIVED_DIR, "results_model~{model}_data~titanic_fold~{fold}.csv")

rule test:
    input:
        test_file = TEST_FILR
        model_file = MODEL_FILE
    output:
        output_file = TEST_PRED_FILE
    script:
        "workflow/test.py"

rule evaluation:
    input:
        test_file = TEST_FILE,
        test_pred_file = TEST_PRED_FILE
    output:
        output_file = EVALUATION_RESULTS_FILE
    params:
        model_name = lambda wildcards.model
    script:
        "workflow/evaluation.py"

# =============================
# Plot
# =============================

FIG_FILE = j(FIG_DIR, "aucroc_curve~data~titanic.png")
rule plot:
    input:
        input_files = expand(EVALUATION_RESULTS_FILE, model=MODEL_LIST, fold=range(N_FOLDS))
    output:
        output_file = FIG_FILE
    script:
        "workflow/plot.py"

# =============================
# All
# =============================
rule all:
    input:
        FIG_FILE
