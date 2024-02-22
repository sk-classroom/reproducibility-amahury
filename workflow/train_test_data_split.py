import numpy as np
import pandas as pd
from sklearn.model_selection import KFold

input_files = snakemake.input["input_file"]
test_files = snakemake.output["test_files"]
train_files = snakemake.output["train_files"]
n_folds = snakemake.params["n_folds"]

data_table = pd.read_csv(input_files)

# Initialize KFold
kf = KFold(n_splits=n_folds, shuffle=True, random_state=42)

# Splitting the data into K folds
for fold, (train_idx, test_idx) in enumerate(kf.split(data_table)):
    train_data = data_table.iloc[train_idx]
    test_data = data_table.iloc[test_idx]

    # Save the kth fold
    train_data.to_csv(train_files[fold], index=False)
    test_data.to_csv(test_files[fold], index=False)
