import numpy as np
import pandas as pd
import pickle

data_file = snakemake.input["test_file"]
model_file = snakemake.input["model_file"]
output_file = snakemake.output["output_file"]

# Load the data
data_table = pd.read_csv(data_file)
X_test = data_table.drop(columns=["Survived"])

# Load the model
with open(model_file, "rb") as f:
    model = pickle.load(f)

# Predict using the model
predictions = model.predict_proba(X_test)[:, 1]

# Save the predictions
pd.DataFrame(predictions).to_csv(output_file, index=False, header=False)
