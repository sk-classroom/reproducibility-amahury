from sklearn.metrics import roc_auc_score
import pandas as pd
import pickle

test_file = snakemake.input["test_file"]
test_pred_file = snakemake.input["test_pred_file"]
model_name = snakemake.params["model_name"]
output_file = snakemake.output["output_file"]

# Load the data
data_table = pd.read_csv(test_file)
y_test = data_table["Survived"]

y_pred = pd.read_csv(test_pred_file, header=None).values


# Compute ROC AUC score
aucroc_score = roc_auc_score(y_test, y_pred)

# Save the AUC ROC score
df = pd.DataFrame([aucroc_score], columns=["AUCROC"])
df["model"] = model_name
df.to_csv(output_file, index=False)
