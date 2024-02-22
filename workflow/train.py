import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier

data_file = snakemake.input["train_file"]
model_file = snakemake.output["model_file"]
model_name = snakemake.params["model_name"]

# Load the data
data_table = pd.read_csv(data_file)
X = data_table.drop(columns=["Survived"])
y = data_table["Survived"]

# Model
if model_name == "Logistic":
    clf = LogisticRegression()
elif model_name == "RandomForest":
    clf = RandomForestClassifier()
elif model_name == "SVC":
    clf = SVC(probability=True)
elif model_name == "KNeighborsClassifier":
    clf = KNeighborsClassifier()
else:
    raise ValueError(f"Unknown model name: {model_name}")

# Train a logistic regression model
estimators = [("reduce_dim", StandardScaler()), ("clf", clf)]
model = Pipeline(estimators)
model.fit(X, y)

# Save the model
import pickle

with open(model_file, "wb") as f:
    pickle.dump(model, f)
