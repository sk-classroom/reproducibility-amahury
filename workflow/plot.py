import pandas as pd

# Load the data
input_files = snakemake.input["input_files"]
data_table = pd.concat([pd.read_csv(f) for f in input_files])
output_file = snakemake.output["output_file"]

# Plot the data
import matplotlib.pyplot as plt
import seaborn as sns

fig, ax = plt.subplots(figsize=(7, 5))

ax = sns.pointplot(data_table, x="model", y="AUCROC", hue="model", ax=ax)

ax.set_xlabel("Model")
ax.set_ylabel("AUC ROC score")
# ax.set_ylim(0.5, 1)

plt.savefig(output_file, dpi=300, bbox_inches="tight")
