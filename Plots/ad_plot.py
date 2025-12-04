import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Load data
df = pd.read_csv("alpha_diversity.csv")

# Filter for 16S only
df = df[df["label"] == "23s"].copy()

print(df)


# Classify section with temperature for 81RABR
def classify(row):
    if row["section"] == "81RABR":
        if row["temp[C]"] == 10.0:
            return "81RABR_10"
        elif row["temp[C]"] == 20.0:
            return "81RABR_20"
        elif row["temp[C]"] == 30.0:
            return "81RABR_30"
        else:
            return "81RABR_unknown"
    return row["section"]


df["group"] = df.apply(classify, axis=1)

group_means = df.groupby("group")["invsimpson"].mean().sort_index()
print("Mean Inverse Simpson Index per Group:")
print(group_means)

# Set order and rename
group_order = ["control", "81RABR_10", "81RABR_20", "81RABR_30", "pilot", "TF"]
group_labels = {
    "control": "Control",
    "pilot": "Pilot",
    "TF": "Trickling Filter",
    "81RABR_10": "Lab-Scale (10°C)",
    "81RABR_20": "Lab-Scale (20°C)",
    "81RABR_30": "Lab-Scale (25°C)",
}

# Color palette (Set1)
palette = dict(zip(group_order, sns.color_palette("Set1", len(group_order))))

# Plot
plt.figure(figsize=(10, 6))
sns.boxplot(data=df, x="group", y="invsimpson", order=group_order, palette=palette)
sns.stripplot(
    data=df,
    x="group",
    y="invsimpson",
    order=group_order,
    color="black",
    size=5,
    jitter=True,
    alpha=0.7,
)

# Customize
plt.xticks(rotation=45, ha="right", fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel("")
plt.ylabel("Inverse Simpson Index", fontsize=14)
# plt.title("Alpha Diversity of Bacterial Communities (16S)", fontsize=16)
# plt.title("Alpha Diversity of Bacterial Communities (16S DADA2)", fontsize=16)
# plt.title("Alpha Diversity of Eukaryotic Communities (18S)", fontsize=16)
plt.title("Alpha Diversity of Phototrophic Communities (23S)", fontsize=16)
plt.xticks(ticks=range(len(group_order)), labels=[group_labels[g] for g in group_order])
sns.despine()
plt.tight_layout()
# plt.savefig("23S_AD.pdf", format="pdf", bbox_inches="tight")
plt.show()
