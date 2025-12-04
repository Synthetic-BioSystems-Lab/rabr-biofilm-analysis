import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Load the data, label pipelines, and combine
df_18s = pd.read_csv("18s_data.csv")

# Filter Data (RABR and Phylum) and exclude same RABR at different time
df_18s = df_18s[
    (df_18s["section"] == "81RABR")
    & (df_18s["level"] == "Phylum")
    & (~df_18s["sample_id"].isin(["R45_11_15_21", "16_R45_11_15_21_16S"]))
].copy()

print(df_18s)

# Write csv
df_18s.to_csv("18s_data_81RABR.csv", index=False)

# Filter temperature and exclude abundunce less than 1%
df_10 = df_18s[df_18s["temp[C]"] == 10.0].copy()
df_10 = df_10[df_10["rel_abund"] >= 1].copy()

df_20 = df_18s[df_18s["temp[C]"] == 20.0].copy()
df_20 = df_20[df_20["rel_abund"] >= 1].copy()

df_30 = df_18s[df_18s["temp[C]"] == 30.0].copy()
df_30 = df_30[df_30["rel_abund"] >= 1].copy()

# Select temp (10, 20, 30)
stacked_df = df_30

# Pivot table: sample_id vs taxon
pivot = stacked_df.pivot_table(
    index="sample_id", columns="taxon", values="rel_abund", aggfunc="sum", fill_value=0
)

# Reorder taxa by total abundance
taxa_order = pivot.sum().sort_values(ascending=False).index.tolist()
pivot = pivot[taxa_order]

label_order = ["R45_10_18_21", "R36_10_28_21", "R58_10_28_21", "R46_11_05_21"]

# labels = label_order

# Plot setup
colors = (
    sns.color_palette("Set1", 9)
    + sns.color_palette("Set2", 8)
    + sns.color_palette("Set3", 10)
)

colors = colors[:27]  # Trim if needed

manual_colors = {
    "Chlorophyta": colors[0],
    "Alveolata": colors[1],
    "Stramenopiles": colors[2],
    "Opisthokonta": colors[3],
    "Rhizaria": colors[4],
    "Proteobacteria": colors[5],
    "Bacteroidetes": colors[6],
    "Firmicutes": colors[7],
    "Actinobacteria": colors[8],
    "Cyanobacteria": colors[9],
    "Chloroflexi": colors[10],
    "Unclassified Eukaryota": colors[11],
    "Epsilonbacteraeota": colors[12],
    "Nitrospirae": colors[13],
    "Thaumarchaeota": colors[14],
    "Euryarchaeota": colors[15],
    "Euglenozoa": colors[16],
    "Acidobacteria": colors[17],
    "Ochrophyta": colors[18],
    "Unclassified": colors[19],
    "Rhodophyta": colors[20],
    "Unclassified Plantae": colors[21],
    "Discosea": colors[22],
    "Streptophyta": colors[23],
    "Gemmatimonadetes": colors[24],
    "Cloacimonetes": colors[25],
}

taxon_colors = {taxon: manual_colors.get(taxon, "#cccccc") for taxon in taxa_order}

custom_order = [
    "Acidobacteria",
    "Actinobacteria",
    "Alveolata",
    "Bacteroidetes",
    "Chloroflexi",
    "Chlorophyta",
    "Cloacimonetes",
    "Cyanobacteria",
    "Discosea",
    "Epsilonbacteraeota",
    "Euglenozoa",
    "Euryarchaeota",
    "Firmicutes",
    "Gemmatimonadetes",
    "Nitrospirae",
    "Ochrophyta",
    "Opisthokonta",
    "Proteobacteria",
    "Rhizaria",
    "Rhodophyta",
    "Streptophyta",
    "Stramenopiles",
    "Thaumarchaeota",
    "Unclassified",
    "Unclassified Eukaryota",
    "Unclassified Plantae",
]

# Keep only taxa in the custom order
present_taxa = [taxon for taxon in custom_order if taxon in pivot.columns]
pivot = pivot[present_taxa]

# Plot single stacked bars per sample_id
fig, ax = plt.subplots(figsize=(12, 6))
bottom = np.zeros(len(pivot))

for taxon in present_taxa:
    heights = pivot[taxon].values
    ax.bar(pivot.index, heights, bottom=bottom, label=taxon, color=taxon_colors[taxon])
    bottom += heights

# Formatting
ax.set_title("Relative Abundance of 18S Eukaryotic Taxa by Sample (25°C)", fontsize=14)
ax.set_ylabel("Relative Abundance (%)", fontsize=16)
ax.set_xticks(np.arange(len(pivot.index)))
ax.set_xticklabels(pivot.index, rotation=45, ha="right", fontsize=12)
ax.set_xlabel("")
ax.set_ylim(0, 100)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(title="Taxon", bbox_to_anchor=(1.05, 1), fontsize=10, title_fontsize=11)
plt.tight_layout()
plt.savefig("18S_Lab_25.pdf", format="pdf", bbox_inches="tight")
plt.show()
