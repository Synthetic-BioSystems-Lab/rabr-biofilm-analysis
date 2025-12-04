import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load the data, label pipelines, and combine
df_18s = pd.read_csv("18s_data.csv")

# Filter Data (RABR and Phylum) and exclude same RABR at different time
df_18s = df_18s[
    (df_18s["section"] == "control") &
    (df_18s["level"] == "Phylum")
].copy()

print(df_18s)

# Write csv
df_18s.to_csv("18s_data_CO.csv", index=False)

#Exclude abundunce less than 1%
stacked_df = df_18s[df_18s["rel_abund"] >= 1].copy()

print(stacked_df)

# Pivot table: sample_id vs taxon
pivot = stacked_df.pivot_table(index="sample_id", columns="taxon", values="rel_abund", aggfunc="sum", fill_value=0)

# Reorder taxa by total abundance
taxa_order = pivot.sum().sort_values(ascending=False).index.tolist()
pivot = pivot[taxa_order]

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
    "Cloacimonetes": colors[25]
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
ax.set_ylabel("Relative Abundance (%)", fontsize=16)
ax.set_xticks(np.arange(len(pivot.index)))
ax.set_xticklabels(pivot.index, rotation=45, ha="right", fontsize=12)
ax.set_title("Relative Abundance of 18S Eukaryotic Taxa by Sample", fontsize=14)
ax.set_xlabel("")
ax.set_ylim(0, 100)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(title="Taxon", bbox_to_anchor=(1.05, 1), fontsize=10, title_fontsize=11)
plt.tight_layout()
plt.show()
