import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Load the data, label pipelines, and combine
mothur_df = pd.read_csv("16s_data_mothur/16s_data_mothur.csv")
dada2_df = pd.read_csv("16s_data_DADA2/16s_data_DADA2.csv")
mothur_df["Pipeline"] = "Mothur"
dada2_df["Pipeline"] = "DADA2"
combined_df = pd.concat([mothur_df, dada2_df], ignore_index=True)

# Filter Data (RABR and Phylum) and exclude same RABR at different time
combined_df = combined_df[
    (combined_df["section"] == "81RABR")
    & (combined_df["level"] == "Phylum")
    & (~combined_df["sample_id"].isin(["12_R60_11_15_21_16S", "16_R45_11_15_21_16S"]))
].copy()

print(combined_df)

# Write csv
combined_df.to_csv("16s_data_81RABR.csv", index=False)

# Filter temperature and exclude abundunce less than 1%
df_10 = combined_df[combined_df["temp[C]"] == 10.0].copy()
df_10 = df_10[df_10["rel_abund"] >= 1].copy()

df_20 = combined_df[combined_df["temp[C]"] == 20.0].copy()
df_20 = df_20[df_20["rel_abund"] >= 1].copy()

df_30 = combined_df[combined_df["temp[C]"] == 30.0].copy()
df_30 = df_30[df_30["rel_abund"] >= 1].copy()

# Select temp (10, 20, 30)
stacked_df = df_30

# Filter for Mothur and DADA2 only
mothur_df = stacked_df[stacked_df["Pipeline"] == "Mothur"].copy()
dada2_df = stacked_df[stacked_df["Pipeline"] == "DADA2"].copy()

# Pivot to get relative abundances by sample
mothur_pivot = mothur_df.pivot_table(
    index="label", columns="taxon", values="rel_abund", aggfunc="sum", fill_value=0
)
dada2_pivot = dada2_df.pivot_table(
    index="label", columns="taxon", values="rel_abund", aggfunc="sum", fill_value=0
)

# Ensure consistent taxon columns
common_taxa = sorted(set(mothur_pivot.columns).union(dada2_pivot.columns))
mothur_pivot = mothur_pivot.reindex(columns=common_taxa, fill_value=0)
dada2_pivot = dada2_pivot.reindex(columns=common_taxa, fill_value=0)

# Reorder taxa: largest total abundance first (for stacking order)
taxa_order = (
    pd.concat([mothur_pivot, dada2_pivot])
    .sum()
    .sort_values(ascending=False)
    .index.tolist()
)
common_taxa = taxa_order

# Build new DataFrame for grouped bar plotting
labels = mothur_pivot.index.tolist()
bar_width = 0.35
spacing = 0.1  # small spacing between Mothur and DADA2

# label_order = ["R60_11_01_21", "R43_11_15_21", "R72_11_15_21"]

# label_order = ["R45_10_18_21", "R36_10_28_21", "R58_10_28_21", "R46_11_05_21"]

# labels = label_order


# Calculate positions
x = np.arange(len(labels)) * 2  # Leave space between each sample group
positions_mothur = x - bar_width / 2 - spacing / 2
positions_dada2 = x + bar_width / 2 + spacing / 2

# Plot setup
fig, ax = plt.subplots(figsize=(12, 6))
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

taxon_colors = {taxon: manual_colors.get(taxon, "#cccccc") for taxon in common_taxa}

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
present_taxa = [taxon for taxon in custom_order if taxon in mothur_pivot.columns]

# Plot Mothur (stacked)
bottom = np.zeros(len(labels))
for taxon in present_taxa:
    heights = mothur_pivot[taxon].values
    ax.bar(
        positions_mothur,
        heights,
        bar_width,
        bottom=bottom,
        label=taxon if taxon not in bottom else "",
        color=taxon_colors[taxon],
    )
    bottom += heights

# Plot DADA2 (stacked)
bottom = np.zeros(len(labels))
for taxon in present_taxa:
    heights = dada2_pivot[taxon].values
    ax.bar(
        positions_dada2, heights, bar_width, bottom=bottom, color=taxon_colors[taxon]
    )
    bottom += heights

# Ticks and labels
xticks_combined = [(m + d) / 2 for m, d in zip(positions_mothur, positions_dada2)]
ax.set_xticks(xticks_combined)
ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=12)
ax.set_ylabel("Relative Abundance (%)")
ax.set_xlabel("")
ax.set_title(
    "Relative Abundance of 16S Bacterial Taxa by Sample (Mothur vs. DADA2, 25°C)",
    fontsize=14,
)

# Aesthetics
ax.set_ylim(0, 100)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(title="Taxon", bbox_to_anchor=(1.05, 1))
ax.set_xlabel("")
ax.set_ylabel("Relative Abundance (%)", fontsize=16)
plt.tight_layout()
plt.savefig("16S_Lab_25.pdf", format="pdf", bbox_inches="tight")
plt.show()
