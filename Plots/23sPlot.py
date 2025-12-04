import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load 23s relative abundance CSVs
phylum_df = pd.read_csv("23S_abund_phylum.csv")
class_df = pd.read_csv("23S_abund_class.csv")
order_df = pd.read_csv("23S_abund_order.csv")
family_df = pd.read_csv("23S_abund_family.csv")
genus_df = pd.read_csv("23S_abund_genus.csv")

# Add 'level' column and combine into one DataFrame
for df, level in zip(
    [phylum_df, class_df, order_df, family_df, genus_df],
    ["Phylum", "Class", "Order", "Family", "Genus"]
):
    df["level"] = level
df = pd.concat([phylum_df, class_df, order_df, family_df, genus_df], ignore_index=True)

print(df)

# Keep only rows from section '81RABR' and only look at one taxonomy
df = df[df["section"] == "81RABR"].copy()
df = df[df["level"] == "Phylum"].copy()

# Extract all classifications that do not appear in all samples (misclassifications, contaminations)
taxon_counts = df.groupby("taxon")["sample_id"].nunique()                       # Step 1: Count how many unique samples each taxon appears in
total_samples = df["sample_id"].nunique()                                       # Step 2: Find total number of unique samples
incomplete_taxa = taxon_counts[taxon_counts < total_samples].index.tolist()     # Step 3: Identify taxa that do NOT appear in every sample
filtered_df = df[df["taxon"].isin(incomplete_taxa)]                             # Step 4: Filter original DataFrame for only those taxa
clean_df = df[~df["taxon"].isin(incomplete_taxa)]                               # Remove rows with taxa not found in every sample

# Group by section, taxon, and level, then calculate mean relative abundance
summary_df = (
    clean_df.groupby(["section", "taxon", "level"])["rel_abund"]
    .agg(mean_abundance="mean", std_abundance="std")
    .reset_index()
)

summary_df = summary_df[summary_df["mean_abundance"] >= 1].copy()

# Remove *
summary_df["taxon"] = summary_df["taxon"].str.replace("*", "", regex=False)

# Define plotting order by average abundance
ordered_taxa = (
    summary_df.groupby("taxon")["mean_abundance"]
    .mean()
    .sort_values(ascending=False)
    .index.tolist()
)

# Ensure categorical ordering
summary_df["taxon"] = pd.Categorical(summary_df["taxon"], categories=ordered_taxa, ordered=True)

# Replace long taxon name with line break version
summary_df["taxon"] = summary_df["taxon"].cat.rename_categories(
    lambda x: "Unclassified\nChlorophyta" if x == "Unclassified Chlorophyta" else x
)
summary_df["taxon"] = summary_df["taxon"].cat.rename_categories(
    lambda x: "Unclassified\nCyanophyceae" if x == "Unclassified Cyanophyceae" else x
)
summary_df["taxon"] = summary_df["taxon"].cat.rename_categories(
    lambda x: "Unclassified\nPlantae" if x == "Unclassified Plantae" else x
)

# Plot
green = sns.color_palette("Set2")[0]  # Index 1 is the green in Set2

plt.figure(figsize=(12, 8))
ax = sns.barplot(
    data=summary_df,
    x="taxon",
    y="mean_abundance",
    errorbar=None,
    color=green,
    width=0.35
)

# Add error bars using actual bar positions
for patch, (_, row) in zip(ax.patches, summary_df.iterrows()):
    mean = row["mean_abundance"]
    std = row["std_abundance"]

    # Compute bounds and clip at zero
    lower = max(0, mean - std)
    upper = mean + std
    yerr_lower = max(0, mean - lower)
    yerr_upper = max(0, upper - mean)

    # Bar center
    bar_x = patch.get_x() + patch.get_width() / 2
    bar_height = patch.get_height()

    ax.errorbar(
        x=bar_x,
        y=bar_height,
        yerr=[[yerr_lower], [yerr_upper]],
        fmt='none',
        capsize=5,
        color='black',
        lw=1.5
    )

# Ensure y-axis starts at 0
ax.set_ylim(bottom=0,top=100)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Format plot
plt.xticks(rotation=45, ha="right", fontsize=18)
plt.yticks(fontsize=18)
plt.ylabel("Mean Relative Abundance (%)", fontsize=18)
plt.xlabel("")
#plt.title("Mean Relative Abundance by Taxon (Mothur vs DADA2) with Error Bars")
plt.tight_layout()
plt.show()
