import matplotlib.pyplot as plt
import pandas as pd


# Load the dataset
dataset = pd.read_csv("../../data/dataset/dataset_merck_all.csv")

dataset = dataset[dataset["bde"].notna()]

primary_som_df = dataset[dataset["primary_som"] == 1]
secondary_som_df = dataset[dataset["secondary_som"] == 1]
tertiary_som_df = dataset[dataset["tertiary_som"] == 1]
som_df = pd.concat([primary_som_df, secondary_som_df, tertiary_som_df])
non_som_df = dataset[
    (dataset["primary_som"] == 0)
    & (dataset["secondary_som"] == 0)
    & (dataset["tertiary_som"] == 0)
]

print(f"number of SOM sites: {len(som_df)}")
print(f"number of non SOM sites: {len(non_som_df)}")

sasa_maestro_min = 0
sasa_maestro_max = 40
bins = [
    sasa_maestro_min + i * (sasa_maestro_max - sasa_maestro_min) / 20 for i in range(21)
]

som_sasa = (
    pd.cut(som_df["sasa_hydrogen_maestro"], bins=bins, include_lowest=True)
    .value_counts()
    .sort_index()
)
non_som_sasa = (
    pd.cut(non_som_df["sasa_hydrogen_maestro"], bins=bins, include_lowest=True)
    .value_counts()
    .sort_index()
)
total_sasa = som_sasa + non_som_sasa

# Create the histogram for reactivity
plt.rcParams["font.size"] = 14
fig, ax = plt.subplots(figsize=(14, 6))

# Plot stacked histograms for SOM and Non-SOM
n, bins, patches = ax.hist(
    som_df["sasa_hydrogen_maestro"],
    bins=bins,
    label=["Experimental SOM"],
    color=["#CC5F5A"],
    edgecolor="black",
    stacked=True,
)

# Annotate each bar with the count
for patch in patches:
    height = patch.get_height()
    if height > 0:  # Only annotate if the height is greater than 0
        ax.annotate(
            f"{int(height)}",
            xy=(
                patch.get_x() + patch.get_width() / 2,
                height + 1 / 25 * 14,
            ),
            xytext=(0, 0),  # No vertical offset
            textcoords="offset points",
            ha="center",
            va="center",
        )
    
# Set labels and legend
ax.set_xlabel("SASA")
ax.set_ylabel("Count")
ax.legend(frameon=False, loc="upper right")

# Add a vertical dash line at x=18
ax.axvline(x=18, color="black", linestyle="--", linewidth=2)
ax.set_xlim(2, 36)
ax.set_ylim(0, 14)

plt.tight_layout()
plt.savefig("./sasa_hydrogen_distribution_all.png")
plt.close()
