import matplotlib.pyplot as plt
import pandas as pd


# Load the dataset
dataset = pd.read_csv("../../data/dataset/dataset_merck_bde.csv")

dataset = dataset[dataset["bde"].notna()]

som_df = dataset[dataset["primary_som"] == 1]
non_som_df = dataset[dataset["primary_som"] == 0]

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
    color=["#08312a"],
    edgecolor="white",
    stacked=True,
)

# Set the last bins to be light green
for patch in patches[5:]:
    patch.set_facecolor('#00E47C')

# Annotate each bar with the count
for patch in patches:
    height = patch.get_height()
    if height > 0: 
        ax.annotate(
            f"{int(height)}",
            xy=(
                patch.get_x() + patch.get_width() / 2,
                height + 1 / 30 * 14,
            ),
            xytext=(0, 0),  
            textcoords="offset points",
            ha="center",
            va="center",
            fontsize=16,
        )



# Set x-axis ticks from 68 to 106 with a step of 4
ax.set_xticks(range(0, 36, 4))

# Set labels and legend
ax.set_xlabel("SASA (Å²)", fontsize=18, fontweight="bold")
ax.set_ylabel("Count", fontsize=18, fontweight="bold")
# ax.legend(frameon=False, loc="upper left")

# Add a vertical dash line at x=10
ax.axvline(x=10, color="black", linestyle="--", linewidth=3)
ax.set_xlim(2, 36)
ax.set_ylim(0, 14)

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

plt.tight_layout()
plt.savefig("./sasa_hydrogen_distribution_all.png")
plt.close()

# Determine the last bin range
last_bin_range = pd.Interval(bins[4], bins[5], closed='right')
print(last_bin_range)

# Filter the som_df for entries in the last bin
last_bin_entries = som_df[som_df["sasa_hydrogen_maestro"].apply(lambda x: x in last_bin_range)]

# Print the Zaretzki index of the entries in the last bin
print("Zaretzki index of the last bin:")
print(last_bin_entries["molecule_title"])
print(last_bin_entries["atomic_index"])  

print(last_bin_entries["bde"])
print(last_bin_entries["sasa_hydrogen_maestro"])