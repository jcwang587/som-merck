import matplotlib.pyplot as plt
import pandas as pd


# Load the dataset
dataset = pd.read_csv("../../data/dataset/dataset_merck_bde.csv")

dataset = dataset[dataset["bde"].notna()]
som_df = dataset[dataset["primary_som"] == 1]

non_som_df = dataset[dataset["primary_som"] == 0]

print(f"number of SOM sites: {len(som_df)}")
print(f"number of non SOM sites: {len(non_som_df)}")

reactivity_min = 0
reactivity_max = 1
bins = [reactivity_min + i * (reactivity_max - reactivity_min) / 20 for i in range(21)]

som_ir = (
    pd.cut(som_df["relative_ir"], bins=bins, include_lowest=True)
    .value_counts()
    .sort_index()
)
non_som_ir = (
    pd.cut(non_som_df["relative_ir"], bins=bins, include_lowest=True)
    .value_counts()
    .sort_index()
)
total_ir = som_ir + non_som_ir

# Create the stacked bar chart for reactivity
plt.rcParams["font.size"] = 14
fig, ax = plt.subplots(figsize=(14, 6))

# Plot stacked histograms for SOM and Non-SOM
n, bins, patches = ax.hist(
    som_df["relative_ir"],
    bins=bins,
    label=["Experimental SOM"],
    color=["#CC5F5A"],
    edgecolor="black",
    stacked=True,
)

# Annotate the bars with counts
for patch in patches:
    height = patch.get_height()
    if height > 0:  # Only annotate if the height is greater than 0
        ax.annotate(
            f"{int(height)}",
            xy=(
                patch.get_x() + patch.get_width() / 2,
                height + 1,
            ),
            xytext=(0, 0),  # No vertical offset
            textcoords="offset points",
            ha="center",
            va="center",
        )

# Set labels and legend
ax.set_xlabel("Relative Intrinsic Reactivity")
ax.set_ylabel("Count")

# Add a vertical dash line at x=0.5
ax.axvline(x=0.5, color="black", linestyle="--", linewidth=2)

ax.set_xlim(-0.05, 1.05)
ax.set_ylim(0, 90)
ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.25, 1))

plt.tight_layout()
plt.savefig("./relative_ir_distribution_som.png")
plt.close()

# Print the Zaretzki index of the entries in the 5th bin
fifth_bin_min = bins[4]
fifth_bin_max = bins[5]
print(som_df[(som_df["relative_ir"] >= fifth_bin_min) & (som_df["relative_ir"] < fifth_bin_max)])

ninth_bin_min = bins[8]
ninth_bin_max = bins[9]
print(som_df[(som_df["relative_ir"] >= ninth_bin_min) & (som_df["relative_ir"] < ninth_bin_max)])

eleventh_bin_min = bins[10]
eleventh_bin_max = bins[11]
print(som_df[(som_df["relative_ir"] >= eleventh_bin_min) & (som_df["relative_ir"] < eleventh_bin_max)])
