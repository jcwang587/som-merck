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

bde_min = 64
bde_max = 108
bins = [bde_min + i * (bde_max - bde_min) / 22 for i in range(23)]

som_bde = (
    pd.cut(som_df["bde"], bins=bins, include_lowest=True).value_counts().sort_index()
)
non_som_bde = (
    pd.cut(non_som_df["bde"], bins=bins, include_lowest=True)
    .value_counts()
    .sort_index()
)
total_bde = som_bde + non_som_bde

# Create the stacked bar chart for reactivity
plt.rcParams["font.size"] = 14
fig, ax = plt.subplots(figsize=(14, 6))

# Plot histogram for SOM only
n, bins, patches = ax.hist(
    som_df["bde"],
    bins=bins,
    label=["Experimental SOM"],
    color=["#08312a"],
    edgecolor="white",
    stacked=True,
)

# Annotate each bar with the count
for patch in patches:
    height = patch.get_height()
    if height > 0: 
        ax.annotate(
            f"{int(height)}",
            xy=(
                patch.get_x() + patch.get_width() / 2,
                height + 1,
            ),
            xytext=(0, 0),  
            textcoords="offset points",
            ha="center",
            va="center",
            fontsize=16,
        )

# rotate the x-axis labels
ax.set_xlabel("BDE (kcal/mol)", fontsize=18, fontweight="bold")
ax.set_ylabel("Count", fontsize=18, fontweight="bold")

# Increase the font size of the tick labels
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

# Add a vertical dash line at x=94
ax.axvline(x=94, color="black", linestyle="--", linewidth=3)
ax.set_xlim(68, 107)
ax.set_ylim(0, 25)
# ax.legend(frameon=False, loc="upper left")

plt.tight_layout()
plt.savefig("./bde_distribution_som.png")
plt.close()
