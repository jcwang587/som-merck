import matplotlib.pyplot as plt
import pandas as pd


# Load the dataset
dataset = pd.read_csv("../../data/dataset/dataset_merck_bde.csv")

dataset = dataset[dataset["bde"].notna()]
som_df = dataset[dataset["primary_som"] == 1]

non_som_df = dataset[dataset["primary_som"] == 0]

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

# Set the last bins to be light green
for patch in patches[:-7]:
    patch.set_facecolor('#00E47C')

# Annotate each bar with the count
for patch in patches:
    height = patch.get_height()
    if height > 0: 
        ax.annotate(
            f"{int(height)}",
            xy=(
                patch.get_x() + patch.get_width() / 2,
                height + 1 / 30 * 18,
            ),
            xytext=(0, 0),  
            textcoords="offset points",
            ha="center",
            va="center",
            fontsize=16,
        )

# Set x-axis ticks from 68 to 106 with a step of 4
ax.set_xticks(range(70, 106, 4))

# rotate the x-axis labels
ax.set_xlabel("BDE (kcal/mol)", fontsize=18, fontweight="bold")
ax.set_ylabel("Count", fontsize=18, fontweight="bold")

# Increase the font size of the tick labels
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

# Add a vertical dash line at x=94
ax.axvline(x=94, color="black", linestyle="--", linewidth=3)
ax.set_xlim(68, 106)
ax.set_ylim(0, 18)
# ax.legend(frameon=False, loc="upper left")

plt.tight_layout()
plt.savefig("./bde_distribution_som.png")
plt.close()

# Determine the last bin range
last_bin_range = pd.Interval(bins[-8], bins[-7], closed='right')
print(last_bin_range)

# Filter the som_df for entries in the last bin
last_bin_entries = som_df[som_df["bde"].apply(lambda x: x in last_bin_range)]

# Print the Zaretzki index of the entries in the last bin
print("Zaretzki index of the last bin:")
print(last_bin_entries["zaretzki_index"])
print(last_bin_entries["atomic_index"])  
