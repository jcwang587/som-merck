import matplotlib.pyplot as plt
import pandas as pd


# Load the dataset
dataset = pd.read_csv("../../data/dataset/dataset_merck_all.csv")

dataset = dataset[dataset["bde"].notna()]
som_df = dataset[dataset["primary_som"] == 1]

non_som_df = dataset[dataset["primary_som"] == 0]

print(f"number of SOM sites: {len(som_df)}")
print(f"number of non SOM sites: {len(non_som_df)}")

bde_min = 70
bde_max = 115
bins = [bde_min + i * (bde_max - bde_min) / 15 for i in range(16)]

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
fig, ax = plt.subplots(figsize=(8, 6))

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
    patch.set_facecolor("#00E47C")

patches[-1].set_facecolor("black")
patches[-3].set_facecolor("grey")

# Annotate each bar with the count
for patch in patches:
    height = patch.get_height()
    if height > 0:
        ax.annotate(
            f"{int(height)}",
            xy=(
                patch.get_x() + patch.get_width() / 2,
                height + 0.5,
            ),
            xytext=(0, 0),
            textcoords="offset points",
            ha="center",
            va="center",
            fontsize=16,
        )

# Set x-axis ticks from 68 to 106 with a step of 4
ax.set_xticks(range(70, 115, 6))

# rotate the x-axis labels
ax.set_xlabel("BDE (kcal/mol)", fontsize=18, fontweight="bold")
ax.set_ylabel("Number of Sites", fontsize=18, fontweight="bold")


# Add a vertical dash line at x=94
ax.axvline(x=94, color="black", linestyle="--", linewidth=3)
# ax.axvline(x=106, color="black", linestyle="--", linewidth=3)
# ax.axvline(x=112, color="black", linestyle="--", linewidth=3)
ax.set_xlim(67, 118)
ax.set_ylim(0, 25)


# Remove the top and right spines
ax.spines["left"].set_visible(False)

# Move the y-axis to the right
ax.yaxis.set_ticks_position("right")
ax.yaxis.set_label_position("right")

# Set the axis line width
ax.spines["bottom"].set_linewidth(1.5)
ax.spines["top"].set_linewidth(1.5)
ax.spines["right"].set_linewidth(1.5)

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

plt.tight_layout()
plt.savefig("./bde_distribution_som.png", transparent=False)
plt.close()

# Determine the last bin range
last_bin_range = pd.Interval(bins[-4], bins[-3], closed="right")
print(last_bin_range)

# Filter the som_df for entries in the last bin
last_bin_entries = som_df[som_df["bde"].apply(lambda x: x in last_bin_range)]

# Print the Zaretzki index of the entries in the last bin
print("Zaretzki index of the last bin:")
print(last_bin_entries["molecule_title"])
print(last_bin_entries["atomic_index"])

print(last_bin_entries["bde"])
print(last_bin_entries["sasa_hydrogen_maestro"])


df_low_bde_som = pd.DataFrame(
    columns=[
        "zaretzki_index",
        "atomic_index",
        "zaretzki_atomic_index",
        "bde",
        "som_level",
    ]
)
for i in range(0, 8):
    df_low_bde_som = pd.concat(
        [
            df_low_bde_som,
            som_df[som_df["bde"].between(bins[i], bins[i + 1])][
                [
                    "zaretzki_index",
                    "atomic_index",
                    "zaretzki_atomic_index",
                    "bde",
                    "som_level",
                ]
            ],
        ],
        ignore_index=True,
    )

df_low_bde_som.to_csv("../../data/bins/bde_som/low_bde_som.csv", index=False)