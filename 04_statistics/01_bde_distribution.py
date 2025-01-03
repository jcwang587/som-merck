import matplotlib.pyplot as plt
import pandas as pd
import sys

sys.path.append("../src/")
from utils import renew_folder


# Load the dataset
dataset = pd.read_csv("../data/dataset/dataset_merck_all.csv")

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

# Plot stacked histograms for SOM and Non-SOM
n, bins, patches = ax.hist(
    [som_df["bde"], non_som_df["bde"]],
    bins=bins,
    label=["SOM", "Non SOM"],
    color=["#CC5F5A", "#82ABA3"],
    edgecolor="black",
    stacked=True,
)

# Annotate each bar with the count
for i in range(len(patches)):
    for j in range(len(patches[i])):
        height = patches[i][j].get_height()
        if height > 0:  # Only annotate if the height is greater than 0
            if i == 0:  # SOM part
                ax.annotate(
                    f"{int(height)}",
                    xy=(
                        patches[i][j].get_x() + patches[i][j].get_width() / 2,
                        height / 2,
                    ),
                    xytext=(0, 0),  # No vertical offset
                    textcoords="offset points",
                    ha="center",
                    va="center",
                )
            else:  # Non-SOM part
                som_height = patches[0][j].get_height()
                ax.annotate(
                    f"{int(height)}",
                    xy=(
                        patches[i][j].get_x() + patches[i][j].get_width() / 2,
                        som_height + height,
                    ),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    color="#82ABA3",
                )

# rotate the x-axis labels
ax.set_xlabel("BDE")
ax.set_ylabel("Count")

# Add a vertical dash line at x=94
ax.axvline(x=94, color="black", linestyle="--", linewidth=2)
ax.set_xlim(62, 110)
ax.set_ylim(0, 100)
ax.legend(frameon=False, loc="upper left")

plt.tight_layout()
plt.savefig("./bde_distribution_som.png")
plt.close()

# Initialize the dataframe first column is the zaretzki index and the second column is the atomic index
df_low_bde_som = pd.DataFrame(
    columns=[
        "zaretzki_index",
        "atomic_index",
        "zaretzki_atomic_index",
        "bde",
        "som_level",
    ]
)
df_low_bde_non_som = pd.DataFrame(
    columns=[
        "zaretzki_index",
        "atomic_index",
        "zaretzki_atomic_index",
        "bde",
        "som_level",
    ]
)

# For loop from 1 to 15
for i in range(1, 15):
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
    df_low_bde_non_som = pd.concat(
        [
            df_low_bde_non_som,
            non_som_df[non_som_df["bde"].between(bins[i], bins[i + 1])][
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


# Save the dataframe to a csv file
df_low_bde_som = df_low_bde_som.drop_duplicates(subset="zaretzki_atomic_index")
df_low_bde_non_som = df_low_bde_non_som.drop_duplicates(subset="zaretzki_atomic_index")
df_low_bde_som.to_csv("../data/bins/bde_som/low_bde_som.csv", index=False)
df_low_bde_non_som.to_csv("../data/bins/bde_som/low_bde_non_som.csv", index=False)
