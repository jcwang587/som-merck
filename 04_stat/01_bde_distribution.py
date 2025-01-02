import matplotlib.pyplot as plt
import pandas as pd
import sys

sys.path.append("../src/")
from utils import renew_folder


# Load the dataset
dataset = pd.read_csv("../../data/dataset/dataset_all_sites.csv")

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
ax.set_ylim(0, 400)
ax.legend(frameon=False, loc="upper left")

plt.tight_layout()
plt.savefig("./bde_distribution_som.png")
plt.close()

# Print the zaretzki index of the samples within the first label class
renew_folder("../../data/bins/bde_som/")

# show all when printing
pd.set_option("display.max_rows", None)
print(
    primary_som_df[primary_som_df["bde"].between(bins[8], bins[9])][
        ["zaretzki_index", "atomic_index", "bde"]
    ]
)
print(
    primary_som_df[primary_som_df["bde"].between(bins[8], bins[9])][
        "zaretzki_index"
    ].nunique()
)

# generate the dataframe for the high bde
for i in range(len(bins) - 1):
    primary_som_df_export = pd.DataFrame()
    primary_som_df_export = pd.concat(
        [
            primary_som_df_export,
            primary_som_df[primary_som_df["bde"].between(bins[i], bins[i + 1])][
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
    primary_som_df_export.to_csv(
        f"../../data/bins/bde_som/primary_som_bde_{bins[i]}_{bins[i+1]}.csv",
        index=False,
    )

for i in range(len(bins) - 1):
    secondary_som_df_export = pd.DataFrame()
    secondary_som_df_export = pd.concat(
        [
            secondary_som_df_export,
            secondary_som_df[secondary_som_df["bde"].between(bins[i], bins[i + 1])][
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
    secondary_som_df_export.to_csv(
        f"../../data/bins/bde_som/secondary_som_bde_{bins[i]}_{bins[i+1]}.csv",
        index=False,
    )

for i in range(len(bins) - 1):
    tertiary_som_df_export = pd.DataFrame()
    tertiary_som_df_export = pd.concat(
        [
            tertiary_som_df_export,
            tertiary_som_df[tertiary_som_df["bde"].between(bins[i], bins[i + 1])][
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
    tertiary_som_df_export.to_csv(
        f"../../data/bins/bde_som/tertiary_som_bde_{bins[i]}_{bins[i+1]}.csv",
        index=False,
    )

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

# Add the zaretzki index and the atomic index to the dataframe
for i in range(15):
    # append the zaretzki index and the atomic index to the dataframe
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

df_low_bde_som.to_csv("../../data/bins/bde_som/low_bde_som.csv", index=False)
df_low_bde_non_som.to_csv("../../data/bins/bde_som/low_bde_non_som.csv", index=False)
print(len(df_low_bde_som))
print(len(df_low_bde_non_som))


# Get the average bde
print(f"average bde of all: {dataset['bde'].mean()}")
