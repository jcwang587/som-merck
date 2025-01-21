import matplotlib.pyplot as plt
import pandas as pd


# Load the dataset
dataset = pd.read_csv("../../data/dataset/dataset_merck_bde.csv")

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
    [som_df["sasa_hydrogen_maestro"], non_som_df["sasa_hydrogen_maestro"]],
    bins=bins,
    label=["Experimental SOM", "Experimental Non SOM"],
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

# Set labels and legend
ax.set_xlabel("SASA")
ax.set_ylabel("Count")
ax.legend(frameon=False, loc="upper right")

# Add a vertical dash line at x=18
ax.axvline(x=18, color="black", linestyle="--", linewidth=2)
ax.set_xlim(-2, 42)
ax.set_ylim(0, 60)

plt.tight_layout()
plt.savefig("./sasa_hydrogen_distribution_all.png")
plt.close()


# Initialize the dataframe first column is the zaretzki index and the second column is the atomic index
df_high_sasa_maestro_som = pd.DataFrame(
    columns=[
        "zaretzki_index",
        "atomic_index",
        "zaretzki_atomic_index",
        "sasa_hydrogen_maestro",
        "som_level",
    ]
)
df_high_sasa_maestro_non_som = pd.DataFrame(
    columns=[
        "zaretzki_index",
        "atomic_index",
        "zaretzki_atomic_index",
        "sasa_hydrogen_maestro",
        "som_level",
    ]
)

for i in range(9, 20):
    df_high_sasa_maestro_som = pd.concat(
        [
            df_high_sasa_maestro_som,
            som_df[som_df["sasa_hydrogen_maestro"].between(bins[i], bins[i + 1])][
                [
                    "zaretzki_index",
                    "atomic_index",
                    "zaretzki_atomic_index",
                    "sasa_hydrogen_maestro",
                    "som_level",
                ]
            ],
        ],
        ignore_index=True,
    )

    df_high_sasa_maestro_non_som = pd.concat(
        [
            df_high_sasa_maestro_non_som,
            non_som_df[
                non_som_df["sasa_hydrogen_maestro"].between(bins[i], bins[i + 1])
            ][
                [
                    "zaretzki_index",
                    "atomic_index",
                    "zaretzki_atomic_index",
                    "sasa_hydrogen_maestro",
                    "som_level",
                ]
            ],
        ],
        ignore_index=True,
    )


# Save the dataframe to a csv file
df_high_sasa_maestro_som = df_high_sasa_maestro_som.drop_duplicates(
    subset="zaretzki_atomic_index"
)
df_high_sasa_maestro_non_som = df_high_sasa_maestro_non_som.drop_duplicates(
    subset="zaretzki_atomic_index"
)
df_high_sasa_maestro_som.to_csv(
    "../../data/bins/sasa_som/high_sasa_maestro_som.csv", index=False
)
df_high_sasa_maestro_non_som.to_csv(
    "../../data/bins/sasa_som/high_sasa_maestro_non_som.csv", index=False
)

# Get the average sasa
print(f"average sasa of all: {dataset['sasa_hydrogen_maestro'].mean()}")
