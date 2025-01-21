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
    [som_df["relative_ir"], non_som_df["relative_ir"]],
    bins=bins,
    label=["Experimental SOM", "Experimental Non SOM"],
    color=["#CC5F5A", "#82ABA3"],
    edgecolor="black",
    stacked=True,
)

# Annotate the bars with counts
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


# Initialize the dataframe first column is the zaretzki index and the second column is the atomic index
df_high_relative_ir_som = pd.DataFrame(
    columns=[
        "zaretzki_index",
        "atomic_index",
        "zaretzki_atomic_index",
        "relative_ir",
        "som_level",
    ]
)
df_high_relative_ir_non_som = pd.DataFrame(
    columns=[
        "zaretzki_index",
        "atomic_index",
        "zaretzki_atomic_index",
        "relative_ir",
        "som_level",
    ]
)

# For loop from 10 to 20
for i in range(10, 20):
    df_high_relative_ir_som = pd.concat(
        [
            df_high_relative_ir_som,
            som_df[som_df["relative_ir"].between(bins[i], bins[i + 1])][
                [
                    "zaretzki_index",
                    "atomic_index",
                    "zaretzki_atomic_index",
                    "relative_ir",
                    "som_level",
                ]
            ],
        ],
        ignore_index=True,
    )
    df_high_relative_ir_non_som = pd.concat(
        [
            df_high_relative_ir_non_som,
            non_som_df[non_som_df["relative_ir"].between(bins[i], bins[i + 1])][
                [
                    "zaretzki_index",
                    "atomic_index",
                    "zaretzki_atomic_index",
                    "relative_ir",
                    "som_level",
                ]
            ],
        ],
        ignore_index=True,
    )


# Save the dataframe to a csv file
df_high_relative_ir_som = df_high_relative_ir_som.drop_duplicates(
    subset="zaretzki_atomic_index"
)
df_high_relative_ir_non_som = df_high_relative_ir_non_som.drop_duplicates(
    subset="zaretzki_atomic_index"
)
df_high_relative_ir_som.to_csv(
    "../../data/bins/relative_ir_som/high_relative_ir_som.csv", index=False
)
df_high_relative_ir_non_som.to_csv(
    "../../data/bins/relative_ir_som/high_relative_ir_non_som.csv", index=False
)

