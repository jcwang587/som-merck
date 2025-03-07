import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
import matplotlib_venn
import pandas as pd


# Load the data
all_data = pd.read_csv("../data/dataset/dataset_merck_all.csv")
bde_non_som = pd.read_csv("../data/bins/bde_som/low_bde_non_som.csv")
relative_ir_non_som = pd.read_csv(
    "../data/bins/relative_ir_som/high_relative_ir_non_som.csv"
)
sasa_non_som = pd.read_csv("../data/bins/sasa_som/high_sasa_maestro_non_som.csv")

print(len(bde_non_som))
print(len(relative_ir_non_som))
print(len(sasa_non_som))

# Plot the venn diagram
plt.figure(figsize=(10, 8))

v = venn3(
    [
        set(bde_non_som["zaretzki_atomic_index"]),
        set(relative_ir_non_som["zaretzki_atomic_index"]),
        set(sasa_non_som["zaretzki_atomic_index"]),
    ],
    ["Low BDE", "High Relative IR", "High SASA"],
    layout_algorithm=matplotlib_venn.layout.venn3.DefaultLayoutAlgorithm(
        fixed_subset_sizes=(1, 1, 1, 1, 1, 1, 1)
    ),
)

# Set the colors for the labels
label_colors = ["#F89898", "#9ECC97", "#9B9AFF"]
for label, color in zip(v.set_labels, label_colors):
    label.set_color(color)
    label.set_fontsize(24)
    label.set_fontweight("bold")

c = venn3_circles(
    [
        set(bde_non_som["zaretzki_atomic_index"]),
        set(relative_ir_non_som["zaretzki_atomic_index"]),
        set(sasa_non_som["zaretzki_atomic_index"]),
    ],
    linestyle="dashed",
    linewidth=2,
    layout_algorithm=matplotlib_venn.layout.venn3.DefaultLayoutAlgorithm(
        fixed_subset_sizes=(1, 1, 1, 1, 1, 1, 1)
    ),
)

# Increase the font size of the numbers inside the circles
for text in v.subset_labels:
    if text is not None:
        text.set_fontsize(24)
        text.set_fontweight("bold")

# Find the unique zaretzki_atomic_index
unique_index = (
    set(bde_non_som["zaretzki_atomic_index"])
    | set(relative_ir_non_som["zaretzki_atomic_index"])
    | set(sasa_non_som["zaretzki_atomic_index"])
)


# Add text at the right bottom of the plot
plt.text(
    0.7,
    0.2,
    f"Others: 1",
    fontsize=24,
    fontweight="bold",
    transform=plt.gcf().transFigure,
)


plt.savefig("./venn_plot_non_som.png", dpi=300, bbox_inches="tight")

# print the zaretzki_atomic_index of the atoms in low bde only
low_bde_only = (
    set(bde_non_som["zaretzki_atomic_index"])
    - set(relative_ir_non_som["zaretzki_atomic_index"])
    - set(sasa_non_som["zaretzki_atomic_index"])
)
print(low_bde_only)
print(len(low_bde_only))

# print the zaretzki_atomic_index of the atoms in high relative ir only
high_relative_ir_only = (
    set(relative_ir_non_som["zaretzki_atomic_index"])
    - set(bde_non_som["zaretzki_atomic_index"])
    - set(sasa_non_som["zaretzki_atomic_index"])
)
print(high_relative_ir_only)
print(len(high_relative_ir_only))

# print the zaretzki_atomic_index of the atoms in high sasa only
high_sasa_only = (
    set(sasa_non_som["zaretzki_atomic_index"])
    - set(bde_non_som["zaretzki_atomic_index"])
    - set(relative_ir_non_som["zaretzki_atomic_index"])
)
print(high_sasa_only)
print(len(high_sasa_only))

# print the zaretzki_atomic_index of the atoms
# remove na rows
all_data = all_data.dropna()
all_non_som_index = all_data[all_data["som"] == 0]["zaretzki_atomic_index"].tolist()

others = set(all_non_som_index) - unique_index
print(others)
print(len(others))
