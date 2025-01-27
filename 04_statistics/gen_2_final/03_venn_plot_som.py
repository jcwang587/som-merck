import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles
import matplotlib_venn
import pandas as pd


# Load the data
all_data = pd.read_csv("../../data/dataset/dataset_merck_all.csv")
bde_som = pd.read_csv("../../data/bins/bde_som/low_bde_som.csv")
sasa_som = pd.read_csv("../../data/bins/sasa_som/high_sasa_maestro_som.csv")

print(len(bde_som))
print(len(sasa_som))

# Plot the venn diagram
plt.figure(figsize=(10, 8))

v = venn2(
    [
        set(sasa_som["zaretzki_atomic_index"]),
        set(bde_som["zaretzki_atomic_index"]),
    ],
    set_labels=None  # Hide the labels
)

# Set custom colors
v.get_patch_by_id('10').set_color('#08312A')
v.get_patch_by_id('01').set_color('#08312A') 
v.get_patch_by_id('11').set_color('#00E47C')  

v.get_patch_by_id('10').set_alpha(1)
v.get_patch_by_id('01').set_alpha(1)
v.get_patch_by_id('11').set_alpha(1)


c = venn2_circles(
    [
        set(sasa_som["zaretzki_atomic_index"]),
        set(bde_som["zaretzki_atomic_index"]),
    ],
    linestyle="dashed",
    linewidth=2
)

# Increase the font size of the numbers inside the circles
for idx, text in enumerate(v.subset_labels):
    if text is not None:
        text.set_fontsize(28)
        text.set_fontweight("bold")
    if idx in [0, 1]:  
        text.set_color('white')

# Find the unique zaretzki_atomic_index
unique_index = (
    set(bde_som["zaretzki_atomic_index"])
    | set(sasa_som["zaretzki_atomic_index"])
)

plt.savefig("./venn_plot_som.png", dpi=300, bbox_inches="tight", transparent=True)

# print the zaretzki_atomic_index of the atoms in low bde only
low_bde_only = set(bde_som["zaretzki_atomic_index"]) - set(sasa_som["zaretzki_atomic_index"])
print(f"low_bde_only: {low_bde_only}")

# print the zaretzki_atomic_index of the atoms in high relative ir only
high_sasa_only = set(sasa_som["zaretzki_atomic_index"]) - set(bde_som["zaretzki_atomic_index"])
print(f"high_sasa_only: {high_sasa_only}")

# print the zaretzki_atomic_index of the atoms
# remove na rows
all_data = all_data.dropna()
all_som_index = all_data[all_data["som"] == 1]["zaretzki_atomic_index"].tolist()

others = set(all_som_index) - unique_index
print(f"others: {others}")

