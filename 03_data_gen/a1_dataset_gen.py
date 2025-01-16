from schrodinger.structure import StructureReader
from schrodinger import adapter
from rdkit import Chem
import pandas as pd
import os
import sys

sys.path.append("../src")
from utils import get_primary_som, get_secondary_som, get_tertiary_som


def compute_relative_values(data_dict, min_val, max_val):
    relative_dict = {}
    for atom_index, value in data_dict.items():
        if pd.isna(value):
            relative_dict[atom_index] = pd.NA
        elif max_val == min_val:
            relative_dict[atom_index] = 0.5
        else:
            relative_dict[atom_index] = (value - min_val) / (max_val - min_val)
    return relative_dict


# Initialize all possible data points
all_c_mae_dir = "../data/merck_zaretzki"
all_c_mae_files = [f for f in os.listdir(all_c_mae_dir) if f.endswith(".mae")]
all_c_mae_files.sort()

# Feature 1: Intrinsic Reactivity
# absolute intrinsic reactivity calculated by epik cyp function
epik_cyp_mae_dir = "../data/ir/mae"
epik_cyp_mae_files = [f for f in os.listdir(epik_cyp_mae_dir) if f.endswith(".mae")]
epik_cyp_mae_files.sort()
ir_dict = {}
for mae_file in epik_cyp_mae_files:
    structure = StructureReader.read(os.path.join(epik_cyp_mae_dir, mae_file))
    mae_file_index = mae_file.split("-")[-1].split("_")[0]
    ir_dict[mae_file_index] = {}
    for atom in structure.atom:
        if "r_cyp_CYP_reactivity" in atom.property:
            ir_dict[mae_file_index][atom.index] = -atom.property["r_cyp_CYP_reactivity"]

# Feature 2: Relative Intrinsic Reactivity
# relative_intrinsic_reactivity within each molecule
relative_ir_dict = {}
for mae_file in epik_cyp_mae_files:
    structure = StructureReader.read(os.path.join(epik_cyp_mae_dir, mae_file))
    mae_file_index = mae_file.split("-")[-1].split("_")[0]
    temp_ir = []
    for atom in structure.atom:
        if "r_cyp_CYP_reactivity" in atom.property:
            temp_ir.append(-atom.property["r_cyp_CYP_reactivity"])

    min_ir = min(temp_ir)
    max_ir = max(temp_ir)

    relative_ir_dict[mae_file_index] = {}
    if max_ir == min_ir:
        for atom in structure.atom:
            if "r_cyp_CYP_reactivity" in atom.property:
                relative_ir_dict[mae_file_index][atom.index] = 0.5
    else:
        for atom in structure.atom:
            if "r_cyp_CYP_reactivity" in atom.property:
                relative_ir_dict[mae_file_index][atom.index] = (
                    -atom.property["r_cyp_CYP_reactivity"] - min_ir
                ) / (max_ir - min_ir)


# Feature 3: BDE
# BDE calculated by qmox
qmox_mae_dir = "../data/qmox_mae"
qmox_mae_files = [f for f in os.listdir(qmox_mae_dir) if f.endswith(".mae")]
qmox_mae_files.sort()

qmox_csv_dir = "../data/qmox_csv"
qmox_csv_files = [f for f in os.listdir(qmox_csv_dir) if f.endswith(".csv")]
qmox_csv_files.sort()

bde_dict = {}
for mae_file in qmox_mae_files:
    structure = StructureReader.read(os.path.join(qmox_mae_dir, mae_file))
    mae_file_index = mae_file.split("_")[0]
    bde_dict[mae_file_index] = {}
    for atom in structure.atom:
        if "r_user_CH-BDE" in atom.property:
            bde_dict[mae_file_index][atom.index] = atom.property["r_user_CH-BDE"]
for csv_file in qmox_csv_files:
    csv_file_index = csv_file.split("_")[0]
    bde_dict[csv_file_index] = {}
    # get the first column as the atomic index
    df = pd.read_csv(os.path.join(qmox_csv_dir, csv_file))
    index_list = df.iloc[:, 0].tolist()
    element_list = [index[0] for index in index_list]
    atomic_index_list = [int(index[1:]) for index in index_list]
    bde_list = df.iloc[:, 1].tolist()
    for i in range(len(atomic_index_list)):
        bde_dict[csv_file_index][atomic_index_list[i]] = bde_list[i]


# Feature 4: SASA Hydrogen Maestro
# SASA of average hydrogen atoms calculated by maestro function
sasa_hydrogen_maestro_dir = "../data/sasa"
sasa_hydrogen_maestro_files = [
    f for f in os.listdir(sasa_hydrogen_maestro_dir) if f.endswith(".mae")
]
sasa_hydrogen_maestro_files.sort()

sasa_hydrogen_maestro_dict = {}
for mae_file in sasa_hydrogen_maestro_files:
    structure = StructureReader.read(os.path.join(sasa_hydrogen_maestro_dir, mae_file))
    mae_file_index = mae_file.split("_")[1]
    sasa_hydrogen_maestro_dict[mae_file_index] = {}

    for atom in structure.atom:
        if atom.element != "H":
            hydrogen_atoms = []
            # Find the hydrogen atoms attached to the carbon atoms
            for bond in atom.bond:
                if bond.atom1.element == "H" or bond.atom2.element == "H":
                    hydrogen_atoms.append(
                        bond.atom1 if bond.atom1.element == "H" else bond.atom2
                    )

            hydrogen_sasa = 0
            # Get the SASA of the hydrogen atoms
            for hydrogen_atom in hydrogen_atoms:
                if "r_user_sasa" in hydrogen_atom.property:
                    hydrogen_sasa += hydrogen_atom.property["r_user_sasa"]

            hydrogen_sasa = (
                hydrogen_sasa / len(hydrogen_atoms)
                if len(hydrogen_atoms) > 0
                else pd.NA
            )
            sasa_hydrogen_maestro_dict[mae_file_index][atom.index] = hydrogen_sasa


# Property 1: Atomic number
# atomic number calculated by schrodinger atomic property "i_m_atomic_number"
atomic_number_dict = {}
for mae_file in all_c_mae_files:
    structure = StructureReader.read(os.path.join(all_c_mae_dir, mae_file))
    mae_file_index = mae_file.split("_")[0]
    atomic_number_dict[mae_file_index] = {}
    for atom in structure.atom:
        atomic_number_dict[mae_file_index][atom.index] = atom.property.get(
            "i_m_atomic_number"
        )


# Property 2: Hydrogen Neighbor
# hydrogen neighbor calculated by rdkit atom.GetTotalNumHs(includeNeighbors=True) function
hydrogen_neighbor_dict = {}
for mae_file in all_c_mae_files:
    structure = StructureReader.read(os.path.join(all_c_mae_dir, mae_file))
    mae_file_index = mae_file.split("_")[0]
    rdkit_mol = adapter.to_rdkit(structure)
    hydrogen_neighbor_dict[mae_file_index] = {}
    for rdkit_atom in rdkit_mol.GetAtoms():
        hydrogen_neighbor_dict[mae_file_index][rdkit_atom.GetIdx() + 1] = (
            rdkit_atom.GetTotalNumHs(includeNeighbors=True)
        )


# Property 3: Is aromatic
# aromatic information calculated by rdkit atom.GetIsAromatic() function
is_aromatic_dict = {}
for mae_file in all_c_mae_files:
    structure = StructureReader.read(os.path.join(all_c_mae_dir, mae_file))
    mae_file_index = mae_file.split("_")[0]
    rdkit_mol = adapter.to_rdkit(structure)
    is_aromatic_dict[mae_file_index] = {}
    for rdkit_atom in rdkit_mol.GetAtoms():
        if rdkit_atom.GetIsAromatic():
            is_aromatic_dict[mae_file_index][rdkit_atom.GetIdx() + 1] = 1
        else:
            is_aromatic_dict[mae_file_index][rdkit_atom.GetIdx() + 1] = 0


# Property 4: Number of heavy atoms
# number of heavy atoms
heavy_atoms_dict = {}
for mae_file in all_c_mae_files:
    structure = StructureReader.read(os.path.join(all_c_mae_dir, mae_file))
    structure.deleteAtoms([atom for atom in structure.atom if atom.element == "H"])
    mae_file_index = mae_file.split("_")[0]
    heavy_atoms_dict[mae_file_index] = len(structure.atom)


# Property 5: Double bonded
# double bonded information calculated by rdkit atom.GetBondWithAtomIdx(idx) function
double_bonded_dict = {}
for mae_file in all_c_mae_files:
    structure = StructureReader.read(os.path.join(all_c_mae_dir, mae_file))
    mae_file_index = mae_file.split("_")[0]
    rdkit_mol = adapter.to_rdkit(structure)
    double_bonded_dict[mae_file_index] = {}
    for rdkit_atom in rdkit_mol.GetAtoms():
        if any(
            bond.GetBondType() == Chem.BondType.DOUBLE for bond in rdkit_atom.GetBonds()
        ):
            double_bonded_dict[mae_file_index][rdkit_atom.GetIdx() + 1] = 1
        else:
            double_bonded_dict[mae_file_index][rdkit_atom.GetIdx() + 1] = 0


# Target: SOM status
molecule_title_dict = {}
primary_som_dict = {}
secondary_som_dict = {}
tertiary_som_dict = {}
som_dict = {}
som_level_dict = {}
for mae_file in all_c_mae_files:
    structure = StructureReader.read(os.path.join(all_c_mae_dir, mae_file))
    mae_file_index = mae_file.split("_")[0]
    molecule_title_dict[mae_file_index] = structure.property.get("s_m_title", "N/A")
    primary_som_dict[mae_file_index] = {}
    secondary_som_dict[mae_file_index] = {}
    tertiary_som_dict[mae_file_index] = {}
    som_dict[mae_file_index] = {}
    som_level_dict[mae_file_index] = {}
    primary_som_all = get_primary_som(structure)
    secondary_som_all = get_secondary_som(structure)
    tertiary_som_all = get_tertiary_som(structure)
    som_all = primary_som_all + secondary_som_all + tertiary_som_all
    for atom in structure.atom:
        primary_som_dict[mae_file_index][atom.index] = (
            1 if atom.index in primary_som_all else 0
        )
        secondary_som_dict[mae_file_index][atom.index] = (
            1 if atom.index in secondary_som_all else 0
        )
        tertiary_som_dict[mae_file_index][atom.index] = (
            1 if atom.index in tertiary_som_all else 0
        )
        som_dict[mae_file_index][atom.index] = 1 if atom.index in som_all else 0
        som_level_dict[mae_file_index][atom.index] = (
            3
            if atom.index in primary_som_all
            else (
                2
                if atom.index in secondary_som_all
                else 1 if atom.index in tertiary_som_all else 0
            )
        )


# Create the dataset csv file
dataset = pd.DataFrame(
    columns=[
        "zaretzki_index",
        "molecule_title",
        "atomic_index",
        "zaretzki_atomic_index",
        "atomic_number",
        "hydrogen_neighbor",
        "is_aromatic",
        "heavy_atoms",
        "double_bonded",
        "ir",
        "relative_ir",
        "bde",
        "sasa_hydrogen_maestro",
        "primary_som",
        "secondary_som",
        "tertiary_som",
        "som",
        "som_level",
    ]
)

for zaretzki_index, atomic_indices in atomic_number_dict.items():
    for atomic_index, atomic_number_value in atomic_indices.items():
        atomic_index = int(atomic_index)
        zaretzki_atomic_index = int(zaretzki_index) * 1000 + atomic_index

        # Check if zaretzki_index exists in reactivity_dict
        bde_value = bde_dict.get(zaretzki_index, {}).get(atomic_index, pd.NA)
        ir_value = ir_dict.get(zaretzki_index, {}).get(atomic_index, pd.NA)
        relative_ir_value = relative_ir_dict.get(zaretzki_index, {}).get(
            atomic_index, pd.NA
        )
        sasa_hydrogen_maestro_value = sasa_hydrogen_maestro_dict.get(
            zaretzki_index, {}
        ).get(atomic_index, pd.NA)

        dataset = dataset.append(
            {
                "zaretzki_index": zaretzki_index,
                "molecule_title": molecule_title_dict[zaretzki_index],
                "atomic_index": atomic_index,
                "zaretzki_atomic_index": zaretzki_atomic_index,
                "atomic_number": atomic_number_dict[zaretzki_index][atomic_index],
                "hydrogen_neighbor": hydrogen_neighbor_dict[zaretzki_index][
                    atomic_index
                ],
                "is_aromatic": is_aromatic_dict[zaretzki_index][atomic_index],
                "heavy_atoms": heavy_atoms_dict[zaretzki_index],
                "double_bonded": double_bonded_dict[zaretzki_index][atomic_index],
                "ir": ir_value,
                "relative_ir": relative_ir_value,
                "bde": bde_value,
                "sasa_hydrogen_maestro": sasa_hydrogen_maestro_value,
                "primary_som": primary_som_dict[zaretzki_index][atomic_index],
                "secondary_som": secondary_som_dict[zaretzki_index][atomic_index],
                "tertiary_som": tertiary_som_dict[zaretzki_index][atomic_index],
                "som": som_dict[zaretzki_index][atomic_index],
                "som_level": som_level_dict[zaretzki_index][atomic_index],
            },
            ignore_index=True,
        )

raw_zaretzki_index = dataset["zaretzki_index"].unique()

dataset.to_csv("../data/dataset/dataset_merck_all.csv", index=False)


dataset = dataset[dataset["atomic_number"] == 6]
dataset = dataset[dataset["is_aromatic"] == 0]
dataset = dataset[dataset["hydrogen_neighbor"] > 0]
dataset = dataset[dataset["double_bonded"] == 0]

# Export the dataset to csv
dataset.to_csv("../data/dataset/dataset_merck_all.csv", index=False)

# Print the ratio of SOM and non-SOM
primary_som_count = dataset["primary_som"].sum()
secondary_som_count = dataset["secondary_som"].sum()
tertiary_som_count = dataset["tertiary_som"].sum()
som_count = dataset["som"].sum()

print(f"primary som ratio: {primary_som_count / len(dataset)}")
print(f"secondary som ratio: {secondary_som_count / len(dataset)}")
print(f"tertiary som ratio: {tertiary_som_count / len(dataset)}")
print(f"som ratio: {som_count / len(dataset)}")

print(f"number of molecules: {len(dataset['zaretzki_index'].unique())}")
print(f"number of sites: {len(dataset)}")

clean_zaretzki_index = dataset["zaretzki_index"].unique()
# print the molecule in the raw_zaretzki_index but not in the clean_zaretzki_index
for molecule in raw_zaretzki_index:
    if molecule not in clean_zaretzki_index:
        print(molecule)
