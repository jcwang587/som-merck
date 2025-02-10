# Script to combine BDE and SASA data

import os
import shutil
from schrodinger.structure import StructureReader


def combine_bde_mae_sasa(bde_file, sasa_file, output_file):
    structure_bde = StructureReader.read(bde_file)
    structure_sasa = StructureReader.read(sasa_file)

    bde_dict = {}
    # Get the SASA property r_user_sasa
    for atom in structure_sasa.atom:
        if "r_user_sasa" in atom.property:
            bde_dict[atom.index] = atom.property["r_user_sasa"]

    # Add the SASA property r_user_sasa to the BDE structure
    for atom in structure_bde.atom:
        if atom.index in bde_dict:
            atom.property["r_user_sasa"] = bde_dict[atom.index]

    # Export the combined structure
    structure_bde.write(output_file)



bde_dir = "./bde_mae"
bde_files = os.listdir(bde_dir)

sasa_dir = "./sasa"
sasa_files = os.listdir(sasa_dir)

output_dir = "./output"
shutil.rmtree(output_dir, ignore_errors=True)
os.makedirs(output_dir, exist_ok=True)

for bde_file, sasa_file in zip(bde_files, sasa_files):
    bde_path = f"{bde_dir}/{bde_file}"
    sasa_path = f"{sasa_dir}/{sasa_file}"

    bde_file = os.path.basename(bde_path)
    sasa_file = os.path.basename(sasa_path)

    structure = StructureReader.read(bde_path)
    output_file = f"{output_dir}/{structure.title}_bde_sasa.mae"
    combine_bde_mae_sasa(bde_path, sasa_path, output_file)
    print(f"Combined {bde_file} and {sasa_file} into {output_file}")
