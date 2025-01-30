# Script to combine BDE and SASA data

from schrodinger.structure import StructureReader


def combine_bde_sasa(bde_file, sasa_file):
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

    # Get the name of th

    # Export the combined structure
    structure_bde.write(f"{structure_bde.title}_bde_sasa.mae")


bde_file = "./bde_sasa/015_CARBAMAZEPINE.mae"
sasa_file = "./bde_sasa/structure_015_sasa.mae"

combine_bde_sasa(bde_file, sasa_file)