from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import cairosvg
import shutil

from schrodinger.structure import StructureReader
from schrodinger.rdkit.rdkit_adapter import to_rdkit

import sys
import glob
import os

sys.path.append("../src")
from utils import renew_folder

renew_folder("../data/svg_merck_merck")
renew_folder("../data/png_merck_merck")

mae_folder = "../data/merck_merck"
mae_files = glob.glob(os.path.join(mae_folder, "*.mae"))

# Iterate over each molecule in the SDF file
for mae_file in mae_files:
    file_name = os.path.basename(mae_file)
    structure = StructureReader.read(mae_file)
    mol = to_rdkit(structure)

    # Remove hydrogens
    mol = Chem.RemoveHs(mol)

    # Compute 2D coordinates
    AllChem.Compute2DCoords(mol)

    # Initialize lists for SOM atom indices
    primary_som_atoms = []
    secondary_som_atoms = []
    tertiary_som_atoms = []

    # Get PRIMARY_SOM atom indices if the value is a integer
    if mol.HasProp("i_sd_PRIMARY\_SOM"):
        primary_som = mol.GetProp("i_sd_PRIMARY\_SOM")
        primary_som_atoms = [int(i) - 1 for i in primary_som.strip().split()]

    # Get PRIMARY_SOM atom indices if the value is a string
    if mol.HasProp("s_sd_PRIMARY\_SOM"):
        primary_som = mol.GetProp("s_sd_PRIMARY\_SOM")
        primary_som_atoms = [int(i) - 1 for i in primary_som.strip().split()]

    # Combine all SOM atom indices
    highlight_atoms = primary_som_atoms + secondary_som_atoms + tertiary_som_atoms

    # Define colors for each SOM type
    primary_color = (0.7, 0, 0)
    secondary_color = (0, 0.7, 0)
    tertiary_color = (0, 0.7, 0.7)  # cyan

    # Create a color dictionary mapping atom indices to colors
    highlight_colors = {}
    for idx in primary_som_atoms:
        highlight_colors[idx] = primary_color
    for idx in secondary_som_atoms:
        highlight_colors[idx] = secondary_color
    for idx in tertiary_som_atoms:
        highlight_colors[idx] = tertiary_color

    # Draw the molecule with highlighted atoms without bonds
    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(500, 300)
    opts = drawer.drawOptions()
    for i in range(mol.GetNumAtoms()):
        opts.atomLabels[i] = f"{mol.GetAtomWithIdx(i).GetSymbol()}{i+1}"
    drawer.DrawMolecule(
        mol,
        highlightAtoms=highlight_atoms,
        highlightAtomColors=highlight_colors,
        highlightBonds=[],
    )
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    # Save the SVG image
    with open(f"../data/svg_merck_merck/{file_name.replace('.mae', '.svg')}", "w") as f:
        f.write(svg)

    cairosvg.svg2png(
        url=f"../data/svg_merck_merck/{file_name.replace('.mae', '.svg')}",
        write_to=f"../data/png_merck_merck/{file_name.replace('.mae', '.png')}",
    )


# remove the svg folder
shutil.rmtree("../data/svg_merck_merck")
