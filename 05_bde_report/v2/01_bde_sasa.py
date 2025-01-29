# Script to combine BDE and SASA data

from schrodinger.structure import StructureReader


def combine_bde_sasa(bde_file, sasa_file):
    structure_bde = StructureReader(bde_file)
    structure_sasa = StructureReader(sasa_file)



bde_file = "05_bde_report/v2/bde_data.sdf"
sasa_file = "05_bde_report/v2/sasa_data.sdf"

combine_bde_sasa(bde_file, sasa_file)