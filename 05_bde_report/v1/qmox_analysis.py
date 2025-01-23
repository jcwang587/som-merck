from pathlib import Path
import os, warnings

warnings.filterwarnings("ignore")


from schrodinger.structure import StructureReader
from schrodinger import adapter

from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Image
from reportlab.platypus import Spacer, Table, TableStyle
from reportlab.lib import colors

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from PIL import Image as ImagePIL
from PIL import ImageChops
import cairosvg

from rdkit import Chem
from rdkit.Chem import rdDepictor

rdDepictor.SetPreferCoordGen(True)
from rdkit.Chem.Draw import rdMolDraw2D

from utils import draw_molecule


def create_image(st, out_dir=Path().resolve(), name="labeled"):
    """
    Create image function compatible with the Schrodinger python
    """

    draw_molecule(st, name)


    out_png = out_dir / f"{name}.png"


    return out_png


def check_missing_sites(structure: StructureReader):
    """Check if the structure has missing sites"""
    for atom in structure.atom:
        if atom.property.get("r_user_CH-BDE") is None:
            print(f"Missing site: {atom.element} {atom.index}")


def autocrop(image: str, bgcolor: str = "white"):
    im = ImagePIL.open(image)
    if im.mode != "RGB":
        im = im.convert("RGB")
    bg = ImagePIL.new("RGB", im.size, bgcolor)
    diff = ImageChops.difference(im, bg)
    bbox = diff.getbbox()
    if bbox:
        image_crop = im.crop(bbox)
        image_crop.save(image)
        width, length = image_crop.size
        return width, length


def create_risk_scale_png(
    medium: float, high: float, filename: Path = Path().resolve() / "risk_scale.png"
):
    """Create risk figure based on medium and high cut-offs"""

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rect1 = Rectangle((50, 50), 50, 40, color="yellow")

    if high > medium:
        rect = Rectangle((0, 50), 50, 40, color="green")
        rect2 = Rectangle((100, 50), 50, 40, color="red")
        ax.add_patch(rect)
        ax.add_patch(rect1)
        ax.add_patch(rect2)
        ax.annotate(str(medium), xy=(44, 94), xytext=(44, 94), fontsize=14)
        ax.annotate(str(high), xy=(94, 94), xytext=(94, 94), fontsize=14)
        ax.annotate("Low", xy=(18, 38), xytext=(18, 38), fontsize=14)
        ax.annotate("Medium", xy=(58, 38), xytext=(58, 38), fontsize=14)
        ax.annotate("High", xy=(115, 38), xytext=(115, 38), fontsize=14)

    else:
        rect = Rectangle((100, 50), 50, 40, color="green")
        rect2 = Rectangle((0, 50), 50, 40, color="red")
        ax.add_patch(rect2)
        ax.add_patch(rect1)
        ax.add_patch(rect)
        ax.annotate(str(high), xy=(44, 94), xytext=(44, 94), fontsize=14)
        ax.annotate(str(medium), xy=(94, 94), xytext=(94, 94), fontsize=14)
        ax.annotate("High", xy=(18, 38), xytext=(18, 38), fontsize=14)
        ax.annotate("Medium", xy=(58, 38), xytext=(58, 38), fontsize=14)
        ax.annotate("Low", xy=(115, 38), xytext=(115, 38), fontsize=14)

    ax.annotate("kcal/mol", xy=(-40, 94), xytext=(-40, 94), fontsize=14)
    plt.xlim(-40, 175)
    plt.ylim(-10, 175)
    plt.axis("off")

    plt.savefig(str(filename), dpi=300)

    autocrop(str(filename))


class CoxidesAnalysis:

    def __init__(
        self,
        dir_report: Path,
        path_structure: Path,
        title: str,
        medium_coff: float,
        high_coff: float,
    ):
        self.dir_report = dir_report
        self.path_structure = path_structure
        self.title = title
        self.medium_coff = medium_coff
        self.high_coff = high_coff
        self.structure = StructureReader.read(self.path_structure)
        self.dft = "B3LYP"
        self.basis = "6-31G(d,p)"

        self.risk_scale = self.dir_report / "C-oxidation_risk_scale.png"
        self.img = create_image(
            self.structure, name=f"{self.dir_report}/{self.title}"
        )

        # Give a fake data
        self.data = self.GenerateDataList()

    def GenerateDataList(self):
        """
        Take information from self.structure to create a Data list for report
        generation
        """

        print("Creating Reports...")
        print("Generating BDE table data...")
        data = [["Atom", "BDE (kcal/mol)", "BDE Risk"]]

        for atom in self.structure.atom:
            try:
                bde = float(atom.property["r_user_CH-BDE"])
                if bde < self.high_coff:
                    propensity = "High"
                if bde <= self.medium_coff and bde >= self.high_coff:
                    propensity = "Moderate"
                if bde > self.medium_coff:
                    propensity = "Low"
                row = [
                    str(atom.element) + str(atom.index),
                    str(atom.property["r_user_CH-BDE"]),
                    propensity,
                ]
                data.append(row)
            except KeyError:
                pass
        
        # Check if SASA property is available
        if "r_user_CH-SASA" in self.structure.atom[1].property:
            data = [["Atom", "BDE (kcal/mol)", "BDE Risk", "SASA (Å²)"]]
            # add SASA for each row in data
            for atom in self.structure.atom:
                try:
                    bde = float(atom.property["r_user_CH-BDE"])
                    if bde < self.high_coff:
                        propensity = "High"
                    if bde <= self.medium_coff and bde >= self.high_coff:
                        propensity = "Moderate"
                    if bde > self.medium_coff:
                        propensity = "Low"
                    row = [
                        str(atom.element) + str(atom.index),
                        str(atom.property["r_user_CH-BDE"]),
                        propensity,
                        str(atom.property["r_user_CH-SASA"]),
                    ]
                    data.append(row)
                except KeyError:
                    pass

        return data

    def CreatePDFReport(self):
        """Create pdf report from a data dictionary and structural image"""
        # Define style

        style = TableStyle(
            [
                ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                ("FONTNAME", (0, 0), (-1, 0), "Courier-Bold"),
                ("FONTSIZE", (0, 0), (-1, 0), 10),
                ("BACKGROUND", (0, 1), (-1, -1), colors.beige),
                ("BOX", (0, 0), (-1, -1), 2, colors.black),
                ("GRID", (0, 1), (-1, -1), 2, colors.black),
            ]
        )
        # Get other styles

        styles = getSampleStyleSheet()
        styleNormal = styles["Normal"]
        styleHeading = styles["Heading1"]
        styleHeading2 = styles["Heading2"]
        styleHeading.alignment = 1

        # Process data
        pil_st = Image(str(self.img), width=4.0 * inch, height=3.6 * inch)
        pil_risk = Image(str(self.risk_scale), width=3.0 * inch, height=0.8 * inch)

        # Create PDF
        story = []

        # Append HEAD of document
        story.append(
            Paragraph("C-oxidation BDE Energy Report for: " + self.title, styleHeading)
        )
        story.append(Spacer(inch, 0.25 * inch))
        story.append(
            Paragraph(
                "This report covers the results for BDE calculations"
                f" performed for: {self.title}. Oxidation propensity is established using "
                "C-H Bond Dissociation Enthalpies (BDE). The lower the C-H BDE values the "
                "higher the propensity for C-oxidation. Details for the DFT calculations and"
                " overall workflow are explained at the end of this document",
                styleNormal,
            )
        )
        story.append(Spacer(inch, 0.05 * inch))

        # Append molecule picture

        story.append(Paragraph("Bond Dissociation Energies (kcal/mol)", styleHeading2))
        story.append(pil_st)
        story.append(Spacer(inch, 0.25 * inch))

        # Append BDE table
        main_table = Table(self.data, colWidths=85, rowHeights=18)
        main_table.setStyle(style)
        story.append(main_table)
        story.append(Spacer(inch, 0.15 * inch))

        # Risk scale
        story.append(Paragraph("Risk Scale:", styleNormal))
        story.append(pil_risk)
        story.append(Spacer(inch, 0.25 * inch))

        # Append final details
        story.append(Paragraph("Calculation details and output files", styleHeading2))
        story.append(Spacer(inch, 0.15 * inch))

        # story.append(Paragraph(f'Report files are available in: {str(self.dir_report)}',
        # styleNormal))
        story.append(Spacer(inch, 0.1 * inch))
        story.append(
            Paragraph(
                "Conformational search calculations were performed only"
                " for the base ground state molecule. The lowest energy conformer was selected"
                " to generate radicals and run optimization DFT calculations. DFT calculations"
                f" were performed using Gaussian with {self.dft} level of theory and {self.basis}"
                " basis set. The BDE protocol was adapted from: <i>Lienard, P., Gavartin, J.,"
                " Boccardi, G., & Meunier, M. (2015). Predicting drug substances autoxidation."
                " Pharmaceutical research, 32, 300-310.</i>",
                styleNormal,
            )
        )

        # Save PDF
        pdf_file = self.dir_report / f"{self.title}_CH-BDE_report.pdf"
        doc = SimpleDocTemplate(
            str(pdf_file),
            pagesize=letter,
            title=f"{self.title} BDE Report",
            author="BDE 2.0",
        )
        doc.build(story)

        # Remove the image file
        os.remove(self.img)

        return pdf_file


if __name__ == "__main__":

    C_HIGH_COFF = 88
    C_MEDIUM_COFF = 94

    mae_dir = Path("./test")
    mae_files = mae_dir.glob("*.mae")

    # Generate the risk scale image
    create_risk_scale_png(
        medium=C_MEDIUM_COFF,
        high=C_HIGH_COFF,
        filename=mae_dir / "C-oxidation_risk_scale.png",
    )

    # Generate the pdf report using the CoxidesAnalysis
    for mae_file in mae_files:
        mae_structure = StructureReader.read(mae_file)

        qmox_analysis = CoxidesAnalysis(
            dir_report=mae_dir,
            path_structure=Path(mae_file),
            title=mae_structure.title,
            medium_coff=C_MEDIUM_COFF,
            high_coff=C_HIGH_COFF,
        )
        qmox_analysis.CreatePDFReport()

    # Remove the risk scale image
    os.remove(mae_dir / "C-oxidation_risk_scale.png")
