from pathlib import Path
import os, subprocess, argparse, sys, getpass, re, shutil, warnings, numpy, json

warnings.filterwarnings('ignore')

from schrodinger.structutils import analyze, rmsd

from contextlib import contextmanager

from schrodinger.utils import log, cmdline
from schrodinger.job import launchapi, jobcontrol
from schrodinger.structure import StructureReader
from schrodinger import structure, adapter

from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Image 
from reportlab.platypus import Spacer, Table, TableStyle 
from reportlab.lib import colors

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

from PIL import Image as ImagePIL
from PIL import ImageChops

from rdkit import Chem
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)
from rdkit.Chem.Draw import rdMolDraw2D


logger = log.get_output_logger("qmox_analysis.py")
log.default_logging_config()


O3_ENE = -75.22035293792803
H_RAD_ENE = -0.49742380015
C_HIGH_COFF = 85
C_MEDIUM_COFF = 96
N_HIGH_CUTOFF = 180
N_MEDIUM_CUTOFF = 203
S_HIGH_CUTOFF = 75
S_MEDIUM_CUTOFF = 69.5
COFF_DIC = {
        'c_medium_coff': C_MEDIUM_COFF,
        'c_high_coff': C_HIGH_COFF,
        'n_medium_coff': N_MEDIUM_CUTOFF,
        'n_high_coff': N_HIGH_CUTOFF,
        's_medium_coff': S_MEDIUM_CUTOFF,
        's_high_coff': S_HIGH_CUTOFF
    }
DEFAULT_JOBNAME = 'qmoxanalysis'
PROGRAM_NAME = 'qmox_analysis'

NMR_ISOTROPIC_PROP = 'r_user_Gaussian_NMR_Isotropic'
CHARGE_PROP = 'r_user_Gaussian_Charge'
DIPOLE_MOMENT_PROP_BASE = 'r_user_Gaussian_Dipole_Moment_'
BASIS_PROP = 's_user_Gaussian_Basis'
ENERGY_PROP = 'r_user_Gaussian_Energy'
THERMAL_ELECTRONIC_SUM_PROP_LEGACY = 'r_user_Gaussian_Sum_of_Thermal_and_Electronic_Enthalpies'   # DK: 2020-06-05: added Matt's fix
THERMAL_ELECTRONIC_SUM_PROP = 'r_user_Gaussian_Sum_of_Electronic_And_Thermal_Free_Energies_(kJ/mol)' # DK: 2020-09-15: corrected label base on Nina's input
METHOD_PROP = 's_user_Gaussian_QM_Method'
SOLVENT_PROP = 's_user_Gaussian_QM_Solvent'
EXIT_STATUS_PROP = 's_user_Gaussian_Exit_Status'
IMAGINARY_FREQ_PROP = 'b_user_Gaussian_Imaginary_Frequencies'

HARTREE_TO_KJMOL = 2625.5  # 1 AU = 2625.50 kJ/mol

dft_scripts = Path('/home/valdivig/Documents/Oxidation_pred/Maestro_script/qmox_analysis')

MWFN_SH = Path(__file__).resolve().parent / 'Multiwfn'
ALIE_SH = Path(__file__).resolve().parent / 'alie_script.sh'
INKSCAPE_DIR = Path(__file__).resolve().parent / 'Inkscape-9c6d41e-x86_64.AppImage'
IE_CMDS =  Path(__file__).resolve().parent / 'ie_input_cmds.txt'
EA_CMDS = Path(__file__).resolve().parent / 'ea_input_cmds.txt'
AC_CMDS = Path(__file__).resolve().parent / 'ac_input_cmds.txt'
# INKSCAPE_DIR = '/apps/prod/DEVELOPMENT/MAS/Formulation/Inkscape-9c6d41e-x86_64.AppImage'
# MWFN_SH = dft_scripts / 'process_gaussian_fchk_1.2.sh'
# CHK2FCHK = dft_scripts / 'formchk_batch.sh'
# '/home/valdivig/Documents/Oxidation_pred/Maestro_script/qmox_analysis/process_gaussian_fchk_1.2.sh'

def check_normal_termination(last_line:str):
    if 'Normal termination' not in last_line:
        return False
    return True
     
def send_email(user_email, title:str, message:str, attachs:list=[]):
    '''Send email to user with title, message, and attachs (list of paths)'''
    
    attach_str = ''
    for attach in attachs:
        attach_str += f' -a {str(attach)}'

    mail_cmd = f'mail -s "{title}"{attach_str} {user_email}<<< "{message}"'
    mail_me = mail_cmd.replace(user_email, 'gabriel.valdivia-berroeta@boehringer-ingelheim.com')
    print(mail_cmd)

    os.system(mail_cmd)
    os.system(mail_me)

def autocrop(image:str, bgcolor:str='white'):
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

def create_risk_scale_png(medium:float, high:float, 
                        filename:Path=Path().resolve() / 'risk_scale.png'):
    '''Create risk figure based on medium and high cut-offs'''

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rect1 = Rectangle((50, 50),50, 40, color='yellow')

    
    if high > medium:
        rect = Rectangle((0, 50),50, 40, color='green')
        rect2 = Rectangle((100, 50),50, 40, color='red')
        ax.add_patch(rect)
        ax.add_patch(rect1)
        ax.add_patch(rect2)
        ax.annotate(str(medium), xy= (44,94),xytext=(44,94), fontsize=14)
        ax.annotate(str(high), xy= (94,94),xytext=(94,94), fontsize=14)
        ax.annotate('Low', xy= (18,38),xytext=(18,38), fontsize=14)
        ax.annotate('Medium', xy= (58,38),xytext=(58,38), fontsize=14)
        ax.annotate('High', xy= (115,38),xytext=(115,38), fontsize=14)

    else:
        rect = Rectangle((100, 50),50, 40, color='green')
        rect2 = Rectangle((0, 50),50, 40, color='red')
        ax.add_patch(rect2)
        ax.add_patch(rect1)
        ax.add_patch(rect)
        ax.annotate(str(high), xy= (44,94),xytext=(44,94), fontsize=14)
        ax.annotate(str(medium), xy= (94,94),xytext=(94,94), fontsize=14)
        ax.annotate('High', xy= (18,38),xytext=(18,38), fontsize=14)
        ax.annotate('Medium', xy= (58,38),xytext=(58,38), fontsize=14)
        ax.annotate('Low', xy= (115,38),xytext=(115,38), fontsize=14)


    ax.annotate('kcal/mol', xy= (-40,94),xytext=(-40,94), fontsize=14)
    plt.xlim(-40,175)
    plt.ylim(-10,175)
    plt.axis('off')

    plt.savefig(str(filename), dpi=300)

    autocrop(str(filename))

def table_to_csv(table:list, name:str):
    '''Convert a table to a csv file'''

    fpath = f'{name}.csv'
    df = pd.DataFrame(table[1:], columns=table[0])
    df.to_csv(fpath, index=False)

    return fpath

def create_image(st, mode='C', 
                out_dir=Path().resolve(), name='labeled',
                inkscape_path:Path=None):
    '''
    Create image function compatible with the Schrodinger python
    Modes: C, N, and S
    SVG format
    '''

    
    molecule = adapter.to_rdkit(st)
     
    for i, atom in enumerate(molecule.GetAtoms()):
        if mode == atom.GetSymbol():
            atom.SetProp('_displayLabel', f'{mode}<sub>{str(i+1)}</sub>')
    
    print('Deleting Hydrogens 2')
    molecule = Chem.RemoveHs(molecule)

    rdDepictor.Compute2DCoords(molecule)
    d = rdMolDraw2D.MolDraw2DSVG(500,500)
    rdMolDraw2D.PrepareAndDrawMolecule(d, molecule)
    d.DrawMolecule(molecule)
    d.FinishDrawing()
    svg_txt = d.GetDrawingText()

    out_svg = out_dir / f'{name}.svg'
    out_png = out_dir / f'{name}.png'
    with open(out_svg, 'w+') as f:
        f.write(svg_txt)
    
    # Convert to png
    # os.system(f'{INKSCAPE_DIR} inkscape {str(out_svg)} ' 
    # f'-o {str(out_png)}')
    os.system(f'convert -density 1200 -resize 500x500 {str(out_svg)} {str(out_png)}')


    return out_png

@contextmanager
def change_working_directory(folder):
    """
    A context manager to temporarily change the working directory to folder
    :param folder: the folder that becomes the working directory
    :type folder: str
    """
    old_folder = os.getcwd()
    os.chdir(folder)
    yield
    os.chdir(old_folder)


def parse_gaussian_log(fname, use_energysum=False):
    """
    Parse a Gaussian log file and return a dictionary of properties loaded
    from the file.

    :param fname: Filename of log to parse
    :type fname: str

    :return: Tuple of dictionary of structure and atom
             property names to property values.
    :rtype: tuple(dict, dict)
    """
    props = {}
    atom_props = {}
    with open(fname, 'r') as fh:
        gauss_log = fh.readlines()
    props[EXIT_STATUS_PROP] = gauss_log[-1]
    for line_idx, line in enumerate(gauss_log):
        if ' Will use up to' in line:
            for idx_2 in range(line_idx, line_idx + 10):
                line_2 = gauss_log[idx_2]
                if line_2[:2] == ' #':
                    fields = line_2.split()
                    method_idx = 2 if fields[0] == '#' else 1
                    qm_method = fields[method_idx].replace('#', '')
                    if '/' in qm_method:
                        qm_method = qm_method.split('/')[0]
                    props[METHOD_PROP] = qm_method
        elif line[:9] == ' SCF Done' and not use_energysum:
        # DK: 2020-09-15: Fix based on Nina's input
            props[ENERGY_PROP] = float(line.split()[4])                
        # DK: 2020-06-05: Added Matt's fix here as additional block
        elif "Sum of electronic and thermal Free Energies" in line and use_energysum:
            props[THERMAL_ELECTRONIC_SUM_PROP] = float(
                line.split()[-1]) * HARTREE_TO_KJMOL
            props[ENERGY_PROP] = float(line.split()[-1])
        elif line[:15] == ' Standard basis':
            basis = line.split('basis: ')[-1].split()[0]
            props[BASIS_PROP] = basis
        elif line[:14] == ' Dipole moment':
            dipole_line = gauss_log[line_idx + 1].split()
            base_prop = DIPOLE_MOMENT_PROP_BASE
            props[base_prop + 'x'] = float(dipole_line[1])
            props[base_prop + 'y'] = float(dipole_line[3])
            props[base_prop + 'z'] = float(dipole_line[5])
            props[base_prop + 'Total'] = float(dipole_line[7])
        elif line[:8] == ' Charge=':
            charge = float(line.split()[1])
            props[CHARGE_PROP] = charge
        elif 'imaginary frequencies' in line:
            props[IMAGINARY_FREQ_PROP] = True
        elif 'Isotropic' in line and 'Anisotropy' in line:
            # NMR shielding constants
            nmr_info = line.split()
            anum = int(nmr_info[0])
            shielding_const = nmr_info[4]
            nmr_prop = (NMR_ISOTROPIC_PROP, float(shielding_const))
            if anum in atom_props:
                atom_props[anum].append(nmr_prop)
            else:
                atom_props[anum] = [nmr_prop]
    return props, atom_props

def find_energy(gaussian_output):
    '''Find the energy of the optimized structure in an ouput file. If the optimization was not completed it returns None'''
    
    props, _ = parse_gaussian_log(gaussian_output)
    
    return props['r_user_Gaussian_Energy']

def calculate_BDE(e_neutral, e_radical, e_h_radical=-0.49742380015):
    '''Calculate BDE form energies in Hartrees. It returns the result in kcal/mol'''
    bde = 627.5*e_radical + 627.5*e_h_radical - 627.5*e_neutral
    return round(bde,2)

def getBoltzmannAverage(values, energies):
    """
    Return the Boltzmann averaged value of a list of values each given a
    certain energy.

    :param values: List of values to average
    :type values: list(float)

    :param energies: List of energies corresponding to each value
    :type energies: list(float)

    :return: Boltzmann averaged value
    :rtype: float
    """
    T = 298  # Kelvin
    R = .0083145  # kJ/molK
    val_arr = numpy.asarray(values)
    min_e = min(energies)
    e_diffs = [e - min_e for e in energies]
    e_arr = numpy.asarray(e_diffs)
    bfs = numpy.exp(-1 * e_arr / (R * T))
    pfxn = sum(bfs)
    probs = bfs / pfxn
    weighted_val = numpy.dot(val_arr, probs)
    return weighted_val

class ProcessMwfn:
    def __init__(self, out_folder:Path, area_coff:float=0) -> None:
        self.out_folder = out_folder
        self.area_coff = area_coff
        self.property_list = ['a_charge','ea_min', 'ea_max', 'area', 
        'ie_min', 'ie_max', 'ie_avg']
        self.column_names = ['index', 'atom'] + self.property_list
        self.conformer_csvs = []
        
        

    def create_initial_csv(self, csv_name:str, dfs:list):
        '''Create csv'''

        for i in range(len(dfs)-1):
            if i == 0:
                df = dfs[i].merge(dfs[i+1], on='index')
            else:
                df = df.merge(dfs[i+1], on='index')

        df_final = df.drop('index', axis=1)
        

        # df_1 = df_ac.merge(df_ea, on='index')       
        # df_2 = df_1.merge(df_ie, on='index')
        # df_final = df_2.drop('index', axis=1)
        # df_final = df_final.dropna().reset_index(drop=True)
        

        path_csv = self.out_folder / f'{csv_name}.csv'
        self.conformer_csvs.append(path_csv)
        df_final.to_csv(path_csv, index=False)
        
        return df_final


    def find_conformer_energy(self, raw_file:Path):
        with open(raw_file, 'r+') as f:
            for line in f.readlines():
                if 'Total energy' in line:
                    return float(line.split()[2])
        
        return False

    def process_txt_files(self):
        self.output_d = {}
        #Find energies for every conformer
        for fchk_file in self.out_folder.glob('*.fchk'):
            self.sname = fchk_file.stem.split('-out_')[0]
            raw_file = self.out_folder / f'{fchk_file.stem}-ac_raw.txt'
            check = self.check_files(fchk_file.stem)
                
            if raw_file.exists() and fchk_file.exists() and check:
                self.output_d[fchk_file.name] = [self.create_initial_csv(fchk_file.stem, check), 
                                            self.find_conformer_energy(raw_file)]
            else:
                print('Mwfn results or gaussian output (.fchk) files were not found or were' 
                f' incorrectly processed for {fchk_file.name}')
        
        # Get sample name
                
    def check_files(self, fname:str):
        '''Check files for a conformer to see if it should be analyzed'''

        ie_txt = self.out_folder / f'{fname}-ie.txt'
        ac_txt = self.out_folder / f'{fname}-ac.txt'
        ea_txt = self.out_folder / f'{fname}-ea.txt'

        if not ac_txt.exists():
            print('Atomic charge info was not found. Please check the bash code.')
            return False    
        if not ie_txt.exists():
            print('Ionization energy info was not found and it will not be included in the summary file')
            return False
        if not ea_txt.exists():
            print('Atomic charge info was not found and it will not be included in the summary file')
            return False
        ac_col_names = ['index', 'atom', 'a_charge']
        ie_col_names = ['index', 'area', 'ie_min', 'ie_max', 'ie_avg']
        ea_col_names = ['index', 'ea_min', 'ea_max']
        
        dfs = list()
        df_ac = pd.read_csv(ac_txt, sep='\t', names=ac_col_names)
        df_ie = pd.read_csv(ie_txt, sep='\t', names=ie_col_names) 
        df_ea = pd.read_csv(ea_txt, sep='\t', names=ea_col_names)
        dfs.append(df_ac)
        dfs.append(df_ea)
        dfs.append(df_ie)
        
        # Check anomalities in df_ac:
        for col in ac_col_names:
            col_type = str(df_ac[col].dtypes)
            if col_type == 'object' and col != 'atom':
                print(f'Problems were found with atomic charge for: {fname}')
                return False
            
        # Check anomalities in df_ie
        for col in ie_col_names:
            col_type = str(df_ie[col].dtypes)
            if col_type == 'object':
                print(f'Problems were found with ionization energy for: {fname}')
                return False
        
        # Check anomalities in df_ea
        for col in ea_col_names:
            col_type = str(df_ea[col].dtypes)
            if col_type == 'object':
                print(f'Problems were found with electron affinity for: {fname}')
                return False
        
        return dfs
    def apply_area_filter(self, row):
        ''' Apply area filter, returns True if area passed, false otherwise'''
        if 'N' in row.atom:
            print(f'Applying filter of {self.area_coff} to {row.atom}')
            print(f'{row.atom} area: {row.area}')
            if float(row.area) < self.area_coff:
                print('Rejected')
                return False
        print('Accepted')
        return True

    def create_unified_dic(self):
        self.boltzmann_d = {}
        self.df_d = {'atom':[]}
        # print(self.output_d)
        for i, value in enumerate(self.output_d.values()):
            print(f'Create unified dic step: {i}')
            for n in range(len(value[0])):
                new_atom = None
                row = value[0].iloc[n]
                check_area = self.apply_area_filter(row)
                if not check_area:
                    print(f'Conformer skipped due to low accesible area: {i}')
                    continue
                keys = [key for key in self.boltzmann_d.keys()]
                if row.atom not in keys:
                    self.boltzmann_d[row.atom] = {}
                    new_atom = row.atom
                for prop in self.property_list:
                    self.df_d[prop] = [] 
                    if row.atom == new_atom:
                        self.boltzmann_d[row.atom][prop] = [row[prop]]
                    else:
                        self.boltzmann_d[row.atom][prop].append(row[prop])                        
                if row.atom == new_atom:
                    self.boltzmann_d[row.atom]['conf_e'] = [value[1]]
                else:
                    self.boltzmann_d[row.atom]['conf_e'].append(value[1])
    
    def create_boltzmann_avg_df(self):
        # print(self.boltzmann_d)
        for key, value in self.boltzmann_d.items():
            self.df_d['atom'].append(key)
            for prop in self.property_list:
                prop_values = value[prop]
                energies = value['conf_e']
                # print(prop_values)
                # print(energies)
                boltzmann_average = getBoltzmannAverage(prop_values, energies)
                self.df_d[prop].append(boltzmann_average)
        
        return pd.DataFrame(self.df_d)
    
    def create_summary_excel(self):
        '''Take final summary excel and conformer excels and put them in one document'''
        fname = self.sname if len(self.sname) < 21 else self.sname[:20] 
        summary_excel_path = self.out_folder / f'{fname}_summary.xlsx'   
        self.boltzmann_df.to_excel(summary_excel_path, index=False, 
                                    sheet_name=summary_excel_path.stem)

        with pd.ExcelWriter(summary_excel_path, mode='a', engine='openpyxl') as writer:
            for csv_file in self.conformer_csvs:
                df = pd.read_csv(csv_file)
                df.to_excel(writer, sheet_name=csv_file.stem, index=False)
        # Create a csv file using the boltzmann
        self.boltzmann_df.to_csv(summary_excel_path.with_suffix('.csv'), index=False) 
        

    def process_mwfn(self):
        '''Main function to execute the class'''
        
        self.process_txt_files()
        self.create_unified_dic()
        self.boltzmann_df = self.create_boltzmann_avg_df()
        self.create_summary_excel()


class CoxidesAnalysis:
    ''' Class to analyze jaguar job calculations'''
    def __init__(self, dir_i:Path, title:str, medium_coff:float, 
                high_coff:float, h_radical_ene:float=-0.49742380015,
                inkscape:Path=None, email:str='', dir_report:Path=Path().resolve):
          
        self.dir = dir_i
        # self.dir_report = self.create_report_dir()
        self.dir_report = dir_report
        
        self.title = title
        self.h_radical_ene = h_radical_ene
        self.username = getpass.getuser()

        self.indexes = list()
        self.gaussian_results = dict()
        self.files_to_send = list()
        
        self.medium_coff = medium_coff
        self.high_coff = high_coff
        self.inkscape = inkscape

        self.email = email

    def analyze(self):
        '''Main function to analyze data in self.dir'''
        
        # This function sets self.title, self.out_files and self.neutral_e
        # self.dft, self.basis
        result_1 = self.analyze_files()
        if not result_1:
            return
        
        result_2 = self.analyze_out_files()        
        if not result_2:
            return

        # self.attach_bde_atom_property()
        self.data = self.GenerateDataList()
        csv = table_to_csv(self.data, f'{self.title}_cox')
        self.files_to_send.append(Path(csv))
        self.img = create_image(self.structure, mode='C', inkscape_path=self.inkscape)
        self.risk_scale = self.dir_report / 'C-oxidation_risk.png'
        create_risk_scale_png(self.medium_coff, self.high_coff, self.risk_scale)
        self.files_to_send.append(self.img)
        print('Creating PDF report')
        pdf_file = self.CreatePDFReport()
        self.files_to_send.append(pdf_file)
        # print('Creating Excel Report')
        # self.CreateExcelReport()
        self.send_email(type='all-good')
        if self.failed:
            self.send_email(type='radical-failed') 
  
    
    def create_report_dir(self):
        '''Create report dir'''
        try:
            os.mkdir(self.dir / 'Reports')
        except:
            pass

        return self.dir / 'Reports'
    
    def analyze_files(self):
        '''Find base compound .mae file to perform analysis '''
        outs = [file for file in self.dir.glob('*.out')]
        if not outs:
            self.send_email(type='all-failed')
            return False
        
        #Find base compound based on name length
        min_length = 999999999 
        for out in outs:
            if len(out.stem) < min_length:
                min_length = len(out.stem)
                self.parent_out_file = out
        
        outs.remove(self.parent_out_file)
        self.out_files = outs
        
        base_props, _ = parse_gaussian_log(self.parent_out_file)
        self.neutral_e = find_H_tot(self.parent_out_file)
    
        # self.neutral_e = base_props['r_user_Gaussian_Energy']
        self.dft = base_props['s_user_Gaussian_QM_Method']
        self.basis = base_props['s_user_Gaussian_Basis']
        parent_status = base_props['s_user_Gaussian_Exit_Status']

        if not 'Normal' in parent_status:
            self.send_email(type='base-failed')
            return False

        self.structure = StructureReader.read(
            str(self.parent_out_file.with_suffix('.maegz')))

        return True

    def analyze_out_files(self):
        '''Analyze file to see if jobs have terminated correctly'''
        
        self.failed = []

        for out_file in self.out_files:
            props, _ = parse_gaussian_log(out_file)
            self.gaussian_results[out_file] = props
            if not 'Normal' in props['s_user_Gaussian_Exit_Status']:
                self.failed.append(out_file)
            else:
                self.attach_bde_atom_property(out_file)
        
        # Return false if all the radical calculations failed

        if len(self.failed) == len(self.out_files):
            self.send_email(type='radical-failed')
            return False
        
        return True
        

    def attach_bde_atom_property(self, out_file:Path):
        '''
        Use the list of out files to attach bond dissoaciation energies to atom 
        properties
        '''
        
        fname = out_file.stem
        index = fname.split('_')[-1] 
        # r_energy = self.gaussian_results[out_file]['r_user_Gaussian_Energy']
        r_energy = find_H_tot(out_file)
        energy = calculate_BDE(self.neutral_e, r_energy, e_h_radical=self.h_radical_ene)
        print(f'Energy for {out_file}: {r_energy} hartrees')
        self.indexes.append(int(index))
        atom = self.structure.atom[int(index)]
        print(f'Label assigned: {index} : {energy}')
        atom.property['r_user_BDE'] = energy

    def GenerateDataList(self):
        ''' 
        Take information from self.structure to create a Data list for report 
        generation
        '''
        
        print('Creating Reports...')
        print('Generating BDE table data...')
        data = [['Atom', 'BDE (kcal/mol)', 'Propensity']]
        
        for atom in self.structure.atom:
            try:
                bde = float(atom.property['r_user_BDE']) 
                if bde < self.high_coff:
                    propensity = 'High'
                if bde <= self.medium_coff and bde >= self.high_coff:
                    propensity = 'Moderate'                    
                if bde > self.medium_coff:
                    propensity = 'Low'
                row = [str(atom.element) + str(atom.index), 
                str(atom.property['r_user_BDE']), propensity]
                data.append(row)
            except KeyError:
                pass
        
        return data
    
    def CreatePDFReport(self):
        '''Create pdf report from a data dictionary and structural image'''
        #Define style
        
        style = TableStyle([
            ('ALIGN',(0,0),(-1,-1),'CENTER'),

            ('FONTNAME', (0,0), (-1,0), 'Courier-Bold'),
            ('FONTSIZE', (0,0), (-1,0), 10),

            ('BACKGROUND',(0,1),(-1,-1),colors.beige),
            ('BOX',(0,0),(-1,-1),2,colors.black),

            ('GRID',(0,1),(-1,-1),2,colors.black),
        ])
        # Get other styles

        styles = getSampleStyleSheet()
        styleNormal = styles['Normal']
        styleHeading = styles['Heading1']
        styleHeading2 = styles['Heading2']
        styleHeading.alignment = 1
        
        # Process data 

        pil_st = Image(str(self.img), width=2.0*inch, height=2.0*inch)
        pil_risk = Image(str(self.risk_scale), width=3.0*inch, height=0.8*inch)
        # Create PDF
        story = []
        
        # Append HEAD of document
        story.append(Paragraph('C-oxidation BDE Energy Report for: ' 
        + self.title, styleHeading))
        story.append(Spacer(inch, .25*inch))
        story.append(Paragraph('This report covers the results for BDE calculations'
        f' performed for: {self.title}. Oxudation propensity is established using '
        'C-H Bond Dissociation Enthalpies (BDE). The lower the C-H BDE values the '
        'higher the propensity for C-oxidation. Details for the DFT calculations and'
        ' overall workflow are explained at the end of this document', styleNormal))
        story.append(Spacer(inch, .05*inch))

        # Append molecule picture

        story.append(Paragraph('Bond Dissociation Energies (kcal/mol)', styleHeading2))
        story.append(pil_st)
        story.append(Spacer(inch, .25*inch))

        # Append BDE table
        main_table = Table(self.data, colWidths=85, rowHeights=18)
        main_table.setStyle(style)
        story.append(main_table)
        story.append(Spacer(inch, .15*inch))

        # Risk scale
        story.append(Paragraph('Risk Scale:', styleNormal))
        story.append(pil_risk)
        story.append(Spacer(inch, .25*inch))


        # Append final details
        story.append(Paragraph('Calculation details and output files', styleHeading2))
        story.append(Spacer(inch, .15*inch))
        
        # story.append(Paragraph(f'Report files are available in: {str(self.dir_report)}', 
        # styleNormal))
        story.append(Spacer(inch, .1*inch))
        story.append(Paragraph('Conformational search calculations were performed only'
        ' for the base ground state molecule. The lowest energy conformer was selected'
        ' to generate radicals and run optimization DFT calculations. DFT calculations'
        f' were performed using Gaussian with {self.dft} level of theory and {self.basis}'
        ' basis set. The BDE protocol was adapted from: <i>Lienard, P., Gavartin, J.,'
        ' Boccardi, G., & Meunier, M. (2015). Predicting drug substances autoxidation.'
        ' Pharmaceutical research, 32, 300-310.</i>', styleNormal))
        story.append(Spacer(inch, .1*inch))

        story.append(Paragraph(f'Output files:\n {str(self.dir_report)}')) 

        
        # Save PDF
        pdf_file = self.dir_report / f'{self.title}_CH-BDE_report.pdf'
        doc = SimpleDocTemplate(str(pdf_file), pagesize = letter, 
        title = f'{self.title} BDE Report', author = "BDE 2.0") 
        doc.build(story)

        return pdf_file

    def CreateExcelReport(self):
        '''Create Excel file with structure image embedded'''
        
        df = pd.DataFrame(self.data[1:], columns = self.data[0])
        df['BDE (kcal/mol)'] = pd.to_numeric(df['BDE (kcal/mol)'])
        
        output_xlsx = self.dir_report / f'{self.title}_BDE_report.xlsx'
        self.files_to_send.append(output_xlsx)
        df.to_excel(output_xlsx, index=False)
    
    def send_email(self, type='all-good'):
        '''Send emails. The accepted types are:
        all-good: All the calculations finished succesfuly
        base-failed: The base calculation failed
        radical-failed: One or more radical calculations failed.
        all-failed: All the calculations failed
        '''

        attach_string = '-a '
        attach_string += ' -a '.join([str(file.as_posix()) for file in self.files_to_send])
        print(self.files_to_send)
        contact_line = 'Please contact the molecular structure group for additional help.'

        if type == 'base-failed':
            mail_str = (f'mail -s "{self.title} C-oxidation calculations failed" '
            f'{self.email}<<<"The base BDE job submitted for {self.title} failed.'
            f' {contact_line}"')
        elif type == 'radical-failed':
            failed_str = ', '.join([file.stem for file in self.failed])
            mail_str = (f'mail -s "{self.title} C-oxidation calculations failed" '
            f'{self.email} <<<"The following BDE jobs submitted for {self.title}'
            f' failed: {failed_str}. {contact_line}"')
        elif type == 'all-failed':
            mail_str = (f'mail -s "{self.title} C-oxidation calculations failed or were'
            f' not found" {self.email}<<<"All the jobs submitted for {self.title} '
            'failed or the data was not found. This is indication of a major problem.'
            f' {contact_line}"')
        elif type == 'all-good':
            mail_str = (f'mail -s "{self.title} C-oxidation calculations" {attach_string}'
            f' {self.email}<<<"Please attached find the result files generated for '
            f'{self.title}. {contact_line}"')

        print(mail_str)
        mail_me = mail_str.replace(self.email, 'gabriel.valdivia-berroeta@boehringer-ingelheim.com')    
        os.system(mail_str)
        os.system(mail_me)


def find_Gibbs(log:Path):
   
    with open(log, 'r') as f:
        log_string = f.read()
    
    Hv_Hartree = find_H_tot(log)
    
    
    temp_string = re.search(r'Vibrational temperatures:\s*([\s\S]*?)(?:\n\n|\n\s*\n)', 
    log_string).group(1)
    numbers_only = re.findall(r'\d+\.?\d*', temp_string)
    #print(temp_string)
    numbers_only = [float(n) for n in numbers_only]
    Vib_temp1=np.array(numbers_only)
    Scale=0.9750 # M. Muzomwe et al. / Natural Science 4 (2012) 286-297
    Vib_temp=Scale*Vib_temp1
    
    Vib_temp_by_temp = Vib_temp/298.15
    partition = Vib_temp_by_temp/((np.exp(Vib_temp_by_temp))-1)
    partition = partition - np.log(1-np.exp(-Vib_temp_by_temp))

    S_vib_kcal_mol_K = 1.987204259*0.001*np.sum(partition)

#     print(log_string)
    lines = log_string.splitlines()
#     print(lines)
    for i in range(len(lines)):
#         print(i)
#         print(lines[i])
        if 'E (Thermal)' in lines[i]:
            start_index = i
        elif 'Vibration' in lines[i]:
            end_index = i-2
    #print(start_index)

    if start_index is not None and end_index is not None:
        table_lines = lines[start_index:end_index+1]

        table = []
        for row in table_lines:
            values = row.split()
            table.append(values)
   # print(table)
    
    for line in table:
        if line[0] == 'Electronic':
            S_electronic_kcal_mol_K = float(line[-1]) * 0.001
        if line[0] == 'Translational':
            S_translational_kcal_mol_K = float(line[-1]) * 0.001
        if line[0] == 'Rotational':
            S_rotational_kcal_mol_K = float(line[-1]) * 0.001
    
    S_kcal_mol_K = (S_rotational_kcal_mol_K + S_translational_kcal_mol_K + 
    S_electronic_kcal_mol_K)
    TS_kcal_mol = 298.15 * S_kcal_mol_K
    TS_Hartree = TS_kcal_mol/627.5095
    
    Gibbs = Hv_Hartree - TS_Hartree
    
    return Gibbs

def find_H_tot(log:Path):
   
    Hv_Hartree = 0
    
    with open(log, 'r') as f:
        log_string = f.read()
    
    try:
        temp_string = re.search(r'Vibrational temperatures:\s*([\s\S]*?)(?:\n\n|\n\s*\n)', 
        log_string).group(1)
    except AttributeError:
        print(f'No data was found for {log}')
        return
    numbers_only = re.findall(r'\d+\.?\d*', temp_string)
    #print(temp_string)
    numbers_only = [float(n) for n in numbers_only]
    Vib_temp1=np.array(numbers_only)
    Scale=0.9750 # M. Muzomwe et al. / Natural Science 4 (2012) 286-297
    Vib_temp=Scale*Vib_temp1
      
    half_temp = Vib_temp/2
    Vib_temp_by_temp = Vib_temp/298.15
    partition = Vib_temp/((np.exp(Vib_temp_by_temp))-1)
    ZPE_kcal_mol = 1.987204259*0.001*np.sum(half_temp)
    HvT_kcal_mol = 1.987204259*0.001*np.sum(partition)
    Hv_kcal_mol = ZPE_kcal_mol+HvT_kcal_mol
    Hv_Hartree = Hv_kcal_mol/627.5095


    log_string = log_string.replace('\n', '')
    log_string = re.sub(r'\s', '', log_string)
    
    if 'MP2' in log_string:
        MP2_values = re.findall(r"MP2=(-?\d+\.\d+)", log_string)  # MP2 energy
        E_electronic = float(MP2_values[-1])
                    
    else:
        hf_values = re.findall(r"HF=(-?\d+\.\d+)", log_string)  # SCF enegry 
        E_electronic = float(hf_values[-1])

    H_tot_Hartree = E_electronic + Hv_Hartree + 4*1.987204259*0.001*298.15/627.5095
    
    
    if H_tot_Hartree != 0:
        
        return H_tot_Hartree
        
    return None

class SoxidesAnalysis:
    def __init__(self, analysis_path:Path, medium_coff:float, high_coff:float,  
                oxygen_h:float=-75.22035293792803, inkscape:Path=None, 
                sample:str='sample', email:str='', dir_report:Path=Path().resolve()):
        
        self.analysis_path = analysis_path
        self.bases = self._find_bases()
        print(self.bases)  
        self.base_h = self._calcuate_boltz_energy(self.bases)
        self.oxygen_h = oxygen_h
        self.title = sample
        self.medium_coff = medium_coff
        self.high_coff = high_coff
        # self.structure = read_mae(str(base.with_suffix('.mae')))
        # self.dir_report = Path(self.analysis_path / 'Reports')
        self.dir_report = dir_report
        self.inkscape = inkscape

        # self.dir_report.mkdir(exist_ok=True)
        self._find_structure()
        self.already_oxidized = self._find_sulfoxides_on_st(self.structure)
        self.email = email

    def analyze(self):
        '''Main function to call'''

        check = self._get_method_basis(self.bases)
        print('Check status:')
        print(check)
        if not check:
            self.send_failed_email()
            print('All calculations for the S molecules failed. Abort')
            return
        
        self.find_sox_data()
        self.calculate_bde()
        self.create_data_list()
        self.csv = table_to_csv(self.data, f'{self.title}_sox')
        self.image_dir = create_image(self.structure, mode='S', 
        inkscape_path=self.inkscape)
        self.risk_scale = self.dir_report / 'S-oxidation_risk.png'
        create_risk_scale_png(self.medium_coff, self.high_coff, self.risk_scale)
        self.pdf_file = self.create_pdf_report()
        self.send_email()
    
    def send_failed_email(self):
        email_msg = (f'Sulfur oxidation for the {self.title} parent molecule failed. Pl'
        'ease check the molecule and the gaussian output files. Resubmission can fix th'
        'e problem. For additional questions please contact the molecular structure group')

        send_email(self.email, f'S-oxidation results for {self.title}',
                    message=email_msg) 
    
    def _find_sulfoxides_on_st(self, st):

        already_oxidized = []
        for atom in st.atom:
            if self.check_sulfoxide(atom):
                already_oxidized.append(atom.index)

        return already_oxidized

    def _find_bases(self):
        
        print(self.analysis_path)
        return [file for file in self.analysis_path.glob('*base*.out')]
    
    def _calcuate_boltz_energy(self, logs:list):
        '''
        logs: list of log_files (Path)
        '''
        enthalpies = []
        gibbs = []
        print(logs)

        for base_log in logs:
            enthalpy = find_H_tot(base_log)
            if not enthalpy:
                continue
            gibb = find_Gibbs(base_log)
            enthalpies.append(enthalpy)
            gibbs.append(gibb)
            print(enthalpies)
            print(gibbs)

        return getBoltzmannAverage(enthalpies, gibbs)


    def _find_structure(self):
        for file in self.analysis_path.glob('*base*.mae'):
            self.structure = structure.StructureReader.read(str(file))

    

    def _get_method_basis(self, bases:list):

        self.opt_funct = None
        self.opt_basis = None
        self.ene_funct = None
        self.ene_basis = None

        for out in bases:
            with open(out, 'r+') as f:
                lines = f.readlines()
            print(out)
            if not check_normal_termination(lines[-1]):
                continue
            for line in lines:
                if '# freq' in line:
                    method = line.split()[2]
                    self.opt_funct = method.split('/')[0]
                    self.opt_basis = method.split('/')[1]
                    break    
            for line in lines:
                if 'Geom=AllCheck Guess=Read' in line:
                    print(line)
                    method = line.split()[1]
                    self.ene_funct = method.split('/')[0]
                    self.ene_basis = method.split('/')[1]
                    break

        if self.opt_funct and self.opt_basis and self.ene_funct and self.ene_basis:
            return True

        return False


    
    def chiral_inverter(self, smiles):
        '''Invert the S=O chirality in an isomeric smiles string'''
        if '[S@@]' in smiles:
            return smiles.replace('[S@@]', '[S@]')
        elif '[S@]' in smiles:
            return smiles.replace('[S@]', '[S@@]')

    # def create_labeled_image(self):
    #     sdf_file = self.dir_report / f'{self.title}.sdf'
    #     self.structure.write(str(sdf_file))
        

    def find_sox_data(self):
        '''
        Find gaussian out files with sox- in it. It does not capture base compound
        Populates smiles:[gaussian_logs], smiles:[sts]          
        '''

        self.smiles_logs = dict()
        self.smiles_sts = dict()
        for file in self.analysis_path.glob('*sox-*.out'):
            st = structure.StructureReader.read(str(file.with_suffix('.mae')))
            smiles = adapter.to_smiles(st)
            if smiles in [i for i in self.smiles_logs.keys()]:
                self.smiles_logs[smiles].append(file)
                self.smiles_sts[smiles].append(st)
            else:
                self.smiles_logs[smiles] = [file]
                self.smiles_sts[smiles] = [st]


    
    def find_oxygens_attached(self, atom):
        '''
        Return index of oxygen attached to sulfur atom and total number of carbons
        '''
        carbons = 0
        oxygens = []
        for atom1 in atom.bonded_atoms:
            if atom1.atomic_number == 8:
                oxygens.append(atom1.index)
            elif atom1.atomic_number == 6:
                carbons += 1
        
        return oxygens, carbons
    
    def check_sulfoxide(self, atom):
        '''Check if the atom is a sulfur sulfoxide'''  

        if not atom.element == 'S':
            return False
        
        oxs_attached, carbons = self.find_oxygens_attached(atom)
        
        if len(oxs_attached) == 1 and carbons == 2 and atom.bond_total == 3:
            return True
        
        return False

    def find_atom_by_coors(self, st, coors):

        for atom in st.atom:
            if atom.xyz == coors:
                return str(atom.index)


    def find_oxidized_label(self, st):
        '''Find chirality and return index'''
        # Find oxidized sulfur. Only one is expected.
        for atom in st.atom:
            
            if not self.check_sulfoxide(atom):
                continue
            
            coors = atom.xyz
            oxygens, _ = self.find_oxygens_attached(atom)
            st.deleteAtoms(oxygens)
            st_renum = rmsd.renumber_conformer(self.structure, st)
            str_index = self.find_atom_by_coors(st_renum, coors)

            if str_index not in self.already_oxidized:
                return str_index

                
            


    def calculate_bde(self):
        ''' Calculate bde using self.enantions'''

        self.atom_bde = dict()
        print(self.smiles_logs)
        for smiles, logs in self.smiles_logs.items():
            so_h = self._calcuate_boltz_energy(logs)
            bde = 627.5095*(self.base_h + self.oxygen_h - so_h)
            atom_label = self.find_oxidized_label(self.smiles_sts[smiles][0])
            self.atom_bde[atom_label] = bde
             
        print(self.atom_bde)

    def create_data_list(self):
        '''Create self.data based on self.atom_bde'''

        self.data = [['Atom', 'S-O BDE (kcal/mol)', 'Propensity']]

        for atom, bde in self.atom_bde.items():
            if bde > self.high_coff:
                propensity = 'High'
            elif bde < self.high_coff and bde > self.medium_coff:
                propensity = 'Medium'
            elif bde < self.medium_coff:
                propensity = 'Low'
            self.data.append([f'S{atom}', str(round(bde, 1)), propensity])

    def send_email(self):
        '''Send email'''
        email_msg = ('Please find attached S-oxidation propensity results for:' 
        f' {self.title}.\nFor additional questions please contact the '
        'molecular structure group')

        attachs = [self.pdf_file, self.image_dir, self.csv]
        send_email(self.email, f'{self.title}: S-oxidation results',
                    message=email_msg, attachs=attachs) 

    def create_pdf_report(self):
        '''Create pdf report from a data dictionary and structural image'''
        #Define style
        
        style = TableStyle([
            ('ALIGN',(0,0),(-1,-1),'CENTER'),

            ('FONTNAME', (0,0), (-1,0), 'Courier-Bold'),
            ('FONTSIZE', (0,0), (-1,0), 10),

            ('BACKGROUND',(0,1),(-1,-1),colors.beige),
            ('BOX',(0,0),(-1,-1),2,colors.black),

            ('GRID',(0,1),(-1,-1),2,colors.black),
        ])
        # Get other styles

        styles = getSampleStyleSheet()
        styleNormal = styles['Normal']
        styleHeading = styles['Heading1']
        styleHeading2 = styles['Heading2']
        styleHeading.alignment = 1
        
        # Process data 

        pil_st = Image(str(self.image_dir), width=2.0*inch, height=2.0*inch)
        pil_risk = Image(str(self.risk_scale), width=3.0*inch, height=0.8*inch)
        # Create PDF
        story = []
        
        # Append HEAD of document
        story.append(Paragraph('S-O Bond dissociation Energy Report for: ' 
        + self.title, styleHeading))
        story.append(Spacer(inch, .25*inch))
        story.append(Paragraph('This report covers the results for S-O BDE calculations '
        f'performed for: {self.title}. S-oxidation propensity is estimated using S-O BDE '
        'calculations. The higher the S-O BDE values the higher the S-oxidation propensity. '
        'Details for the DFT calculations are found at the end of this document.', 
        styleNormal))
        story.append(Spacer(inch, .05*inch))

        # Append molecule picture

        story.append(Paragraph('Bond Dissociation Energies (kcal/mol)', styleHeading2))
        story.append(pil_st)
        story.append(Spacer(inch, .25*inch))

        # Append BDE table
        main_table = Table(self.data, colWidths=85, rowHeights=18)
        main_table.setStyle(style)
        story.append(main_table)
        story.append(Spacer(inch, .15*inch))

        # Risk scale
        story.append(Paragraph('Risk Scale:', styleNormal))
        story.append(pil_risk)
        story.append(Spacer(inch, .25*inch))


        # Append final details
        story.append(Paragraph('Calculation details and output files', styleHeading2))
        story.append(Spacer(inch, .15*inch))
        
        # story.append(Paragraph(f'Report files are available in: {str(self.dir_report)}', 
        # styleNormal))
        # story.append(Spacer(inch, .1*inch))
        story.append(Paragraph(f'DFT calculations were performed using Gaussian with ' 
        f'{self.opt_funct}/{self.opt_basis} method for optimization and frequency '
        f'calculations. For single point energy calculations, the {self.ene_funct}/'
        f'{self.ene_basis} method was employed. The S-O BDE protocol was implemented '
        'from: <i>Bose, A., Valdivia-Berroeta, G. A., & Gonnella, N. C. (2023). '
        'Predicting Autoxidation of Sulfides in APIs and Drug-like Molecules using '
        'QM/DFT Methods. Journal of Chemical Information and Modeling, 64(1), 128-137'
        '</i>', styleNormal))

        story.append(Spacer(inch, .1*inch))

        # Save PDF
        pdf_file = self.dir_report / f'{self.title}_SO-BDE_report.pdf'
        story.append(Paragraph(f'Output files:\n {str(pdf_file.parent)}')) 
        doc = SimpleDocTemplate(str(pdf_file), pagesize = letter, 
        title = f'{self.title} BDE Report', author = "BDE 2.0") 
        doc.build(story)

        return pdf_file


def get_tertiary_nitrogens(st) -> list:
    '''Get tertiary nitrogens'''

    tertiary_Ns = []
    for atom in st.atom:        
        if atom.formal_charge == 0:
            if atom.element == 'N':
                bonded_atoms = [b_atom.element for b_atom in atom.bonded_atoms]
                if 'H' not in bonded_atoms and len(bonded_atoms) > 1:
                    tertiary_Ns.append(atom.element + str(atom.index))
                
    return tertiary_Ns

class NoxidesAnalysis:

    def __init__(self, maegz_file:Path, alie_loc:Path, 
                medium_coff:float, high_coff:float,
                inkspace:Path, out_dir:Path=Path().resolve(), 
                sname:str='sample', email:str=''):
        
        self.alie_csv = self.find_alie_summary(alie_loc)
        self.out_dir = out_dir
        self.sname = sname
        self.medium_coff = medium_coff
        self.high_coff = high_coff
        self.basis, self.theory = self.get_basis_theory()
        self.structure = StructureReader.read(str(maegz_file))
        self.inkspace = inkspace
        self.email = email

    def analyze(self):
        '''Main function for the class'''
        
        df = self.get_nitrogen_df()
        if df is None:
            return
        self.process_N_df(df)
        # self.create_labeled_image()
        self.image_dir = create_image(self.structure, mode='N', 
        inkscape_path=self.inkspace)
        self.risk_scale = self.out_dir / 'N-oxidation_risk.png'
        create_risk_scale_png(self.medium_coff, self.high_coff, self.risk_scale)
        self.create_pdf_report()
        self.send_email()

    def find_alie_summary(self, alie_loc:Path):
        summary_csvs = [file for file in alie_loc.glob('*_summary.csv')]
        return summary_csvs[0]

    
    def get_basis_theory(self):
        '''Get basis and theory'''
        return '6-31G**', 'B3LYP'

    def get_nitrogen_df(self):
        '''Get nitrogen df. Only tertiary nitrogens are considered'''
        
        df_i = pd.read_csv(self.alie_csv)
        tert_Ns = get_tertiary_nitrogens(self.structure)
        #Obtain initial N df
        df_N = df_i[df_i['atom'].str.contains("N")]
        if df_N.empty:
            self.send_failed_email()
            return None
        print(df_N)
        Ns_switch_labels = [i[1:] + i[0] for i in tert_Ns]
        print(Ns_switch_labels)
        df_N = df_N[df_N['atom'].isin(Ns_switch_labels)]
        if df_N.empty:
            self.send_failed_email()
            return None
        print(df_N)
        df_N.to_csv(self.out_dir / './Ns_data.csv', index=False)
        
        return df_N

    def process_N_df(self, df_N:pd.DataFrame):
        '''Process N_df to add risk and change values to kcal/mol'''

        def add_risk(row):
            if row['alie_min'] < self.high_coff:
                return 'High'
            elif row['alie_min'] > self.high_coff and row['alie_min'] < self.medium_coff:
                return 'Medium'
            else:
                return 'Low'

        print(df_N)

        df_N['alie_min'] = df_N.apply(lambda row: row.ie_min*23.0609, axis=1)
        df_N['risk'] = df_N.apply(lambda row: add_risk(row), axis=1) 
        
        print(df_N)
        self.df_final = df_N[['atom','alie_min','risk', 'area']]
        self.final_csv = self.out_dir / f'{self.sname}_N-oxidation_report.csv'
        self.df_final.to_csv(self.final_csv, index=False)
        
        return self.df_final

    
    
    def create_pdf_report(self):
        '''Create pdf report from a data dictionary and structural image'''
        #Define style
        
        style = TableStyle([
            ('ALIGN',(0,0),(-1,-1),'CENTER'),

            ('FONTNAME', (0,0), (-1,0), 'Courier-Bold'),
            ('FONTSIZE', (0,0), (-1,0), 10),

            ('BACKGROUND',(0,1),(-1,-1),colors.beige),
            ('BOX',(0,0),(-1,-1),2,colors.black),

            ('GRID',(0,1),(-1,-1),2,colors.black),
        ])
        # Get other styles

        styles = getSampleStyleSheet()
        styleNormal = styles['Normal']
        styleHeading = styles['Heading1']
        styleHeading2 = styles['Heading2']
        styleHeading.alignment = 1
        
        # Process data
        self.data = self.df_final.values.tolist()
        for row in self.data:
            # Invert N label
            row[0] = row[0][-1] +  row[0][:-1]
            # Round up ALIEmin and area
            row[1] = round(row[1], 2)
            row[3] = round(row[3], 2)
        self.data.insert(0, ['Atom', 'ALIEmin (kcal/mol)', 'Ox. Risk', 'Area (A^2)'])

        pil_st = Image(str(self.image_dir), width=2.0*inch, height=2.0*inch)
        pil_risk = Image(str(self.risk_scale), width=3.0*inch, height=0.8*inch)

        
        # Create PDF
        story = []
        
        # Append HEAD of document
        story.append(Paragraph('N-oxidation Propensity Report for: ' + self.sname, 
        styleHeading))
        story.append(Spacer(inch, .25*inch))
        story.append(Paragraph(f'This report covers the results for N-oxide propensity'
        f' calculations performed for: {self.sname}. N-oxidation propensity is established ' 
        'using Average Local Ionization Energy (ALIE). Lower ALIE values indicate higher '
        'propensity for N-oxidation.'
        , styleNormal))
        story.append(Spacer(inch, .05*inch))

        # Append molecule picture

        story.append(Paragraph('Average Local Ionization Energies (kcal/mol)', 
        styleHeading2))
        story.append(pil_st)
        story.append(Spacer(inch, .25*inch))

        # Append ALIE table
        main_table = Table(self.data, colWidths=105, rowHeights=18)
        main_table.setStyle(style)
        story.append(main_table)
        story.append(Spacer(inch, .15*inch))

        # Risk scale
        story.append(Paragraph('Risk Scale:', styleNormal))
        story.append(pil_risk)
        story.append(Spacer(inch, .25*inch))

        story.append(Paragraph(f'Nitrogens showing Area values less than 1 A should be '
        f' investigated for potential hydrogen bond formation. If the hydrogen bond form'
        'ation needs to be avoided, please run a calculation using Advanced Options > C'
        'onformational search parameters > Restrained torsions.', styleNormal))

        # Append final details
        story.append(Paragraph('Calculation details and output files', styleHeading2))
        story.append(Spacer(inch, .15*inch))
        
        story.append(Paragraph(f'DFT calculations were performed using Gaussian with' 
        f' {self.theory} level of theory and {self.basis} basis set.', styleNormal))
        story.append(Spacer(inch, .1*inch))

        story.append(Paragraph(f'N-oxidation ALIE protocol was incorporated from: ' 
        f'<i>Valdivia-Berroeta, G. A., & Gonnella, N. C. (2023). N-oxidation Regiosel'
        'ectivity and Risk Prediction Using DFT-ALIE Calculations. Pharmaceutical'
        ' Research, 1-11</i>', styleNormal))
        story.append(Spacer(inch, .1*inch))
        story.append(Paragraph(f'Output files:\n {str(self.out_dir)}')) 

            
        # Save PDF
        pdf_file = self.out_dir / f'{self.sname}_N-oxidation_report.pdf'
        doc = SimpleDocTemplate(str(pdf_file), pagesize = letter, 
        title = f'{self.sname} N-oxidation Report', 
        author = "N-oxidation Propensity 1.0") 
        doc.build(story)
        self.pdf_file = pdf_file

        return pdf_file

    # def find_dft_information(self):

    

    def send_failed_email(self):
        email_msg = ('Nitrogen atoms propense to N-oxidation were not found in:' 
        f' {self.sname}.\nFor additional questions please contact the '
        'molecular structure group')

        send_email(self.email, f'N-oxidation results for {self.sname}',
                    message=email_msg) 

    def send_email(self):
        '''Send email'''
        email_msg = ('Please find attached N-oxidation propensity results for:' 
        f' {self.sname}.\nFor additional questions please contact the '
        'molecular structure group')
        xlsx_path = self.final_csv.parent / f'{self.sname}_summary.xlsx'
        shutil.copy(str(self.alie_csv.with_suffix('.xlsx')), str(xlsx_path))
        attachs = [
                    self.final_csv, 
                    self.pdf_file, 
                    self.image_dir, 
                    xlsx_path,
                ]
        send_email(self.email, f'{self.sname}: N-oxidation results',
                    message=email_msg, attachs=attachs) 


if __name__ == '__main__':

    logger.info('Runtime command: %s' % subprocess.list2cmdline(sys.argv))



 