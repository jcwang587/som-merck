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


if __name__ == '__main__':

    logger.info('Runtime command: %s' % subprocess.list2cmdline(sys.argv))

    # Generate the pdf report using the CoxidesAnalysis

 
