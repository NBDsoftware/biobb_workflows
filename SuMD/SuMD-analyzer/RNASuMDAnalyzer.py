import os
import sys
import tqdm
import yaml
import logging
from glob import glob
from collections import Counter
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.colors as mplcolors
import seaborn as sns
import MDAnalysis as mda
import MDAnalysis.analysis.rms
import MDAnalysis.analysis.align as align
from scipy.stats import iqr
from scipy.interpolate import make_interp_spline
import barnaba as bb
import multiprocessing
import argparse

##################################################

HEADER = '''
#      ____  _   _____   _____       __  _______ 
#     / __ \/ | / /   | / ___/__  __/  |/  / __ \ 
#    / /_/ /  |/ / /| | \__ \/ / / / /|_/ / / / / 
#   / _, _/ /|  / ___ |___/ / /_/ / /  / / /_/ /  
#  /_/ |_/_/ |_/_/  |_/____/\__,_/_/  /_/_____/   
#     /   |  ____  ____ _/ __  ______ ___  _____  
#    / /| | / __ \/ __ `/ / / / /_  // _ \/ ___/  
#   / ___ |/ / / / /_/ / / /_/ / / //  __/ /      
#  /_/  |_/_/ /_/\__,_/_/\__, / /___\___/_/       
#                       /____/                    

               M.Pavan 13/01/2022
              P.Navarro 28/06/2023
'''

HELP = '''
Available protocols:
    geometry
    mmgbsa
    intEnergy
    perResRec
    perResLig
    matrix
    video 
'''

LOG_FILENAME = "sumd_analysis.log"

# Set up logging
# Create logger
logger = logging.getLogger("SuMD Analyzer")

# Create file handler
file_handler = logging.FileHandler(LOG_FILENAME, mode = "w")
file_handler.setLevel(logging.INFO)

# Create console handler
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)

# Create formatter and add it to file handler 
file_format = logging.Formatter('%(name)s [%(levelname)s]: %(message)s')
file_handler.setFormatter(file_format)

# Create formatter and add it to console handler
console_format = logging.Formatter('%(name)s [%(levelname)s]: %(message)s')
console_handler.setFormatter(console_format)

# Set logger level
logger.setLevel(logging.INFO)

# Add file handler to the logger
logger.addHandler(file_handler)

# Add console handler to the logger
logger.addHandler(console_handler)

logger.info(HEADER)
logger.info("=================================== \n")
logger.info(HEADER)
logger.info("=================================== \n")

##################################################
### MULTIPROCESSING FUNCTION AND MODULES

class parallelizer(object):
    def __init__(self, args, func, num_procs, desc):
        self.num_procs = num_procs
        self.args = args
        self.func = func
        self.desc = desc

    def start(self):
        pass

    def end(self):
        pass

    def run(args, func, num_procs, desc):
        return MultiThreading(args, func, num_procs, desc)

def MultiThreading(args, func, num_procs, desc):
    results = []
    tasks = []
    for index,item in enumerate(args):
        task = (index, (func, item))
        tasks.append(task)
    
    results = start_processes(tasks, num_procs, desc)
    return results

def worker(input, output):
    for seq, job in iter(input.get, "STOP"):
        func, args = job
        result = func(*args)
        ret_val = (seq, result)
        output.put(ret_val)

def start_processes(inputs, num_procs, desc):
    task_queue = multiprocessing.Queue()
    done_queue = multiprocessing.Queue()

    for item in inputs:
        task_queue.put(item)

    pbar = tqdm.tqdm(total=len(inputs), desc=desc)
    
    for i in range(num_procs):
        multiprocessing.Process(target=worker, args=(task_queue, done_queue)).start()

    results = []
    for i in range(len(inputs)):
        results.append(done_queue.get())
        pbar.update(1)

    for i in range(num_procs):
        task_queue.put("STOP")

    results.sort(key=lambda tup: tup[0])
    return [item[1] for item in map(list, results)]

##################################################

def plotRMSDvsReference(top, trj, ref_pdb, receptor_backbone_selection, ligand_selection):
    ### this function calculate and plots the ligand RMSD vs reference
    logger.info('\nCalculating ligand RMSD vs reference...')
    
    # Load trajectory and topology
    u = mda.Universe(top,trj)

    # Create time array
    time = np.array([i for i in range(0, len(u.trajectory))]) * add / 1000000 
    
    # Load reference structure
    ref = mda.Universe(ref_pdb)
    
    # Create RMSD object
    R = MDAnalysis.analysis.rms.RMSD(u, ref,
        select=receptor_backbone_selection,             #  selection for alignment/superposition - backbone of receptor
        groupselections=[ligand_selection])     #  selections for RMSD calculation after alignment
    
    R.run()

    rmsd = R.rmsd.T     # transpose makes it easier for plotting

    logger.info(rmsd)

    with open('Geometry/time_and_rmsd.txt','w') as f:
        for ts,el in zip(time,rmsd[3]):
            output = str(ts) + ',' + str(el) + '\n'
            f.write(output) 
   
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(time, rmsd[3], 'k-', linewidth=0.75)
    # ax.legend(loc="best")
    ax.set_xlabel("Time (ns)", fontsize=14)
    ax.set_ylabel(r"RMSD ($\AA$)", fontsize=14)
    ax.set_title('Ligand RMSD to reference', fontsize=16)
    plt.tight_layout()
    fig.savefig("Geometry/rmsd_lig_vs_reference.png", dpi=300)
    plt.close()

def plotDynamicRMSDvsReference():
    #### this function dynamically plots the ligand RMSD vs reference for video purpose 
    logger.info('\nDynamically plotting RMSD to reference...')

    df = pd.read_csv('Geometry/time_and_rmsd.txt', header=None, names=['Time','RMSD'])
    time = list(df['Time'])
    space = list(df['RMSD'])
    count = 0

    if not os.path.exists('Geometry/provideo'):
        os.makedirs('Geometry/provideo')

    for el in time:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(time[:count+1], space[:count+1], 'k-')
        ax.set_xlabel("Time (ns)", fontsize=14)
        ax.set_ylabel(r"RMSD ($\AA$)", fontsize=14)
        ax.set_title('Ligand RMSD$_{backbone}$ to reference', fontsize=16)
        plt.xlim(min(time)-1, max(time)+1)
        plt.ylim(min(space)-1, max(space)+1)
        plt.tight_layout()
        fig.savefig(f"Geometry/provideo/rmsd_lig_vs_reference_{count}.png", dpi=300)   
        plt.clf()
        count += 1
        plt.close()

def plotCDMDist(top, trj, ligand_selection, binding_site_selection):
    #### this function calculate and plots the distance between the center of mass of the ligand and the binding site at each frame
    logger.info('\nCalculating center of mass distance...')   

    u = mda.Universe(top,trj)
    time = np.array([i for i in range(0, len(u.trajectory))]) * add / 1000000

    cdmDistance = []
    for ts in u.trajectory:
        ligandCDM = u.select_atoms(ligand_selection).center_of_mass()
        binding_siteCDM = u.select_atoms(binding_site_selection).center_of_mass()
        distance = np.linalg.norm(ligandCDM-binding_siteCDM)
        cdmDistance.append(distance)

    with open('Geometry/time_and_space.txt','w') as f:
        for ts,el in zip(time,cdmDistance):
            output = str(ts) + ',' + str(el) + '\n'
            f.write(output)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(time, cdmDistance, 'k-')
    ax.set_xlabel("Time (ns)", fontsize=14)
    ax.set_ylabel(r"cdm$_{bind site-lig}$ distance ($\AA$)", fontsize=14)
    ax.set_title('cdm$_{bind site-lig}$ distance', fontsize=16)
    plt.tight_layout()
    fig.savefig("Geometry/cdm_distance.png", dpi=300)
    plt.close()
    
def plotDynamicCDMDist():
    #### this function dynamically plots the distance between the center of mass of the ligand and the binding site for video purpose 
    logger.info('\nDynamically plotting center of mass distance...')

    df = pd.read_csv('Geometry/time_and_space.txt', header=None, names=['Time','Space'])
    time = list(df['Time'])
    space = list(df['Space'])
    count = 0

    if not os.path.exists('Geometry/provideo'):
        os.makedirs('Geometry/provideo')

    for el in time:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(time[:count+1], space[:count+1], 'k-')
        ax.set_xlabel("Time (ns)", fontsize=14)
        ax.set_ylabel(r"cdm$_{rec-lig}$ distance ($\AA$)", fontsize=14)
        ax.set_title('cdm$_{rec-lig}$ distance', fontsize=16)
        plt.xlim(min(time)-1, max(time)+1)
        plt.ylim(min(space)-1, max(space)+1)
        plt.tight_layout()
        fig.savefig(f"Geometry/provideo/cdm_distance_{count}.png", dpi=300)   
        plt.clf()
        plt.close()

        count += 1
    
def plotRMSDRec(top, trj, receptor_backbone_selection):
    #### this function calculate and plots the receptor RMSD
    logger.info('\nCalculating receptor RMSD...')

    u = mda.Universe(top,trj)
    time = np.array([i for i in range(0, len(u.trajectory))]) * add / 1000000 
    ref = mda.Universe(top,trj)

    # Create RMSD object
    R = MDAnalysis.analysis.rms.RMSD(u, ref,
        select=receptor_backbone_selection,                                           #  selection for alignment/superposition
        groupselections=[receptor_backbone_selection])     #  selections for RMSD calculation after alignment
    
    R.run()
    rmsd = R.rmsd.T   # transpose makes it easier for plotting

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(time, rmsd[3], 'k-', linewidth=0.75, color='purple')
    # ax.legend(loc="best")
    ax.set_xlabel("Time (ns)", fontsize=14)
    ax.set_ylabel(r"RMSD ($\AA$)", fontsize=14)
    ax.set_title('Receptor RMSD$_{backbone}$', fontsize=16)
    plt.tight_layout()
    fig.savefig("Geometry/rmsd_rec.png", dpi=300)
    plt.close()
    
def plotRMSDLig(top, trj, receptor_backbone_selection, ligand_selection):
    #### this function calculate and plots the ligand RMSD
    logger.info('\nCalculating ligand RMSD...')
    
    u = mda.Universe(top,trj)
    time = np.array([i for i in range(0, len(u.trajectory))]) * add / 1000000 
    ref = mda.Universe(top,trj)

    R = MDAnalysis.analysis.rms.RMSD(u, ref,
        select=receptor_backbone_selection,
        groupselections=[ligand_selection])
    
    R.run()
    rmsd = R.rmsd.T   # transpose makes it easier for plotting

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(time, rmsd[3], 'k-', linewidth=0.75, color='green')
    # ax.legend(loc="best")
    ax.set_xlabel("Time (ns)", fontsize=14)
    ax.set_ylabel(r"RMSD ($\AA$)", fontsize=14)
    ax.set_title('Ligand RMSD$_{backbone}$', fontsize=16)
    plt.tight_layout()
    fig.savefig("Geometry/rmsd_lig.png", dpi=300)
    plt.close()
  
def plotRMSFRec(top, trj, receptor_backbone_selection):
    #### this function calculate and plots the receptor RMSF
    logger.info('\nCalculating receptor RMSF...')

    u = mda.Universe(top,trj)
    
    # Align trajectory to the average structure
    average = align.AverageStructure(u, u, select=receptor_backbone_selection, ref_frame=0).run()
    ref = average.results.universe
    align.AlignTraj(u, ref, select=receptor_backbone_selection,in_memory=True).run()

    # Select the receptor backbone atoms
    backbone_atoms = u.select_atoms(receptor_backbone_selection)  

    # Get the residue number for each backbone atom
    residue_per_backbone_atom = backbone_atoms.resnums

    # Get the unique residue numbers
    residues = list(set(residue_per_backbone_atom))

    # Calculate the RMSF per backbone atom
    RMSF = MDAnalysis.analysis.rms.RMSF(backbone_atoms).run()

    plt.figure()

    rmsf_per_residue = []

    # Obtain the RMSF per residue
    for i in range (0, len(residues)):
        rmsf_per_residue.append(RMSF.rmsf[residue_per_backbone_atom == residues[i]].mean())

    plt.plot(residues, rmsf_per_residue, linewidth=0.75, color='purple')
    plt.title('Receptor RMSF$_{backbone}$', fontsize=16)
    plt.xlabel('Residue number', fontsize=14)
    plt.ylabel('RMSF (Å)', fontsize=14)
    plt.xlim(min(residues), max(residues))
    plt.ylim(0, max(rmsf_per_residue) + 0.1*max(rmsf_per_residue))
    plt.tight_layout()
    plt.savefig('Geometry/rmsf_rec.png', dpi=300)
    plt.close()

    #### RMSF through time per residue
    rmsf_timeline = []
    count= 1

    for ts in u.trajectory:
        if count <= len(u.trajectory):
            # Calculate the RMSF per backbone atom
            RMSF = MDAnalysis.analysis.rms.RMSF(backbone_atoms).run(start=0, stop=count)

            rmsf_per_residue = []

            # Obtain the RMSF per residue
            for i in range (0, len(residues)):
                rmsf_per_residue.append(RMSF.rmsf[residue_per_backbone_atom == residues[i]].mean())

            rmsf_timeline.append(rmsf_per_residue)
            count +=1
        else:
            break

    rmsf_timeline = np.array(rmsf_timeline)

    if os.path.exists('Geometry/rmsf_timeline_rec.tsv'):
        os.system('rm Geometry/rmsf_timeline_rec.tsv')

    for frame, rmsf_array in enumerate(rmsf_timeline):
        for resid, rmsf_val in zip(residues, rmsf_array):
            with open('Geometry/rmsf_timeline_rec.tsv','a') as f:
                line = str((frame)*add/1000000) + '\t' + str(resid) + '\t' + str(rmsf_val) + '\n'
                f.write(line)

    plt.figure()

    sod='Geometry/rmsf_timeline_rec.tsv'
    # cm = plt.cm.get_cmap('RdYlBu_r')
    cm = plt.cm.get_cmap('Purples')
    x = np.genfromtxt(sod, usecols = [0], delimiter='\t')
    y = np.genfromtxt(sod, usecols = [1], delimiter='\t')
    colors = np.genfromtxt(sod, usecols = [2], delimiter='\t')
    colors_list = list(colors)
    colors_list.sort()

    #### use scipy.stats to calculate IQR (Inter Quartile Range)
    IQR = iqr(colors_list, axis=0, nan_policy='omit')
    #### use numpy to calculate first quartile and third quartile
    Q1 = np.nanquantile(colors_list, 0.02)
    mask_min = Q1-1.5*IQR
    Q3 = np.nanquantile(colors_list, 0.98)
    mask_max = Q3+1.5*IQR

    plt.scatter(x,y,c=colors, s=15, cmap=cm, marker='s', vmin=Q1, vmax=Q3, linewidths= 0, edgecolors='none')
    cbar = plt.colorbar()
    cbar.set_label('RMSF (Å)', rotation=270, labelpad=15)
    plt.title('Receptor RMSF$_{backbone}$', fontsize=16)
    plt.ylabel('Residue number', fontsize=14)
    plt.xlabel('Time (ns)', fontsize=14)
    plt.xlim(min(x), max(x))
    plt.ylim(min(y), max(y))
    plt.tight_layout()
    plt.savefig('Geometry/rmsf_timeline_rec.png', dpi=300)
    plt.close()

"""
def plotRMSFLig():
    #### this function calculate and plots the ligand RMSF
    logger.info('\nCalculating ligand RMSF...')

    u = mda.Universe(top,trj)
    ref = mda.Universe(top,trj)
    align.AlignTraj(u,ref,
                select=,  
                match_atoms=True,  # whether to match atoms based on mass
                in_memory=True
               ).run()
    
    nonSTDresidues = 'G5 C3'    
    calphas = u.select_atoms(f'''({ligandSel} or resname {nonSTDresidues}) and name C5' ''')
    resids = calphas.resnums + numShiftLig
    rmsfer = MDAnalysis.analysis.rms.RMSF(calphas).run()

    plt.figure()
    plt.plot(rmsfer.rmsf, resids, linewidth=0.75, color='green')
    plt.title('Ligand RMSF$_{backbone}$', fontsize=16)
    plt.ylabel('Residue number', fontsize=14)
    plt.xlabel('RMSF (Å)', fontsize=14)
    plt.ylim(min(resids), max(resids))
    plt.tight_layout()
    plt.savefig('Geometry/rmsf_lig.png', dpi=300)
    plt.close()

    #### RMSF through time per residue
    rmsf_timeline = []
    count= 1

    for ts in u.trajectory:
        if count <= len(u.trajectory):
            rmsfer = MDAnalysis.analysis.rms.RMSF(calphas).run(start=0, stop=count)
            rmsf_timeline.append(rmsfer.rmsf)
            count +=1
        else:
            break

    rmsf_timeline = np.array(rmsf_timeline)

    if os.path.exists('Geometry/rmsf_timeline_lig.tsv'):
        os.system('rm Geometry/rmsf_timeline_lig.tsv')

    for frame, rmsf_array in enumerate(rmsf_timeline):
        for resid, rmsf_val in zip(resids, rmsf_array):
            with open('Geometry/rmsf_timeline_lig.tsv','a') as f:
                line = str((frame)*add/1000000) + '\t' + str(resid) + '\t' + str(rmsf_val) + '\n'
                f.write(line)

    plt.figure()

    sod='Geometry/rmsf_timeline_lig.tsv'
    # cm = plt.cm.get_cmap('RdYlBu_r')
    cm = plt.cm.get_cmap('Greens')
    x = np.genfromtxt(sod, usecols = [0], delimiter='\t')
    y = np.genfromtxt(sod, usecols = [1], delimiter='\t')
    colors = np.genfromtxt(sod, usecols = [2], delimiter='\t')
    colors_list = list(colors)
    colors_list.sort()

    #### use scipy.stats to calculate IQR (Inter Quartile Range)
    IQR = iqr(colors_list, axis=0, nan_policy='omit')
    #### use numpy to calculate first quartile and third quartile
    Q1 = np.nanquantile(colors_list, 0.02)
    mask_min = Q1-1.5*IQR
    Q3 = np.nanquantile(colors_list, 0.98)
    mask_max = Q3+1.5*IQR

    plt.scatter(x,y,c=colors, s=15, cmap=cm, marker='s', vmin=Q1, vmax=Q3, linewidths= 0, edgecolors='none')
    cbar = plt.colorbar()
    cbar.set_label('RMSF (Å)', rotation=270, labelpad=15)
    plt.title('Ligand RMSF$_{backbone}$', fontsize=16)
    plt.ylabel('Residue number', fontsize=14)
    plt.xlabel('Time (ns)', fontsize=14)
    plt.xlim(min(x), max(x))
    plt.ylim(min(y), max(y))
    plt.tight_layout()
    plt.savefig('Geometry/rmsf_timeline_lig.png', dpi=300)
    plt.close()
"""

def plotRgRec(top, trj, receptorSel, add):
    #### this function calculate and plots the receptor Rg
    logger.info('\nCalculating receptor radius of gyration...')

    u = mda.Universe(top,trj)
    Rgyr = []
    receptor = u.select_atoms(receptorSel)

    for i, ts in enumerate(u.trajectory):
        ns = i*add/1000000
        Rgyr.append((ns, receptor.radius_of_gyration()))

    Rgyr = np.array(Rgyr)

    plt.figure()
    ax = plt.subplot(111)
    ns = Rgyr[:,0]
    radius = Rgyr[:,1]
    xnew = np.linspace(ns.min(), ns.max(), 300) 
    spl1 = make_interp_spline(ns, radius, k=3)
    power_smooth1 = spl1(xnew)
    ax.plot(xnew, power_smooth1, linewidth=0.75, color='purple')
    ax.set_xlabel("Time (ns)", fontsize=14)
    ax.set_ylabel(r"radius of gyration $R_G$ ($\AA$)", fontsize=14)
    plt.title("Receptor Radius of Gyration", fontsize=16)
    plt.tight_layout()
    ax.figure.savefig("Geometry/Rgyr_receptor.png",dpi=300)
    plt.close()

def plotRgLig(top, trj, ligandSel, add):
    #### this function calculate and plots the ligand Rg
    logger.info('\nCalculating ligand radius of gyration...')

    u = mda.Universe(top,trj)
    Rgyr = []
    protein = u.select_atoms(ligandSel)

    for i, ts in enumerate(u.trajectory):
        ns = i*add/1000000
        Rgyr.append((ns, protein.radius_of_gyration()))

    Rgyr = np.array(Rgyr)

    plt.figure()
    ax = plt.subplot(111)
    ns = Rgyr[:,0]
    radius = Rgyr[:,1]
    xnew = np.linspace(ns.min(), ns.max(), 300) 
    spl1 = make_interp_spline(ns, radius, k=3)
    power_smooth1 = spl1(xnew)
    ax.plot(xnew, power_smooth1, linewidth=0.75, color='green')
    ax.set_xlabel("Time (ns)", fontsize=14)
    ax.set_ylabel(r"radius of gyration $R_G$ ($\AA$)", fontsize=14)
    plt.title("Ligand Radius of Gyration", fontsize=16)
    plt.tight_layout()
    ax.figure.savefig("Geometry/Rgyr_ligand.png",dpi=300)
    plt.close()
    
def plotERMSD(top, trj, ref_pdb, ref_bool, add):
    #### this function exploits the barnaba python package to calculate ERMSD
    logger.info('\nCalculating ERMSD using barnaba...')

    if ref_bool:
        native = ref_pdb
    else:
        native = top

    # calculate eRMSD between native and all frames in trajectory
    ermsd = bb.ermsd(native,trj,topology=top)
    time = np.array([i for i in range(0, len(ermsd))]) * add / 1000000

    # plot time series
    plt.figure()
    plt.plot(time, ermsd, linewidth=0.75, color='orange')
    plt.ylabel("eRMSD from native")
    plt.xlabel('Time (ns)')
    plt.tight_layout()
    plt.savefig("Geometry/ERMSD_vs_time.png",dpi=300)
    plt.clf()

    # make histogram
    plt.figure()
    plt.hist(ermsd,density=True,bins=50, stacked=True, alpha=0.5, color='orange')
    plt.xlabel("eRMSD from native")
    plt.ylabel('Probability density')
    plt.tight_layout()
    plt.savefig("Geometry/ERMSD_distribution.png",dpi=300)
    plt.clf()
    
    # calculate RMSD
    rmsd = bb.rmsd(native,trj,topology=top,heavy_atom=False)
    # plot time series
    plt.figure()
    plt.plot(time, rmsd, linewidth=0.75)
    plt.ylabel("RMSD from native (nm)")
    plt.xlabel('Time (ns)')
    plt.tight_layout()
    plt.savefig("Geometry/RMSD_vs_time.png",dpi=300)
    plt.clf()

    # make histogram
    plt.hist(rmsd,density=True,bins=50, stacked=True, alpha=0.5)
    plt.xlabel("RMSD from native (nm)")
    plt.ylabel('Probability density')
    plt.tight_layout()
    plt.savefig("Geometry/RMSD_distribution.png",dpi=300)
    plt.clf()
    
    # combined plot
    plt.xlabel("eRMSD from native")
    plt.ylabel("RMSD from native (nm)")
    plt.axhline(0.4,ls = "--", c= 'k')
    plt.axvline(0.7,ls = "--", c= 'k')
    plt.scatter(ermsd,rmsd,s=2.5)
    plt.tight_layout()
    plt.savefig("Geometry/barnaba.png",dpi=300)
    plt.clf()

    plt.close()
    
    a = 'Geometry/ERMSD_vs_time.png'
    b = 'Geometry/ERMSD_distribution.png'
    c = 'Geometry/RMSD_vs_time.png'
    d = 'Geometry/RMSD_distribution.png'

    os.system(f'montage -tile 2x2 -geometry 1920x1080 {a} {b} {c} {d} Geometry/merged_barnaba.png') 

def mountPanelRec():
    #### this function mounts a 4 tile panel for receptor geometric analysis
    logger.info('\nMounting receptor 4 tile panel...')   

    a = 'Geometry/rmsd_rec.png'
    b = 'Geometry/Rgyr_receptor.png'
    c = 'Geometry/rmsf_rec.png'
    d = 'Geometry/rmsf_timeline_rec.png'

    os.system(f'montage -tile 2x2 -geometry 1920x1080 {a} {b} {c} {d} Geometry/merged_rec.png')    
    
def mountPanelLig():
    #### this function mounts a 4 tile panel for ligand geometric analysis
    logger.info('\nMounting ligand 4 tile panel...')  

    a = 'Geometry/rmsd_lig.png'
    b = 'Geometry/Rgyr_ligand.png'
    c = 'Geometry/rmsf_lig.png'
    d = 'Geometry/rmsf_timeline_lig.png'

    os.system(f'montage -tile 2x2 -geometry 1920x1080 {a} {b} {c} {d} Geometry/merged_lig.png') 
    
##################################################################################################################################################
##################################################################################################################################################

def calcMMGBSA(trj, receptor_topology, ligand_topology, dry_complex_topology):
    #### this function calculates MMGBSA interaction energy between receptor and ligand
    logger.info('\nCalculating MMGBSA interaction energy...')

    mmgbsa = '''Input file for running PB and GB
&general
verbose=1,
/
&gb
igb=5, saltcon=0.100
/
'''

    with open('mmgbsa.in','w') as f:
        f.write(mmgbsa)

    os.system(f''' python $AMBERHOME/bin/MMPBSA.py -i mmgbsa.in -o MMGBSA/gb_TOT.dat -eo MMGBSA/gb_frame.dat -cp {dry_complex_topology} -rp {receptor_topology} -lp {ligand_topology} -y {trj} > /dev/null 2>&1 ''')

    # Remove all the temporal files named _MMPBSA*
    os.system('rm _MMPBSA*')

    os.system('rm mmgbsa.in')
    
def processMMGBSA():
    #### this function parses MMGBSA output files in order to extract meaningful values for plotting purpose
    logger.info('\nProcessing MMGBSA output files...')

    with open('MMGBSA/gb_frame.dat','r') as f, open('MMGBSA/gb_DELTA.dat','w') as n:
        for line in f:
            if line.startswith('DELTA Energy Terms'):
                while line != '':
                    n.write(line)
                    line = f.readline()
            else:
                continue
    
    best_frame = 0
    best_mmgbsa = 0

    with open('MMGBSA/gb_DELTA.dat','r') as f:
        for line in f:
            if line.startswith('DELTA Energy Terms') or line.startswith('Frame'):
                continue
            elif not line.startswith('') or not line.startswith('\n'):
                frame = line.split(',')[0]
                gb = line.split(',')[7]
                if float(gb) <= best_mmgbsa:
                    best_mmgbsa = float(gb)
                    best_frame = frame   

    with open('MMGBSA/gb_TOT.dat','r') as f:
        for line in f:
            if line.startswith('DELTA TOTAL'):
                mmgbsa = str("{:.2f}".format(float(line.split()[2].strip())))
            else:
                pass                 
    
def plotMMGBSAvsDistance(ref_bool):
    #### this function plots MMGBSA interaction energy as a function of the cdm distance
    logger.info('\nPlotting MMGBSA vs cdm distance...')

    if ref_bool:
        txt = 'time_and_rmsd.txt'
    else:
        txt = 'time_and_space.txt'

    df = pd.read_csv(f'Geometry/{txt}', header=None, names=['Time','Distance'])
    gb_list = []

    with open('MMGBSA/gb_DELTA.dat','r') as f:
        for line in f:
            if line.startswith('DELTA Energy Terms') or line.startswith('Frame'):
                continue
            elif not line.startswith('') or not line.startswith('\n'):
                gb = float(line.split(',')[7])
                gb_list.append(gb)

    x = list(df['Distance'])
    y = gb_list
    colors = gb_list    
    colors_list = list(colors)
    colors_list.sort()

    #### use scipy.stats to calculate IQR (Inter Quartile Range)
    IQR = iqr(colors_list, axis=0, nan_policy='omit')
    #### use numpy to calculate first quartile and third quartile
    Q1 = np.nanquantile(colors_list, 0.25)
    mask_min = Q1-1.5*IQR
    Q3 = np.nanquantile(colors_list, 0.75)
    mask_max = Q3+1.5*IQR    

    plt.figure()
    plt.scatter(x, y, c=colors, cmap='RdYlBu_r')
    plt.colorbar()
    plt.title('Interaction Energy Landscape', fontsize=16)
    plt.ylabel('MMGBSA (kcal/mol)', fontsize=14)
    plt.xlabel('dcm$_{bs-lig}$ ($\AA$)', fontsize=14)
    plt.tight_layout()
    plt.savefig('MMGBSA/mmgbsa_vs_distance.png',dpi=300)
    plt.close()

def plotMMGBSAvsTime(ref_bool):
    #### this function plots MMGBSA interaction energy as a function of the simulation time
    logger.info('\nPlotting MMGBSA vs simulation time...')

    if ref_bool:
        txt = 'time_and_rmsd.txt'
    else:
        txt = 'time_and_space.txt'

    df = pd.read_csv(f'Geometry/{txt}', header=None, names=['Time','Distance'])
    gb_list = []

    with open('MMGBSA/gb_DELTA.dat','r') as f:
        for line in f:
            if line.startswith('DELTA Energy Terms') or line.startswith('Frame'):
                continue
            elif not line.startswith('') or not line.startswith('\n'):
                gb = float(line.split(',')[7])
                gb_list.append(gb)

    x = list(df['Time'])
    y = gb_list
    colors = gb_list    
    colors_list = list(colors)
    colors_list.sort()

    #### use scipy.stats to calculate IQR (Inter Quartile Range)
    IQR = iqr(colors_list, axis=0, nan_policy='omit')
    #### use numpy to calculate first quartile and third quartile
    Q1 = np.nanquantile(colors_list, 0.25)
    mask_min = Q1-1.5*IQR
    Q3 = np.nanquantile(colors_list, 0.75)
    mask_max = Q3+1.5*IQR

    plt.figure()
    plt.scatter(x, y, c=colors, cmap='RdYlBu_r')
    plt.colorbar()
    plt.title('Interaction Energy Profile', fontsize=16)
    plt.ylabel('MMGBSA (kcal/mol)', fontsize=14)
    plt.xlabel('Time (ns)', fontsize=14)
    plt.tight_layout()
    plt.savefig('MMGBSA/mmgbsa_vs_time.png',dpi=300)
    plt.close()
    
##################################################################################################################################################
##################################################################################################################################################

def InteractionEnergyNAMD(trj, traj_type, timestep, add, receptorSel, ligandSel, dry_complex_topology, namdPATH, stride):
    #### this function calculates interaction energy between receptor and ligand with NAMD
    logger.info('\nCalculating interaction energy (NAMD)...')

    output_basename = 'interactionEnergy' 
    vmdFile = 'interactionEnergy.tcl'

    with open(vmdFile, 'w') as f:
        f.write(f'''mol new {dry_complex_topology}
mol addfile {trj} type {traj_type}  first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
set prot [atomselect top "{receptorSel}"]
set ligand [atomselect top "{ligandSel}"]
global env
set Arch [vmdinfo arch]
set vmdEnv $env(VMDDIR)
puts $vmdEnv
source $vmdEnv/plugins/noarch/tcl/namdenergy1.4/namdenergy.tcl
namdenergy -exe {namdPATH} -elec -vdw -sel $ligand $prot -ofile "{output_basename}.dat" -tempname "{output_basename}_temp" -ts {add} -timemult {timestep} -stride {stride} -switch  7.5 -cutoff 9 -par {parameters}
quit''')

    os.system(f'vmd -dispdev text -e {vmdFile} > /dev/null 2>&1')

    os.system('mv interactionEnergy.dat InteractionEnergy/interactionEnergy.dat')

    os.remove(vmdFile)

def plotIntEnergyvsDistance(ref_bool):
    #### this function plots NAMD interaction energy as a function of the cdm distance
    logger.info('\nPlotting interaction energy vs cdm distance...')

    if ref_bool:
        txt = 'time_and_rmsd.txt'
    else:
        txt = 'time_and_space.txt'

    df = pd.read_csv(f'Geometry/{txt}', header=None, names=['Time','Distance'])
    df2 = pd.read_csv('InteractionEnergy/interactionEnergy.dat', sep='\s+')

    gb_list = list(df2['Total'])
    x = list(df['Distance'])
    y = gb_list
    colors = gb_list    
    colors_list = list(colors)
    colors_list.sort()

    #### use scipy.stats to calculate IQR (Inter Quartile Range)
    IQR = iqr(colors_list, axis=0, nan_policy='omit')
    #### use numpy to calculate first quartile and third quartile
    Q1 = np.nanquantile(colors_list, 0.02)
    mask_min = Q1-1.5*IQR
    Q3 = np.nanquantile(colors_list, 0.98)
    mask_max = Q3+1.5*IQR

    plt.figure()
    divnorm = mplcolors.TwoSlopeNorm(vmin=mask_min, vcenter=0, vmax=mask_max)
    plt.scatter(x, y, c=colors, cmap='RdYlBu_r', norm=divnorm)
    cbar = plt.colorbar()
    cbar.set_label('Interaction Energy (Kcal/mol)', rotation=270, labelpad=15)
    plt.title('Distance-dependent\nInteraction Energy Profile')
    plt.ylabel('Interaction Energy (kcal/mol)')

    if ref_bool:
        plt.xlabel('Ligand RMSD$_{backbone}$ to reference ($\AA$)')
    else:
        plt.xlabel('dcm$_{bs-lig}$ ($\AA$)')

    plt.tight_layout()
    plt.savefig('InteractionEnergy/intEnergy_vs_distance.png',dpi=300)
    plt.close()

def plotIntEnergyvsTime(ref_bool):
    #### this function plots NAMD interaction energy as a function of the cdm distance
    logger.info('\nPlotting interaction energy vs simulation time...')

    if ref_bool:
        txt = 'time_and_rmsd.txt'
    else:
        txt = 'time_and_space.txt'

    df = pd.read_csv(f'Geometry/{txt}', header=None, names=['Time','Distance'])
    df2 = pd.read_csv('InteractionEnergy/interactionEnergy.dat', sep='\s+')

    gb_list = list(df2['Total'])
    x = list(df['Time'])
    y = gb_list
    colors = gb_list    
    colors_list = list(colors)
    colors_list.sort()

    #### use scipy.stats to calculate IQR (Inter Quartile Range)
    IQR = iqr(colors_list, axis=0, nan_policy='omit')
    #### use numpy to calculate first quartile and third quartile
    Q1 = np.nanquantile(colors_list, 0.02)
    mask_min = Q1-1.5*IQR
    Q3 = np.nanquantile(colors_list, 0.98)
    mask_max = Q3+1.5*IQR

    plt.figure()
    divnorm = mplcolors.TwoSlopeNorm(vmin=mask_min, vcenter=0, vmax=mask_max)
    plt.scatter(x, y, c=colors, cmap='RdYlBu_r', norm=divnorm)
    cbar = plt.colorbar()
    cbar.set_label('Interaction Energy (Kcal/mol)', rotation=270, labelpad=15)
    plt.title('Time-dependent\nInteraction Energy Profile')
    plt.ylabel('Interaction Energy (kcal/mol)')
    plt.xlabel('Time (ns)')
    plt.tight_layout()
    plt.savefig('InteractionEnergy/intEnergy_vs_time.png',dpi=300)
    plt.close()

def plotInteractionEnergyLandscape(ref_bool):
    #### this function plots NAMD interaction energy as a function of the cdm distance/RMSD and of the simulation time
    logger.info('\nPlotting interaction energy landscape...')

    if ref_bool:
        txt = 'time_and_rmsd.txt'
    else:
        txt = 'time_and_space.txt'

    df = pd.read_csv(f'Geometry/{txt}', header=None, names=['Time','Distance'])
    df2 = pd.read_csv('InteractionEnergy/interactionEnergy.dat', sep='\s+')

    gb_list = list(df2['Total'])
    x = list(df['Time'])
    y = list(df['Distance'])
    colors = gb_list    
    colors_list = list(colors)
    colors_list.sort()

    #### use scipy.stats to calculate IQR (Inter Quartile Range)
    IQR = iqr(colors_list, axis=0, nan_policy='omit')
    #### use numpy to calculate first quartile and third quartile
    Q1 = np.nanquantile(colors_list, 0.02)
    mask_min = Q1-1.5*IQR
    Q3 = np.nanquantile(colors_list, 0.98)
    mask_max = Q3+1.5*IQR  

    plt.figure()
    divnorm = mplcolors.TwoSlopeNorm(vmin=mask_min, vcenter=0, vmax=mask_max)
    plt.scatter(x, y, c=colors, cmap='RdYlBu_r', norm=divnorm)
    cbar = plt.colorbar()
    cbar.set_label('Interaction Energy (Kcal/mol)', rotation=270, labelpad=15)
    plt.title('Interaction Energy Landscape')

    if ref_bool:
        plt.ylabel('Ligand RMSD$_{backbone}$ to reference ($\AA$)')
    else:
        plt.ylabel('dcm$_{bs-lig}$ ($\AA$)')

    plt.xlabel('Time (ns)')
    plt.tight_layout()
    plt.savefig('InteractionEnergy/intEnergyLandscape.png',dpi=300)
    plt.close()

def plotDynamicInteractionEnergyLandscape(ref_bool):
    #### this function plots NAMD interaction energy as a function of the cdm distance and of the simulation time
    logger.info('\nPlotting dynamic interaction energy landscape...')

    if not os.path.exists('InteractionEnergy/provideo'):
        os.makedirs('InteractionEnergy/provideo')
    if ref_bool:
        txt = 'time_and_rmsd.txt'
    else:
        txt = 'time_and_space.txt'

    df = pd.read_csv(f'Geometry/{txt}', header=None, names=['Time','Distance'])
    df2 = pd.read_csv('InteractionEnergy/interactionEnergy.dat', sep='\s+')

    gb_list = list(df2['Total'])
    x = list(df['Time'])
    y = list(df['Distance'])
    colors = gb_list    
    colors_list = list(colors)
    colors_list.sort()

    #### use scipy.stats to calculate IQR (Inter Quartile Range)
    IQR = iqr(colors_list, axis=0, nan_policy='omit')
    #### use numpy to calculate first quartile and third quartile
    Q1 = np.nanquantile(colors_list, 0.02)
    mask_min = Q1-1.5*IQR
    Q3 = np.nanquantile(colors_list, 0.98)
    mask_max = Q3+1.5*IQR

    divnorm = mplcolors.TwoSlopeNorm(vmin=mask_min, vcenter=0, vmax=mask_max)
    count = 0      

    for el in range(0, len(x)):
        plt.figure()
        plt.scatter(x[:count+1], y[:count+1], c=colors[:count+1], cmap='RdYlBu_r', norm=divnorm)
        cbar = plt.colorbar()
        cbar.set_label('Interaction Energy (Kcal/mol)', rotation=270, labelpad=15)
        plt.title('Interaction Energy Landscape')

        if ref_bool:
            plt.ylabel('Ligand RMSD$_{backbone}$ to reference ($\AA$)')
        else:
            plt.ylabel('dcm$_{bs-lig}$ ($\AA$)')

        plt.xlabel('Time (ns)')
        plt.xlim(min(x)-1, max(x)+1)
        plt.ylim(min(y)-1, max(y)+1)
        plt.tight_layout()
        plt.savefig(f'InteractionEnergy/provideo/intEnergyLandscape_{count}.png', dpi=300)
        plt.clf()
        plt.close()
        count += 1
        
##################################################################################################################################################
##################################################################################################################################################

#############################            
#######    receptor   #######
#############################

def getContactsRec(top, trj, receptorSel, ligandSel, distanceCutoff):
    #### this function parses the MD trajectory and returns the number of contacts between the ligand and each receptor residue
    logger.info('\nCalculating receptor-ligand contacts...')

    u = mda.Universe(top,trj)
    contactSel = f"({receptorSel}) and same residue as around {distanceCutoff} ({ligandSel})"
    contactsList = []

    #### iterate through each trajectory frame
    for ts in u.trajectory:
        #### create a ResidueGroup containing residues that are in contact with the ligand
        contacts = u.select_atoms(contactSel).residues
        contactsResidsList = [rg.resid for rg in contacts]
        contactsList.extend(contactsResidsList)

    #### create a sorted list of all residues in contact with the associated number of contacts
    countList = sorted(Counter(contactsList).items(), key = lambda kv: kv[1], reverse=True)

    #### write the output file and remove temporary one
    with open('per_residue_receptor/contactedResidues.txt','w') as f:
        for tuple in countList:
            out = str(tuple[0]) + ' ' + str(tuple[1]) +'\n'
            f.write(out)

def topContactsRec(numResidRec):
    #### this function defines the residues to be considered for the analysis
    logger.info('\nDefining most contacted receptor residues...')

    #### check if the number of residues to consider for the analysis is greater than the total number of contacted residues
    with open('per_residue_receptor/contactedResidues.txt','r') as f:
        count = 0
        for line in f:
            count +=1

    if count < numResidRec:
        numResidRec = count

    with open('per_residue_receptor/contactedResidues.txt','r') as r, open('per_residue_receptor/topContacts.txt','w') as a:
        for i in range(0, numResidRec):
            resid = r.readline().split()[0].strip()
            a.write(resid + ' ')
            
def getResnamesRec(top, trj, receptorSel):
    #### this function, required by Giovanni Bolcato, associate a residue name to each resid for the most contacted residues along MD trajectory
    logger.info('\nExtracting receptor residue names...')

    u = mda.Universe(top,trj)

    with open('per_residue_receptor/topContacts.txt', 'r') as f:
        resids = f.read().split()
        resnames = []

        for resid in resids:
            residue = u.select_atoms(f'{receptorSel} and resid {resid}').residues
            resname = str([rg.resname for rg in residue]).lstrip("\'\[").rstrip("\'\]")
            resnames.append(resname)

    with open('per_residue_receptor/residue_labels.txt','w') as f:
        for resname, resid in zip(resnames, resids):
            output = resname + ' ' + resid + '\n'
            f.write(output)

"""          
def interactionEnergyNAMDRec():
     #### this function calculates elec, vdw and total energy between ligand and each protein residue defined in the previous steps
     logger.info('\nCalculating receptor per-residue interaction energy with NAMD...\n')

     #### determine which residue to consider for the interaction energy analysis
     with open('per_residue_receptor/topContacts.txt','r') as f:
         resids = f.read().split()

     #### iter through each residue and calculate interaction energy with namdenergy plugin from vmd (namd2 exe required)
     for resid in resids:
         logger.info(f'resid {resid}')
         receptor_selection = f'{receptorSel} and resid {resid}'
         ligand_selection = f'{ligandSel}'
         output_basename = f'interactionEnergy_{resid}'         vmdFile = f'interactionEnergy_{resid}.tcl'

         with open(vmdFile, 'w') as f:
             f.write(f'''mol new {dry_complex_topology}
 mol addfile {trj} type {traj_type}  first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
 set prot [atomselect top "{receptor_selection}"]
 set ligand [atomselect top "{ligand_selection}"]
 global env
 set Arch [vmdinfo arch]
 set vmdEnv $env(VMDDIR)
 puts $vmdEnv
 source $vmdEnv/plugins/noarch/tcl/namdenergy1.4/namdenergy.tcl
 namdenergy -exe {namdPATH} -elec -vdw -sel $ligand $prot -ofile "{output_basename}.dat" -tempname "{output_basename}_temp" -ts {add} -timemult {timestep} -stride {stride} -switch  7.5 -cutoff 9 -par {parameters}
 quit''')

         os.system(f'vmd -dispdev text -e {vmdFile} > /dev/null 2>&1')

         os.system(f'mv interactionEnergy_{resid}.dat per_residue_receptor/')

         os.remove(vmdFile)
""" 

def interactionEnergyNAMDRec(n_procs):

    #### this function calculates elec, vdw and total energy between ligand and each protein residue defined in the previous steps
    logger.info('\nCalculating per-residue interaction energy with NAMD...\n')

    #### determine which residue to consider for the interaction energy analysis
    with open('per_residue_receptor/topContacts.txt','r') as f:
        resids = f.read().split()

    #### iter through each residue and calculate interaction energy with namdenergy plugin from vmd (namd2 exe required)
    tuple_resids = []
    for resid in resids:
        resid = (int(resid),)
        tuple_resids.append(resid)

    parallelizer.run(tuple_resids, calc_namd, n_procs, 'Calculating Interaction Energy Ligand')

def calc_namd(resid):
    receptor_selection = f'{ligandSel} and resid {resid}'
    ligand_selection = f'{receptorSel}'
    output_basename = f'interactionEnergy_{resid}'
    vmdFile = f'interactionEnergy_{resid}.tcl'

    with open(vmdFile, 'w') as f:
        f.write(f'''mol new {dry_complex_topology}
mol addfile {trj} type {traj_type}  first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
set prot [atomselect top "{receptor_selection}"]
set ligand [atomselect top "{ligand_selection}"]
global env
set Arch [vmdinfo arch]
set vmdEnv $env(VMDDIR)
puts $vmdEnv
source $vmdEnv/plugins/noarch/tcl/namdenergy1.4/namdenergy.tcl
namdenergy -exe {namdPATH} -elec -vdw -sel $ligand $prot -ofile "{output_basename}.dat" -tempname "{output_basename}_temp" -ts {add} -timemult {timestep} -stride {stride} -switch  7.5 -cutoff 9 -par {parameters} -updatesel
quit''')

    os.system(f'vmd -dispdev text -e {vmdFile} > /dev/null 2>&1')
    os.system(f'mv interactionEnergy_{resid}.dat per_residue_ligand/')
    os.remove(vmdFile)

""" 
def calc_namd(resid):
    receptor_selection = f'{receptorSel} and resid {resid}'
    ligand_selection = f'{ligandSel}'
    output_basename = f'interactionEnergy_{resid}'
    vmdFile = f'interactionEnergy_{resid}.tcl'

    with open(vmdFile, 'w') as f:
        f.write(f'''mol new {dry_complex_topology}
mol addfile {trj} type {traj_type}  first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
set prot [atomselect top "{receptor_selection}"]
set ligand [atomselect top "{ligand_selection}"]
global env
set Arch [vmdinfo arch]
set vmdEnv $env(VMDDIR)
puts $vmdEnv
source $vmdEnv/plugins/noarch/tcl/namdenergy1.4/namdenergy.tcl
namdenergy -exe {namdPATH} -elec -vdw -sel $ligand $prot -ofile "{output_basename}.dat" -tempname "{output_basename}_temp" -ts {add} -timemult {timestep} -stride {stride} -switch  7.5 -cutoff 9 -par {parameters} -updatesel
quit''')
    os.system(f'vmd -dispdev text -e {vmdFile} > vmd_log.txt') # With output from vmd
    # os.system(f'vmd -dispdev text -e {vmdFile} > /dev/null 2>&1') # No output from vmd
    os.system(f'mv interactionEnergy_{resid}.dat per_residue_ligand/')
    os.remove(vmdFile)
"""      

def HeatmapProcessingRec():
    #### this function combines the single namd output files into single csv that can be manipulated for plotting purpose
    logger.info('\nProcessing receptor output files...')

    #### create a list of all interaction energy files ordered by ascending resid number
    intFile = sorted(glob('per_residue_receptor/interactionEnergy_*.dat'), key=lambda x: int(os.path.basename(x).split('_')[1].rstrip('.dat')))

    #### create a pandas dataframe to manipulate data from the interaction energy file
    df = pd.read_csv(intFile[0], sep='\s+')
    resid = os.path.basename(intFile[0]).split('_')[1].rstrip('.dat')

    elec = df['Elec'].rename(resid)
    vdw = df['VdW'].rename(resid)
    total = df['Total'].rename(resid)

    #### iterate through each interaction energy file in order to concatenate it to the first one and create a single csv file for each replica
    for i in range(1,len(intFile)):
        resid1 = os.path.basename(intFile[i]).split('_')[1].rstrip('.dat')
        df1 = pd.read_csv(intFile[i], sep='\s+')

        elec1 = df1['Elec'].rename(resid1)
        vdw1 = df1['VdW'].rename(resid1)
        total1 = df1['Total'].rename(resid1)
        elec = pd.concat([elec, elec1], axis=1)
        vdw = pd.concat([vdw, vdw1], axis=1)
        total = pd.concat([total, total1], axis=1)

#    elec.to_csv('per_residue_receptor/elecPerResidue.csv')
#    vdw.to_csv('per_residue_receptor/vdwPerResidue.csv')
    total.to_csv('per_residue_receptor/totalPerResidue.csv')
    
def determineHeatmapMaskValuesRec():
    #### this function defines the max and min values to use as mask values for the heatmap generation
    logger.info('\nDetermining mask values for receptor heatmap generation...')
    
    #### initialize the list that will store each value from each trajectory
    arrList = []
    df = pd.read_csv('per_residue_receptor/totalPerResidue.csv', index_col = 0).T

    for arr in df.values:
    #### convert each numpy array (each dataframe row) into a list and append it to the arrList
        arr2list = list(arr)
        arrList.extend(arr2list)    

    #### sort arrList in order to define quartile and IQR values            
    arrList.sort()

    #### use scipy.stats to calculate IQR (Inter Quartile Range)
    IQR = iqr(arrList, axis=0, nan_policy='omit')
    #### use numpy to calculate first quartile and third quartile
    Q1 = np.nanquantile(arrList, 0.25)
    Q3 = np.nanquantile(arrList, 0.75)
    #### use IQR to determine mask values for the heatmap
    mask_min = Q1-1.5*IQR
    mask_max = Q3+1.5*IQR

    return mask_min, mask_max
    
def plotHeatmapRec(mask_min, mask_max, numShiftRec, timestep, stride):
    #### this function creates the heatmap for each replica using the mask values previously calculated
    logger.info('\nPlotting Receptor Per-Residue Interaction Energy Heatmap...')

    basename = os.path.basename(os.getcwd())
    df = pd.read_csv('per_residue_receptor/totalPerResidue.csv', index_col = 0).T

    #### use the temp df to create a list of ns and correct resids to use as columns and indexes for the definitive df
    ns = []
    resid_correct = []

    for i in df.columns:
        ns.append((i)*add/1000000)

    with open('per_residue_receptor/residue_labels.txt', 'r') as f:
        lines = f.readlines()
        for i in df.index:
            for line in lines:
                if line.endswith(' ' + str(i) + '\n'):
                    resname = line.rstrip('\n').split(' ')[0]
            new_label = resname + ' ' + str(int(i)+numShiftRec)
            resid_correct.append(new_label)

    #### create definitive df
    df2 = pd.DataFrame(df.values, columns = ns, index = resid_correct)

    #### create heatmap
    fig = plt.figure(facecolor='white')
    xtl = int(5 * 1000000 / (timestep*stride))
    ax = sns.heatmap(df2, cmap='RdBu_r', center = 0, vmin = mask_min, vmax = mask_max, yticklabels=1, xticklabels=xtl)
    cbar_axes = ax.figure.axes[-1]
    cbar_axes.set_ylabel('Interaction Energy (Kcal/mol)', rotation=270, labelpad=15)
    ax.tick_params(axis='both', which='minor', labelsize=6)
    ax.set_title('Per-Residue Interaction Energy (Receptor)')
    ax.set_ylabel('Residue')
    ax.set_xlabel('Time (ns)')
    plt.tight_layout()
    fig.savefig('per_residue_receptor/receptor_heatmap.png', dpi=300)
    plt.close()

    return df2
    
def plotDynamicHeatmapRec(df2, mask_min, mask_max, timestep, stride):
    #### this function creates the heatmap for each replica using the mask values previously calculated
    logger.info('\nPlotting Dynamic Receptor Per-Residue Interaction Energy Heatmap...')

    if not os.path.exists('per_residue_receptor/provideo'):    
        os.makedirs('per_residue_receptor/provideo')

    columns = df2.columns
    count = 0

    for col in columns:
        columns_list = columns[0:count+1]
        mask = df2.isin(df2[columns_list])

        #### create heatmap
        basename = os.path.basename(os.getcwd())
        fig = plt.figure(facecolor='white')
        xtl = int(5 * 1000000 / (timestep*stride))
        ax = sns.heatmap(df2, cmap='RdBu_r', center = 0, mask = ~mask, vmin = mask_min, vmax = mask_max, yticklabels=1, xticklabels=xtl)
        cbar_axes = ax.figure.axes[-1]
        cbar_axes.set_ylabel('Interaction Energy (Kcal/mol)', rotation=270, labelpad=15)
        ax.tick_params(axis='both', which='minor', labelsize=6)
        ax.set_title('Per-Residue Interaction Energy (Receptor)')
        ax.set_ylabel('Residue')
        ax.set_xlabel('Time (ns)')
        plt.tight_layout()
        fig.savefig(f'per_residue_receptor/provideo/receptor_heatmap_{count}.png', dpi=300)
        plt.clf()
        plt.close()

        count +=1
        
#############################            
#######     ligand    #######
#############################
            
def getContactsLig(top, trj, ligandSel, receptorSel, distanceCutoff):

    #### this function parses the MD trajectory and returns the number of contacts between the receptor and each ligand residue
    logger.info('\nCalculating ligand-protein contacts...')

    u = mda.Universe(top,trj)
    contactSel = f"({ligandSel}) and same residue as around {distanceCutoff} ({receptorSel})"
    contactsList = []

    #### iterate through each trajectory frame
    for ts in u.trajectory:
        #### create a ResidueGroup containing residues that are in contact with the ligand
        contacts = u.select_atoms(contactSel).residues
        contactsResidsList = [rg.resid for rg in contacts]
        contactsList.extend(contactsResidsList)

    #### create a sorted list of all residues in contact with the associated number of contacts
    countList = sorted(Counter(contactsList).items(), key = lambda kv: kv[1], reverse=True)

    #### write the output file and remove temporary one
    with open('per_residue_ligand/contactedResidues.txt','w') as f:
        for tuple in countList:
            out = str(tuple[0]) + ' ' + str(tuple[1]) +'\n'
            f.write(out)
            
def topContactsLig(numResidLig):
    #### this function defines the residues to be considered for the analysis
    logger.info('\nDefining most contacted ligand residues...')

    #### check if the number of residues to consider for the analysis is greater than the total number of contacted residues
    with open('per_residue_ligand/contactedResidues.txt','r') as f:
        count = 0
        for line in f:
            count +=1

    if count < numResidLig:
        numResidLig = count

    with open('per_residue_ligand/contactedResidues.txt','r') as r, open('per_residue_ligand/topContacts.txt','w') as a:
        for i in range(0, numResidLig):
            resid = r.readline().split()[0].strip()
            a.write(resid + ' ')
            
def getResnamesLig(top, trj):
    #### this function, required by Giovanni Bolcato, associate a residue name to each resid for the most contacted residues along MD trajectory
    logger.info('\nExtracting ligand residue names...')

    u = mda.Universe(top,trj)
    with open('per_residue_ligand/topContacts.txt', 'r') as f:
        resids = f.read().split()
        resnames = []

        for resid in resids:
            residue = u.select_atoms(f'resid {resid}').residues
            resname = str([rg.resname for rg in residue]).lstrip("\'\[").rstrip("\'\]")
            resnames.append(resname)

    with open('per_residue_ligand/residue_labels.txt','w') as f:
        for resname, resid in zip(resnames, resids):
            output = resname + ' ' + resid + '\n'
            f.write(output)

"""       
def interactionEnergyNAMDLig():
    #### this function calculates elec, vdw and total energy between ligand and each protein residue defined in the previous steps
   logger.info('\nCalculating ligand per-residue interaction energy with NAMD...\n')

    #### determine which residue to consider for the interaction energy analysis
    with open('per_residue_ligand/topContacts.txt','r') as f:
        resids = f.read().split()

    #### iter through each residue and calculate interaction energy with namdenergy plugin from vmd (namd2 exe required)
    for resid in resids:
        logger.info(f'resid {resid}')
        
        receptor_selection = f'{ligandSel} and resid {resid}'
        ligand_selection = f'{receptorSel}'
        output_basename = f'interactionEnergy_{resid}'
        vmdFile = f'interactionEnergy_{resid}.tcl'

        with open(vmdFile, 'w') as f:
            f.write(f'''mol new {dry_complex_topology}
mol addfile {trj} type {traj_type}  first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
set prot [atomselect top "{receptor_selection}"]
set ligand [atomselect top "{ligand_selection}"]
global env
set Arch [vmdinfo arch]
set vmdEnv $env(VMDDIR)
puts $vmdEnv
source $vmdEnv/plugins/noarch/tcl/namdenergy1.4/namdenergy.tcl
namdenergy -exe {namdPATH} -elec -vdw -sel $ligand $prot -ofile "{output_basename}.dat" -tempname "{output_basename}_temp" -ts {add} -timemult {timestep} -stride {stride} -switch  7.5 -cutoff 9 -par {parameters}
quit''')

        os.system(f'vmd -dispdev text -e {vmdFile} > /dev/null 2>&1')

        os.system(f'mv interactionEnergy_{resid}.dat per_residue_ligand/')

        os.remove(vmdFile)    
""" 

def interactionEnergyNAMDLig(n_procs):

    #### this function calculates elec, vdw and total energy between ligand and each protein residue defined in the previous steps
    logger.info('\nCalculating per-residue interaction energy with NAMD...\n')

    #### determine which residue to consider for the interaction energy analysis
    with open('per_residue_receptor/topContacts.txt','r') as f:
        resids = f.read().split()

    #### iter through each residue and calculate interaction energy with namdenergy plugin from vmd (namd2 exe required)
    tuple_resids = []
    for resid in resids:
        resid = (int(resid),)
        tuple_resids.append(resid)

    parallelizer.run(tuple_resids, calc_namd, n_procs, 'Calculating Interaction Energy Ligand')

def HeatmapProcessingLig():
    #### this function combines the single namd output files into single csv that can be manipulated for plotting purpose
    logger.info('\nProcessing ligand output files...')

    #### create a list of all interaction energy files ordered by ascending resid number
    intFile = sorted(glob('per_residue_ligand/interactionEnergy_*.dat'), key=lambda x: int(os.path.basename(x).split('_')[1].rstrip('.dat')))

    #### create a pandas dataframe to manipulate data from the interaction energy file
    df = pd.read_csv(intFile[0], sep='\s+')
    resid = os.path.basename(intFile[0]).split('_')[1].rstrip('.dat')
    elec = df['Elec'].rename(resid)
    vdw = df['VdW'].rename(resid)
    total = df['Total'].rename(resid)

    #### iterate through each interaction energy file in order to concatenate it to the first one and create a single csv file for each replica
    for i in range(1,len(intFile)):
        resid1 = os.path.basename(intFile[i]).split('_')[1].rstrip('.dat')
        df1 = pd.read_csv(intFile[i], sep='\s+')

        elec1 = df1['Elec'].rename(resid1)
        vdw1 = df1['VdW'].rename(resid1)
        total1 = df1['Total'].rename(resid1)
        elec = pd.concat([elec, elec1], axis=1)
        vdw = pd.concat([vdw, vdw1], axis=1)
        total = pd.concat([total, total1], axis=1)
        
    elec.to_csv('per_residue_ligand/elecPerResidue.csv')
    vdw.to_csv('per_residue_ligand/vdwPerResidue.csv')
    total.to_csv('per_residue_ligand/totalPerResidue.csv')
    
def determineHeatmapMaskValuesLig():
    #### this function defines the max and min values to use as mask values for the heatmap generation
    logger.info('\nDetermining mask values for ligand heatmap generation...')

    #### initialize the list that will store each value from each trajectory
    arrList = []
    df = pd.read_csv('per_residue_ligand/totalPerResidue.csv', index_col = 0).T
    for arr in df.values:

    #### convert each numpy array (each dataframe row) into a list and append it to the arrList
        arr2list = list(arr)
        arrList.extend(arr2list) 

    #### sort arrList in order to define quartile and IQR values            
    arrList.sort()

    #### use scipy.stats to calculate IQR (Inter Quartile Range)
    IQR = iqr(arrList, axis=0, nan_policy='omit')
    #### use numpy to calculate first quartile and third quartile
    Q1 = np.nanquantile(arrList, 0.25)
    Q3 = np.nanquantile(arrList, 0.75)
    #### use IQR to determine mask values for the heatmap
    mask_min = Q1-1.5*IQR
    mask_max = Q3+1.5*IQR

    return mask_min, mask_max

def plotHeatmapLig(mask_min, mask_max, numShiftLig, add):
    #### this function creates the heatmap for each replica using the mask values previously calculated
    logger.info('\nPlotting Ligand Per-Residue Interaction Energy Heatmap...')

    basename = os.path.basename(os.getcwd())
    df = pd.read_csv('per_residue_ligand/totalPerResidue.csv', index_col = 0).T

    #### use the temp df to create a list of ns and correct resids to use as columns and indexes for the definitive df
    ns = []
    resid_correct = []

    for i in df.columns:
        ns.append((i)*add/1000000)

    with open('per_residue_ligand/residue_labels.txt', 'r') as f:
        lines = f.readlines()
        for i in df.index:
            for line in lines:
                if line.endswith(' ' + str(i) + '\n'):
                    resname = line.rstrip('\n').split(' ')[0]
            new_label = resname + ' ' + str(int(i)+numShiftLig)
            resid_correct.append(new_label)

    #### create definitive df
    df2 = pd.DataFrame(df.values, columns = ns, index = resid_correct)

    #### create heatmap
    fig = plt.figure(facecolor='white')
    xtl = int(5 * 1000000 / (timestep*stride))
    ax = sns.heatmap(df2, cmap='RdBu_r', center = 0, vmin = mask_min, vmax = mask_max, yticklabels=1, xticklabels=xtl)
    cbar_axes = ax.figure.axes[-1]
    cbar_axes.set_ylabel('Interaction Energy (Kcal/mol)', rotation=270, labelpad=15)
    ax.tick_params(axis='both', which='minor', labelsize=6)
    ax.set_title('Per-Residue Interaction Energy (Ligand)')
    ax.set_ylabel('Residue')
    ax.set_xlabel('Time (ns)')
    plt.tight_layout()
    fig.savefig('per_residue_ligand/ligand_heatmap.png', dpi=300)
    plt.close()

    return df2

def plotDynamicHeatmapLig(df2, mask_min, mask_max):
    #### this function creates the heatmap for each replica using the mask values previously calculated
    logger.info('\nPlotting Dynamic Ligand Per-Residue Interaction Energy Heatmap...')

    if not os.path.exists('per_residue_ligand/provideo'):    
        os.makedirs('per_residue_ligand/provideo')

    columns = df2.columns
    count = 0

    for col in columns:
        columns_list = columns[0:count+1]
        mask = df2.isin(df2[columns_list])

        #### create heatmap
        fig = plt.figure(facecolor='white')
        xtl = int(5 * 1000000 / (timestep*stride))
        ax = sns.heatmap(df2, cmap='RdBu_r', center = 0, mask=~mask, vmin = mask_min, vmax = mask_max, yticklabels=1, xticklabels=xtl)
        cbar_axes = ax.figure.axes[-1]
        cbar_axes.set_ylabel('Interaction Energy (Kcal/mol)', rotation=270, labelpad=15)
        ax.tick_params(axis='both', which='minor', labelsize=6)
        ax.set_title('Per-Residue Interaction Energy (Ligand)')
        ax.set_ylabel('Residue')
        ax.set_xlabel('Time (ns)')
        plt.tight_layout()
        fig.savefig(f'per_residue_ligand/provideo/ligand_heatmap_{count}.png', dpi=300)
        plt.clf()
        plt.close()
        count +=1

#############################            
#### interaction matrix  ####
#############################

def interactionEnergyNAMDPairwiseMatrix(top, dry_complex_topology, ligandSel, receptorSel, traj_type, namdPATH, parameters):
    #### this function calculates a pairwise per-residue elec, vdw and total interaction energy value between receptor and ligand with NAMD
    logger.info('\nCalculating pairwise per-residue interaction energy with NAMD...\n')

    u = mda.Universe(top)

    with open('per_residue_ligand/topContacts.txt','r') as f:
        ligand_resids = f.read().split()

    with open('per_residue_receptor/topContacts.txt','r') as f:
        resids = f.read().split()

    for lr in ligand_resids:
        logger.info('#### ' + lr)
        ligSel = f'({ligandSel}) and resid {lr}'

        for r in resids:
            logger.info(r)
            output_basename = f'interactionEnergy_{lr}_{r}' %(lr, r)
            protSel = f'({receptorSel}) and resid {r}'
            vmdFile = f'interactionEnergy_{lr}_{r}.tcl'

            with open(vmdFile, 'w') as f:
                f.write(f'''mol new {dry_complex_topology}
mol addfile {trj} type {traj_type}  first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
set prot [atomselect top "{protSel}"]
set ligand [atomselect top "{ligSel}"]
global env
set Arch [vmdinfo arch]
set vmdEnv $env(VMDDIR)
puts $vmdEnv
source $vmdEnv/plugins/noarch/tcl/namdenergy1.4/namdenergy.tcl
namdenergy -exe {namdPATH} -elec -vdw -sel $ligand $prot -ofile "{output_basename}.dat" -tempname "{output_basename}_temp" -ts {add} -timemult {timestep} -stride {stride} -switch  7.5 -cutoff 9 -par {parameters}
quit''')

            os.system(f'vmd -dispdev text -e {vmdFile} > /dev/null 2>&1')
            os.system(f'mv interactionEnergy_{lr}_{r}.dat per_residue_matrix/')
            os.remove(vmdFile)

def HeatmapProcessingPairwiseMatrix():
    #### this function combines the single namd output files into single csv that can be manipulated for plotting purpose
    logger.info('\nProcessing pairwise matrix output files...')

    resids = []
    intFileList = []

    #### create a 2d list containing each atom name-protein residue interaction energy file
    with open('per_residue_ligand/topContacts.txt','r') as f:
        ligand_resids = sorted(f.read().split(), key=lambda x: int(x))

    for lr in ligand_resids:
        intFile = sorted(glob(f'per_residue_matrix/interactionEnergy_{lr}_*.dat'), key=lambda x: int(os.path.basename(x).split('_')[2].rstrip('.dat')))
        intFileList.append(intFile)

    #### start with the first ligand residue
    df = pd.read_csv(intFileList[0][0], sep='\s+')

    resid = os.path.basename(intFileList[0][0]).split('_')[2].rstrip('.dat')
    lr = os.path.basename(intFileList[0][0]).split('_')[1]

    resids.append(resid)

    elec = df['Elec'].mean()
    vdw = df['VdW'].mean()
    total = df['Total'].mean()

    elecList = [elec]
    vdwList = [vdw]
    totalList = [total]

    for i in range(1, len(intFile)):
        df = pd.read_csv(intFileList[0][i], sep='\s+')
        resid = os.path.basename(intFileList[0][i]).split('_')[2].rstrip('.dat')

        resids.append(resid)
        elec = df['Elec'].mean()
        vdw = df['VdW'].mean()
        total = df['Total'].mean()
        elecList.append(elec)
        vdwList.append(vdw)
        totalList.append(total)

    elecSeries = pd.Series(elecList, index=resids, name=lr)
    elecDf = pd.DataFrame(elecSeries)
    vdwSeries = pd.Series(vdwList, index=resids, name=lr)
    vdwDf = pd.DataFrame(vdwSeries)
    totalSeries = pd.Series(totalList, index=resids, name=lr)
    totalDf = pd.DataFrame(totalSeries)

    for row in range(1, len(intFileList)):
        df = pd.read_csv(intFileList[row][0], sep='\s+')
        resid = os.path.basename(intFileList[row][0]).split('_')[2].rstrip('.dat')
        lr = os.path.basename(intFileList[row][0]).split('_')[1]
        elec = df['Elec'].mean()
        vdw = df['VdW'].mean()
        total = df['Total'].mean()
        elecList = [elec]
        vdwList = [vdw]
        totalList = [total]

        for col in range(1, len(intFileList[row])):
            df = pd.read_csv(intFileList[row][col], sep='\s+')
            elec = df['Elec'].mean()
            vdw = df['VdW'].mean()
            total = df['Total'].mean()
            elecList.append(elec)
            vdwList.append(vdw)
            totalList.append(total)

        elecSeries1 = pd.Series(elecList, index=resids, name=lr)
        elecDf1 = pd.DataFrame(elecSeries1)
        vdwSeries1 = pd.Series(vdwList, index=resids, name=lr)
        vdwDf1 = pd.DataFrame(vdwSeries1)
        totalSeries1 = pd.Series(totalList, index=resids, name=lr)
        totalDf1 = pd.DataFrame(totalSeries1)
        elecDf = pd.concat([elecDf, elecDf1], axis=1)
        vdwDf = pd.concat([vdwDf, vdwDf1], axis=1)
        totalDf = pd.concat([totalDf, totalDf1], axis=1)

#    elecDf.to_csv('per_residue_matrix/elec.csv')
#    vdwDf.to_csv('per_residue_matrix/vdw.csv')
    totalDf.to_csv('per_residue_matrix/total.csv')

def determineHeatmapMaskValuesPairwiseMatrix():
    #### this function defines the max and min values to use as mask values for the heatmap generation
    logger.info('\nDetermining mask values for pairwise matrix heatmap generation...')

    df = pd.read_csv('per_residue_matrix/total.csv', index_col = 0).T
    arrList = []

    for arr in df.values:
        #### convert each numpy array (each dataframe row) into a list and append it to the arrList
        arr2list = list(arr)
        arrList.extend(arr2list)

    #### sort arrList in order to define quartile and IQR values            
    arrList.sort()

    #### use scipy.stats to calculate IQR (Inter Quartile Range)
    IQR = iqr(arrList, axis=0)
    #### use numpy to calculate first quartile and third quartile
    Q1 = np.quantile(arrList, 0.02)
    Q3 = np.quantile(arrList, 0.98)
    #### use IQR to determine mask values for the heatmap
    mask_min = Q1-1.5*IQR
    mask_max = Q3+1.5*IQR

    return mask_min, mask_max

def plotHeatmapPairwiseMatrix(mask_min, mask_max, numShiftRec, numShiftLig):
    #### this function creates the heatmap for each replica using the mask values previously calculated
    logger.info('\nPlotting Pairwise Per-Residue Interaction Energy Heatmap...')

    df = pd.read_csv('per_residue_matrix/total.csv', index_col = 0).T

    #### use the temp df to create a list of ns and correct resids to use as columns and indexes for the definitive df
    lr_correct = []
    resid_correct = []

    with open('per_residue_receptor/residue_labels.txt', 'r') as f:
        lines = f.readlines()
        for i in df.columns:
            for line in lines:
                if line.endswith(' ' + str(i) + '\n'):
                    resname = line.rstrip('\n').split(' ')[0]
            new_label = resname + ' ' + str(int(i)+numShiftRec)
            resid_correct.append(new_label)

    with open('per_residue_ligand/residue_labels.txt', 'r') as f:
        lines = f.readlines()
        for i in df.index:
            for line in lines:
                if line.endswith(' ' + str(i) + '\n'):
                    resname = line.rstrip('\n').split(' ')[0]
            new_label = resname + ' ' + str(int(i)+numShiftLig)
            lr_correct.append(new_label)

    #### create definitive df
    df2 = pd.DataFrame(df.values, columns = resid_correct, index = lr_correct)  

    #### create heatmap
    # build the figure instance with the desired height
    fig, ax = plt.subplots(facecolor='white')
    ax = sns.heatmap(df2, cmap='RdBu_r', center = 0, vmin=mask_min, vmax=mask_max, yticklabels=1, xticklabels=1)
    cbar_axes = ax.figure.axes[-1]
    cbar_axes.set_ylabel('Interaction Energy (Kcal/mol)', rotation=270, labelpad=15)
    ax.tick_params(labelsize=8, axis='x', rotation=45)
    ax.set_title('Pairwise Interaction Matrix')
    ax.set_ylabel('Ligand Residue')
    ax.set_xlabel('Receptor Residue')
    plt.tight_layout()
    fig.savefig('per_residue_matrix/matrix.png', dpi=300)

def DynamicHeatmapProcessingPairwiseMatrix():
    #### this function reprocess the single namd output files into single csv that can be manipulated for plotting purpose
    logger.info('\nRe-processing pairwise matrix output files...')

    if not os.path.exists('per_residue_matrix/provideo'):
        os.makedirs('per_residue_matrix/provideo')

    intFileList = []
    #### create a 2d list containing each atom name-protein residue interaction energy file
    with open('per_residue_ligand/topContacts.txt','r') as f:
        ligand_resids = sorted(f.read().split(), key=lambda x: int(x))

    for lr in ligand_resids:
        intFile = sorted(glob(f'per_residue_matrix/interactionEnergy_{lr}_*.dat'), key=lambda x: int(os.path.basename(x).split('_')[2].rstrip('.dat')))
        intFileList.append(intFile)

    u = mda.Universe(top,trj)
    frame_count = 0
    trj_lenght = len(u.trajectory)

    for i in range(0, trj_lenght):
        resids = []

        #### start with the first ligand residue
        df = pd.read_csv(intFileList[0][0], sep='\s+')
        resid = os.path.basename(intFileList[0][0]).split('_')[2].rstrip('.dat')
        lr = os.path.basename(intFileList[0][0]).split('_')[1]
        resids.append(resid)
        elec = df['Elec'][frame_count]
        vdw = df['VdW'][frame_count]
        total = df['Total'][frame_count]
        elecList = [elec]
        vdwList = [vdw]
        totalList = [total]

        for i in range(1, len(intFile)):
            df = pd.read_csv(intFileList[0][i], sep='\s+')
            resid = os.path.basename(intFileList[0][i]).split('_')[2].rstrip('.dat')
            resids.append(resid)
            elec = df['Elec'][frame_count]
            vdw = df['VdW'][frame_count]
            total = df['Total'][frame_count]
            elecList.append(elec)
            vdwList.append(vdw)
            totalList.append(total)

        elecSeries = pd.Series(elecList, index=resids, name=lr)
        elecDf = pd.DataFrame(elecSeries)
        vdwSeries = pd.Series(vdwList, index=resids, name=lr)
        vdwDf = pd.DataFrame(vdwSeries)
        totalSeries = pd.Series(totalList, index=resids, name=lr)
        totalDf = pd.DataFrame(totalSeries)

        for row in range(1, len(intFileList)):
            df = pd.read_csv(intFileList[row][0], sep='\s+')
            resid = os.path.basename(intFileList[row][0]).split('_')[2].rstrip('.dat')
            lr = os.path.basename(intFileList[row][0]).split('_')[1]
            elec = df['Elec'][frame_count]
            vdw = df['VdW'][frame_count]
            total = df['Total'][frame_count]
            elecList = [elec]
            vdwList = [vdw]
            totalList = [total]

            for col in range(1, len(intFileList[row])):
                df = pd.read_csv(intFileList[row][col], sep='\s+')
                elec = df['Elec'][frame_count]
                vdw = df['VdW'][frame_count]
                total = df['Total'][frame_count]
                elecList.append(elec)
                vdwList.append(vdw)
                totalList.append(total)

            elecSeries1 = pd.Series(elecList, index=resids, name=lr)
            elecDf1 = pd.DataFrame(elecSeries1)
            vdwSeries1 = pd.Series(vdwList, index=resids, name=lr)
            vdwDf1 = pd.DataFrame(vdwSeries1)
            totalSeries1 = pd.Series(totalList, index=resids, name=lr)
            totalDf1 = pd.DataFrame(totalSeries1)
            elecDf = pd.concat([elecDf, elecDf1], axis=1)
            vdwDf = pd.concat([vdwDf, vdwDf1], axis=1)
            totalDf = pd.concat([totalDf, totalDf1], axis=1)

#         elecDf.to_csv('per_residue_matrix/elec.csv')
#         vdwDf.to_csv('per_residue_matrix/vdw.csv')
        totalDf.to_csv(f'per_residue_matrix/provideo/total_{frame_count}.csv')
        frame_count +=1

def plotDynamicHeatmapPairwiseMatrix(mask_min, mask_max, numShiftRec, numShiftLig):
    #### this function creates the heatmap for each replica using the mask values previously calculated
    logger.info('\nPlotting Dynamic Pairwise Per-Residue Interaction Energy Heatmap...')

    #### use the temp df to create a list of ns and correct resids to use as columns and indexes for the definitive df
    ctrl = True

    for file in sorted(glob('per_residue_matrix/provideo/total_*.csv'),key=lambda x: int(os.path.basename(x).split('_')[1].rstrip('.csv'))):
        df = pd.read_csv(file, index_col = 0).T

        if ctrl == True:
            lr_correct = []
            resid_correct = []

            with open('per_residue_receptor/residue_labels.txt', 'r') as f:
                lines = f.readlines()
                for i in df.columns:
                    for line in lines:
                        if line.endswith(' ' + str(i) + '\n'):
                            resname = line.rstrip('\n').split(' ')[0]
                    new_label = resname + ' ' + str(int(i)+numShiftRec)
                    resid_correct.append(new_label)

            with open('per_residue_ligand/residue_labels.txt', 'r') as f:
                lines = f.readlines()
                for i in df.index:
                    for line in lines:
                        if line.endswith(' ' + str(i) + '\n'):
                            resname = line.rstrip('\n').split(' ')[0]
                    new_label = resname + ' ' + str(int(i)+numShiftLig)
                    lr_correct.append(new_label)

            ctrl = False

        numFrame = os.path.basename(file).split('_')[1].rstrip('.csv')
        simTime = "{:.2f}".format(int(numFrame) * add / 1000000)

        #### create definitive df
        df2 = pd.DataFrame(df.values, columns = resid_correct, index = lr_correct) 

        #### create heatmap
        fig, ax = plt.subplots(facecolor='white')
        ax = sns.heatmap(df2, cmap='RdBu_r', center = 0, vmin=mask_min, vmax=mask_max, yticklabels=1, xticklabels=1)
        cbar_axes = ax.figure.axes[-1]
        cbar_axes.set_ylabel('Interaction Energy (Kcal/mol)', rotation=270, labelpad=15)
        ax.tick_params(labelsize=8, axis='x', rotation=45)
        ax.set_title(f'Pairwise Interaction Matrix\n(Time = {simTime} ns)')
        ax.set_ylabel('Ligand Residue')
        ax.set_xlabel('Receptor Residue')
        plt.tight_layout()
        fig.savefig(f'per_residue_matrix/provideo/pairwise_matrix_{numFrame}.png', dpi=300) 
        ax.clear()
        plt.clf()
        
def mountEnergyPanel():
    #### this function mounts a 4 tile panel for energetic analysis
    logger.info('\nMounting energetic 4 tile panel...')

    a = 'InteractionEnergy/intEnergyLandscape.png'
    b = 'per_residue_matrix/matrix.png'
    c = 'per_residue_receptor/receptor_heatmap.png'
    d = 'per_residue_ligand/ligand_heatmap.png'

    os.system(f'montage -tile 2x2 -geometry 1920x1080 {a} {b} {c} {d} InteractionEnergy/merged_energy.png') 

def mountVideo():
    #### this function mount the final video for suMD analysis
    logger.info('\nMounting final video...')
    
    frames = os.path.abspath("frame/")
    iel = os.path.abspath("InteractionEnergy/provideo")
    perresidue = os.path.abspath("per_residue_receptor/provideo")
    matrix = os.path.abspath("per_residue_matrix/provideo")

    frames_list = sorted(glob(os.path.join(frames, "*.tga")), key=lambda x: int(x.split('.')[1].rstrip('.tga')))
    iel_list = sorted(glob(os.path.join(iel,"*.png")), key=lambda x: int(os.path.basename(x).split('_')[1].rstrip('.png')))
    perresidue_list = sorted(glob(os.path.join(perresidue,"*.png")), key=lambda x: int(os.path.basename(x).split('_')[2].rstrip('.png')))
    matrix_list = sorted(glob(os.path.join(matrix,"*.png")), key=lambda x: int(os.path.basename(x).split('_')[2].rstrip('.png')))

#    if not os.path.exists('video'):
#        os.makedirs('video')

    for a,b,c,d in zip(frames_list, iel_list, perresidue_list, matrix_list):
        output_num = os.path.basename(a).split('.')[1].rstrip('.tga')
        logger.info(output_num)
        os.system(f'montage -tile 2x2 -geometry 1920x1080 {a} {b} {c} {d} video/merged_{output_num}.tga')

    os.system('avconv -f image2 -i video/merged_%04d.tga -vb 10M -vcodec h264 video/output.mp4')
    #os.system('avconv -y -i video/output.mp4 -i ~/Pictures/Logo_MMS.png -filter_complex "overlay=x=(main_w-overlay_w)/2:y=(main_h-overlay_h)/2" video/video_logo.mp4')

def get_parameters(configuration_path):
    """
    Function to read global parameters from configuration file. Exits if configuration file is not found.
    """
    
    # Read global parameters from configuration file
    if os.path.isfile(configuration_path):
        with open(configuration_path) as config_file:
            global_parameters = yaml.load(config_file, Loader = yaml.FullLoader)
    else:
        logger.error("Configuration file not found")
        sys.exit()
    
    return global_parameters

def check_input_files(*files):
    """
    Function that checks the existence of input files one by one. Exits if any of the files is not found 
    and prints the name of the missing file.
    """
    
    for file in files:

        if file is None:
            logger.error(f"File {file} is not given, added to input.yml!")

        if not os.path.isfile(file):
            logger.error(f"File {file} not found!")
            sys.exit()


if __name__ == '__main__':

    parser = argparse.ArgumentParser("RNASuMDAnalyzer")

    # Add task argument
    parser.add_argument("task", help="Task to be performed", choices=['geometry', 'mmgbsa', 'intEnergy', 'perResRec', 'perResLig', 'matrix', 'video'])

    # Add YML input file argument
    parser.add_argument("input", help="Input YML file")

    # Parse arguments
    args = parser.parse_args()

    global_parameters = get_parameters(args.input)

    top = global_parameters.get('top', None)
    trj = global_parameters.get('trj', None)
    traj_type = global_parameters.get('traj_type', None)

    dry_complex_topology = global_parameters.get('dry_complex_topology', None)
    receptor_topology = global_parameters.get('receptor_topology', None)
    ligand_topology = global_parameters.get('ligand_topology', None)
    parameters = global_parameters.get('parameters', None)

    # Check if the input files exist
    check_input_files(top, trj, dry_complex_topology, receptor_topology, ligand_topology, parameters)

    receptorSel = global_parameters['receptorSel']
    ligandSel = global_parameters['ligandSel']

    binding_site_resids = global_parameters['binding_site_resids']
    ligand_resids = global_parameters['ligand_resids']
    receptor_resids = global_parameters['receptor_resids']

    receptor_chain = global_parameters['receptor_chain']
    ligand_chain = global_parameters['ligand_chain']

    ligand_selection = f"resid {ligand_resids} and chainID {ligand_chain}"
    binding_site_selection = f"resid {binding_site_resids} and chainID {receptor_chain}"

    receptor_backbone_selection = f"{global_parameters['receptor_backbone']} and resid {receptor_resids} and chainID {receptor_chain}"
    receptor_calphas_selection = f"{global_parameters['receptor_calphas']} and resid {receptor_resids} and chainID {receptor_chain}"

    distanceCutoff = global_parameters['distanceCutoff']

    timestep = global_parameters['timestep']
    stride = global_parameters['stride']
    
    numResidRec = global_parameters['numResidRec']
    numResidLig = global_parameters['numResidLig']

    namdPATH = global_parameters['namdPATH']

    check_input_files(namdPATH)
    
    n_procs = global_parameters['n_procs']

    ref_bool = global_parameters['ref_bool']
    ref_pdb = global_parameters['ref_pdb']

    workdir = os.getcwd()
    add = timestep*stride

    if args.task == 'geometry':
        if not os.path.exists('Geometry'):
            os.makedirs(os.path.join(workdir, 'Geometry'))
        if ref_bool:
            plotRMSDvsReference(top, trj, ref_pdb, receptor_backbone_selection, ligand_selection)
#           plotDynamicRMSDvsReference()
        else:
            plotCDMDist(top, trj, ligand_selection, binding_site_selection)
#           plotDynamicCDMDist()
        plotRMSDRec(top, trj, receptor_backbone_selection)
        plotRMSDLig(top, trj, receptor_backbone_selection, ligand_selection)
        plotRMSFRec(top, trj, receptor_backbone_selection)
        # plotRMSFLig()
        plotRgRec(top, trj, receptorSel, add)
        plotRgLig(top, trj, ligandSel, add)
        plotERMSD(top, trj, ref_pdb, ref_bool, add)
        mountPanelRec()
        mountPanelLig()
    elif args.task == 'mmgbsa':
        if not os.path.exists('MMGBSA'):
            os.makedirs(os.path.join(workdir, 'MMGBSA'))
        calcMMGBSA(trj, receptor_topology, ligand_topology, dry_complex_topology)
        processMMGBSA()    
        plotMMGBSAvsDistance(ref_bool)
        plotMMGBSAvsTime(ref_bool)
    elif args.task == 'intEnergy':
        if not os.path.exists('InteractionEnergy'):
            os.makedirs(os.path.join(workdir, 'InteractionEnergy'))
        InteractionEnergyNAMD(trj, traj_type, timestep, add, receptorSel, ligandSel, dry_complex_topology, namdPATH, stride)
        plotIntEnergyvsDistance(ref_bool)
        plotIntEnergyvsTime(ref_bool)
        plotInteractionEnergyLandscape(ref_bool)
        plotDynamicInteractionEnergyLandscape(ref_bool)
    elif args.task == 'perResRec':
        if not os.path.exists('per_residue_receptor'):    
            os.makedirs('per_residue_receptor')
        getContactsRec(top, trj, receptorSel, ligandSel, distanceCutoff)
        topContactsRec(parameters['numResidRec'])
        getResnamesRec(top, trj, receptorSel)
        interactionEnergyNAMDRec(n_procs)
        HeatmapProcessingRec()
        mask_min, mask_max = determineHeatmapMaskValuesRec()
        df2 = plotHeatmapRec(mask_min, mask_max)
        plotDynamicHeatmapRec(df2, mask_min, mask_max)
    elif args.task == 'perResLig':
        if not os.path.exists('per_residue_ligand'):
            os.makedirs('per_residue_ligand') 
        getContactsLig()
        topContactsLig(parameters['numResidLig'])
        getResnamesLig()
        interactionEnergyNAMDLig()
        HeatmapProcessingLig()
        mask_min, mask_max = determineHeatmapMaskValuesLig()
        df2 = plotHeatmapLig(mask_min, mask_max)
        plotDynamicHeatmapLig(df2, mask_min, mask_max)
    elif args.task == 'matrix':
        if not os.path.exists('per_residue_matrix'):
            os.makedirs('per_residue_matrix')
        interactionEnergyNAMDPairwiseMatrix()
        HeatmapProcessingPairwiseMatrix()
        mask_min, mask_max = determineHeatmapMaskValuesPairwiseMatrix()
        plotHeatmapPairwiseMatrix(mask_min, mask_max)
        DynamicHeatmapProcessingPairwiseMatrix()
        plotDynamicHeatmapPairwiseMatrix(mask_min, mask_max, numShiftRec, numShiftLig)
        mountEnergyPanel()
    elif args.task == 'video':
        if not os.path.exists('video'):
            os.makedirs('video')
        mountVideo()    
    else:
        logger.info(HELP)


    
