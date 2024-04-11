# Script to visualize the ensembles:
#   1. The centroid of each ensemble with the RMSF 
#   2. The trajectory of each ensemble aligned to the centroid

import os
import sys
import pymol 
import glob
import yaml
import zipfile
import logging
from pathlib import Path


#############
# Functions #
#############

def load_model(model_name: str, model_path: str, models_cartoon_color: str, residues_to_highlight: dict) -> None:
    """
    Function to load a model in pymol and highlight some residues 

    Inputs:
    -------
        
        model_name             (str): name of the model to load (object name in pymol)
        model_path             (str): path to the model to load (pdb file)
        models_cartoon_color   (str): color of the model
        residues_to_highlight (dict): dictionary with the special residues to highlight, e.g. {"i. 114-124": "cyan"}
    """

    # Load model
    pymol.cmd.load(model_path, model_name)

    # Change color of the model
    pymol.cmd.color(models_cartoon_color, model_name)

    # Change color of the special residues
    for residues in residues_to_highlight.keys():
        pymol.cmd.color(residues_to_highlight[residues], f"{model_name} and {residues}")

    return

def find_models(pockets_summary: dict, output_folder: str, show_all: bool) -> dict:
    """
    Function to create a dictionary with information for each model.

    E.g. {"model_name": {"model_path": "path/to/model.pdb", "pockets_zip_path": "path/to/filtered_pockets.zip", pockets: ["pocket1", "pocket2"]}

    Inputs:
    -------
        
        pockets_summary (dict): dictionary with the summary of the pockets for each model
        output_folder    (str): path to the output folder
        show_all        (bool): whether to show all models or only the models with pockets
    """

    # Initialize dictionary
    Models = {}

    # Counter for the number of models with pockets
    models_with_pockets = 0

    # For each model
    for model_name in pockets_summary.keys():

        # Find the model path
        model_path = os.path.join(output_folder, f"{model_name}/model.pdb")

        # Find filtered pockets path
        pockets_zip_path = os.path.join(output_folder, f"{model_name}/step8_filter_residue_com/filtered_pockets.zip")

        # Find pockets list
        pockets_list = pockets_summary[model_name]["pockets"]

        # If the model path exists
        if os.path.exists(model_path):

            # If the model has pockets
            if os.path.exists(pockets_zip_path):
                # Save the model information
                Models[model_name] = {"model_path": model_path, "pockets_zip_path": pockets_zip_path, "pockets": pockets_list}
                models_with_pockets += 1

            elif show_all:
                # Save the model information
                Models[model_name] = {"model_path": model_path, "pockets_zip_path": None}
    
    # Log
    if models_with_pockets == 0:
        logger.warning("No models with filtered pockets found")
    else:
        logger.info(f"Found {models_with_pockets} models with filtered pockets")

    return Models

def find_pockets(Models: dict) -> None:
    """
    For each model with a pockets zip file, extract the pockets in their corresponding folder and add the paths to the dictionary.

    E.g. {"model_name": {"model_path": "path/to/model.pdb", "pockets_zip_path": "path/to/filtered_pockets.zip", pockets: ["pocket1", "pocket2"], "pockets_paths": ["path/to/pocket1_vert.pqr", "path/to/pocket2_vert.pqr"]}

    Inputs:
    -------
        
        Models (dict): dictionary with the information of each model
    """

    # For each model
    for model_name in Models.keys():

        # Find the pockets zip path
        pockets_zip_path = Models[model_name]["pockets_zip_path"]

        # Find the pockets list (names of the pockets)
        pockets_list = Models[model_name]["pockets"]

        # If the pockets zip file is not None
        if pockets_zip_path is not None:

            # Check the existence of the pockets zip file
            if not os.path.exists(pockets_zip_path):
                logger.error(f"Pockets zip file {pockets_zip_path} for model {model_name} does not exist")
                continue

            # Find the folder where the pockets will be extracted
            pocket_folder = Path(pockets_zip_path).parent

            # Unzip the pockets
            with zipfile.ZipFile(pockets_zip_path, 'r') as zip_ref:
                zip_ref.extractall(pocket_folder)

            # Add list with the paths of the filtered pockets
            Models[model_name]["pockets_paths"] = [os.path.join(pocket_folder, f"{pocket}_vert.pqr") for pocket in pockets_list]

    return

#############
# Constants #
#############

# Define constants
LOG_FILENAME = "visualize_results.log"

# Set up logging
# Create logger
logger = logging.getLogger("PyMOL_pockets")

# Create file handler
file_handler = logging.FileHandler(LOG_FILENAME, mode = "w")
file_handler.setLevel(logging.INFO)

# Create formatter and add it to file handler 
file_format = logging.Formatter('%(name)s [%(levelname)s]: %(message)s')
file_handler.setFormatter(file_format)

# Set logger level
logger.setLevel(logging.INFO)

# Add file handler to the logger
logger.addHandler(file_handler)

logger.info("Cavity analysis visualization with PyMOL")
logger.info("========================================\n")

####################
# Important Inputs #
####################

# Output folder of the cavity analysis workflow
output_folder = "/home/pnavarro/2023_IT/tests/BiobbWorkflows/Single_protein/Cavity_analysis/output_traj_b"

# Experimental structures to add to the visualization
experimental_folder = "/home/pnavarro/2023_ENSEM/P53/Experimental_structures/Reference_structures/prepared"

# Wether to show all the representative structures or just the representative structures with pockets
show_all = False

# Distance to the pocket to show the atoms explicitly
licorice_distance_threshold = 6


#################
# Other options #
#################

# Summary file with all the pockets
pockets_summary_path = f"{output_folder}/summary_by_drug_score.yml"

# Colors for the experimental structures
experimental_colors = ["tv_red", "salmon", "orange", "yellow", "tv_green", "cyan", "marine", 
                "purple", "white", "chocolate", "forest", "magenta", "teal", "deepblue", 
                "limon", "splitpea", "blue", "gray", "red", "pink" ]
experimental_colors.extend(experimental_colors)

# Color of the representative structures
models_cartoon_color = "white"

# Dictionary with special residues to highlight, e.g. {"i. 114-124": "cyan"}. See the end of the script for more examples. 
residues_to_highlight = {
}

##################
# Initialization #
##################

# Check existence of output folder
if not os.path.exists(output_folder):
    logger.error(f"Output folder {output_folder} does not exist")
    sys.exit()

# Check existence of summary file
if not os.path.exists(pockets_summary_path):
    logger.error(f"Summary file {pockets_summary_path} does not exist")
    sys.exit()

# Read summary yaml file as dictionary
with open(pockets_summary_path, 'r') as file:
    pockets_summary = yaml.load(file, Loader=yaml.FullLoader)

# If pockets summary is empty, exit
if len(pockets_summary) == 0:
    logger.error("No pockets found in the summary file")
    sys.exit()

# Create dictionary with information per model
Models = find_models(pockets_summary, output_folder, show_all)

# Find and extract pockets for each model
find_pockets(Models)

###########################
# Load models and pockets #
###########################

# Find all model names
model_names = list(Models.keys())

# Debug 
logger.info(f"Models: {model_names}")

# Load the first model as reference
load_model("reference", Models[model_names[0]]["model_path"], "green", residues_to_highlight)

for name in model_names:

    # Log 
    logger.info(f"Loading model {name}")

    # Load the model to be aligned (the one that will be shown in the visualization)
    load_model(name, Models[name]["model_path"], models_cartoon_color, residues_to_highlight)

    # Align model to reference
    pymol.cmd.align(name, "reference")

    # Load the pockets if any
    if Models[name].get("pockets_paths") is not None:

        # Load the model to be used to align pockets (this object will be used temporarily and then deleted)
        pymol.cmd.load(Models[name]["model_path"], f"{name}_temp")
        
        # For each pocket
        for pocket_name, pocket_path in zip(Models[name]["pockets"], Models[name]["pockets_paths"]):

            ## Load and align the pocket to the reference

            # Load the pocket
            pymol.cmd.load(pocket_path, pocket_name)

            # Log 
            logger.info(f"Loading pocket {pocket_name}")

            # Create a temporal object with the pocket and the temporal model
            pymol.cmd.copy_to(f"{name}_and_pocket", f"{name}_temp or {pocket_name}")

            # Delete the original pocket
            pymol.cmd.delete(pocket_name)

            # Align temporal object to reference 
            pymol.cmd.align(f"{name}_and_pocket", "reference")

            # Name for the aligned pocket 
            aligned_pocket = f"{pocket_name}_{name}"

            # Extract the pocket from the temporal object
            pymol.cmd.extract(aligned_pocket,  f"{name}_and_pocket and resname STP")

            # Delete the temporal object
            pymol.cmd.delete(f"{name}_and_pocket")

            ## Configure the visualization of the pocket and its surroundings

            # Hide licorice representation
            pymol.cmd.hide("licorice", aligned_pocket)

            # Color pocket by atom type using util.cbag
            pymol.cmd.util.cbag(aligned_pocket)

            # Show mesh representation
            pymol.cmd.show("mesh", aligned_pocket)

            # Show licorice representation of the atoms close to the pocket (show licorice, (all within 5 of pocket_name) and not pocket_name)
            pymol.cmd.show("licorice", f"({name} within {licorice_distance_threshold} of {aligned_pocket}) and not {aligned_pocket}")

            # Color the side-chain atoms close to the pocket by atom type using util.cbac
            pymol.cmd.util.cbaw(f"({name} within {licorice_distance_threshold} of {aligned_pocket}) and not {aligned_pocket}")

        # Delete the temporal model
        pymol.cmd.delete(f"{name}_temp")

        # List all the aligned pockets for this model
        aligned_pockets = [f"{pocket_name}_{name}" for pocket_name in Models[name]["pockets"]]

        # Group model and pockets
        pymol.cmd.group(f"{name}_group", f"{name} {'or'.join(aligned_pockets)}")

#################################
# Load experimental structures  #
#################################

if experimental_folder is not None:

    # Find the experimental paths
    experimental_paths = sorted(glob.glob(f"{experimental_folder}/*.pdb"))

    # Load experimental structures
    for model_index, model_path in enumerate(experimental_paths):

        # Find model name
        model_name = Path(model_path).stem

        # Load model
        pymol.cmd.load(model_path, model_name)

        # Align to the reference structure
        result = pymol.cmd.align(model_name, "reference", cycles = 20)

        # Set a random color to the model
        pymol.cmd.color(experimental_colors[model_index], model_name)


#####################
# General settings  #
#####################

# Hide all hydrogens (hide licorice, hydrogens)
pymol.cmd.hide("licorice", "hydrogens")

# Hide reference model
pymol.cmd.hide("everything", "reference")

# Set min_mesh_spacing to 0.4
pymol.cmd.set("min_mesh_spacing", 0.4)

#########
# Notes #
#########

# Example of residues_to_highlight dictionary:

#residues_to_highlight = {
#    "i. 114-124": "cyan",
#    "i. 125-127": "blue",
#    "i. 131-136": "magenta",
#    "i. 139-146": "yellow",
#    "i. 277-289": "red",
#    "i. 221-230": "orange",
#    "i. 231-236": "green"
# }

# residues_to_highlight = {
#     "i. 21-31": "cyan",
#     "i. 32-34": "blue", 
#     "i. 38-43": "magenta",
#     "i. 46-53": "yellow",
#     "i. 184-196": "red",
#     "i. 128-137": "orange",
#     "i. 138-143": "green"
# }