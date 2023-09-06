# Script to visualize the ensembles:
#   1. The centroid of each ensemble with the RMSF 
#   2. The trajectory of each ensemble aligned to the centroid

import os
import pymol 
import glob
import yaml
import zipfile
import logging
import regex as re
from pathlib import Path


#############
# Functions #
#############

def load_model(model_name, model_path, models_cartoon_color, residues_to_highlight):

    """
    Function to load a model and highlights the special residues 

    Inputs:
    -------
        
        model_name (str): name of the model to load (object name in pymol)
        model_path (str): path to the model to load (pdb file)
        models_cartoon_color (str): color of the model
        residues_to_highlight (dict): dictionary with the special residues to highlight, e.g. {"i. 114-124": "cyan"}
    """

    # Log
    logger.info(f"Loading {model_name}")

    # Load model
    pymol.cmd.load(model_path, model_name)

    # Change color of the model
    pymol.cmd.color(models_cartoon_color, model_name)

    # Change color of the special residues
    for residues in residues_to_highlight.keys():
        pymol.cmd.color(residues_to_highlight[residues], f"{model_name} and {residues}")

    return

def sort_pockets(pockets_summary_path: str, output_folder: str):
    """
    Function that reads the YAML summary file with all models and pockets and sorts the models by:

        1. Volume of the largest pocket in that model
        2. Druggability score of the best pocket in that model
        3. Score of the best pocket in that model
    
    It then writes new YAML files with the sorted models and their pockets.
    """

    # Find the YAML summary file name
    pockets_summary_filename = Path(pockets_summary_path).stem

    # Read the YAML summary file
    with open(pockets_summary_path, "r") as f:
        pockets_summary = yaml.load(f, Loader = yaml.FullLoader)

    # Sort the pockets by volume
    sorted_pockets = dict(sorted(pockets_summary.items(), key = lambda x: largest_volume(x[1]), reverse = True))

    # Write the sorted pockets to a new YAML file
    with open(f"{output_folder}/{pockets_summary_filename}_by_volume.yml", "w") as f:
        yaml.dump(sorted_pockets, f, sort_keys = False)

    # Sort the pockets by druggability score
    sorted_pockets = dict(sorted(pockets_summary.items(), key = lambda x: highest_drug_score(x[1]), reverse = True))

    # Write the sorted pockets to a new YAML file 
    with open(f"{output_folder}/{pockets_summary_filename}_by_drug_score.yml", "w") as f:
        yaml.dump(sorted_pockets, f, sort_keys = False)

    # Sort the pockets by score
    sorted_pockets = dict(sorted(pockets_summary.items(), key = lambda x: highest_score(x[1]), reverse = True))

    # Write the sorted pockets to a new YAML file
    with open(f"{output_folder}/{pockets_summary_filename}_by_score.yml", "w") as f:
        yaml.dump(sorted_pockets, f, sort_keys = False)
    
    return

def largest_volume(model: dict):
    """
    Function to sort the pockets by volume
    """

    # Find the largest volume
    largest_volume = 0
    for pocket in model["pockets"]:
        if model[pocket]["volume"] > largest_volume:
            largest_volume = model[pocket]["volume"]

    return largest_volume

def highest_drug_score(model: dict):
    """
    Function to sort the pockets by druggability score
    """

    # Find the highest druggability score
    highest_drug_score = 0
    for pocket in model["pockets"]:
        if model[pocket]["druggability_score"] > highest_drug_score:
            highest_drug_score = model[pocket]["druggability_score"]

    return highest_drug_score

def highest_score(model: dict):
    """
    Function to sort the pockets by score
    """

    # Find the highest score
    highest_score = 0
    for pocket in model["pockets"]:
        if model[pocket]["score"] > highest_score:
            highest_score = model[pocket]["score"]

    return highest_score

#############
# Constants #
#############

# Define constants
LOG_FILENAME = "visualize_results.log"

# Set up logging
# Create logger
logger = logging.getLogger("PYMOL_POCKETS")

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

logger.info("PYMOL_POCKETS")
logger.info("==============\n")

##########
# Inputs #
##########

# Change to the output folder of the cavity analysis
#output_folder = "/home/nbd-pablo/Documents/2023_ENSEM/P53/round3_new/clustering/extended_pocket/states/pocket_analysis_extendedRdown"
output_folder = "/home/nbd-pablo/Documents/2023_ENSEM/P53/round0/trajectories/clustering/extended_pocket/states/pocket_analysis_allstate9_vol300_ds3_sc3"

# Change to the input folder of the cavity analysis or path were the models are stored
# input_folder = "/home/nbd-pablo/Documents/2023_ENSEM/P53/round3_new/clustering/extended_pocket/states/input/extended_R_down"
input_folder = "/home/nbd-pablo/Documents/2023_ENSEM/P53/round0/trajectories/clustering/extended_pocket/states/input/all_state9"

# Summary file with all the pockets
pockets_summary_path = f"{output_folder}/pocket_analysis.yml"

# Experimental structure to compare the models to
experimental_folder = "/home/nbd-pablo/Documents/2023_ENSEM/P53/Crystals_and_models/prepared"
experimental_colors = ["tv_red", "salmon", "orange", "yellow", "tv_green", "cyan", "marine", 
                "purple", "white", "chocolate", "forest", "magenta", "teal", "deepblue", 
                "limon", "splitpea", "blue", "gray", "red", "pink" ]
experimental_colors.extend(experimental_colors)

# Wether to show all the models or just the ones with pockets
show_all_models = False

# Wether to sort the pockets by volume, druggability score and score or not
sort_pockets_summary = True

# Distance to the pocket to show the atoms explicitly
distance_to_pocket = 8

# Color of the models
models_cartoon_color = "white"

# Should be ok by default
filtered_pockets_folder = "step5_filter_residue_com"
filtered_pockets_filename = "filtered_pockets.zip"

# Default dictionary with special residues to highlight
residues_to_highlight = {
}
    
# Dictionary with special residues to highlight, e.g. {"i. 114-124": "cyan"} 
#residues_to_highlight = {
#    "i. 114-124": "cyan",
#    "i. 125-127": "blue",
#    "i. 131-136": "magenta",
#    "i. 139-146": "yellow",
#    "i. 277-289": "red",
#    "i. 221-230": "orange",
#    "i. 231-236": "green"
# }

# Different version of the dictionary (for simulation structures)
# residues_to_highlight = {
#     "i. 21-31": "cyan",
#     "i. 32-34": "blue", 
#     "i. 38-43": "magenta",
#     "i. 46-53": "yellow",
#     "i. 184-196": "red",
#     "i. 128-137": "orange",
#     "i. 138-143": "green"
# }

# List pdb files inside the input folder
model_paths = glob.glob(f"{input_folder}/*.pdb")

# List directories inside the output folder
pockets_folders = glob.glob(f"{output_folder}/*/")

# Add filtered pockets folder to each path 
pockets_folders = [f"{pocket_folder}{filtered_pockets_folder}/" for pocket_folder in pockets_folders]

# Add filtered pockets filename to each path
filtered_pockets_paths = [f"{pocket_folder}{filtered_pockets_filename}" for pocket_folder in pockets_folders]

# Sort pockets
if sort_pockets_summary:
    sort_pockets(pockets_summary_path, output_folder)

##########################
# Unzip filtered pockets #
##########################

# Unzip each filtered pocket inside the corresponding folder
for filtered_pockets_path in filtered_pockets_paths:

    # Get the folder where the filtered pockets are
    pocket_folder = Path(filtered_pockets_path).parent

    # If the file exists
    if os.path.exists(filtered_pockets_path):

        # Unzip the filtered pockets
        with zipfile.ZipFile(filtered_pockets_path, 'r') as zip_ref:
            zip_ref.extractall(pocket_folder)

###########################
# Load models and pockets #
###########################

# Load the first model as reference
reference_model = Path(model_paths[0]).stem
load_model(reference_model, model_paths[0], models_cartoon_color, residues_to_highlight)

# Load all models and corresponding filtered pockets if any
for model_path in model_paths:

    # Find model name
    model_name = Path(model_path).stem

    # Load all models
    if show_all_models:
        load_model(model_name, model_path, models_cartoon_color, residues_to_highlight)

    # Create corresponding filtered pockets path
    filtered_pockets_path = os.path.join(output_folder, f"{model_name}/{filtered_pockets_folder}/")

    # Check if the filtered pockets path exists
    if os.path.exists(filtered_pockets_path):

        # Find all the pockets inside the filtered pockets path (vertices as pqr files)
        pocket_paths = glob.glob(f"{filtered_pockets_path}/*.pqr")

        # If pockets were found, load them
        if len(pocket_paths) > 0:
            
            # Load only the models with pockets
            if not show_all_models:
                load_model(model_name, model_path, models_cartoon_color, residues_to_highlight)

            pocket_objects = []

            # For each pocket
            for pocket_path in pocket_paths:

                # Find pocket name
                pocket_name = Path(pocket_path).stem

                # Remove "_vert" from the pocket name
                pocket_name = pocket_name.replace("_vert", f"_{model_name}")

                # Log
                logger.info(f"Loading {pocket_name}")

                # Save pocket name
                pocket_objects.append(pocket_name)

                # Load pocket
                pymol.cmd.load(pocket_path, pocket_name)

                # Log 
                logger.info(f"Loaded {pocket_path}")

                # Hide licorice representation (hide licorice, pocket_name)
                pymol.cmd.hide("licorice", pocket_name)

                # Color pocket by atom type using util.cbag
                pymol.cmd.util.cbac(pocket_name)

                # Show mesh representation (show mesh, pocket_name)
                pymol.cmd.show("mesh", pocket_name)

                # Show licorice representation of the atoms close to the pocket (show licorice, (all within 5 of pocket_name) and not pocket_name)
                pymol.cmd.show("licorice", f"({model_name} within {distance_to_pocket} of {pocket_name}) and not {pocket_name}")

                # Color the side-chain atoms close to the pocket by atom type using util.cbac
                pymol.cmd.util.cbaw(f"({model_name} within {distance_to_pocket} of {pocket_name}) and not {pocket_name}")
            
            # Copy the last pocket into temporal new object before aligning to reference model
            pymol.cmd.copy_to(f"{model_name}_and_pockets", f"{model_name} or {pocket_objects[-1]}")

            # Delete the original objects
            pymol.cmd.delete(pocket_objects[-1])
            pymol.cmd.delete(model_name)
            
            # Align temporal object to reference 
            if model_name != reference_model:
                pymol.cmd.align(f"{model_name}_and_pockets", reference_model)

            # Extract the pocket  
            pymol.cmd.extract(f"{pocket_objects[-1]}", f"{model_name}_and_pockets and resname STP")

            # Extract the model
            pymol.cmd.extract(f"{model_name}", f"{model_name}_and_pockets and not resname STP")

            # Delete temporal object
            pymol.cmd.delete(f"{model_name}_and_pockets")

            # Group model and pockets together
            pymol.cmd.group(f"{model_name}_group", f"{model_name} {' '.join(pocket_objects)}")

        else:
            logger.info(f"No pockets found for {model_name}")

            # If we are showing all models, align it to the reference
            if show_all_models and model_name != reference_model:
                pymol.cmd.align(model_name, reference_model)
    else:
        logger.info(f"No {filtered_pockets_folder} folder found for {model_name}")

        # If we are showing all models, align it to the reference
        if show_all_models and model_name != reference_model:
            pymol.cmd.align(model_name, reference_model)


#################################
# Load experimental structures  #
#################################

# Find the experimental paths
experimental_paths = sorted(glob.glob(f"{experimental_folder}/*.pdb"))

# Load experimental structures
for model_index, model_path in enumerate(experimental_paths):

    # Find model name
    model_name = Path(model_path).stem

    # Load model
    pymol.cmd.load(model_path, model_name)

    # Align to the reference structure. NOTE: "and not {flexible_loop_selection} and name CA"
    result = pymol.cmd.align(model_name , reference_model, cycles = 20)

    # Set a random color to the model
    pymol.cmd.color(experimental_colors[model_index], model_name)


#####################
# General settings  #
#####################

# Hide all hydrogens (hide licorice, hydrogens)
pymol.cmd.hide("licorice", "hydrogens")

# Set min_mesh_spacing to 0.4
pymol.cmd.set("min_mesh_spacing", 0.4)