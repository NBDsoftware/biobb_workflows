# Script to visualize the ensembles:
#   1. The centroid of each ensemble with the RMSF 
#   2. The trajectory of each ensemble aligned to the centroid

import os
import pymol 
import glob
import zipfile
import logging
import regex as re
from pathlib import Path

def load_model(model_name, model_path, residues_to_highlight):

    """
    Function to load a model and highlight the special residues 
    """

    # Log
    logger.info(f"Loading {model_name}")

    # Load model
    pymol.cmd.load(model_path, model_name)

    # Change color of the model
    pymol.cmd.color("gray", model_name)

    # Change color of the special residues
    pymol.cmd.color("cyan", f"{model_name} and {residues_to_highlight}")

    return

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
output_folder = "/home/nbd-pablo/Documents/2023_ENSEM/P53/Crystals_and_models/prepared/pocket_analysis"

# Change to the input folder of the cavity analysis or path were the models are stored
input_folder = "/home/nbd-pablo/Documents/2023_ENSEM/P53/Crystals_and_models/prepared"

# Wether to show all the models or just the ones with pockets
show_all_models = False

# Should be ok by default
filtered_pockets_folder = "step5_filter_residue_com"
filtered_pockets_filename = "filtered_pockets.zip"

# Special residues to highlight
residues_to_highlight = f"i. 114-124" # f"i. 21-31"

# List pdb files inside the input folder
model_paths = glob.glob(f"{input_folder}/*.pdb")

# List directories inside the output folder
pockets_folders = glob.glob(f"{output_folder}/*/")

# Add filtered pockets folder to each path 
pockets_folders = [f"{pocket_folder}{filtered_pockets_folder}/" for pocket_folder in pockets_folders]

# Add filtered pockets filename to each path
filtered_pockets_paths = [f"{pocket_folder}{filtered_pockets_filename}" for pocket_folder in pockets_folders]

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

# Load all models and corresponding filtered pockets if any
for model_path in model_paths:

    # Find model name
    model_name = Path(model_path).stem

    # Load all models
    if show_all_models:
        load_model(model_name, model_path, residues_to_highlight)

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
                load_model(model_name, model_path, residues_to_highlight)

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

                # Show mesh representation (show mesh, pocket_name)
                pymol.cmd.show("mesh", pocket_name)

                # Show licorice representation of the atoms close to the pocket (show licorice, (all within 5 of pocket_name) and not pocket_name)
                pymol.cmd.show("licorice", f"({model_name} within 6 of {pocket_name}) and not {pocket_name}")
            
            # Group models and pockets together
            pymol.cmd.group(f"{model_name}_group", f"{model_name} {' '.join(pocket_objects)}")

        else:
            logger.info(f"No pockets found for {model_name}")

    else:
        logger.info(f"No {filtered_pockets_folder} folder found for {model_name}")

# Hide all hydrogens (hide licorice, hydrogens)
pymol.cmd.hide("licorice", "hydrogens")