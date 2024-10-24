#!/bin/bash
##### Number of tasks
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH -t 0-01:00:00
#####Output and error log files
#SBATCH -e gmx_%j.e
#SBATCH -o gmx_%j.o


# Purge loaded modules
module purge 

# Load GROMACS module
ml GROMACS

# Activate conda environment
module load Miniconda3
source activate /shared/work/BiobbWorkflows/envs/biobb_md

# Path to the workflow
REPO_PATH=/shared/work/BiobbWorkflows/src/biobb_workflows
WF_PATH=$REPO_PATH/Single_protein/1.MD_with_gromacs/

# Input arguments
INPUT_PDB=$1
PDB_CHAIN=$2
OUTPUT_DIR=$3
LIGAND_PARAMETERS=$4
NUM_PARTS=$5
NSTEPS=$6
DO_FINAL_ANALYSIS=$7

# Read the condition to do the final analysis
if [ "$DO_FINAL_ANALYSIS" == "true" ]; then
    FINAL_ANALYSIS="--final_analysis"
else
    FINAL_ANALYSIS=""
fi

# Create folder for the HIS protonation
mkdir -p $OUTPUT_DIR/his_protonation
HIS_PDB=$OUTPUT_DIR/his_protonation/model.pdb

# Remove water molecules, keep only protein residues and run Reduce to add hydrogens
pdb4amber -i $INPUT_PDB -o $HIS_PDB -d -p --reduce

# Check the HIS_PDB file exists
if [ ! -f $HIS_PDB ]; then
    echo "Error: pdb4amber failed to generate the HIS_PDB file."
    exit 1
fi

# Move log file if it exists
LOG_FILE=reduce_info.log
if [ -f $LOG_FILE ]; then
    mv $LOG_FILE $OUTPUT_DIR/his_protonation
fi

# Extract column 4 from lines containing "CA HI" into a bash array
histidine_residues=($(grep "CA  HI" $HIS_PDB | awk '{print $4}'))

# Check if any histidine residues were found
if [ ${#histidine_residues[@]} -eq 0 ]; then

    echo "No histidine residues found. Launching the workflow."

    python $WF_PATH/biobb_md_setup_mutation.py --config input_HPC.yml --num_parts $NUM_PARTS --input_pdb $INPUT_PDB --output $OUTPUT_DIR --pdb_chains $PDB_CHAIN --nsteps $NSTEPS $FINAL_ANALYSIS --ligand_parameters $LIGAND_PARAMETERS

elif [ ${#histidine_residues[@]} -ne 0 ]; then

    echo "Histidine residues found, launching the workflow with the following protonation states: ${histidine_residues[@]}"

    # Initialize an array to store the numbers
    histidine_numbers=""

    # Loop through the array and replace histidine names with the corresponding values
    for residue in "${histidine_residues[@]}"; do
        case $residue in
            HID)
                histidine_numbers+="0 "  # Store 0 for HISD
                ;;
            HIE)
                histidine_numbers+="1 "  # Store 1 for HISE
                ;;
            HIP)
                histidine_numbers+="2 "  # Store 2 for HISH
                ;;
            *)
                echo "Unknown residue: $residue"
                ;;
        esac
    done

    python $WF_PATH/biobb_md_setup_mutation.py --config input_HPC.yml --num_parts $NUM_PARTS --input_pdb $INPUT_PDB --output $OUTPUT_DIR --pdb_chains $PDB_CHAIN --nsteps $NSTEPS $FINAL_ANALYSIS --his "$histidine_numbers" --ligand_parameters $LIGAND_PARAMETERS
fi
