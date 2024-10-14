#!/bin/bash
##### Number of tasks
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=52
#SBATCH -t 0-00:30:00
#####Output and errorlog files
#SBATCH -J 3FDL_1
#SBATCH -e gmx_%j.e
#SBATCH -o gmx_%j.o


# Purge loaded modules
module purge 

# Load conda / miniconda module
module load anaconda

# Load slurm to call srun in the computation node
module load slurm

# module load oneapi/2024.1 plumed/2.9.0
# Load GROMACS module

ml intel
ml openmpi
ml mkl
module load gromacs/2024.1-ompi

# Activate previously created conda environment from environment.yml
source activate /home/nost/nost674589/.conda/envs/biobb_sp_md
SECONDS=0
# Path to the workflow
REPO_PATH=/gpfs/scratch/cns14/nost674589
WF_PATH=$REPO_PATH/1.MD_with_mutation

# Input files
INPUT_PDB=/gpfs/scratch/cns14/nost674589/input_structures/apo
PDB_CODE=$1
NUM_PARTS=$2
NSTEPS=$3
ANALYSIS=$4
OUTPUT_DIR=/gpfs/scratch/cns14/nost674589/Gromacs_wf_2/output_traj
PDB_IN=$INPUT_PDB/$PDB_CODE.pdb
PDB_OUT=$INPUT_PDB/his_info/$PDB_CODE.pdb


pdb4amber -i $PDB_IN -o $PDB_OUT -d -p --reduce

# Check if the file exists
if [ ! -f $PDB_IN ]; then
    echo "File not found!"
    exit 1
fi

# Extract column 4 from lines containing "CA HI" into a bash array
histidine_residues=($(grep "CA  HI" $PDB_OUT | awk '{print $4}'))

# Check if any histidine residues were found
if [ ${#histidine_residues[@]} -eq 0 ]; then
    echo "No histidine residues found."
    exit 0
fi

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

# Print the final array of histidine numbers
echo "Histidine protonation states: ${histidine_numbers[@]}"

echo  $PYTHONPATH
which python

# Crear el directorio de salida si no existe
mkdir -p $OUTPUT_DIR
#--his $histidine_numbers
python $WF_PATH/new_wf_biobb_3.py --config proves.yml --num_parts $NUM_PARTS --input_pdb $PDB_IN --output $OUTPUT_DIR/$PDB_CODE --pdb_chains A --nsteps $NSTEPS --analysis $ANALYSIS --his "$histidine_numbers"

echo $SECONDS
