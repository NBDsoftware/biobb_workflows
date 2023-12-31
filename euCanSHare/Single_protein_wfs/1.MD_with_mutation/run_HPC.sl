#!/bin/bash
#SBATCH --job-name=single_protein_wf1
#SBATCH --nodes=1                 # node count
#SBATCH --ntasks=8                                  # total number of tasks across all nodes
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=1000
# #SBATCH --constraint=gpu
# #SBATCH --gres=gpu:1
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err
# #SBATCH --mail-type=begin                              # send email when job begins
# #SBATCH --mail-type=end                                # send email when job ends
# #SBATCH --mail-user=your@email.com

# Purge loaded modules
module purge 

# Load conda / miniconda module
module load Miniconda3/4.9.2

# Load GROMACS module
# module load GROMACS/2022.3-intel-2021b-CUDA.11.6.0
module load GROMACS/2022.3-intel-2021b

# Activate previously created conda environment from environment.yml
source activate /path/to/conda/envs/biobb_sp_md

# Launch example 1
# INPUT_PDB=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/euCanShare/Single_protein_wfs/1.MD_with_mutation/input/P38alpha4LOO.pdb
# python biobb_md_setup_mutation.py --config input_HPC.yml --num_trajs 2 --input_pdb $INPUT_PDB

# Launch example 2
# INPUT_PDB=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/euCanShare/Single_protein_wfs/1.MD_with_mutation/input/P38alpha4LOO.pdb
# OUTPUT_PATH=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/euCanShare/Single_protein_wfs/1.MD_with_mutation/output2
# python biobb_md_setup_mutation.py --config input_HPC.yml --num_trajs 1 --input_pdb $INPUT_PDB --output $OUTPUT_PATH --setup_only

# Launch example 3
# INPUT_PDB=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/euCanShare/Single_protein_wfs/1.MD_with_mutation/input/P38alpha4LOO.pdb
# OUTPUT_PATH=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/euCanShare/Single_protein_wfs/1.MD_with_mutation/output3
# python biobb_md_setup_mutation.py --config input_HPC.yml --num_trajs 1 --input_pdb $INPUT_PDB --output $OUTPUT_PATH --setup_only --fix_ss --fix_backbone

# Launch example 4
# INPUT_GRO=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/euCanShare/Single_protein_wfs/1.MD_with_mutation/output3/step7_genion/genion.gro
# INPUT_TOP=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/euCanShare/Single_protein_wfs/1.MD_with_mutation/output3/step7_genion/genion_top.zip
# OUTPUT_PATH=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/euCanShare/Single_protein_wfs/1.MD_with_mutation/output4
# python biobb_md_setup_mutation.py --config input_HPC.yml --num_trajs 1 --input_gro $INPUT_GRO --input_top $INPUT_TOP --output $OUTPUT_PATH 

# Launch example 5
# INPUT_PDB=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/euCanShare/Single_protein_wfs/1.MD_with_mutation/input/P38alpha4LOO.pdb
# OUTPUT_PATH=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/euCanShare/Single_protein_wfs/1.MD_with_mutation/output5
# python biobb_md_setup_mutation.py --config input_HPC.yml --num_trajs 1 --input_pdb $INPUT_PDB --output $OUTPUT_PATH --setup_only --pdb_chains A B 

# Launch example 6
# INPUT_PDB=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/euCanShare/Single_protein_wfs/1.MD_with_mutation/input/P38alpha4LOO.pdb
# OUTPUT_PATH=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/euCanShare/Single_protein_wfs/1.MD_with_mutation/output6
# python biobb_md_setup_mutation.py --config input_HPC.yml --num_trajs 1 --input_pdb $INPUT_PDB --output $OUTPUT_PATH --setup_only --pdb_chains A --mutation_list A:Arg5Ala A:Pro6Ala

# Launch example 8
INPUT_PDB=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/euCanShare/Single_protein_wfs/1.MD_with_mutation/input/P38alpha4LOO.pdb
OUTPUT_PATH=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/euCanShare/Single_protein_wfs/1.MD_with_mutation/output8
python biobb_md_setup_mutation.py --config input_HPC.yml --num_trajs 1 --input_pdb $INPUT_PDB --output $OUTPUT_PATH --pdb_chains A 
