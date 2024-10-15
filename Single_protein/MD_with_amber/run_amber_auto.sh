#!/bin/bash -l
#SBATCH -A p200581
#SBATCH --job-name="AMBER_prod"
#SBATCH -p gpu
#SBATCH -q default
#SBATCH --time=2-00:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-task=1
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err

# Purge loaded modules
module purge 

ml CUDA/11.7.0
ml GCCcore/10.2.0 

# Load conda / miniconda module
module load Miniconda3

# Activate previously created conda environment from environment.yml
#source activate /home/users/u102183/.conda/envs/amber_setup_env
conda activate amber_setup_env

output_path=/project/scratch/p200581/output_traj
pdb_code=$1
md_type=$2
input_path=/project/scratch/p200581/dynamic_pdb_90_lt500/non_gapped

echo  $PYTHONPATH
which python

# Launch workflow
if [ "$md_type" = "free" ]; then
    pdb_code_type=$pdb_code"_free"
    python /home/users/u102183/just_run/workflow_free.py --config /home/users/u102183/just_run/workflow_free.yml --output $output_path/$pdb_code_type --input_pdb $input_path/$pdb_code.pdb
elif [ "$md_type" = "gamd" ]; then
    pdb_code_type=$pdb_code"_gamd"
    python /home/users/u102183/just_run/workflow_gamd.py --config /home/users/u102183/just_run/workflow_gamd.yml --output $output_path/$pdb_code_type --input_pdb $input_path/$pdb_code.pdb
else
    echo "Incorrect md_type value. Must be 'gamd' or 'free'."
fi