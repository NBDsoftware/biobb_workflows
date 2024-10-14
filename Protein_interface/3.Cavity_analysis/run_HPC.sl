#!/bin/bash
#SBATCH --job-name=PP_cavity_analysis_wf
#SBATCH --ntasks=1                                      # total number of tasks across all nodes
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=2000
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err
# #SBATCH --mail-type=begin                              # send email when job begins
# #SBATCH --mail-type=end                                # send email when job ends
# #SBATCH --mail-user=your@email.com

# Purge loaded modules
module purge 

# Load conda / miniconda module
module load Miniconda3/4.9.2

# Load SLURM
module load slurm/slurm/21.08.6 

# Activate previously created conda environment from environment.yml
source activate /home/pnavarro/.conda/envs/protein_protein_wf3

# Launch workflow
# python biobb_cavity_analysis.py --config input_HPC.yml 

# Launch workflow 2
INPUT_ZIP=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/euCanShare/Protein_interface_wfs/3.Cavity_analysis/input/top_refined_poses.zip
OUTPUT_PATH=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/euCanShare/Protein_interface_wfs/3.Cavity_analysis/output2
OUTPUT_SUMMARY=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/euCanShare/Protein_interface_wfs/3.Cavity_analysis/output2/summary_cavities.yml
python biobb_cavity_analysis.py --config input_HPC.yml --input_zip $INPUT_ZIP --output $OUTPUT_PATH --output_summary $OUTPUT_SUMMARY
