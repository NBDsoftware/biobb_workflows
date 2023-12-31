#!/bin/bash
#SBATCH --job-name=SuMD_analysis
#SBATCH --ntasks=5                                      # total number of tasks across all nodes
#SBATCH --nodes=1
#SBATCH --time=00:40:00
#SBATCH --mem-per-cpu=2000
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err

# Purge loaded modules
module purge 

# Load conda
module load Miniconda3/4.9.2

# Load vmd
module load VMD/1.9.4a57-intel-2021b

# Activate previously created conda environment
source activate /home/pnavarro/.conda/envs/sumd_analyzer

# Export the right PYTHONPATH
export PYTHONPATH=/home/pnavarro/.conda/envs/sumd_analyzer/lib/python3.10/site-packages/
export PATH=/home/pnavarro/.conda/envs/sumd_analyzer/bin:$PATH

# Paths
INPUT_YML=input_analysis.yml
PATH_TO_SUMDANALYZER=/shared/work/BiobbWorkflows/src/biobb_workflows/Other/SuMD-analyzer

# Launch Geometry analysis
python $PATH_TO_SUMDANALYZER/RNASuMDAnalyzer.py geometry $INPUT_YML

# Launch MMGBSA analysis
python $PATH_TO_SUMDANALYZER/RNASuMDAnalyzer.py mmgbsa $INPUT_YML

# Launch Interaction energy (NAMD) analysis
python $PATH_TO_SUMDANALYZER/RNASuMDAnalyzer.py intEnergy $INPUT_YML