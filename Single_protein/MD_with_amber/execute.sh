#!/bin/bash -l
#SBATCH -A p200581
#SBATCH --job-name="execution"
#SBATCH -p gpu
#SBATCH -q default
#SBATCH --time=2-00:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-task=1
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err

pdb_code=$1

sbatch run_amber_auto.sh $pdb_code free

sleep 30

sbatch run_amber_auto.sh $pdb_code free

sleep 30

sbatch run_amber_auto.sh $pdb_code gamd

sleep 30

sbatch run_amber_auto.sh $pdb_code gamd

sleep 30

sbatch run_amber_auto.sh $pdb_code gamd

sleep 30

sbatch run_amber_auto.sh $pdb_code gamd