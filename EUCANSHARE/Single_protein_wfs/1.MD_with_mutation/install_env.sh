#!/bin/bash

# To install in HPC (Tirant)

module purge

. /storage/apps/MINICONDA/3/etc/profile.d/conda.sh

# conda deactivate

# conda env remove -n eucanshare_wf1

export KEY_MODELLER=HERE YOUR MODELLER KEY

conda env create -f environment.yml

