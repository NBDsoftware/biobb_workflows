#!/bin/bash

rm -rf ./output
python biobb_docking_htvs.py --config input.yml --lig-lib ligand_lib.txt 