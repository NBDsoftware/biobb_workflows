#!/bin/bash

python biobb_md_setup_mutation.py -i input/P38alpha4LOO.pdb --n_trajs 2  --mut_list 'A:Y9V' --config input.yml  

# python biobb_md_setup_mutation.py -i pdb:1AKI --op free --mut_list 'A:F38C,A:N39W,A:T40G' --config input.yml  
