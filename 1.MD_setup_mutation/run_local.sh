#!/bin/bash

rm -rf output

python biobb_md_setup_mutation.py -i pdb:4LOO -o output --op free --mut_list 'A:Y9V' --config input.yml  

# python biobb_md_setup_mutation.py -i pdb:1AKI -o output --op free --mut_list 'A:F38C,A:N39W,A:T40G' --config input.yml  
