#!/bin/bash

INPUT_FOLDER=/home/nbd-pablo/repos/biobb_workflows/OTHER/5.SuMD/input
STRUCTURE=$INPUT_FOLDER/complex.gro
TOPOLOGY=$INPUT_FOLDER/topology.zip

python biobb_SuMD.py -input $STRUCTURE -topology $TOPOLOGY -config input.yml

