#!/bin/bash

simname=$1
filename=$2
z=$3

python quasarscan/create_qso_endpoints.py $simname $filename
srun -n 32 -u python -u quasarscan/get_coldens.py quasarscan/output/${simname}coldensinfo/0_of_448-agoraions_z$z.txt 96 p
#python quasarscan/get_coldens.py quasarscan/output/${simname}coldensinfo/0_of_448-agoraions_z$z.txt 96 p

echo "remember to rm $filename"
