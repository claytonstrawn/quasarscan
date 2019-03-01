#!/bin/bash

#from VELA01 to VELA07, then VELA16 to VELA35
#10MpcBox_csf512_a0.020.d 
#from a0.330 to a0.500

#SOURCEDIR1="/home/c/ceverino/VELA"
#SOURCEDIR2="/project/projectdirs/mp363/VELA_ANALYSIS"
simname=$1
filename=$2
z=$3

python quasarscan/create_qso_endpoints.py /global/cscratch1/sd/cstrawn/$filename $simname
srun -n 32 -u python -u quasarscan/get_coldens.py quasarscan/output/${simname}coldensinfo/0_of_448-agoraions_z$z.txt 96 p

echo "remember to rm $filename"
