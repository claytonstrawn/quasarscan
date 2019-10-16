#!/bin/bash

simname=$1
filename=$2
z=$3
n=$4

if (($n == '0')); then
    python quasarscan/create_qso_endpoints.py $simname $filename
fi
srun -n 128 -u python -u quasarscan/get_coldens.py quasarscan/output/${simname}coldensinfo/${n}_of_384-agoraions_z${z}.txt 128 p

#in case you want to run one line interactively:
#srun -n 32 -u python -u quasarscan/get_coldens.py quasarscan/output/AGORA_v1_art_01coldensinfo/96_of_448-agoraions_z15.0.txt 96 p
#python quasarscan/get_coldens.py quasarscan/output/${simname}coldensinfo/0_of_448-agoraions_z$z.txt 96 p
#not in parallel:
#python quasarscan/get_coldens.py quasarscan/output/NIHAO_v1_tipsy_02coldensinfo/0_of_448-agoraions_z0.31.txt 2 n
#python quasarscan/get_coldens.py quasarscan/output/${simname}coldensinfo/${n}_of_448-agoraions_z$z.txt 96 np

echo "remember to rm $filename"
