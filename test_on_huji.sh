#!/bin/bash

srun -n 32 -u python -u quasarscan/get_coldens.py quasarscan/output/NIHAO_v1_tipsy_02coldensinfo/82_of_448-agoraions_z0.43.txt 96 p
#in case you want to run one line interactively:
#srun -n 32 -u python -u quasarscan/get_coldens.py quasarscan/output/AGORA_v1_art_01coldensinfo/96_of_448-agoraions_z15.0.txt 96 p
#python quasarscan/get_coldens.py quasarscan/output/${simname}coldensinfo/0_of_448-agoraions_z$z.txt 96 p
#not in parallel:
#python quasarscan/get_coldens.py quasarscan/output/NIHAO_v1_tipsy_02coldensinfo/0_of_448-agoraions_z0.31.txt 2 n
#python quasarscan/get_coldens.py quasarscan/output/${simname}coldensinfo/0_of_448-agoraions_z$z.txt 96 p

echo "remember to rm $filename"
