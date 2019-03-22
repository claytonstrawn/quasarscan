#!/bin/bash

source myenv/bin/activate.csh
python quasarscan/create_qso_endpoints.py NIHAO_v1_tipsy_01 files_to_process/g2.63e10.00144
srun -n 10 -u python -u quasarscan/get_coldens.py quasarscan/output/NIHAO_v1_tipsy_01coldensinfo/0_of_448-agoraions_z3.0.txt 96 p
python quasarscan/create_qso_endpoints.py NIHAO_v1_tipsy_01 files_to_process/g2.63e10.00224
srun -n 10 -u python -u quasarscan/get_coldens.py quasarscan/output/NIHAO_v1_tipsy_01coldensinfo/0_of_448-agoraions_z2.0.txt 96 p
python quasarscan/create_qso_endpoints.py NIHAO_v1_tipsy_01 files_to_process/g2.63e10.00304
srun -n 10 -u python -u quasarscan/get_coldens.py quasarscan/output/NIHAO_v1_tipsy_01coldensinfo/0_of_448-agoraions_z1.5.txt 96 p
