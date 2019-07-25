#!/bin/bash

echo """'g2.63e10': 'NIHAO_v1_tipsy_01',
'g2.19e11': 'NIHAO_v1_tipsy_02',
'g5.02e11': 'NIHAO_v1_tipsy_03',
'g5.31e11': 'NIHAO_v1_tipsy_04',
'g5.36e11': 'NIHAO_v1_tipsy_05',
'g5.38e11': 'NIHAO_v1_tipsy_06',
'g6.96e11': 'NIHAO_v1_tipsy_07',
'g7.08e11': 'NIHAO_v1_tipsy_08',
'g7.44e11': 'NIHAO_v1_tipsy_09',
'g7.55e11': 'NIHAO_v1_tipsy_10',
'g7.66e11': 'NIHAO_v1_tipsy_11',
'g8.06e11': 'NIHAO_v1_tipsy_12',
'g8.13e11': 'NIHAO_v1_tipsy_13',
'g8.26e11': 'NIHAO_v1_tipsy_14',
'g1.12e12': 'NIHAO_v1_tipsy_15',
'g1.77e12': 'NIHAO_v1_tipsy_16',
'g1.92e12': 'NIHAO_v1_tipsy_17',
'g2.79e12': 'NIHAO_v1_tipsy_18'"""

d=$1
num=$2
if [ -z "$num" ]
then
        exit 1
else
        n=NIHAO_v1_tipsy_$num
fi



mkdir files_to_process/$n
cp /vol/sci/astro/cosmo/nas2/Data/nihao/$d/$d.00112.* files_to_process/$n
cp /vol/sci/astro/cosmo/nas2/Data/nihao/$d/$d.00160.* files_to_process/$n
cp /vol/sci/astro/cosmo/nas2/Data/nihao/$d/$d.00192.* files_to_process/$n
cp /vol/sci/astro/cosmo/nas2/Data/nihao/$d/$d.00240.* files_to_process/$n
cp /vol/sci/astro/cosmo/nas2/Data/nihao/$d/$d.00320.* files_to_process/$n
cp /vol/sci/astro/cosmo/nas2/Data/nihao/$d/$d.00432.* files_to_process/$n
cp /vol/sci/astro/cosmo/nas2/Data/nihao/$d/$d.00640.* files_to_process/$n
python quasarscan/create_qso_endpoints.py $n files_to_process/$n/$d.00112
python quasarscan/get_coldens.py quasarscan/output/${n}coldensinfo/0_of_448-agoraions_z4.0.txt 96 n
python quasarscan/create_qso_endpoints.py $n files_to_process/$n/$d.00160
python quasarscan/get_coldens.py quasarscan/output/${n}coldensinfo/0_of_448-agoraions_z3.0.txt 96 n
python quasarscan/create_qso_endpoints.py $n files_to_process/$n/$d.00192
python quasarscan/get_coldens.py quasarscan/output/${n}coldensinfo/0_of_448-agoraions_z2.0.txt 96 n
python quasarscan/create_qso_endpoints.py $n files_to_process/$n/$d.00240
python quasarscan/get_coldens.py quasarscan/output/${n}coldensinfo/0_of_448-agoraions_z1.0.txt 96 n
python quasarscan/create_qso_endpoints.py $n files_to_process/$n/$d.00320
python quasarscan/get_coldens.py quasarscan/output/${n}coldensinfo/0_of_448-agoraions_z1.5.txt 96 n
python quasarscan/create_qso_endpoints.py $n files_to_process/$n/$d.00432
python quasarscan/get_coldens.py quasarscan/output/${n}coldensinfo/0_of_448-agoraions_z0.5.txt 96 n
python quasarscan/create_qso_endpoints.py $n files_to_process/$n/$d.00640
python quasarscan/get_coldens.py quasarscan/output/${n}coldensinfo/0_of_448-agoraions_z0.0.txt 96 n
