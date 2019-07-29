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
cp /vol/sci/astro/cosmo/nas2/Data/nihao/$d/$d.param files_to_process/$n
cp /vol/sci/astro/cosmo/nas2/Data/nihao/$d/$d.00112.* files_to_process/$n
cp /vol/sci/astro/cosmo/nas2/Data/nihao/$d/$d.00112 files_to_process/$n
cp /vol/sci/astro/cosmo/nas2/Data/nihao/$d/$d.00160.* files_to_process/$n
cp /vol/sci/astro/cosmo/nas2/Data/nihao/$d/$d.00160 files_to_process/$n
cp /vol/sci/astro/cosmo/nas2/Data/nihao/$d/$d.00240.* files_to_process/$n
cp /vol/sci/astro/cosmo/nas2/Data/nihao/$d/$d.00240 files_to_process/$n
cp /vol/sci/astro/cosmo/nas2/Data/nihao/$d/$d.00320.* files_to_process/$n
cp /vol/sci/astro/cosmo/nas2/Data/nihao/$d/$d.00320 files_to_process/$n
cp /vol/sci/astro/cosmo/nas2/Data/nihao/$d/$d.00432.* files_to_process/$n
cp /vol/sci/astro/cosmo/nas2/Data/nihao/$d/$d.00432 files_to_process/$n
cp /vol/sci/astro/cosmo/nas2/Data/nihao/$d/$d.00640.* files_to_process/$n
cp /vol/sci/astro/cosmo/nas2/Data/nihao/$d/$d.00640 files_to_process/$n
cp /vol/sci/astro/cosmo/nas2/Data/nihao/$d/$d.01024.* files_to_process/$n
cp /vol/sci/astro/cosmo/nas2/Data/nihao/$d/$d.01024 files_to_process/$n
