#!/bin/bash

#from VELA01 to VELA07, then VELA16 to VELA35
#10MpcBox_csf512_a0.020.d 
#from a0.330 to a0.500

#SOURCEDIR1="/home/c/ceverino/VELA"
#SOURCEDIR2="/project/projectdirs/mp363/VELA_ANALYSIS"
module load python
source activate myenv
v="_v${1}_"
i=$2
a=$3
z=$4
n=$5
if [[ $n = '0' ]]
then
    python quasarscan/quasar_scan.py n /global/cscratch1/sd/cstrawn/VELA$v$i/10MpcBox_csf512_a0.$a.d VELA$v$i
fi
srun -n 96 -u python -u quasarscan/get_coldens.py -fn quasarscan/output2.0/VELA$v${i}coldensinfo/${n}_of_416-_OI_OII_OIII_OIV_OV_OVI_OVII_OVIII_OIX_HI_HII_z$z.txt -s 96 -p

echo "remember to rm -rf VELA$v$i"
