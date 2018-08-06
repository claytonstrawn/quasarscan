#!/bin/bash

#from VELA01 to VELA07, then VELA16 to VELA35
#10MpcBox_csf512_a0.020.d 
#from a0.330 to a0.500

SOURCEDIR1="/home/c/ceverino/VELA"
SOURCEDIR2="/project/projectdirs/mp363/VELA_ANALYSIS"

i=$2

if [[ $1 = "1" ]]
then
    v=""
    DIRNAME=$SOURCEDIR2
    mkdir VELA$v$i
    cd VELA$v$i
    scp $DIRNAME/VELA$i/*a0.$a* .
else
    v="_v"$1"_"
    DIRNAME=$SOURCEDIR1"_v"$1
    cd $SCRATCH
    mkdir VELA$v$i
    cd VELA$v$i
    #hsi "cd $DIRNAME/VELA$i; get *a0.$a*"
fi
cd

module load python
source activate myenv

quasarscan/quasar_scan.py n /global/cscratch1/sd/cstrawn/VELA$v$i/10MpcBox_csf512_a0.200.d VELA$v$i
srun -n 32 -u python -u quasarscan/get_coldens.py -fn quasarscan/output/VELA$v${i}coldensinfo/0_of_400-allions_z4.0.txt -s 96
quasarscan/quasar_scan.py n /global/cscratch1/sd/cstrawn/VELA$v$i/10MpcBox_csf512_a0.250.d VELA$v$i
srun -n 32 -u python -u quasarscan/get_coldens.py -fn quasarscan/output/VELA$v${i}coldensinfo/0_of_400-allions_z3.0.txt -s 96
quasarscan/quasar_scan.py n /global/cscratch1/sd/cstrawn/VELA$v$i/10MpcBox_csf512_a0.330.d VELA$v$i
srun -n 32 -u python -u quasarscan/get_coldens.py -fn quasarscan/output/VELA$v${i}coldensinfo/0_of_400-allions_z2.0.txt -s 96
quasarscan/quasar_scan.py n /global/cscratch1/sd/cstrawn/VELA$v$i/10MpcBox_csf512_a0.400.d VELA$v$i
srun -n 32 -u python -u quasarscan/get_coldens.py -fn quasarscan/output/VELA$v${i}coldensinfo/0_of_400-allions_z1.5.txt -s 96
quasarscan/quasar_scan.py n /global/cscratch1/sd/cstrawn/VELA$v$i/10MpcBox_csf512_a0.500.d VELA$v$i
srun -n 32 -u python -u quasarscan/get_coldens.py -fn quasarscan/output/VELA$v${i}coldensinfo/0_of_400-allions_z1.0.txt -s 96

echo "remember to rm -rf VELA$v$i"
source deactivate
