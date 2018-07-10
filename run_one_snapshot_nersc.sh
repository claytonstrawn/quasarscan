#!/bin/bash

#from VELA01 to VELA07, then VELA16 to VELA35
#10MpcBox_csf512_a0.020.d 
#from a0.330 to a0.500

SOURCEDIR1="/home/c/ceverino/VELA"
SOURCEDIR2="/project/projectdirs/mp363/VELA_ANALYSIS"

i=$2
a=$3

mkdir "VELA"$v$i
cd "VELA"$v$i

if [[ $1 = "1" ]]
then
    v=""
    DIRNAME=$SOURCEDIR2
    scp $DIRNAME/VELA$i/*a0.$a* .
else
    v="_v"$1"_"
    DIRNAME=$SOURCEDIR"_v"$1
    hsi "cd "$DIRNAME"/VELA"$i"; get *a0.$a*"
fi
cd ..



module load python
source activate myenv

python quasarscan/quasar_scan.py 'n' 'VELA'$v$i'/10MpcBox_csf512_a0.'$a'.d' 'VELA'$v$i '-p'

rm -rf 'VELA'$v$i
source deactivate
