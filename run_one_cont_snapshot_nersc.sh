#!/bin/bash

i=$2
z=$3
srun=$4

v="_v"$1"_"
cd "VELA"$v$i
ls
cd ..

module load python
source activate myenv
echo VELA$v$i
if [[ $srun = "srun" ]]
then
    srun -c 50 python quasarscan/quasar_scan.py 'c' '-sz' VELA$v$i $z '-p' '-s' 36
else
    python quasarscan/quasar_scan.py 'c' '-sz' VELA$v$i $z '-p' '-s' 36
fi

echo "be sure to rm -rf VELA$v$i"
source deactivate

