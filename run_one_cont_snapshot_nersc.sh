#!/bin/bash

i=$2
z=$3

v="_v"$1"_"
cd "VELA"$v$i
ls
cd ..

module load python
source activate myenv
echo VELA$v$i

python quasarscan/quasar_scan.py 'c' '-sz' VELA$v$i $z '-p' '-s' 12

echo "be sure to rm -rf VELA$v$i"
source deactivate
