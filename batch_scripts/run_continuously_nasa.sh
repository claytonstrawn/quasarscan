#!/bin/bash

n=$1

module use -a /swbuild/analytix/tools/modulefiles
module load miniconda3/v4
source activate my_env

python quasarscan/decide_next_onlyVELA.py
until ! [ -s quasarscan/nextfile.sh ]
do
    chmod u+x quasarscan/nextfile.sh
    quasarscan/./nextfile.sh
    python quasarscan/decide_next_onlyVELA.py $n
done